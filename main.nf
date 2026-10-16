#!/usr/bin/env nextflow

// Define pipeline parameters with default values
params.methid_map      = false
params.input_dir       = false
params.sample_sheet    = false
params.ref_genome      = false
params.min_cov         = 2
params.outdir          = 'results'
params.logdir          = "${params.outdir}/logs"

// Sequencing platform, as in smMethID:
//   'ont'   - ONT modBAMs (+ .bai), run through modkit pileup
//   '5base' - 5CmC / Rastair mod_pileup.txt (Illumina 5-base or Watchmaker data)
params.platform        = false
// Optional: restart from an intermediate instead of the platform's raw input:
//   'shortbed'  - <outdir>/shortbeds from an earlier run (either platform)
//   'bedmethyl' - modkit pileup output (ont only)
params.start_from      = false
// Optional override for how input files are matched (defaults set per type below)
params.input_glob      = false

// Which modifications modkit counts as methylated (ont only):
// '5mC_5hmC' combines 5mC + 5hmC to match WGBS-based maps; '5mC' counts 5mC alone
params.mod_bases       = '5mC_5hmC'

// Balanced marker selection, as in smMethID: 0 = every atlas CpG; N = for each cell type, the N CpGs
// where it differs most from its closest other cell type, pooled. marker_ranking 'standardized'
// divides that difference by the replicate SD and needs var_atlas_file (build_atlas.nf output).
params.markers_per_cell_type = 0
params.marker_ranking        = 'gap'
params.var_atlas_file        = false

// Optional metadata sheet (CSV, or TSV if named .tsv/.txt): the sample sheet's columns plus 'label' and
// 'exposure'. When given, results are also summarized per cell type by label and by label x exposure.
params.metadata_sheet  = false

// --- Input type handling ---

if (!(params.platform in ['ont', '5base'])) {
    exit 1, "Set platform to 'ont' or '5base' in params.json (--platform)"
}
if (params.start_from && !(params.start_from in ['shortbed', 'bedmethyl'])) {
    exit 1, "Invalid --start_from '${params.start_from}'. Must be 'shortbed', 'bedmethyl', or false"
}
if (params.start_from == 'bedmethyl' && params.platform != 'ont') {
    exit 1, "start_from 'bedmethyl' is only for platform 'ont' (modkit output)"
}

// What the input files are: the platform's raw input, unless restarting from an intermediate
def input_type = params.start_from ?: (params.platform == 'ont' ? 'bam' : 'mod_pileup')

def valid_mod_bases = ['5mC_5hmC', '5mC']
if (!valid_mod_bases.contains(params.mod_bases)) {
    exit 1, "Invalid --mod_bases '${params.mod_bases}'. Must be one of: ${valid_mod_bases.join(', ')}"
}

// Default file patterns per input type
def default_globs = [
    bam       : '*.bam',
    bedmethyl : '*.bed',
    shortbed  : '*_short.bed',
    // 5CmC writes one <sample>...__5cmc-vX/ directory per sample, each with a mod_pileup.txt
    mod_pileup: '*/mod_pileup.txt'
]
def input_glob = params.input_glob ?: default_globs[input_type]

// Create necessary directories
new File(params.outdir).mkdirs()
new File(params.logdir).mkdirs()

// Log parameters for reproducibility
def summary = """
         M E T H Y L A T I O N   D E C O N V O L U T I O N
         ================================================
         Atlas Map       : ${params.methid_map}
         Sample Sheet    : ${params.sample_sheet}
         Input Dir       : ${params.input_dir}
         Platform        : ${params.platform}
         Start From      : ${params.start_from ?: '(raw input)'}
         Input Type      : ${input_type}
         Input Pattern   : ${params.sample_sheet ? '(from sample sheet)' : input_glob}
         Mod Bases       : ${input_type == 'bam' ? params.mod_bases : '(n/a)'}
         Ref Genome      : ${params.ref_genome}
         Output Dir      : ${params.outdir}
         Min Coverage    : ${params.min_cov}
         Markers/cell type: ${params.markers_per_cell_type ?: 'all CpGs'}
         Marker ranking  : ${params.marker_ranking}
         Metadata Sheet  : ${params.metadata_sheet}
         """
log.info summary
new File("${params.outdir}/manifest.txt").write(summary)

// Input Validation
if (!params.methid_map) { exit 1, "Methylation ID map not specified! --methid_map <file>" }
// Check the metadata sheet up front, so a typo doesn't surface only after every sample is deconvolved
if (params.metadata_sheet) {
    def metadata_file = file(params.metadata_sheet)
    if (!metadata_file.exists()) { exit 1, "Metadata sheet not found: ${params.metadata_sheet}" }
    def metadata_sep = metadata_file.name.toLowerCase() =~ /\.(tsv|txt)$/ ? '\t' : ','
    def metadata_header = metadata_file.readLines()[0].split(metadata_sep).collect { it.trim().replaceAll(/^"|"$/, '') }
    def metadata_missing = ['CGID', 'label', 'exposure'].findAll { !(it in metadata_header) }
    if (metadata_missing) {
        exit 1, "Metadata sheet ${params.metadata_sheet} needs columns CGID, label and exposure. Missing: ${metadata_missing.join(', ')}"
    }
}
if (!params.input_dir && !params.sample_sheet) {
    exit 1, "No inputs specified! Use --input_dir <dir> or --sample_sheet <file>"
}
// Reference genome is only needed when we actually run modkit pileup
if (!(params.marker_ranking in ['gap', 'standardized'])) { exit 1, "marker_ranking must be 'gap' or 'standardized'" }
if (params.marker_ranking == 'standardized' && params.markers_per_cell_type && !params.var_atlas_file) {
    exit 1, "marker_ranking 'standardized' needs var_atlas_file (atlas_beta_binomial_var.tsv from build_atlas.nf)"
}
if (input_type == 'bam' && !params.ref_genome) {
    exit 1, "Reference genome not specified! --ref_genome <file> (required for platform 'ont')"
}


// --- Build the input channel according to the input file type ---

// Derive a sample id from an input file name when no sample sheet is given
def sampleIdFromFile(f, input_type) {
    def sample_id = f.baseName
    // Strip the _short suffix so sample ids match those produced from BAM input
    if (input_type == 'shortbed') {
        sample_id = sample_id.replaceAll(/_short$/, '')
    }
    // 5CmC files are all named mod_pileup.txt, so take the sample id from the parent
    // directory and drop the 5CmC and CG-Align suffixes
    // (e.g. SampleABC_bt2-bwa_umi_true_hg19_wm_idt_spike_v1.0__5cmc-v1.1/mod_pileup.txt -> SampleABC)
    if (input_type == 'mod_pileup') {
        def base = f.name == 'mod_pileup.txt' ? f.parent.name : f.name.replaceAll(/_?mod_pileup\.txt$/, '')
        sample_id = base
            .replaceAll(/__5cmc.*$/, '')
            .replaceAll(/_(bt2-bwa|bt2|bwa)_(umi|std|none)_(true|false)_.*$/, '')
    }
    // Sanitize for use in file names by replacing spaces and slashes
    return sample_id.replaceAll("[ /]", "_")
}

if (params.sample_sheet) {

    // Sample sheet, as in smMethID: CSV (or TSV if named .tsv/.txt) with a CGID column and a path
    // column that depends on the input: bam_path (ont), sample_dir (5base: the 5CmC sample folder
    // holding mod_pileup.txt), or path (start_from: the shortbed / bedmethyl file). Samples are named by CGID.
    def path_col = params.start_from ? 'path' : (params.platform == 'ont' ? 'bam_path' : 'sample_dir')
    def sep = params.sample_sheet.toString().toLowerCase() =~ /\.(tsv|txt)$/ ? '\t' : ','
    Channel
        .fromPath(params.sample_sheet, checkIfExists: true)
        .splitCsv(header: true, sep: sep)
        .map { row ->
            if (!row.containsKey('CGID') || !row.containsKey(path_col)) {
                exit 1, "Sample sheet ${params.sample_sheet} needs columns CGID and ${path_col}. Columns present: ${row.keySet().join(', ')}"
            }
            def sample_id = row.CGID.trim().replaceAll("[ /]", "_")
            def f = file(row[path_col].trim())
            if (input_type == 'mod_pileup') {
                if (!f.isDirectory()) { exit 1, "5CmC sample dir not found for ${sample_id}: ${f}" }
                f = f.resolve('mod_pileup.txt')
            }
            if (!f.exists()) {
                exit 1, "Input file not found for sample ${sample_id}: ${f}"
            }
            return [ sample_id, f ]
        }
        .set { sample_files_ch }

} else {

    Channel
        .fromPath("${params.input_dir}/${input_glob}", checkIfExists: true)
        .map { f -> [ sampleIdFromFile(f, input_type), f ] }
        .set { sample_files_ch }
}

// BAMs travel with their index
if (input_type == 'bam') {
    sample_files_ch
        .map { sample_id, bam ->
            def bai = file("${bam}.bai")
            if (!bai.exists()) {
                exit 1, "BAM index not found for ${bam.name}. Please run 'samtools index ${bam.name}'"
            }
            return [ sample_id, [bam, bai] ]
        }
        .set { input_ch }
} else {
    input_ch = sample_files_ch
}


// --- Processes ---

process MODKIT_PILEUP {
    tag "Modkit pileup on ${sample_id}"
    publishDir "${params.outdir}/bedmethyl", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_and_bai)

    output:
    tuple val(sample_id), path("${sample_id}.bed")

    script:
    def bam_file = bam_and_bai[0]
    def mod_args = params.mod_bases == '5mC' ? '--modified-bases 5mC' : '--modified-bases 5mC 5hmC --combine-mods'
    """
    echo "Modkit version: \$(modkit --version)"
    modkit pileup "${bam_file}" "${sample_id}.bed" \\
        --ref "${params.ref_genome}" \\
        ${mod_args} \\
        --cpg \\
        --combine-strands
    """
}

process BEDMETHYL_TO_SHORTBED {
    tag "Convert to shortbed for ${sample_id}"
    publishDir "${params.outdir}/shortbeds", mode: 'copy'

    input:
    tuple val(sample_id), path(bedmethyl_file)

    output:
    tuple val(sample_id), path("${sample_id}_short.bed")

    script:
    """
    awk '{ print \$1, \$2, \$3, \$12, \$10 }' "${bedmethyl_file}" > "${sample_id}_short.bed"
    """
}

process MOD_PILEUP_TO_SHORTBED {
    tag "Convert 5CmC mod_pileup to shortbed for ${sample_id}"
    publishDir "${params.outdir}/shortbeds", mode: 'copy'

    input:
    tuple val(sample_id), path(pileup_file)

    output:
    tuple val(sample_id), path("${sample_id}_short.bed")

    // Rastair mod_pileup.txt is one row per (CpG, strand) with header
    //   #chr start end name beta_est strand unmod mod no_snp snp coverage genotype ...
    // `start` is 0-based like modkit bedMethyl; the '-' strand row sits at start+1.
    // Drop rows where the sample genotype isn't the reference C (C/C on +, G/G on -), since
    // C>T SNPs look like methylation calls, then shift '-' rows back 1bp and sum both strands
    // to match `modkit pileup --cpg --combine-strands`. Output: chr start end n_mod n_valid.
    script:
    """
    awk 'BEGIN { FS = OFS = "\\t" }
         /^#/ { next }
         (\$6 == "+" && \$12 == "C/C") || (\$6 == "-" && \$12 == "G/G") {
           start = (\$6 == "-") ? \$2 - 1 : \$2
           key = \$1 OFS start
           mod[key] += \$8
           cov[key] += \$7 + \$8
         }
         END { for (k in mod) if (cov[k] > 0) { split(k, a, OFS); print a[1], a[2], a[2] + 1, mod[k], cov[k] } }' \\
      "${pileup_file}" \\
      | sort -k1,1 -k2,2n > "${sample_id}_short.bed"
    """
}

process DECONVOLVE {
    tag "Deconvolve ${sample_id}"
    publishDir "${params.outdir}/deconvolution/per_sample", mode: 'copy'

    input:
    tuple val(sample_id), path(shortbed_file)
    path(map_file)

    output:
    path("${sample_id}_deconvolution.csv"), emit: deconv
    path("${sample_id}_coverage.tsv"),      emit: cov

    script:
    // Read off the shared filesystem rather than staged, like ref_genome
    def var_arg = params.var_atlas_file ? "--var_atlas_file \"${file(params.var_atlas_file)}\"" : ""
    """
    Rscript ${baseDir}/bin/deconvolve.r \\
        --map_file "${map_file}" \\
        --shortbed "${shortbed_file}" \\
        --sample_name "${sample_id}" \\
        --output_dir "." \\
        --threads "${task.cpus}" \\
        --min_cov "${params.min_cov}" \\
        --markers_per_cell_type "${params.markers_per_cell_type}" \\
        --marker_ranking "${params.marker_ranking}" ${var_arg}
    """
}

process COMBINE {
    tag "Combine deconvolution results"
    publishDir "${params.outdir}/deconvolution", mode: 'copy'

    input:
    path(deconv_files)
    path(cov_files)
    path(map_file)

    output:
    path("deconvolved_samples.csv"), emit: deconv
    path("coverage_*.tsv")
    path("manifest_*.txt")

    script:
    """

    Rscript ${baseDir}/bin/combine.r \\
        --input_dir "." \\
        --output_dir "." \\
        --map_file "${map_file}"
    """
}

process SUMMARIZE {
    tag "Summarize by label and exposure"
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path(deconv_file)
    path(metadata_file)

    output:
    path("summary_*.csv")
    path("deconvolved_samples_with_metadata.csv")

    script:
    """
    Rscript ${baseDir}/bin/summarize.r \\
        --deconv_file "${deconv_file}" \\
        --metadata_file "${metadata_file}" \\
        --output_dir "."
    """
}

// Define the workflow logic
workflow {

    map_ch = file(params.methid_map)   // value channel, reused for every sample

    // Stage 1: Run modkit on every input BAM (skipped unless input is BAM)
    if (input_type == 'bam') {
        MODKIT_PILEUP(input_ch)
        bedmethyl_ch = MODKIT_PILEUP.out
    }
    else {
        bedmethyl_ch = input_ch
    }

    // Stage 2: Convert bedmethyl files to shortbed format (skipped if input is already shortbed)
    if (input_type == 'shortbed') {
        shortbed_ch = input_ch
    }
    else if (input_type == 'mod_pileup') {
        MOD_PILEUP_TO_SHORTBED(input_ch)
        shortbed_ch = MOD_PILEUP_TO_SHORTBED.out
    }
    else {
        BEDMETHYL_TO_SHORTBED(bedmethyl_ch)
        shortbed_ch = BEDMETHYL_TO_SHORTBED.out
    }

    // Stage 3: Deconvolve each sample in parallel (one task per sample)
    DECONVOLVE(shortbed_ch, map_ch)

    // Stage 4: Combine all per-sample results into a single table
    COMBINE(
        DECONVOLVE.out.deconv.collect(),
        DECONVOLVE.out.cov.collect(),
        map_ch
    )

    // Stage 5 (optional): summarize by label and by label x exposure
    if (params.metadata_sheet) {
        SUMMARIZE(COMBINE.out.deconv, file(params.metadata_sheet))
    }
}
