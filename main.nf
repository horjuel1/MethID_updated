#!/usr/bin/env nextflow

// Define pipeline parameters with default values
params.methid_map      = false
params.input_dir       = false
params.sample_sheet    = false
params.ref_genome      = false
params.min_cov         = 2
params.outdir          = 'results'
params.logdir          = "${params.outdir}/logs"

// Which stage to start from: 'bam', 'bedmethyl', 'shortbed', or 'mod_pileup'
// ('mod_pileup' = 5CmC / Rastair mod_pileup.txt from Illumina 5-base or Watchmaker data)
params.input_file_type = 'bam'
// Optional override for how input files are matched (defaults set per type below)
params.input_glob      = false

// Which modifications modkit counts as methylated (BAM input only):
// '5mC_5hmC' combines 5mC + 5hmC to match WGBS-based maps; '5mC' counts 5mC alone
params.mod_bases       = '5mC_5hmC'

// Balanced marker selection, as in smMethID: 0 = every atlas CpG; N = for each cell type, the N CpGs
// where it differs most from its closest other cell type, pooled. marker_ranking 'standardized'
// divides that difference by the replicate SD and needs var_atlas_file (build_atlas.nf output).
params.markers_per_cell_type = 0
params.marker_ranking        = 'gap'
params.var_atlas_file        = false

// --- Input type handling ---

def valid_types = ['bam', 'bedmethyl', 'shortbed', 'mod_pileup']
if (!valid_types.contains(params.input_file_type)) {
    exit 1, "Invalid --input_file_type '${params.input_file_type}'. Must be one of: ${valid_types.join(', ')}"
}

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
def input_glob = params.input_glob ?: default_globs[params.input_file_type]

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
         Input File Type : ${params.input_file_type}
         Input Pattern   : ${params.sample_sheet ? '(from sample sheet)' : input_glob}
         Mod Bases       : ${params.input_file_type == 'bam' ? params.mod_bases : '(n/a)'}
         Ref Genome      : ${params.ref_genome}
         Output Dir      : ${params.outdir}
         Min Coverage    : ${params.min_cov}
         Markers/cell type: ${params.markers_per_cell_type ?: 'all CpGs'}
         Marker ranking  : ${params.marker_ranking}
         """
log.info summary
new File("${params.outdir}/manifest.txt").write(summary)

// Input Validation
if (!params.methid_map) { exit 1, "Methylation ID map not specified! --methid_map <file>" }
if (!params.input_dir && !params.sample_sheet) {
    exit 1, "No inputs specified! Use --input_dir <dir> or --sample_sheet <file>"
}
// Reference genome is only needed when we actually run modkit pileup
if (!(params.marker_ranking in ['gap', 'standardized'])) { exit 1, "marker_ranking must be 'gap' or 'standardized'" }
if (params.marker_ranking == 'standardized' && params.markers_per_cell_type && !params.var_atlas_file) {
    exit 1, "marker_ranking 'standardized' needs var_atlas_file (atlas_beta_binomial_var.tsv from build_atlas.nf)"
}
if (params.input_file_type == 'bam' && !params.ref_genome) {
    exit 1, "Reference genome not specified! --ref_genome <file> (required when --input_file_type bam)"
}


// --- Build the input channel according to the input file type ---

// Derive a sample id from an input file name when no sample sheet is given
def sampleIdFromFile(f, input_file_type) {
    def sample_id = f.baseName
    // Strip the _short suffix so sample ids match those produced from BAM input
    if (input_file_type == 'shortbed') {
        sample_id = sample_id.replaceAll(/_short$/, '')
    }
    // 5CmC files are all named mod_pileup.txt, so take the sample id from the parent
    // directory and drop the 5CmC and CG-Align suffixes
    // (e.g. SampleABC_bt2-bwa_umi_true_hg19_wm_idt_spike_v1.0__5cmc-v1.1/mod_pileup.txt -> SampleABC)
    if (input_file_type == 'mod_pileup') {
        def base = f.name == 'mod_pileup.txt' ? f.parent.name : f.name.replaceAll(/_?mod_pileup\.txt$/, '')
        sample_id = base
            .replaceAll(/__5cmc.*$/, '')
            .replaceAll(/_(bt2-bwa|bt2|bwa)_(umi|std|none)_(true|false)_.*$/, '')
    }
    // Sanitize for use in file names by replacing spaces and slashes
    return sample_id.replaceAll("[ /]", "_")
}

if (params.sample_sheet) {

    // Sample sheet (TSV) with 'id', 'label', and the input file path in a 'path' column
    // ('bam_path' is also accepted). Samples are named <id>_<label>.
    Channel
        .fromPath(params.sample_sheet, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def input_path = row.path ?: row.bam_path
            if (!row.id || !row.label || !input_path) {
                exit 1, "Sample sheet must contain 'id', 'label', and 'path' (or 'bam_path') columns. Found: ${row.keySet().join(', ')}"
            }
            def f = file(input_path)
            if (!f.exists()) {
                exit 1, "Input file not found for sample ${row.id}: ${input_path}"
            }
            def sample_id = "${row.id}_${row.label}".replaceAll("[ /]", "_")
            return [ sample_id, f ]
        }
        .set { sample_files_ch }

} else {

    Channel
        .fromPath("${params.input_dir}/${input_glob}", checkIfExists: true)
        .map { f -> [ sampleIdFromFile(f, params.input_file_type), f ] }
        .set { sample_files_ch }
}

// BAMs travel with their index
if (params.input_file_type == 'bam') {
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
    path("deconvolved_samples.csv")
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

// Define the workflow logic
workflow {

    map_ch = file(params.methid_map)   // value channel, reused for every sample

    // Stage 1: Run modkit on every input BAM (skipped unless input is BAM)
    if (params.input_file_type == 'bam') {
        MODKIT_PILEUP(input_ch)
        bedmethyl_ch = MODKIT_PILEUP.out
    }
    else {
        bedmethyl_ch = input_ch
    }

    // Stage 2: Convert bedmethyl files to shortbed format (skipped if input is already shortbed)
    if (params.input_file_type == 'shortbed') {
        shortbed_ch = input_ch
    }
    else if (params.input_file_type == 'mod_pileup') {
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
}
