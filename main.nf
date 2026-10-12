#!/usr/bin/env nextflow

// Define pipeline parameters with default values
params.methid_map   = false
params.sample_sheet = false
params.ref_genome   = false
params.min_cov      = 2
params.threads      = 4
params.outdir       = 'results'
params.logdir       = "${params.outdir}/logs"
// params.keep_bed_files = true

// Create necessary directories
new File(params.outdir).mkdirs()
new File(params.logdir).mkdirs()

// Log parameters for reproducibility
def summary = """
         M E T H Y L A T I O N   D E C O N V O L U T I O N
         ================================================
         Atlas Map    : ${params.methid_map}
         Sample Sheet : ${params.sample_sheet}
         Ref Genome   : ${params.ref_genome}
         Output Dir   : ${params.outdir}
         Min Coverage : ${params.min_cov}
         """
log.info summary
new File("${params.outdir}/manifest.txt").write(summary)

// Input Validation
if (!params.methid_map)   { exit 1, "Methylation ID map not specified! --methid_map <file>" }
if (!params.sample_sheet) { exit 1, "Sample sheet not specified! --sample_sheet <file>" }
if (!params.ref_genome)   { exit 1, "Reference genome not specified! --ref_genome <file>" }


// Create the input channel from the sample sheet
Channel
    .fromPath(params.sample_sheet)
    .splitCsv(header: true, sep: '\t')
    .map { row ->
        // Ensure required columns exist
        if (!row.id || !row.label || !row.bam_path) {
            exit 1, "Sample sheet must contain 'id', 'label', and 'bam_path' columns. Found: ${row.keySet().join(', ')}"
        }
        def bam = file(row.bam_path)
        def bai = file(row.bam_path + ".bai")
        if (!bam.exists()) {
            exit 1, "BAM file not found for sample ${row.id}: ${row.bam_path}"
        }
        if (!bai.exists()) {
            exit 1, "BAM index not found for ${bam.name}. Please run 'samtools index ${bam.name}'"
        }
        // Sanitize id and label for use in file names by replacing spaces and slashes
        def cleaned_id = row.id.replaceAll("[ /]", "_")
        def cleaned_label = row.label.replaceAll("[ /]", "_")
        def sample_id = "${cleaned_id}_${cleaned_label}"
        return [ sample_id, [bam, bai] ]
    }
    .set { bam_files_ch }


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
    """
    echo "Modkit version:"
    modkit --version
    modkit pileup "${bam_file}" "${sample_id}.bed" \\
        --modified-bases 5mC 5hmC \\
        --cpg \\
        --reference "${params.ref_genome}" \\
        --combine-mods \\
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

process DECONVOLVE {
    tag "Deconvolve all samples"
    publishDir "${params.outdir}/deconvolution", mode: 'copy'

    input:
    path(shortbed_files) // Collect to ensure this step runs after all shortbeds are made
    path(map_file)
    path(sample_sheet)

    output:
    path("deconvolved_samples.csv")
    path("coverage_*.tsv")
    path("manifest_*.txt")

    script:
    """
    #Disable fail on unbound variable setting
    set +u
    Rscript ${baseDir}/bin/deconvolve.r \\
        --map_file "${map_file}" \\
        --output_dir "." \\
        --sample_sheet "${sample_sheet}" \\
        --input_dir "${params.outdir}/shortbeds" \\
        --threads "${task.cpus}" \\
        --min_cov "${params.min_cov}"
    """
}

// 5. Define the workflow logic
workflow {
    // Stage 1: Run modkit on all input BAMs defined in the sample sheet
    MODKIT_PILEUP(bam_files_ch)

    // Stage 2: Convert bedmethyl files to shortbed format
    BEDMETHYL_TO_SHORTBED(MODKIT_PILEUP.out)

    // Stage 3: Collect all shortbed files and run deconvolution
    DECONVOLVE(
        BEDMETHYL_TO_SHORTBED.out.map { it[1] }.collect(), // Collect file paths to ensure completion
        file(params.methid_map),
        file(params.sample_sheet)
    )
}
//     // Stage 4: remove bedmethyl and shortbed files if not needed
//     if (!params.keep_bed_files) {
//         def bedmethyl_dir = new File("${params.outdir}/bedmethyl")
//         if (bedmethyl_dir.deleteDir()) {
//             log.info "Deleted directory: ${bedmethyl_dir}"
//         } else {
//             log.warn "Failed to delete directory: ${bedmethyl_dir}"
//         }

//         def shortbed_dir = new File("${params.outdir}/shortbeds")
//         if (shortbed_dir.deleteDir()) {
//             log.info "Deleted directory: ${shortbed_dir}"
//         } else {
//             log.warn "Failed to delete directory: ${shortbed_dir}"
//         }
//     }

// }
