#!/usr/bin/env nextflow

// Atlas construction pipeline -- separate entrypoint from main.nf, since
// this builds the reference atlas rather than deconvolving samples against one.
//
// Three stages: pool replicates into atlas.tsv, fit the beta-binomial
// (mean + variance per CpG), then separability + marker balance.
//
// NOTE on --mapped_bed_dirs: bin/build_atlas.r derives each cell type's
// label from its INPUT DIRECTORY's own path structure (the 2nd/3rd-to-last
// path components, e.g. .../loyfer/Blood-T-CD4/mapped_beds -> label
// "loyfer_Blood-T-CD4"), not from any argument. Nextflow's normal `path`
// staging would symlink each directory into the work dir under its own
// basename, discarding that parent-path context -- so these are passed
// through as plain absolute-path strings (`val`, not `path`) and read
// directly off the shared filesystem, the same way main.nf already reads
// params.ref_genome inline rather than staging it.

params.mapped_bed_dirs   = false   // comma-separated list of per-cell-type mapped.bed directories. Mutually exclusive with --config.
params.config            = false   // JSON manifest: { "<source>": { "base_dir": "...", "cell_types": ["...", ...] } } --
                                    // bin/build_atlas.r has supported this alternative to --mapped_bed_dirs since it was
                                    // written (expects <base_dir>/<cell_type>/<bed_subdir> per source/cell_type), but this
                                    // Nextflow wrapper never actually wired it through until now. Mutually exclusive with
                                    // --mapped_bed_dirs.
params.bed_subdir        = false   // only used with --config; overrides bin/build_atlas.r's own "bed_hg19" default if set
params.samples_to_exclude = false  // optional comma-separated sample (bed file) names to exclude
<<<<<<< HEAD
params.summarize_method  = 'pooled_counts'  // 'pooled_counts' (recommended), 'mean', or 'median' -- see bin/build_atlas.r
params.min_coverage      = 10
// Changed from the original false default: the beta-binomial fit below
// needs the per-replicate raw counts, and it runs by default now, so raw
// atlas output defaults on too. Set false only if you're disabling every
// downstream stage that needs it (beta-binomial fit AND F-stat/GSS marker
// selection) and just want the original pooled atlas.tsv.
params.get_raw_atlas     = true
=======
params.min_coverage      = 1
>>>>>>> 96053d0b95b798260e6b399258f93b9e86044f1f
params.output_name       = 'atlas.tsv'
params.cpus              = 4
params.outdir            = 'results_atlas'

// Stage 2: beta-binomial fit
params.bb_min_df_for_dispersion = 30
params.bb_min_phi               = 1e-4
params.bb_max_phi               = 0.999

// Separability is also recomputed on balanced marker sets (as main.nf's
// markers_per_cell_type would select): comma-separated N values, '' to skip
params.separability_marker_sets  = '5000,20000'

new File(params.outdir).mkdirs()

def summary = """
         A T L A S   C O N S T R U C T I O N
         ====================================
         Mapped bed dirs         : ${params.mapped_bed_dirs}
         Config                  : ${params.config}
         Min coverage            : ${params.min_coverage}
         Output dir              : ${params.outdir}
         """
log.info summary
new File("${params.outdir}/manifest.txt").write(summary)

// Input validation
if (!params.mapped_bed_dirs && !params.config) { exit 1, "Provide either --mapped_bed_dirs (comma-separated list of per-cell-type mapped.bed directories) or --config (a JSON source/cell_type manifest)" }
if (params.mapped_bed_dirs && params.config) { exit 1, "Pass either --mapped_bed_dirs or --config, not both." }

// Two mutually-exclusive input modes: --mapped_bed_dirs resolves to a plain
// list of directories at the Nextflow level (validated here); --config just
// gets handed to bin/build_atlas.r as-is, which resolves <base_dir>/
// <cell_type>/<bed_subdir> itself -- Nextflow doesn't need to know the
// individual directory paths in that mode at all.
//
// Both modes are threaded through as plain val() strings, read directly off
// the shared filesystem inside the task rather than staged in (see the
// --mapped_bed_dirs note at the top of this file) -- which means both need
// to be ABSOLUTE paths: a relative --config path resolves against whatever
// directory `nextflow run` was launched from, not the task's own isolated
// work directory, so a relative path silently fails to be found once inside
// the task ("--config file not found") even though it exists right where
// the user typed it. Resolving here, once, up front, means the user doesn't
// have to remember to always pass --config as an absolute path themselves.
def input_dirs = []
def config_abs_path = null   // resolved absolute --config path, if that mode is used (see note above)
if (params.mapped_bed_dirs) {
    input_dirs = params.mapped_bed_dirs.split(',').collect { it.trim() }
    input_dirs.each { d ->
        if (!file(d).isDirectory()) { exit 1, "--mapped_bed_dirs entry is not a directory: ${d}" }
    }
} else {
    if (!file(params.config).exists()) { exit 1, "--config file not found: ${params.config}" }
    // NOTE: reassigning params.config directly here is silently a no-op --
    // Nextflow DSL2 warns "params.config is defined multiple times --
    // Assignments following the first are ignored" and keeps the original
    // (possibly relative) value. A separate local variable is required.
    config_abs_path = file(params.config).toAbsolutePath().toString()
}

// --- Processes ---

process BUILD_ATLAS {
    tag params.config ? "Build atlas from --config ${params.config}" : "Build atlas from ${input_dirs.size()} cell-type director(ies)"
    publishDir params.outdir, mode: 'copy'
    cpus params.cpus

    input:
    val(dirs)        // used in --mapped_bed_dirs mode; empty list in --config mode
    val(config_path)  // used in --config mode; the string 'NO_CONFIG' in --mapped_bed_dirs mode

    output:
    path("atlas.tsv"),              emit: atlas
    path("atlas_coverage.tsv"),     emit: coverage
    path("atlas_raw_counts"),       emit: raw_counts

    script:
    def exclude_arg = params.samples_to_exclude
        ? "--samples_to_exclude " + params.samples_to_exclude.split(',').collect { "\"${it.trim()}\"" }.join(' ')
        : ''
    def input_arg = config_path != 'NO_CONFIG'
        ? "--config \"${config_path}\"" + (params.bed_subdir ? " --bed_subdir \"${params.bed_subdir}\"" : "")
        : "--input_dirs " + dirs.collect { "\"${it}\"" }.join(' ')
    """
    Rscript ${baseDir}/bin/build_atlas.r \\
        ${input_arg} \\
        --output_file "atlas.tsv" \\
        --min_coverage "${params.min_coverage}" \\
        --cpus "${task.cpus}" \\
        ${exclude_arg}
    """
}

process ATLAS_BETA_BINOMIAL_FIT {
    tag "Beta-binomial hyperparameter fit"
    publishDir params.outdir, mode: 'copy'

    input:
    path(atlas)
    path(raw_counts)

    output:
    path("atlas_beta_binomial_mean.tsv"),       emit: mean_atlas
    path("atlas_beta_binomial_var.tsv"),        emit: var_atlas
    path("atlas_beta_binomial_dispersion.tsv"), emit: dispersion

    script:
    """
    Rscript ${baseDir}/bin/fit_atlas_beta_binomial.r \\
        --atlas_file "${atlas}" \\
        --raw_counts_dir "${raw_counts}" \\
        --output_prefix "atlas_beta_binomial" \\
        --min_coverage "${params.min_coverage}" \\
        --min_df_for_dispersion "${params.bb_min_df_for_dispersion}" \\
        --min_phi "${params.bb_min_phi}" \\
        --max_phi "${params.bb_max_phi}"
    """
}

process ATLAS_SEPARABILITY {
    tag "Pairwise separability + marker balance"
    publishDir params.outdir, mode: 'copy'

    input:
    path(atlas)

    output:
    path("atlas_pair_separability.csv")
    path("atlas_marker_balance.tsv")
    path("atlas_confusability_map*.png")

    script:
    """
    Rscript ${baseDir}/bin/atlas_separability.r \\
        --map_file "${atlas}" \\
        --pair_table_out "atlas_pair_separability.csv" \\
        --confusability_output_file "atlas_confusability_map.png" \\
        --marker_sets "${params.separability_marker_sets}"
    """
}

// --- Workflow ---

workflow {
    def config_path_val = config_abs_path ?: 'NO_CONFIG'
    BUILD_ATLAS(input_dirs, config_path_val)
    ATLAS_BETA_BINOMIAL_FIT(BUILD_ATLAS.out.atlas, BUILD_ATLAS.out.raw_counts)
    ATLAS_SEPARABILITY(ATLAS_BETA_BINOMIAL_FIT.out.mean_atlas)
}
