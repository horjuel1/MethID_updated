#!/usr/bin/env nextflow

// Atlas construction pipeline -- separate entrypoint from main.nf, since
// this builds the reference atlas (bin/build_atlas.r, bin/
// marker_selection_fstat.py) rather than deconvolving samples against one.
//
// Extended per the unified single-molecule EM framework spec (layer A,
// Build Item 1a) with a beta-binomial hyperparameter fit and lineage
// clustering on top of the original coverage-weighted pooled atlas --
// both new stages run by default; --run_beta_binomial_fit/
// --run_gss_marker_selection/--run_lineage_clustering false turns any of
// them off individually if you only want the original atlas.tsv.
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

params.mapped_bed_dirs   = false   // comma-separated list of per-cell-type mapped.bed directories
params.samples_to_exclude = false  // optional comma-separated sample (bed file) names to exclude
params.summarize_method  = 'pooled_counts'  // 'pooled_counts' (recommended), 'mean', or 'median' -- see bin/build_atlas.r
params.min_coverage      = 1
// Changed from the original false default: the beta-binomial fit below
// needs the per-replicate raw counts, and it runs by default now, so raw
// atlas output defaults on too. Set false only if you're disabling every
// downstream stage that needs it (beta-binomial fit AND F-stat/GSS marker
// selection) and just want the original pooled atlas.tsv.
params.get_raw_atlas     = true
params.output_name       = 'atlas.tsv'
params.cpus              = 4
params.outdir            = 'results_atlas'

// Stage 2: beta-binomial hyperparameter fit (layer A)
params.run_beta_binomial_fit    = true
params.bb_min_df_for_dispersion = 30
params.bb_min_phi               = 1e-4
params.bb_max_phi               = 0.999

// Stage 3: GSS marker selection (doc 10), now also the input to lineage
// clustering's marker-restricted distance metric (layer B's repurposing)
params.run_gss_marker_selection = true
params.gss_min_gss              = 0.3
params.gss_n_markers_per_celltype = 50
params.gss_secondary_test_alpha = 0.05
params.gss_require_secondary_pass = false

// Stage 4: pairwise separability + lineage clustering (Build Item 1a)
params.run_lineage_clustering = true
params.n_lineages           = false  // exactly one of n_lineages / lineage_cut_height required if run_lineage_clustering
params.lineage_cut_height   = false

// Legacy: F-statistic marker selection on the raw atlas, kept for
// side-by-side comparison against GSS, not run by default
params.run_marker_selection_fstat = false
params.manifest                    = false  // TSV: group(cell type), sample column name in the raw atlas
params.marker_selection_name       = 'atlas_fstats.tsv'

new File(params.outdir).mkdirs()

def summary = """
         A T L A S   C O N S T R U C T I O N
         ====================================
         Mapped bed dirs         : ${params.mapped_bed_dirs}
         Summarize method        : ${params.summarize_method}
         Min coverage            : ${params.min_coverage}
         Get raw atlas           : ${params.get_raw_atlas}
         Output dir              : ${params.outdir}
         Run beta-binomial fit   : ${params.run_beta_binomial_fit}
         Run GSS marker selection: ${params.run_gss_marker_selection}
         Run lineage clustering  : ${params.run_lineage_clustering}
         N lineages / cut height : ${params.n_lineages} / ${params.lineage_cut_height}
         Run legacy F-stat markers: ${params.run_marker_selection_fstat}
         """
log.info summary
new File("${params.outdir}/manifest.txt").write(summary)

// Input validation
if (!params.mapped_bed_dirs) { exit 1, "Provide --mapped_bed_dirs (comma-separated list of per-cell-type mapped.bed directories)" }
if ((params.run_marker_selection_fstat || params.run_beta_binomial_fit || params.run_gss_marker_selection) && !params.get_raw_atlas) {
    exit 1, "--run_marker_selection_fstat/--run_beta_binomial_fit/--run_gss_marker_selection need per-sample raw beta/counts, not just the pooled atlas -- pass --get_raw_atlas true (the default) or disable all three"
}
if (params.run_marker_selection_fstat && !params.manifest) { exit 1, "--run_marker_selection_fstat requires --manifest" }
if (params.run_lineage_clustering) {
    if (!params.n_lineages && !params.lineage_cut_height) {
        exit 1, "--run_lineage_clustering requires exactly one of --n_lineages / --lineage_cut_height"
    }
    if (params.n_lineages && params.lineage_cut_height) {
        exit 1, "Pass at most one of --n_lineages / --lineage_cut_height, not both"
    }
}

def input_dirs = params.mapped_bed_dirs.split(',').collect { it.trim() }
input_dirs.each { d ->
    if (!file(d).isDirectory()) { exit 1, "--mapped_bed_dirs entry is not a directory: ${d}" }
}

// --- Processes ---

process BUILD_ATLAS {
    tag "Build atlas from ${input_dirs.size()} cell-type director(ies)"
    publishDir params.outdir, mode: 'copy'
    cpus params.cpus

    input:
    val(dirs)

    output:
    path("atlas.tsv"),              emit: atlas
    path("atlas_coverage.tsv"),     emit: coverage
    path("atlas_raw.tsv"),          optional: true, emit: raw_atlas
    path("atlas_raw_counts.tsv"),   optional: true, emit: raw_counts

    script:
    def dirs_arg = dirs.collect { "\"${it}\"" }.join(' ')
    def exclude_arg = params.samples_to_exclude
        ? "--samples_to_exclude " + params.samples_to_exclude.split(',').collect { "\"${it.trim()}\"" }.join(' ')
        : ''
    def raw_arg = params.get_raw_atlas ? "--get_raw_atlas TRUE" : ""
    """
    Rscript ${baseDir}/bin/build_atlas.r \\
        --input_dirs ${dirs_arg} \\
        --output_file "atlas.tsv" \\
        --summarize_method "${params.summarize_method}" \\
        --min_coverage "${params.min_coverage}" \\
        --cpus "${task.cpus}" \\
        ${raw_arg} \\
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
        --raw_counts_file "${raw_counts}" \\
        --output_prefix "atlas_beta_binomial" \\
        --min_coverage "${params.min_coverage}" \\
        --min_df_for_dispersion "${params.bb_min_df_for_dispersion}" \\
        --min_phi "${params.bb_min_phi}" \\
        --max_phi "${params.bb_max_phi}"
    """
}

process GSS_MARKER_SELECTION {
    tag "GSS marker selection"
    publishDir params.outdir, mode: 'copy'

    input:
    path(atlas)      // the POOLED per-cell-type atlas (one column per cell type) --
                      // NOT the raw per-replicate atlas, which has one column per
                      // (cell_type, sample) and would make every replicate look like
                      // its own "cell type" to this script
    path(coverage)

    output:
    path("gss_markers.tsv"),         emit: markers
    path("gss_markers_summary.tsv"), emit: summary

    script:
    def secondary_flag = params.gss_require_secondary_pass ? "--require_secondary_pass" : ""
    """
    Rscript ${baseDir}/bin/marker_selection_gss.r \\
        --map_file "${atlas}" \\
        --coverage_file "${coverage}" \\
        --min_gss "${params.gss_min_gss}" \\
        --n_markers_per_celltype "${params.gss_n_markers_per_celltype}" \\
        --secondary_test_alpha "${params.gss_secondary_test_alpha}" \\
        --output_file "gss_markers.tsv" \\
        --summary_file "gss_markers_summary.tsv" \\
        ${secondary_flag}
    """
    // NOTE: run without --nn_table (stage 4 pairwise gap-filling) by
    // design -- that table comes from ATLAS_SEPARABILITY below, which
    // itself optionally consumes THIS process's markers for its own
    // marker-restricted distance metric. Rather than a two-pass
    // separability run to break that cycle, gap-filling is left as a
    // manual follow-up (rerun this process with --nn_table pointed at
    // atlas_nearest_neighbor.csv) for anyone who wants it.
}

process MAKE_MARKER_MAP {
    tag "Subset atlas to GSS markers"
    publishDir params.outdir, mode: 'copy'

    input:
    path(atlas)     // same atlas GSS selected markers on (beta-binomial mean, or pooled if that fit is off)
    path(markers)

    output:
    path("methid_map_gss.tsv"), emit: methid_map

    script:
    """
    Rscript ${baseDir}/bin/subset_atlas_to_markers.r \\
        --atlas_file "${atlas}" \\
        --marker_file "${markers}" \\
        --output_file "methid_map_gss.tsv"
    """
}

process ATLAS_SEPARABILITY {
    tag "Pairwise separability + lineage clustering"
    publishDir params.outdir, mode: 'copy'

    input:
    path(atlas)
    path(markers)

    output:
    path("atlas_pair_separability.csv"), emit: pair_table
    path("atlas_nearest_neighbor.csv"),  emit: nn_table
    path("atlas_confusability_map.png")
    path("atlas_correlation_map.png")
    path("atlas_density_pairs.png")
    path("lineage_map.tsv"),             optional: true, emit: lineage_map
    path("lineage_map_dendrogram.png"),  optional: true

    script:
    def marker_arg = params.run_gss_marker_selection ? "--marker_file \"${markers}\"" : ""
    def lineage_args = params.run_lineage_clustering
        ? "--lineage_map_out \"lineage_map.tsv\" " + (params.n_lineages ? "--n_lineages ${params.n_lineages}" : "--lineage_cut_height ${params.lineage_cut_height}")
        : ""
    """
    Rscript ${baseDir}/bin/atlas_separability.r \\
        --map_file "${atlas}" \\
        ${marker_arg} \\
        --pair_table_out "atlas_pair_separability.csv" \\
        --nn_table_out "atlas_nearest_neighbor.csv" \\
        --confusability_output_file "atlas_confusability_map.png" \\
        --correlation_output_file "atlas_correlation_map.png" \\
        --density_output_file "atlas_density_pairs.png" \\
        ${lineage_args}
    """
}

process MARKER_SELECTION_FSTAT {
    tag "F-statistic marker selection (legacy)"
    publishDir params.outdir, mode: 'copy'
    cpus params.cpus

    input:
    path(raw_atlas)
    path(manifest)

    output:
    path("${params.marker_selection_name}")

    script:
    """
    python3 ${baseDir}/bin/marker_selection_fstat.py \\
        --map "${raw_atlas}" \\
        --manifest "${manifest}" \\
        --out_file "${params.marker_selection_name}"
    """
}

process WRITE_ATLAS_MANIFEST {
    tag "Write output manifest"
    // Deliberately no publishDir/output: block -- write_output_manifest.r
    // both reads (walks) AND writes directly into params.outdir itself
    // (an absolute path on the shared filesystem, not the task's own work
    // dir), since its whole job is to describe what's already sitting
    // there. Nextflow's normal stage-locally-then-publish flow doesn't fit
    // a process whose output IS a description of the publish directory.

    input:
    val(ready)  // dependency trigger only (a collected list from every upstream stage) -- forces
                // this to run after everything else has published to params.outdir on disk

    script:
    """
    Rscript ${baseDir}/bin/write_output_manifest.r \\
        --outdir "${params.outdir}" \\
        --output_tsv "OUTPUT_MANIFEST.tsv" \\
        --output_md "OUTPUT_MANIFEST.md"
    """
}

// --- Workflow ---

workflow {
    // Stage 1: pool per-sample meth/cov into a coverage-weighted atlas
    // (+ the raw per-sample beta/counts every downstream stage below needs)
    BUILD_ATLAS(input_dirs)

    // Stage 2: beta-binomial hyperparameter fit (layer A)
    def mean_atlas_ch
    def bb_signal
    if (params.run_beta_binomial_fit) {
        ATLAS_BETA_BINOMIAL_FIT(BUILD_ATLAS.out.atlas, BUILD_ATLAS.out.raw_counts)
        mean_atlas_ch = ATLAS_BETA_BINOMIAL_FIT.out.mean_atlas
        bb_signal = ATLAS_BETA_BINOMIAL_FIT.out.dispersion
    } else {
        mean_atlas_ch = BUILD_ATLAS.out.atlas
        bb_signal = Channel.value('beta_binomial_fit_skipped')
    }

    // Stage 3: GSS marker selection -- when skipped, ATLAS_SEPARABILITY
    // still needs SOME file staged for its (mandatory) markers input; a
    // named empty placeholder, same pattern main.nf already uses for
    // NO_SPIKE_IN, and the process script omits --marker_file whenever
    // this placeholder is what got passed.
    def markers_ch
    def gss_signal
    if (params.run_gss_marker_selection) {
        GSS_MARKER_SELECTION(mean_atlas_ch, BUILD_ATLAS.out.coverage)
        markers_ch = GSS_MARKER_SELECTION.out.markers
        // MethID: turn the markers into a map main.nf can use directly
        MAKE_MARKER_MAP(mean_atlas_ch, markers_ch)
        gss_signal = GSS_MARKER_SELECTION.out.summary.mix(MAKE_MARKER_MAP.out.methid_map)
    } else {
        def no_markers = file("${params.outdir}/NO_MARKER_FILE")
        no_markers.text = ''
        markers_ch = Channel.value(no_markers)
        gss_signal = Channel.value('gss_marker_selection_skipped')
    }

    // Stage 4: pairwise separability + lineage clustering (Build Item 1a),
    // on the beta-binomial mean atlas (falls back to the plain pooled
    // atlas if stage 2 was disabled) restricted to GSS markers (falls
    // back to every complete-row CpG if stage 3 was disabled)
    ATLAS_SEPARABILITY(mean_atlas_ch, markers_ch)

    // Legacy: F-statistic marker selection, off by default
    def fstat_signal
    if (params.run_marker_selection_fstat) {
        MARKER_SELECTION_FSTAT(BUILD_ATLAS.out.raw_atlas, file(params.manifest))
        fstat_signal = MARKER_SELECTION_FSTAT.out
    } else {
        fstat_signal = Channel.value('fstat_marker_selection_skipped')
    }

    ready_ch = Channel.empty().mix(
        BUILD_ATLAS.out.atlas, bb_signal, gss_signal, ATLAS_SEPARABILITY.out.pair_table, fstat_signal
    ).collect()
    WRITE_ATLAS_MANIFEST(ready_ch)
}
