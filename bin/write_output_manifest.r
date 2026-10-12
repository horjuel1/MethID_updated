#!/usr/bin/env Rscript

# Generates a real, walked-from-disk output manifest for either pipeline
# (build_atlas.nf or main.nf) -- not a static doc that can silently drift
# out of sync with what actually got produced. Every file under --outdir is
# matched against a fixed pattern table (regex on the path relative to
# --outdir); matches are grouped by which pattern they hit (a per-sample
# file class collapses to one manifest row with a count, not N rows) and
# unmatched files are listed individually and flagged rather than hidden,
# so a genuinely new/unexpected output can't quietly go undocumented.

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(argparser)))

p <- arg_parser("Walk a pipeline output directory and write a real (not static) manifest of every file produced")
p <- add_argument(p, "--outdir", type = "character", help = "Pipeline output directory (params.outdir) to walk")
p <- add_argument(p, "--output_tsv", type = "character", default = "OUTPUT_MANIFEST.tsv", help = "Output TSV path (relative to --outdir unless absolute) [OUTPUT_MANIFEST.tsv]")
p <- add_argument(p, "--output_md",  type = "character", default = "OUTPUT_MANIFEST.md",  help = "Output Markdown path (relative to --outdir unless absolute) [OUTPUT_MANIFEST.md]")

args <- parse_args(p)
resolve <- function(path) if (grepl("^/", path)) path else file.path(args$outdir, path)

# -----------------------------------------------
# Pattern table: regex (matched against the path relative to --outdir,
# forward slashes), description, producing process, scope. Matched
# top-to-bottom, first hit wins -- more specific patterns listed first.
# Shared across both pipelines since a single manifest script serves both.
# -----------------------------------------------

patterns <- list(
  # --- build_atlas.nf ---
  list(rx = "(^|/)methid_map_gss\\.tsv$",              desc = "MethID map: the atlas subset to GSS marker CpGs -- pass to main.nf as --methid_map",                                                     proc = "MAKE_MARKER_MAP", scope = "global"),
  list(rx = "(^|/)atlas_beta_binomial_mean\\.tsv$",    desc = "Beta-binomial posterior mean E[beta_kj] per CpG per cell type -- the atlas deconvolve_em.r consumes as --map_file",                    proc = "ATLAS_BETA_BINOMIAL_FIT", scope = "global"),
  list(rx = "(^|/)atlas_beta_binomial_var\\.tsv$",     desc = "Beta-binomial posterior variance Var[beta_kj] per CpG per cell type",                                                                     proc = "ATLAS_BETA_BINOMIAL_FIT", scope = "global"),
  list(rx = "(^|/)atlas_beta_binomial_dispersion\\.tsv$", desc = "Per-cell-type fitted overdispersion phi_k (intra-class correlation) + estimation diagnostics",                                          proc = "ATLAS_BETA_BINOMIAL_FIT", scope = "global"),
  list(rx = "(^|/)atlas_raw_counts\\.tsv$",            desc = "Per-replicate raw meth/cov counts, one column pair per (cell_type, sample) -- input to the beta-binomial fit",                             proc = "BUILD_ATLAS", scope = "global"),
  list(rx = "(^|/)atlas_raw\\.tsv$",                   desc = "Per-replicate raw beta values, one column per (cell_type, sample)",                                                                       proc = "BUILD_ATLAS", scope = "global"),
  list(rx = "(^|/)atlas_coverage\\.tsv$",               desc = "Per-CpG pooled coverage and replicate-sample count backing the pooled atlas",                                                              proc = "BUILD_ATLAS", scope = "global"),
  list(rx = "(^|/)atlas\\.tsv$",                        desc = "Pooled per-CpG methylation atlas (coverage-weighted mean beta), one column per cell type",                                                  proc = "BUILD_ATLAS", scope = "global"),
  list(rx = "(^|/)gss_markers(_summary)?\\.tsv$",       desc = "GSS-selected marker CpGs per cell type (and per-cell-type selection summary)",                                                              proc = "GSS_MARKER_SELECTION", scope = "global"),
  list(rx = "(^|/)atlas_fstats.*\\.tsv$",               desc = "Legacy F-statistic marker selection output (kept for comparison against GSS)",                                                             proc = "MARKER_SELECTION", scope = "global"),
  list(rx = "(^|/)lineage_map\\.tsv$",                  desc = "Cell-type-to-lineage grouping from the confusability dendrogram -- layer D's lambda(k) prerequisite",                                     proc = "ATLAS_SEPARABILITY", scope = "global"),
  list(rx = "(^|/)lineage_map.*dendrogram\\.png$",      desc = "Confusability dendrogram with the lineage cut marked",                                                                                     proc = "ATLAS_SEPARABILITY", scope = "global"),
  list(rx = "(^|/)atlas_pair_separability\\.csv$",      desc = "Full pairwise cell-type confusability table",                                                                                              proc = "ATLAS_SEPARABILITY", scope = "global"),
  list(rx = "(^|/)atlas_nearest_neighbor\\.csv$",       desc = "Per-cell-type nearest-(most confusable) neighbor summary",                                                                                  proc = "ATLAS_SEPARABILITY", scope = "global"),
  list(rx = "(^|/)atlas_(confusability|correlation)_map\\.png$", desc = "Pairwise cell-type confusability/correlation heatmap",                                                                            proc = "ATLAS_SEPARABILITY", scope = "global"),
  list(rx = "(^|/)atlas_density_pairs\\.png$",          desc = "2D beta-density plot for the most-confusable cell-type pairs",                                                                             proc = "ATLAS_SEPARABILITY", scope = "global"),

  # --- main.nf: calibration sub-workflow ---
  list(rx = "(^|/)synthetic_reads.*molecules\\.bed$",   desc = "Synthetic single-molecule reads drawn from the beta-binomial atlas, known ground truth",                                                   proc = "GENERATE_SYNTHETIC_READS", scope = "global (per atlas)"),
  list(rx = "(^|/)synthetic_reads.*truth\\.tsv$",       desc = "Ground truth (true_cell_type, theta_true, injected_condition) for the synthetic reads",                                                    proc = "GENERATE_SYNTHETIC_READS", scope = "global (per atlas)"),
  list(rx = "(^|/)confusion_matrix(_lineage)?\\.tsv$",  desc = "Cell-type / lineage confusion matrix from the synthetic-read calibration run",                                                              proc = "CALIBRATE_LINEAGE_CONFUSION", scope = "global (per atlas)"),
  list(rx = "(^|/)delta_calibration\\.tsv$",            desc = "Layer D delta calibration curve + chosen value (row 1)",                                                                                    proc = "CALIBRATE_LINEAGE_CONFUSION", scope = "global (per atlas)"),
  list(rx = "(^|/)omega_null_validation\\.tsv$",        desc = "omega_i false-positive-rate / power validation against synthetic ground truth",                                                            proc = "CALIBRATE_LINEAGE_CONFUSION", scope = "global (per atlas)"),
  list(rx = "(^|/)calibration_report\\.tsv$",           desc = "Single-row calibration summary (chosen delta, achieved false-lineage rate, omega_i validation pass/fail)",                                proc = "CALIBRATE_LINEAGE_CONFUSION", scope = "global (per atlas)"),

  # --- main.nf: real per-sample stages ---
  list(rx = "(^|/)calls/.*_calls\\.tsv$",               desc = "Per-read, per-CpG methylation calls from `modkit extract calls`",                                                                          proc = "MODKIT_EXTRACT", scope = "per-sample"),
  list(rx = "(^|/)molecule_beds/.*_molecules\\.bed$",   desc = "Molecule-bed (chrom, start, end, read_id, call), one row per (read, CpG)",                                                                   proc = "CALLS_TO_MOLECULEBED", scope = "per-sample"),
  list(rx = "(^|/)spike_in/.*_spikein\\.tsv$",          desc = "Per-sample spike-in fp/fn/theta_failed_prior calibration",                                                                                  proc = "SPIKE_IN_CALIBRATION", scope = "per-sample"),
  list(rx = "(^|/)spike_in_stats\\.tsv$",                desc = "Combined spike-in calibration table (one row per sample), consumed by DECONVOLVE_EM",                                                       proc = "SPIKE_IN_COMBINE", scope = "global"),
  list(rx = "(^|/).*_em_restarts\\.tsv$",                desc = "Per-restart final log-likelihood from layer E's multi-start EM (one row per restart)",                                                     proc = "DECONVOLVE_EM", scope = "per-sample"),
  list(rx = "(^|/).*_molecule_assignments\\.tsv\\.gz$",  desc = "Per-molecule posteriors, best_cell_type, entropy_H, margin_delta, resolution, omega_z/p/outlier, theta_converted_posterior",               proc = "DECONVOLVE_EM", scope = "per-sample"),
  list(rx = "(^|/).*_deconvolution\\.csv$",              desc = "Per-sample cell-type mixture fractions (pi_k, %)",                                                                                        proc = "DECONVOLVE_EM", scope = "per-sample"),
  list(rx = "(^|/).*_coverage\\.tsv$",                   desc = "Per-sample summary stats: molecule counts, EM convergence, theta_failed_fraction, resolution breakdown, mean entropy_H/margin_delta, omega_outlier_fraction", proc = "DECONVOLVE_EM", scope = "per-sample"),
  list(rx = "(^|/)deconvolved_samples\\.csv$",           desc = "All samples' pi_k combined into one wide table",                                                                                           proc = "COMBINE", scope = "global"),
  list(rx = "(^|/)coverage_.*\\.tsv$",                    desc = "All samples' coverage.tsv rows combined",                                                                                                  proc = "COMBINE", scope = "global"),
  list(rx = "(^|/)combined_confidence_summary\\.tsv$",   desc = "Just the confidence/uncertainty-quantification columns from coverage.tsv, across all samples",                                             proc = "COMBINE", scope = "global"),
  list(rx = "(^|/)manifest_.*\\.txt$",                   desc = "Run provenance (atlas file, output dir, date)",                                                                                            proc = "COMBINE", scope = "global"),

  # --- shared / top-level ---
  list(rx = "(^|/)manifest\\.txt$",                      desc = "Pipeline parameter summary, written at launch",                                                                                            proc = "(workflow entrypoint)", scope = "global"),
  list(rx = "(^|/)OUTPUT_MANIFEST\\.(tsv|md)$",          desc = "This manifest",                                                                                                                            proc = "WRITE_MANIFEST", scope = "global"),
  list(rx = "(^|/)\\.nextflow.*", desc = "Nextflow's own work/log files (not a pipeline output)", proc = "(nextflow)", scope = "internal")
)

all_files <- list.files(args$outdir, recursive = TRUE, full.names = FALSE, all.files = FALSE, no.. = TRUE)
if (length(all_files) == 0) stop("--outdir ", args$outdir, " has no files -- nothing to manifest.")

match_pattern <- function(relpath) {
  for (i in seq_along(patterns)) {
    if (grepl(patterns[[i]]$rx, relpath, perl = TRUE)) return(i)
  }
  NA_integer_
}
pattern_idx <- vapply(all_files, match_pattern, integer(1))

file_sizes <- file.size(file.path(args$outdir, all_files))

grouped <- data.table(path = all_files, pattern_idx = pattern_idx, size = file_sizes)

rows <- list()
for (i in seq_along(patterns)) {
  hits <- grouped[pattern_idx == i]
  if (nrow(hits) == 0) next
  pat <- patterns[[i]]
  rows[[length(rows) + 1]] <- data.table(
    pattern = pat$desc, produced_by = pat$proc, scope = pat$scope,
    n_files = nrow(hits), total_size_bytes = sum(hits$size),
    example_path = hits$path[1]
  )
}
manifest_dt <- if (length(rows) > 0) rbindlist(rows) else data.table()

unmatched <- grouped[is.na(pattern_idx)]
if (nrow(unmatched) > 0) {
  unmatched_dt <- data.table(
    pattern = paste0("[UNRECOGNIZED] ", unmatched$path), produced_by = "(unknown -- check which process wrote this)",
    scope = "?", n_files = 1L, total_size_bytes = unmatched$size, example_path = unmatched$path
  )
  manifest_dt <- rbindlist(list(manifest_dt, unmatched_dt), fill = TRUE)
}

setorder(manifest_dt, scope, produced_by)
fwrite(manifest_dt, resolve(args$output_tsv), sep = "\t")

fmt_size <- function(b) {
  if (is.na(b)) return("?")
  units <- c("B", "KB", "MB", "GB", "TB")
  i <- 1
  while (b >= 1024 && i < length(units)) { b <- b / 1024; i <- i + 1 }
  sprintf("%.1f %s", b, units[i])
}

md_lines <- c(
  "# Output manifest",
  sprintf("Generated from a real walk of `%s` -- %d file(s) total, %s.", args$outdir, length(all_files), fmt_size(sum(file_sizes))),
  ""
)
if (nrow(unmatched) > 0) {
  md_lines <- c(md_lines, sprintf("**%d file(s) did not match any known pattern** -- flagged below, not hidden.", nrow(unmatched)), "")
}
for (scope_val in unique(manifest_dt$scope)) {
  md_lines <- c(md_lines, sprintf("## %s", scope_val), "")
  sub <- manifest_dt[scope == scope_val]
  for (r in seq_len(nrow(sub))) {
    md_lines <- c(md_lines, sprintf(
      "- **%s** (%s) -- %d file(s), %s. Example: `%s`",
      sub$pattern[r], sub$produced_by[r], sub$n_files[r], fmt_size(sub$total_size_bytes[r]), sub$example_path[r]
    ))
  }
  md_lines <- c(md_lines, "")
}
writeLines(md_lines, resolve(args$output_md))

cat(sprintf("Wrote %s and %s: %d pattern group(s), %d unmatched file(s)\n",
            resolve(args$output_tsv), resolve(args$output_md), length(rows), nrow(unmatched)))
