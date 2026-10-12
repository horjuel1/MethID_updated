#!/usr/bin/env Rscript
# -----------------------------
# Load required libraries
# -----------------------------
library(argparser)
library(data.table)
library(jsonlite)

# Coverage-weighted rewrite of the original atlas-build script (doc 8 core
# fix). The original combined replicate samples within a cell type by
# computing each sample's beta = meth/cov, THEN discarding meth/cov and
# averaging the per-sample beta ratios unweighted (rowMeans/median across
# sample columns). That treats a CpG covered by 1 read in one sample
# exactly the same as one covered by 500 reads in another -- for a cell
# type with few replicates or shallow coverage at a site, combined_beta
# can collapse straight to a single noisy 0/1 value. Empirically, on the
# atlas this produced, 3.6-8% of CpGs sit at exactly beta=0 or 1 for most
# cell types, and 47.9% for one shallow-source cell type -- far more than
# plausible as genuine always/never-methylated biology at that rate.
#
# Fix: pool meth and cov across replicate samples FIRST (sum both), then
# divide -- combined_beta = sum(meth) / sum(cov) -- so a high-coverage
# sample properly outweighs a low-coverage one instead of counting equally.
# The pooled coverage itself is carried through as a real output (a
# sibling _coverage.tsv file, not mixed into the main beta table, so every
# existing downstream consumer of the main atlas file keeps working
# unmodified) so future work -- a --min_coverage filter in
# pairwise_discriminating_cpgs.r, a beta-binomial group test replacing the
# marker-selection F-statistic, deconvolve_em.r distinguishing a
# well-supported atlas value from a single-low-coverage one -- has
# something to filter/weight on, which the original atlas gave it no way
# to do.
#
# --summarize_method gains "pooled_counts" (the fix above) as the new
# default; "mean"/"median" (the original unweighted-per-sample-beta
# behavior) are kept as opt-in for direct comparison, not removed.
#
# --get_raw_atlas is also extended: alongside the original one-beta-column-
# per-sample wide table, it now writes a companion <output>_raw_counts.tsv
# with per-sample meth AND cov columns (not just beta). The original raw
# atlas only ever carried beta, so nothing downstream of it (e.g. the
# marker-selection F-statistic script) could be upgraded to a proper
# count-based test even in principle -- this unblocks that without
# changing the original raw beta table's shape.
#
# --config is an alternative to spelling out every --input_dirs entry by
# hand: a JSON file of { "<source>": { "base_dir": "...", "cell_types":
# ["...", ...] } }, one entry per reference source. For each source x
# cell_type, this expects <base_dir>/<cell_type>/<bed_subdir> to exist and
# contain at least one *.bed file -- missing/empty directories are warned
# on and skipped rather than failing the whole run, and the run only
# aborts if a source's cell_types list is empty (impossible to leave a
# source silently unconfigured, unlike a bash associative-array entry
# that's simply never assigned) or if NO valid directory is found at all.
# Absorbs what was previously an external bash "gather all valid bed
# directories" pre-processing step directly into this script, so a single
# `Rscript build_atlas.r --config atlas_sources.json ...` invocation is
# self-contained. --input_dirs (the original flat-list interface) is kept
# as a mutually-exclusive alternative -- build_atlas.nf still drives it
# that way (a Nextflow process assembling the list itself), and it's the
# simpler path when there's only a handful of directories to point at
# directly rather than a whole source/cell-type config to maintain.

# -----------------------------
# Parse command-line arguments
# -----------------------------
p <- arg_parser("Create base methylation atlas (coverage-weighted pooling)")
p <- add_argument(p, "--input_dirs", nargs = Inf, default = NA,
                  help = "One or more directories containing mapped.bed files. Mutually exclusive with --config.")
p <- add_argument(p, "--config", type = "character", default = NA,
                  help = "JSON file: { \"<source>\": { \"base_dir\": \"...\", \"cell_types\": [\"...\"] } }. For each source x cell_type, expects <base_dir>/<cell_type>/<bed_subdir> to exist with >=1 *.bed file. Mutually exclusive with --input_dirs.")
p <- add_argument(p, "--bed_subdir", type = "character", default = "bed_hg19",
                  help = "Subdirectory name appended after cell_type when expanding --config [bed_hg19]")
p <- add_argument(p, "--output_file",
                  help = "Path to write the final combined output (TSV)")
p <- add_argument(p, "--summarize_method", default = "pooled_counts",
                  help = "Summary method: 'pooled_counts' (coverage-weighted, recommended default), 'mean', or 'median' (legacy, unweighted across per-sample beta ratios)")
p <- add_argument(p, "--min_coverage", default = 1,
                  help = "Minimum pooled coverage (sum of cov across replicate samples) required to report a non-NA combined_beta at a CpG for a cell type; below this, combined_beta is set NA rather than reported from too little evidence [1, i.e. only require any data at all]")
p <- add_argument(p, "--cpus", default = 2,
                  help = "Number of CPU cores to use for parallel processing")
p <- add_argument(p, "--samples_to_exclude", nargs = Inf,
                  help = "List of sample names to exclude from processing")
p <- add_argument(p, "--get_raw_atlas", default = FALSE,
                  help = "Flag to indicate if raw atlas (+ raw counts) should be generated")
p <- add_argument(p, "--exclude_chroms", default = "chrX,chrY,chrM",
                  help = "Comma-separated contigs to drop from the map before writing (sex + mitochondrial by default). Pass '' to keep every contig.")

args <- parse_args(p)

exclude_chroms <- trimws(strsplit(args$exclude_chroms, ",")[[1]])
exclude_chroms <- exclude_chroms[nzchar(exclude_chroms)]

# -----------------------------
# Resolve --input_dirs vs --config into one input_dirs vector, plus a
# parallel input_labels vector (source_celltype) aligned with it.
# -----------------------------
input_dirs_given <- !all(is.na(args$input_dirs))
config_given <- !is.na(args$config)
if (input_dirs_given && config_given) stop("Pass either --input_dirs or --config, not both.")
if (!input_dirs_given && !config_given) stop("Pass either --input_dirs (one or more directories) or --config (a JSON source/cell_type manifest).")

if (config_given) {
  if (!file.exists(args$config)) stop("--config file not found: ", args$config)
  sources <- fromJSON(args$config, simplifyVector = FALSE)
  if (length(sources) == 0) stop("--config ", args$config, " has no sources defined.")

  for (src in names(sources)) {
    if (is.null(sources[[src]]$cell_types) || length(sources[[src]]$cell_types) == 0) {
      stop("--config: no cell_types defined for source '", src, "' -- every source must list at least one cell type (an empty/omitted list is treated as a config error, not silently skipped).")
    }
    if (is.null(sources[[src]]$base_dir) || !nzchar(sources[[src]]$base_dir)) {
      stop("--config: no base_dir defined for source '", src, "'.")
    }
  }

  input_dirs <- character(0)
  input_labels <- character(0)
  for (src in names(sources)) {
    base_src_dir <- sources[[src]]$base_dir
    cell_types <- unlist(sources[[src]]$cell_types)
    for (cell in cell_types) {
      bed_dir <- file.path(base_src_dir, cell, args$bed_subdir)
      has_bed <- dir.exists(bed_dir) && length(list.files(bed_dir, pattern = "\\.bed$")) > 0
      if (has_bed) {
        input_dirs <- c(input_dirs, bed_dir)
        input_labels <- c(input_labels, paste(src, cell, sep = "_"))
      } else {
        message("WARNING: Missing or empty directory: ", bed_dir)
      }
    }
  }
  if (length(input_dirs) == 0) stop("--config: no valid bed directories found across any source/cell_type.")
} else {
  input_dirs <- args$input_dirs
  # Original path-derived label (third-from-last / second-from-last path
  # component) -- kept exactly as before for this mode, since existing
  # callers (build_atlas.nf) already produce directory paths shaped for it.
  input_labels <- vapply(input_dirs, function(d) {
    path_parts <- strsplit(d, "/")[[1]]
    paste(path_parts[length(path_parts) - 2], path_parts[length(path_parts) - 1], sep = "_")
  }, character(1), USE.NAMES = FALSE)
}

output_file    <- args$output_file
summarize_method <- match.arg(args$summarize_method, choices = c("pooled_counts", "mean", "median"))
min_coverage   <- as.numeric(args$min_coverage)
cpus           <- as.integer(args$cpus)
files_to_ignore <- args$samples_to_exclude
get_raw_atlas  <- args$get_raw_atlas
setDTthreads(cpus)

cat("===== Script Parameters =====\n")
print(args)
cat("============================\n")

# -----------------------------
# Initialize storage for cell-type data
# -----------------------------
all_cell_data <- list()      # per cell type: pooled table (chr,start,end,total_meth,total_cov,n_samples,n_samples_covered,mean_beta,median_beta,pooled_beta,combined_beta)
all_raw_beta  <- list()      # per cell type: wide, one beta column per sample (original --get_raw_atlas shape)
all_raw_counts <- list()     # per cell type: wide, one meth + one cov column per sample (new)

# -----------------------------
# Process each input directory
# -----------------------------
for (i in seq_along(input_dirs)) {
  input_dir <- input_dirs[i]
  label     <- input_labels[i]

  # List bed files
  files <- list.files(input_dir,
                      pattern   = "*.bed$",
                      full.names = TRUE)

  # Exclude specific files if provided
  if (!is.null(files_to_ignore) && length(files_to_ignore) > 0) {
    base_filenames <- basename(files)
    files <- files[!base_filenames %in% files_to_ignore]
    message("Excluding files based on: ", paste(files_to_ignore, collapse = ", "))
  }

  # -----------------------------
  # Read every replicate sample into one long table: chr,start,end,meth,cov,beta,sample
  # -----------------------------
  samples_long <- rbindlist(lapply(files, function(file) {
    message("Loading file: ", file)
    dt <- fread(file, header = FALSE, sep = "\t")
    setnames(dt, c("chr", "start", "end", "meth", "cov"))
    dt[, end := start + 1]  # 1bp intervals, matching the original script
    dt[, beta := fifelse(cov > 0, meth / cov, NA_real_)]
    dt[, sample := sub(".bed", "", basename(file), fixed = TRUE)]
    dt
  }))

  # -----------------------------
  # Pool meth/cov across replicates FIRST, then derive beta -- the core fix.
  # Also keep the legacy unweighted mean/median of per-sample beta ratios
  # for direct comparison via --summarize_method.
  # -----------------------------
  pooled <- samples_long[, .(
    total_meth         = sum(meth, na.rm = TRUE),
    total_cov          = sum(cov, na.rm = TRUE),
    n_samples          = uniqueN(sample),
    n_samples_covered  = sum(cov > 0, na.rm = TRUE),
    mean_beta          = mean(beta, na.rm = TRUE),
    median_beta         = as.numeric(median(beta, na.rm = TRUE))
  ), by = .(chr, start, end)]

  pooled[, pooled_beta := fifelse(total_cov > 0, total_meth / total_cov, NA_real_)]

  combined_col <- switch(summarize_method,
    pooled_counts = "pooled_beta",
    mean          = "mean_beta",
    median        = "median_beta"
  )
  pooled[, combined_beta := get(combined_col)]
  pooled[total_cov < min_coverage, combined_beta := NA_real_]
  pooled[, combined_beta := round(combined_beta, 4)]

  all_cell_data[[label]] <- pooled

  if (isTRUE(get_raw_atlas) || identical(get_raw_atlas, "TRUE")) {
    # Prefix every per-sample column with this cell type's label before any
    # cross-cell-type merge -- sample filenames are not guaranteed unique
    # across different cell-type directories (e.g. every cell type having
    # its own "sample1.bed"), and merge()'s .x/.y collision suffixing would
    # otherwise silently discard which cell type a raw column belongs to.
    raw_beta_wide <- dcast(samples_long, chr + start + end ~ sample, value.var = "beta")
    sample_cols <- setdiff(names(raw_beta_wide), c("chr", "start", "end"))
    setnames(raw_beta_wide, sample_cols, paste(label, sample_cols, sep = "_"))
    all_raw_beta[[label]] <- raw_beta_wide

    counts_long <- melt(samples_long, id.vars = c("chr", "start", "end", "sample"),
                         measure.vars = c("meth", "cov"), variable.name = "metric")
    counts_long[, col := paste(label, sample, metric, sep = "_")]
    raw_counts_wide <- dcast(counts_long, chr + start + end ~ col, value.var = "value")
    all_raw_counts[[label]] <- raw_counts_wide
  }
}

# -----------------------------
# Combine combined_beta (+ pooled coverage) across all cell types
# -----------------------------
final_combined_data <- NULL
final_coverage_data <- NULL

for (label in names(all_cell_data)) {
  dt <- all_cell_data[[label]][, .(chr, start, end, combined_beta)]
  setnames(dt, "combined_beta", label)

  cov_dt <- all_cell_data[[label]][, .(chr, start, end, total_cov, n_samples_covered)]
  setnames(cov_dt, c("total_cov", "n_samples_covered"), c(paste0(label, "_cov"), paste0(label, "_nsamples")))

  if (is.null(final_combined_data)) {
    final_combined_data <- dt
    final_coverage_data <- cov_dt
  } else {
    final_combined_data <- merge(final_combined_data, dt, by = c("chr", "start", "end"), all = TRUE)
    final_coverage_data <- merge(final_coverage_data, cov_dt, by = c("chr", "start", "end"), all = TRUE)
  }
}

# If raw atlas is requested (without summarizing samples), it is generated here
if (isTRUE(get_raw_atlas) || identical(get_raw_atlas, "TRUE")) {
  raw_atlas <- Reduce(function(a, b) merge(a, b, by = c("chr", "start", "end"), all = TRUE), all_raw_beta)
  raw_counts <- Reduce(function(a, b) merge(a, b, by = c("chr", "start", "end"), all = TRUE), all_raw_counts)

  beta_cols <- setdiff(names(raw_atlas), c("chr", "start", "end"))
  raw_atlas[, (beta_cols) := lapply(.SD, function(x) round(x, 4)), .SDcols = beta_cols]

  raw_atlas <- raw_atlas[!chr %in% exclude_chroms]
  raw_counts <- raw_counts[!chr %in% exclude_chroms]

  fwrite(raw_atlas, file = sub(".tsv$", "_raw.tsv", output_file), sep = "\t", quote = FALSE, row.names = FALSE)
  message("Raw MethID Map written to: ", sub(".tsv$", "_raw.tsv", output_file))

  fwrite(raw_counts, file = sub(".tsv$", "_raw_counts.tsv", output_file), sep = "\t", quote = FALSE, row.names = FALSE)
  message("Raw per-sample meth/cov counts written to: ", sub(".tsv$", "_raw_counts.tsv", output_file))
}

# -----------------------------
# Write the final output
# -----------------------------

final_combined_data <- final_combined_data[!chr %in% exclude_chroms]
final_coverage_data <- final_coverage_data[!chr %in% exclude_chroms]

fwrite(final_combined_data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Done! MethID Map written to: ", output_file)

coverage_file <- sub(".tsv$", "_coverage.tsv", output_file)
fwrite(final_coverage_data, file = coverage_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Pooled per-cell-type coverage written to: ", coverage_file)
