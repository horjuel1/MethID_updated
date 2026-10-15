#!/usr/bin/env Rscript
# -----------------------------
# Load required libraries
# -----------------------------
library(argparser)
library(data.table)
library(jsonlite)

# Builds the reference atlas from per-cell-type mapped.bed directories.
#
# Replicate samples within a cell type are pooled as sum(meth) / sum(cov), so
# a high-coverage sample outweighs a low-coverage one.
#
# Pooled coverage is written to a sibling _coverage.tsv (not mixed into the
# beta table, so atlas consumers are unaffected), and per-sample meth/cov
# counts to <output>_raw_counts/<cell_type>.tsv.gz, one file per cell type,
# for the beta-binomial fit.
#
# Memory: cell types are processed one at a time; each one's replicate data
# is written out and freed before the next is read, and only the pooled beta
# and coverage columns are kept for the final genome-wide merge.
#
# Inputs: either --input_dirs (a flat list; how build_atlas.nf drives it) or
# --config, a JSON file of { "<source>": { "base_dir": "...", "cell_types":
# ["...", ...] } } resolved to <base_dir>/<cell_type>/<bed_subdir>. Missing or
# empty directories are warned on and skipped; the run aborts only if a
# source's cell_types list is empty or no valid directory is found at all.

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
p <- add_argument(p, "--min_coverage", default = 1,
                  help = "Minimum pooled coverage (sum of cov across replicate samples) required to report a non-NA combined_beta at a CpG for a cell type; below this, combined_beta is set NA rather than reported from too little evidence [1, i.e. only require any data at all]")
p <- add_argument(p, "--cpus", default = 2,
                  help = "Number of CPU cores to use for parallel processing")
p <- add_argument(p, "--samples_to_exclude", nargs = Inf,
                  help = "List of sample names to exclude from processing")
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
min_coverage   <- as.numeric(args$min_coverage)
cpus           <- as.integer(args$cpus)
files_to_ignore <- args$samples_to_exclude
setDTthreads(cpus)

cat("===== Script Parameters =====\n")
print(args)
cat("============================\n")

# -----------------------------
# Output locations
# -----------------------------
raw_dir <- sub("\\.tsv$", "_raw_counts", output_file)
unlink(raw_dir, recursive = TRUE)
dir.create(raw_dir, recursive = TRUE)

final_combined_data <- NULL
final_coverage_data <- NULL

# -----------------------------
# Process each cell type
# -----------------------------
for (i in seq_along(input_dirs)) {
  input_dir <- input_dirs[i]
  label     <- input_labels[i]

  files <- list.files(input_dir, pattern = "*.bed$", full.names = TRUE)
  if (!is.null(files_to_ignore) && length(files_to_ignore) > 0) {
    files <- files[!basename(files) %in% files_to_ignore]
    message("Excluding files based on: ", paste(files_to_ignore, collapse = ", "))
  }

  # Every replicate sample in one long table: chr,start,end,meth,cov,sample
  samples_long <- rbindlist(lapply(files, function(file) {
    message("Loading file: ", file)
    dt <- fread(file, header = FALSE, sep = "\t")
    setnames(dt, c("chr", "start", "end", "meth", "cov"))
    dt[, end := start + 1]  # 1bp intervals, matching the original script
    dt[, sample := sub(".bed", "", basename(file), fixed = TRUE)]
    dt
  }))
  samples_long <- samples_long[!chr %in% exclude_chroms]

  # Pool meth/cov across replicates, then derive beta.
  pooled <- samples_long[, .(
    total_meth         = sum(meth, na.rm = TRUE),
    total_cov          = sum(cov, na.rm = TRUE),
    n_samples_covered  = sum(cov > 0, na.rm = TRUE)
  ), by = .(chr, start, end)]
  pooled[, combined_beta := fifelse(total_cov > 0, total_meth / total_cov, NA_real_)]
  pooled[total_cov < min_coverage, combined_beta := NA_real_]
  pooled[, combined_beta := round(combined_beta, 4)]

  # Per-sample counts for this cell type, one meth + one cov column per
  # sample, prefixed with the cell type (sample filenames aren't unique
  # across cell-type directories).
  counts_long <- melt(samples_long, id.vars = c("chr", "start", "end", "sample"),
                      measure.vars = c("meth", "cov"), variable.name = "metric")
  counts_long[, col := paste(label, sample, metric, sep = "_")]
  fwrite(dcast(counts_long, chr + start + end ~ col, value.var = "value"),
         file.path(raw_dir, paste0(label, ".tsv.gz")), sep = "\t", quote = FALSE)
  rm(samples_long, counts_long)

  dt <- pooled[, .(chr, start, end, combined_beta)]
  setnames(dt, "combined_beta", label)
  cov_dt <- pooled[, .(chr, start, end, total_cov, n_samples_covered)]
  setnames(cov_dt, c("total_cov", "n_samples_covered"), c(paste0(label, "_cov"), paste0(label, "_nsamples")))
  rm(pooled)

  if (is.null(final_combined_data)) {
    final_combined_data <- dt
    final_coverage_data <- cov_dt
  } else {
    final_combined_data <- merge(final_combined_data, dt, by = c("chr", "start", "end"), all = TRUE)
    final_coverage_data <- merge(final_coverage_data, cov_dt, by = c("chr", "start", "end"), all = TRUE)
  }
  rm(dt, cov_dt)
  invisible(gc())
  message(sprintf("Done with %s (%d of %d)", label, i, length(input_dirs)))
}
message("Per-sample meth/cov counts written to: ", raw_dir)

# -----------------------------
# Write the final output
# -----------------------------

fwrite(final_combined_data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Done! MethID Map written to: ", output_file)

coverage_file <- sub(".tsv$", "_coverage.tsv", output_file)
fwrite(final_coverage_data, file = coverage_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Pooled per-cell-type coverage written to: ", coverage_file)
