#!/usr/bin/env Rscript

#load libraries
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(argparser)))

#------------------------------
# Receive Arguments
#------------------------------

p <- arg_parser("Summarize deconvolution results by label and by label x exposure")
p <- add_argument(p, "--deconv_file", type = "character", help = "deconvolved_samples.csv from combine.r (cell types x samples)")
p <- add_argument(p, "--metadata_file", type = "character", help = "Metadata sheet (CSV, or TSV if named .tsv/.txt) with CGID, label and exposure columns")
p <- add_argument(p, "--output_dir", type = "character", default = ".", help = "Output directory")

args <- parse_args(p)
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Read deconvolution results and metadata
# -----------------------------

# Cell types as rows (first column), one column per sample; reshape to one row per (sample, cell type)
deconv <- fread(args$deconv_file, header = TRUE)
setnames(deconv, 1, "cell_type")
cell_type_order <- deconv$cell_type
results <- melt(deconv, id.vars = "cell_type", variable.name = "CGID", value.name = "percent", variable.factor = FALSE)

sep <- if (grepl("\\.(tsv|txt)$", tolower(args$metadata_file))) "\t" else ","
metadata <- fread(args$metadata_file, sep = sep, colClasses = "character")
missing_cols <- setdiff(c("CGID", "label", "exposure"), names(metadata))
if (length(missing_cols) > 0) {
  stop("Metadata sheet needs columns CGID, label and exposure. Missing: ", paste(missing_cols, collapse = ", "))
}

# Match the pipeline's sample naming (spaces and slashes become underscores)
metadata[, CGID := gsub("[ /]", "_", trimws(CGID))]
metadata[, label := trimws(label)]
metadata[, exposure := trimws(exposure)]
metadata[is.na(label) | label == "", label := "NA"]
metadata[is.na(exposure) | exposure == "", exposure := "NA"]

duplicated_ids <- unique(metadata$CGID[duplicated(metadata$CGID)])
if (length(duplicated_ids) > 0) stop("CGID appears more than once in the metadata sheet: ", paste(duplicated_ids, collapse = ", "))

# -----------------------------
# Join and report what matched
# -----------------------------

samples <- unique(results$CGID)
no_metadata <- setdiff(samples, metadata$CGID)
no_results <- setdiff(metadata$CGID, samples)
message(sprintf("Samples with deconvolution results: %d; with metadata: %d", length(samples), length(intersect(samples, metadata$CGID))))
if (length(no_metadata) > 0) message("No metadata (left out of the summaries): ", paste(no_metadata, collapse = ", "))
if (length(no_results) > 0) message("In the metadata sheet but not deconvolved: ", paste(no_results, collapse = ", "))
if (length(intersect(samples, metadata$CGID)) == 0) {
  stop("No sample names match a CGID in the metadata sheet. Deconvolved samples: ", paste(head(samples, 10), collapse = ", "))
}

merged <- merge(results, metadata[, .(CGID, label, exposure)], by = "CGID")
merged[, cell_type := factor(cell_type, levels = cell_type_order)]

# -----------------------------
# Summaries
# -----------------------------

summarize_by <- function(dt, by) {
  out <- dt[, .(
    n_samples = .N,
    mean      = mean(percent),
    sd        = sd(percent),
    median    = median(percent),
    q25       = quantile(percent, 0.25, names = FALSE),
    q75       = quantile(percent, 0.75, names = FALSE),
    min       = min(percent),
    max       = max(percent)
  ), by = by]
  stat_cols <- c("mean", "sd", "median", "q25", "q75", "min", "max")
  out[, (stat_cols) := lapply(.SD, round, 2), .SDcols = stat_cols]
  setorderv(out, by)
  out[, cell_type := as.character(cell_type)]
  out
}

by_label <- summarize_by(merged, c("label", "cell_type"))
by_label_exposure <- summarize_by(merged, c("label", "exposure", "cell_type"))

per_sample <- merged[, .(CGID, label, exposure, cell_type = as.character(cell_type), percent)]
setorder(per_sample, label, exposure, CGID)

fwrite(by_label,          file.path(args$output_dir, "summary_by_label.csv"))
fwrite(by_label_exposure, file.path(args$output_dir, "summary_by_label_exposure.csv"))
fwrite(per_sample,        file.path(args$output_dir, "deconvolved_samples_with_metadata.csv"))

# Sample counts per group, plus which samples could not be matched
counts <- unique(merged[, .(CGID, label, exposure)])[, .(n_samples = .N), by = .(label, exposure)]
setorder(counts, label, exposure)
fwrite(counts, file.path(args$output_dir, "summary_sample_counts.csv"))
unmatched <- rbind(
  data.table(CGID = no_metadata, issue = rep("deconvolved but not in metadata sheet", length(no_metadata))),
  data.table(CGID = no_results,  issue = rep("in metadata sheet but not deconvolved", length(no_results)))
)
fwrite(unmatched, file.path(args$output_dir, "summary_unmatched_samples.csv"))

#Quit R session
quit(save = "no")
