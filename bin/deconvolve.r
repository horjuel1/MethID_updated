#!/usr/bin/env Rscript

#load libraries
suppressWarnings(suppressMessages(library(GenomicRanges)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(gridExtra)))
suppressWarnings(suppressMessages(library(limSolve)))
suppressWarnings(suppressMessages(library(argparser)))


#------------------------------
# Receive Arguments
#------------------------------

p <- arg_parser("Deconvolve a single sample against a methylation atlas")
p <- add_argument(p, "--map_file", type = "character", help = "Atlas map file (TSV format)")
p <- add_argument(p, "--shortbed", type = "character", help = "Shortbed file for a single sample")
p <- add_argument(p, "--sample_name", type = "character", help = "Sample name used to label output files")
p <- add_argument(p, "--output_dir", type = "character", default = ".", help = "Output directory")
p <- add_argument(p, "--threads", type = "integer", default = 1, help = "Number of threads for data.table")
p <- add_argument(p, "--allow_feature_overlap", type = "logical", default = TRUE, help = "If FALSE, do not allow features to overlap in the atlas")
p <- add_argument(p, "--min_cov", type = "integer", default = 2, help = "Minimum coverage for a site to be considered [2]")
p <- add_argument(p, "--markers_per_cell_type", type = "integer", default = 0, help = "Keep, for every cell type, the N atlas CpGs where it differs most from its closest other cell type (either direction), so each cell type is represented equally; 0 = all atlas CpGs [0]")
p <- add_argument(p, "--marker_ranking", type = "character", default = "gap", help = "How --markers_per_cell_type ranks CpGs: gap (difference to the closest other cell type) or standardized (that difference divided by the combined replicate SD of the two cell types, from --var_atlas_file) [gap]")
p <- add_argument(p, "--var_atlas_file", type = "character", default = NA, help = "fit_atlas_beta_binomial.r's <prefix>_var.tsv; needed for --marker_ranking standardized")
p <- add_argument(p, "--marker_var_floor", type = "double", default = 0.0025, help = "Standardized ranking: minimum variance per cell type (0.0025 = SD 0.05) [0.0025]")
p <- add_argument(p, "--beta_eps", type = "double", default = 0.01, help = "Atlas betas are clipped to [eps, 1-eps] when ranking markers, as in smMethID [0.01]")

args <- parse_args(p)
if (!args$marker_ranking %in% c("gap", "standardized")) stop("--marker_ranking must be gap or standardized")
if (args$marker_ranking == "standardized" && args$markers_per_cell_type > 0 && is.na(args$var_atlas_file)) stop("--marker_ranking standardized needs --var_atlas_file")
map_file <- args$map_file
sample_file <- args$shortbed
sample_name <- args$sample_name
output_dir <- args$output_dir
threads <- args$threads
allow_feature_overlap <- args$allow_feature_overlap
min_cov <- args$min_cov

setDTthreads(threads = threads)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------
# Generate GRanges, find overlaps, and deconvolve
# -----------------------------------------------

#Read and convert atlas to GRanges
atlases <- fread(map_file, data.table = FALSE)

# Balanced marker selection, same as smMethID's deconvolve_em.r: for each cell type, how far
# it stands apart from its closest other cell type at every CpG (in either direction); the top
# N per cell type are kept and pooled, so a distinctive cell type can't dominate the fit through
# sheer numbers of CpGs. gap: the difference in beta. standardized: each difference divided by
# the combined replicate SD of the two cell types, so a large difference counts only if
# replicates agree. Ranking uses betas clipped to [eps, 1-eps]; the fit uses the original betas.
if (args$markers_per_cell_type > 0) {
  cell_types <- setdiff(colnames(atlases), c("chr", "start", "end"))
  atlases <- atlases[rowSums(is.na(atlases[, cell_types, drop = FALSE])) == 0, ]
  eps <- args$beta_eps
  clipped <- lapply(atlases[, cell_types, drop = FALSE], function(x) pmin(pmax(x, eps), 1 - eps))

  vdt <- NULL
  if (args$marker_ranking == "standardized") {
    var_atlas <- fread(args$var_atlas_file, data.table = FALSE)
    row_idx <- match(paste(atlases$chr, atlases$start), paste(var_atlas$chr, var_atlas$start))
    # no variance estimate: treat as maximally uncertain
    vdt <- lapply(setNames(cell_types, cell_types), function(ct) { v <- var_atlas[[ct]][row_idx]; v[is.na(v)] <- 0.25; v })
    rm(var_atlas)
  }

  marker_score <- function(k) {
    score <- NULL
    for (j in setdiff(cell_types, k)) {
      d <- abs(clipped[[k]] - clipped[[j]])
      if (!is.null(vdt)) d <- d / sqrt(pmax(vdt[[k]], args$marker_var_floor) + pmax(vdt[[j]], args$marker_var_floor))
      score <- if (is.null(score)) d else pmin(score, d)
    }
    score
  }

  keep <- rep(FALSE, nrow(atlases))
  for (k in cell_types) {
    gap <- marker_score(k)
    top <- order(gap, decreasing = TRUE, method = "radix")[seq_len(min(args$markers_per_cell_type, length(gap)))]
    keep[top] <- TRUE
    message(sprintf("Markers for %s: %d, weakest %s to nearest other cell type %.3f",
                    k, length(top), args$marker_ranking, gap[top[length(top)]]))
  }
  message(sprintf("Balanced markers: %d atlas CpGs kept of %d", sum(keep), nrow(atlases)))
  atlases <- atlases[keep, ]
  rm(clipped, vdt)
}

atlas_granges <- makeGRangesFromDataFrame(atlases, seqnames.field = "chr", start.field = "start", end.field = "end")

# Read sample file and convert to GRanges
sample <- fread(sample_file, data.table = FALSE)

#Check if the sample file has 4 or 5 columns
if (ncol(sample) == 5) {
  # If it has 5 columns, assume the last column is coverage
  message(paste("Sample file", sample_file, "has 5 columns. Assuming: seqnames, start, end , meth, cov."))
  colnames(sample) <- c("seqnames", "start", "end", "meth", "cov")
  # Filter out sites with coverage below the minimum threshold
  sample <- sample[sample$cov >= min_cov, ]
  sample$beta_value <- round(sample$meth / sample$cov * 100, 2)  # Convert to percentage
  sample <- sample[, c("seqnames", "start", "end", "beta_value")]
} else if (ncol(sample) == 4) {
  message(paste("Sample file", sample_file, "has 4 columns. Assuming: seqnames, start, end, beta_value."))
  colnames(sample) <- c("seqnames", "start", "end", "beta_value")
} else {
  stop(paste("Sample file", sample_file, "does not have 4 or 5 columns."))
}

sample_sites <- nrow(sample)
sample_granges <- makeGRangesFromDataFrame(sample, keep.extra.columns = TRUE)

# Find overlaps
ol <- findOverlaps(atlas_granges, sample_granges)
ol <- as.data.frame(ol)

if (!allow_feature_overlap) {
  ol <- ol %>%
    distinct(subjectHits, .keep_all = TRUE)
}

# Extract methylation values for overlapping sample regions
overlap_meth <- split(sample_granges$beta_value[ol$subjectHits], ol$queryHits)

#  NA for bins that have no overlaps
sample_binned <- rep(NA, length(atlas_granges))

# Fill in mean beta values for bins with overlaps
sample_binned[as.integer(names(overlap_meth))] <- sapply(overlap_meth, function(x) mean(x, na.rm = TRUE)) / 100
overlap_sites <- sum(!is.na(sample_binned))

# Deconvolution
atlas_matrix <- atlases %>%
  select(-chr, -start, -end) %>%
  as.matrix()
idx <- !is.na(sample_binned) & rowSums(is.na(atlas_matrix)) == 0
fit <- suppressWarnings(lsei(
  atlas_matrix[idx, ],
  sample_binned[idx],
  G = diag(ncol(atlas_matrix)),
  H = rep(0, ncol(atlas_matrix))
))

# -----------------------------
# Write per-sample results & coverage
# -----------------------------

# Save deconvolution result (normalized cell-type fractions)
out_file <- file.path(output_dir, paste0(sample_name, "_deconvolution.csv"))
write.table(fit$X * 100 / sum(fit$X),
            file = out_file,
            sep = ",", col.names = FALSE, quote = FALSE)

# Save per-sample coverage row
cov_df <- data.frame(
  sample        = sample_name,
  sample_sites  = sample_sites,
  overlap_sites = overlap_sites,
  atlas_length  = nrow(atlases)
)
write.table(cov_df,
            file = file.path(output_dir, paste0(sample_name, "_coverage.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

#Quit R session
quit(save = "no")
