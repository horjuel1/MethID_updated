#!/usr/bin/env Rscript

#load libraries
suppressWarnings(suppressMessages(library(GenomicRanges)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(gridExtra)))
suppressWarnings(suppressMessages(library(limSolve)))
suppressWarnings(suppressMessages(library(argparser)))


#------------------------------
# Receive Arguments
#------------------------------

p <- arg_parser("Deconvolve methylation atlas from a sample sheet")
p <- add_argument(p, "--map_file", type = "character", help = "Atlas map file (TSV format)")
p <- add_argument(p, "--output_dir", type = "character", help = "Output directory")
p <- add_argument(p, "--sample_sheet", type = "character", help = "Sample sheet with id, label, and bam_path columns")
p <- add_argument(p, "--input_dir", type = "character", help = "Directory containing the intermediate shortbed files")
p <- add_argument(p, "--threads", type = "integer", default = 1, help = "Number of threads for data.table")
p <- add_argument(p, "--allow_feature_overlap", type = "logical", default = TRUE, help = "If FALSE, do not allow features to overlap in the atlas")
p <- add_argument(p, "--overwrite", default = TRUE, help = "Overwrite existing files in the output directory")
p <- add_argument(p, "--min_cov", type = "integer", default = 2, help = "Minimum coverage for a site to be considered [2]")

args <- parse_args(p)
map_file <- args$map_file
output_dir <- args$output_dir
sample_sheet_file <- args$sample_sheet
input_dir <- args$input_dir
threads <- args$threads
allow_feature_overlap <- args$allow_feature_overlap
overwrite <- args$overwrite
min_cov <- args$min_cov

sample_sheet <- fread(sample_sheet_file)
setDTthreads(threads = threads)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Process sample sheet to locate input files
# -----------------------------

# Replace spaces and slashes in the sample names with underscores to match Nextflow's sanitization
sample_sheet[, cleaned_id := gsub("[[:space:]/]", "_", id)]
sample_sheet[, cleaned_label := gsub("[[:space:]/]", "_", label)]

# Create the keys (unique sample identifiers) by combining cleaned 'id' and cleaned 'label' columns
keys <- paste(sample_sheet$cleaned_id, sample_sheet$cleaned_label, sep = "_")

# Create the values (full paths to shortbed files)
values <- file.path(input_dir, paste0(keys, "_short.bed"))

# Create a named list where names are sample IDs and values are file paths
samples <- setNames(as.list(values), keys)


# -----------------------------------------------
# Generate GRanges, find overlaps, and deconvolve
# -----------------------------------------------

#Read and convert atlas to GRanges
atlases <- fread(map_file, data.table = FALSE)
atlas_granges <- makeGRangesFromDataFrame(atlases, seqnames.field = "chr", start.field = "start", end.field = "end")

#Track lengths for coverage
n <- length(names(samples))
cov_df <- data.frame(
  row.names   = names(samples),
  sample_sites  = rep(NA, n),
  overlap_sites = rep(NA, n),
  atlas_length  = rep(nrow(atlases), n)
)

#Loop over each sample
for (sample_name in names(samples)){
    sample_file <- samples[[sample_name]]
  
  out_file <- file.path(output_dir, paste0(sample_name, "_deconvolution.csv"))

  # Check if the output file already exists and overwrite is FALSE
  if (file.exists(out_file) && !overwrite) {
    message(paste("Output file already exists:", out_file, "- skipping."))
    next
  }

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

  cov_df[sample_name, "sample_sites"] <- nrow(sample)
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
  cov_df[sample_name, "overlap_sites"] <- sum(!is.na(sample_binned))

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
# Write Results & Manifest & Coverage
# -----------------------------

  # Save results
  write.table(fit$X * 100 / sum(fit$X),
              file = out_file,
              sep = ",", col.names = FALSE, quote = FALSE)
}

#Manifest file
manifest <- c(
  "Atlas file" = map_file,
  "Output directory" = output_dir,
  "Sample sheet" =  sample_sheet_file
)
writeLines(paste(names(manifest), manifest, sep = ": "), file.path(output_dir, paste0("manifest_", Sys.Date(), ".txt")))

#Coverage file
write.table(cov_df,
              file = file.path(output_dir, paste0("coverage_", Sys.Date(), ".tsv")),
              sep = "\t",
              col.names = NA, # Use NA for row.names = TRUE
              quote = FALSE)

# -----------------------------
# Write Combined Results
# -----------------------------

#Let's generate a combined table of all the deconvolved samples
deconvolved_files <- list.files(output_dir, pattern = "_deconvolution.csv", full.names = TRUE)

if (length(deconvolved_files) > 0) {
    deconvolved_list <- lapply(deconvolved_files, function(x) {
      sample_name <- gsub("_deconvolution.csv", "", basename(x))
      deconvolved_data <- read.csv(x, header = FALSE)
      # Set row names from the first column
      rownames(deconvolved_data) <- deconvolved_data[, 1]
      deconvolved_data <- deconvolved_data[, -1, drop = FALSE]
      # Rename the remaining column
      colnames(deconvolved_data) <- sample_name
      return(deconvolved_data)
    })
    final_table <- do.call(cbind, deconvolved_list) %>% round(2)

    fwrite(final_table,
          file = file.path(output_dir, paste0("deconvolved_samples", ".csv")),
          sep = ",",
          col.names = TRUE,
          row.names = TRUE
    )

    #Remove all individual deconvolution files
    file.remove(deconvolved_files)
} else {
    message("No deconvolution files were generated to combine.")
}

#Quit R session
quit(save = "no")