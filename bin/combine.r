#!/usr/bin/env Rscript

#load libraries
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(argparser)))

#------------------------------
# Receive Arguments
#------------------------------

p <- arg_parser("Combine per-sample deconvolution and coverage results into summary tables")
p <- add_argument(p, "--input_dir", type = "character", default = ".", help = "Directory containing per-sample *_deconvolution.csv and *_coverage.tsv files")
p <- add_argument(p, "--output_dir", type = "character", default = ".", help = "Output directory")
p <- add_argument(p, "--map_file", type = "character", help = "Atlas map file (recorded in the manifest for provenance)")

args <- parse_args(p)
input_dir <- args$input_dir
output_dir <- args$output_dir
map_file <- args$map_file

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Combine deconvolution results
# -----------------------------

deconvolved_files <- list.files(input_dir, pattern = "_deconvolution\\.csv$", full.names = TRUE)

if (length(deconvolved_files) > 0) {
    deconvolved_list <- lapply(deconvolved_files, function(x) {
      sample_name <- gsub("_deconvolution.csv", "", basename(x))
      deconvolved_data <- read.csv(x, header = FALSE)
      # Set row names from the first column (cell types)
      rownames(deconvolved_data) <- deconvolved_data[, 1]
      deconvolved_data <- deconvolved_data[, -1, drop = FALSE]
      # Rename the remaining column to the sample name
      colnames(deconvolved_data) <- sample_name
      return(deconvolved_data)
    })
    final_table <- do.call(cbind, deconvolved_list) %>% round(2)

    fwrite(final_table,
          file = file.path(output_dir, "deconvolved_samples.csv"),
          sep = ",",
          col.names = TRUE,
          row.names = TRUE
    )
} else {
    message("No deconvolution files were found to combine.")
}

# -----------------------------
# Combine coverage results
# -----------------------------

coverage_files <- list.files(input_dir, pattern = "_coverage\\.tsv$", full.names = TRUE)

if (length(coverage_files) > 0) {
    cov_df <- rbindlist(lapply(coverage_files, fread))
    fwrite(cov_df,
          file = file.path(output_dir, paste0("coverage_", Sys.Date(), ".tsv")),
          sep = "\t")
}

# -----------------------------
# Write manifest
# -----------------------------

manifest <- c(
  "Atlas file" = map_file,
  "Output directory" = output_dir,
  "Date" = as.character(Sys.Date())
)
writeLines(paste(names(manifest), manifest, sep = ": "),
           file.path(output_dir, paste0("manifest_", Sys.Date(), ".txt")))

#Quit R session
quit(save = "no")
