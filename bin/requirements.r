options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  Ncpus = max(1, parallel::detectCores() - 1)
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

install_and_verify_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    install.packages(missing)
  }

  still_missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing) > 0) {
    stop(
      "Failed to install CRAN packages: ",
      paste(still_missing, collapse = ", ")
    )
  }
}

install_and_verify_bioc <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  }

  still_missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing) > 0) {
    stop(
      "Failed to install Bioconductor packages: ",
      paste(still_missing, collapse = ", ")
    )
  }
}

cran_packages <- c(
  "tidyverse",
  "data.table",
  "ggplot2",
  "gridExtra",
  "limSolve",
  "argparser"
)

bioc_packages <- c(
  "GenomicRanges",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19"
)

install_and_verify_cran(cran_packages)
install_and_verify_bioc(bioc_packages)

message("All R dependencies are installed.")
