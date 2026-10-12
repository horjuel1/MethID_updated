options(
  repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/2026-05-01"),
  BioC_mirror = "https://packagemanager.posit.co/bioconductor",
  HTTPUserAgent = sprintf(
    "R/%s R (%s)",
    getRversion(),
    paste(getRversion(), R.version$platform, R.version$arch, R.version$os)
  ),
  Ncpus = min(8, max(1, parallel::detectCores() - 1))
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
  "gridExtra",
  "limSolve",
  "argparser",
  "ggplot2",
  "jsonlite"
)

bioc_packages <- c(
  "GenomicRanges"
)

install_and_verify_cran(cran_packages)
install_and_verify_bioc(bioc_packages)

message("All R dependencies are installed.")
