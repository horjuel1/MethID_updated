#!/usr/bin/env Rscript

# Beta-binomial hyperparameter fit on top of build_atlas.r's pooled atlas.
# Reads build_atlas.r's per-cell-type count files
# (<output>_raw_counts/<cell_type>.tsv.gz: one meth + one cov column per
# replicate, "<cell_type>_<sample>_meth"/"..._cov") one cell type at a time,
# and its plain <output>.tsv for the cell-type list and CpG rows. Output
# tables have the atlas's rows, in its order.
#
# METHOD: per cell type k, a single shared overdispersion parameter phi_k
# (the intra-class correlation of the beta-binomial, i.e. the phi in
# a=mu(1-phi)/phi, b=(1-mu)(1-phi)/phi) is estimated via method of moments
# -- the same style used for beta-binomial dispersion in methylKit/RADMeth
# and equivalent to the classic Fleiss-Cuzick moment estimator for
# intra-cluster correlation on binary/proportion data:
#
#   1. For every CpG j with >=2 covered replicates for cell type k, compute
#      the Pearson-type statistic S_j = sum_s cov_sj*(p_sj-mu_j)^2/[mu_j(1-mu_j)]
#      against the CpG's own pooled proportion mu_j (already computed by
#      build_atlas.r as pooled_beta). Under pure-binomial (no overdispersion)
#      sampling, E[S_j] = n_j - 1.
#   2. Pool across all qualifying CpGs: VIF_k = sum(S_j) / sum(n_j - 1), a
#      variance-inflation-factor-style ratio (>1 under real overdispersion).
#   3. Convert to the intra-class-correlation parameterization via the
#      standard relation VIF ~= 1 + (nbar-1)*rho: rho_k = (VIF_k-1)/(nbar_k-1),
#      clipped to [--min_phi, --max_phi].
#
# Cell types with too few multi-replicate CpGs to estimate their own phi_k
# (e.g. a single-replicate cell type) fall back to the coverage-weighted
# pooled phi across every OTHER cell type that could estimate one, flagged
# via used_global_fallback in atlas_dispersion.tsv rather than silently
# guessing or crashing.
#
# Why not plain a_kj = mu_kj*(1-phi_k)/phi_k: that gives every CpG
# for a cell type the SAME total pseudo-count weight (1-phi_k)/phi_k
# regardless of that CpG's own pooled coverage, so Var[beta_kj] would not
# shrink as coverage grows. This script
# instead treats phi_k's implied pseudo-count weight as an
# empirical-Bayes PRIOR (a0_kj = mu_kj*m_k, b0_kj = (1-mu_kj)*m_k, with
# m_k = (1-phi_k)/phi_k), then forms the per-CpG POSTERIOR by adding that
# CpG's own already-pooled (meth_kj, cov_kj) on top:
#   a_kj = a0_kj + meth_kj,  b_kj = b0_kj + (cov_kj - meth_kj)
# Because the prior mean is set to the CpG's own mu_kj, E[beta_kj] =
# a_kj/(a_kj+b_kj) is unchanged (still exactly mu_kj, the same point
# estimate atlas.tsv already has) -- only Var[beta_kj] moves, and it now
# correctly shrinks as a CpG's own real pooled coverage grows, which is
# what the atlas-completeness / "can this atlas value be trusted" use case
# actually needs.

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(argparser)))

p <- arg_parser("Fit beta-binomial hyperparameters (a_kj, b_kj, E[beta_kj], Var[beta_kj]) on top of build_atlas.r's pooled atlas")
p <- add_argument(p, "--atlas_file",      type = "character", help = "build_atlas.r's pooled <output>.tsv -- authoritative cell-type label list and pooled beta (E[beta_kj])")
p <- add_argument(p, "--raw_counts_dir", type = "character", help = "build_atlas.r's <output>_raw_counts directory (one <cell_type>.tsv.gz per cell type)")
p <- add_argument(p, "--output_prefix",   type = "character", default = "atlas_beta_binomial", help = "Output file prefix [atlas_beta_binomial]")
p <- add_argument(p, "--min_coverage",    type = "double",  default = 1, help = "Minimum pooled coverage for a CpG's mean/variance to be reported (NA below this), matching build_atlas.r's own --min_coverage convention [1]")
p <- add_argument(p, "--min_df_for_dispersion", type = "integer", default = 30, help = "Minimum sum(total_cov_j - n_j) coverage-based degrees of freedom across multi-replicate CpGs required to estimate a cell type's own phi_k; below this, fall back to the global pooled phi [30]")
p <- add_argument(p, "--min_phi", type = "double", default = 1e-4, help = "Floor on phi_k (rho), keeps a_kj/b_kj finite and Var[beta_kj] > 0 [1e-4]")
p <- add_argument(p, "--max_phi", type = "double", default = 0.999, help = "Ceiling on phi_k (rho) [0.999]")

args <- parse_args(p)

atlas <- fread(args$atlas_file, data.table = TRUE, select = 1:3)
cell_types <- setdiff(names(fread(args$atlas_file, nrows = 0)), c("chr", "start", "end"))
if (length(cell_types) == 0) stop("--atlas_file has no cell-type columns beyond chr/start/end.")
pos <- atlas[, .(chr, start, end)]
rm(atlas)

# One cell type's per-replicate counts -> pooled totals plus the matrices the
# dispersion estimate needs. Columns are "<cell_type>_<sample>_meth/_cov".
read_cell_type <- function(ct) {
  f <- file.path(args$raw_counts_dir, paste0(ct, ".tsv.gz"))
  if (!file.exists(f)) stop("No raw counts file for cell type '", ct, "': ", f)
  raw <- fread(f, data.table = TRUE)
  meth_cols <- grep("_meth$", names(raw), value = TRUE)
  cov_cols  <- grep("_cov$",  names(raw), value = TRUE)
  meth_keys <- sub("_meth$", "", meth_cols)
  cov_keys  <- sub("_cov$",  "", cov_cols)
  keys <- intersect(meth_keys, cov_keys)
  if (length(keys) == 0) stop("Cell type '", ct, "' has no matched meth/cov column pairs in ", f)
  meth_mat <- as.matrix(raw[, meth_cols[match(keys, meth_keys)], with = FALSE])
  cov_mat  <- as.matrix(raw[, cov_cols[match(keys, cov_keys)],  with = FALSE])
  total_meth <- rowSums(meth_mat, na.rm = TRUE)
  total_cov  <- rowSums(cov_mat,  na.rm = TRUE)
  list(pos = raw[, .(chr, start, end)], n_replicates = length(keys),
       mu = ifelse(total_cov > 0, total_meth / total_cov, NA_real_),
       total_meth = total_meth, total_cov = total_cov,
       n_covered = rowSums(!is.na(cov_mat) & cov_mat > 0),
       meth_mat = meth_mat, cov_mat = cov_mat)
}

# -----------------------------------------------
# Method-of-moments phi_k per cell type
# -----------------------------------------------

# For observation (replicate s, CpG j) with coverage cov_sj, the
# beta-binomial model gives Var(Y_sj) = cov_sj*mu_j*(1-mu_j)*[1+(cov_sj-1)*phi],
# so the Pearson residual cov_sj*(p_sj-mu_j)^2/[mu_j(1-mu_j)] has expectation
# 1+(cov_sj-1)*phi per observation. Summed over a CpG's covered replicates:
#   E[S_j] = (n_j - 1) + phi * (total_cov_j - n_j)
# so phi_hat = [sum_j(S_j - (n_j-1))] / [sum_j(total_cov_j - n_j)].
estimate_phi <- function(fit) {
  q <- which(fit$n_covered >= 2 & !is.na(fit$mu) & fit$mu > 0 & fit$mu < 1)
  if (length(q) == 0) return(list(sum_S = NA_real_, sum_meandf = 0, sum_covdf = 0))
  cov_q <- fit$cov_mat[q, , drop = FALSE]
  p_q <- fit$meth_mat[q, , drop = FALSE] / cov_q
  m <- fit$mu[q]
  terms <- cov_q * (p_q - m)^2
  terms[is.na(cov_q) | cov_q <= 0] <- 0
  S <- rowSums(terms) / (m * (1 - m))
  list(sum_S = sum(S), sum_meandf = sum(fit$n_covered[q] - 1), sum_covdf = sum(fit$total_cov[q] - fit$n_covered[q]))
}

# Pass 1: dispersion statistics, one cell type at a time
phi_stats <- list()
n_replicates <- integer(0)
for (ct in cell_types) {
  fit <- read_cell_type(ct)
  phi_stats[[ct]] <- estimate_phi(fit)
  n_replicates[ct] <- fit$n_replicates
  rm(fit); invisible(gc())
}

moments_to_phi <- function(s, min_phi, max_phi) {
  if (is.na(s$sum_S) || is.na(s$sum_covdf) || s$sum_covdf <= 0) return(NA_real_)
  raw <- (s$sum_S - s$sum_meandf) / s$sum_covdf
  min(max(raw, min_phi), max_phi)
}

phi_k <- vapply(phi_stats, moments_to_phi, numeric(1), min_phi = args$min_phi, max_phi = args$max_phi)
names(phi_k) <- cell_types

can_estimate <- vapply(phi_stats, function(s) !is.na(s$sum_covdf) && s$sum_covdf >= args$min_df_for_dispersion, logical(1))

if (!any(can_estimate)) {
  stop("No cell type has >= --min_df_for_dispersion (", args$min_df_for_dispersion,
       ") degrees of freedom for its own phi_k estimate -- lower --min_df_for_dispersion or check --raw_counts_file has >=2 replicates for at least one cell type.")
}

# Global fallback: coverage(df)-weighted pooled phi across every cell type
# that COULD estimate its own, for any cell type that couldn't.
df_weights <- vapply(phi_stats, function(s) if (is.na(s$sum_covdf)) 0 else s$sum_covdf, numeric(1))
global_phi <- sum(phi_k[can_estimate] * df_weights[can_estimate]) / sum(df_weights[can_estimate])
global_phi <- min(max(global_phi, args$min_phi), args$max_phi)

used_global_fallback <- !can_estimate
if (any(used_global_fallback)) {
  message(sprintf("Cell type(s) with < %d degrees of freedom for their own phi_k -- using the df-weighted global pooled phi (%.4f) instead: %s",
                   args$min_df_for_dispersion, global_phi, paste(cell_types[used_global_fallback], collapse = ", ")))
  phi_k[used_global_fallback] <- global_phi
}

dispersion_dt <- data.table(
  cell_type = cell_types,
  phi_k = round(phi_k, 6),
  n_replicates = as.integer(n_replicates[cell_types]),
  dispersion_covdf = round(vapply(phi_stats, function(s) s$sum_covdf, numeric(1)), 1),
  used_global_fallback = used_global_fallback
)
fwrite(dispersion_dt, paste0(args$output_prefix, "_dispersion.tsv"), sep = "\t")
cat("Per-cell-type dispersion fit:\n"); print(dispersion_dt)

# -----------------------------------------------
# Per-CpG posterior a_kj, b_kj, E[beta_kj], Var[beta_kj], one column per
# cell type, same row order/shape as --atlas_file
# -----------------------------------------------

# Pass 2: per-CpG posterior a_kj, b_kj -> E[beta_kj], Var[beta_kj], aligned
# to the atlas rows (CpGs this cell type never covered stay NA).
mean_dt <- copy(pos)
var_dt  <- copy(pos)

for (ct in cell_types) {
  f <- read_cell_type(ct)
  m_k <- (1 - phi_k[[ct]]) / phi_k[[ct]]   # prior pseudo-count weight from phi_k
  a_post <- f$mu * m_k + f$total_meth
  b_post <- (1 - f$mu) * m_k + (f$total_cov - f$total_meth)
  e_beta <- a_post / (a_post + b_post)
  var_beta <- (a_post * b_post) / ((a_post + b_post)^2 * (a_post + b_post + 1))
  below_min_cov <- is.na(f$total_cov) | f$total_cov < args$min_coverage
  e_beta[below_min_cov] <- NA_real_
  var_beta[below_min_cov] <- NA_real_

  res <- f$pos[, `:=`(e = round(e_beta, 4), v = signif(var_beta, 4))]
  mean_dt[res, (ct) := i.e, on = .(chr, start, end)]
  var_dt[res,  (ct) := i.v, on = .(chr, start, end)]
  rm(f, res); invisible(gc())
}

fwrite(mean_dt, paste0(args$output_prefix, "_mean.tsv"), sep = "\t", quote = FALSE)
fwrite(var_dt,  paste0(args$output_prefix, "_var.tsv"),  sep = "\t", quote = FALSE)

cat(sprintf("\nWrote %s_mean.tsv, %s_var.tsv, %s_dispersion.tsv (%d CpGs x %d cell types)\n",
            args$output_prefix, args$output_prefix, args$output_prefix, nrow(pos), length(cell_types)))
