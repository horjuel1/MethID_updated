#!/usr/bin/env Rscript

# Beta-binomial hyperparameter fit on top of build_atlas.r's already-pooled
# counts (doc-8-era atlas layer, per the unified EM framework spec's layer
# A). build_atlas.r's --get_raw_atlas already retains everything this
# script needs -- per-replicate meth/cov counts, one column pair per
# (cell_type, sample) -- so this reads its <output>_raw_counts.tsv output
# rather than re-deriving anything from the mapped.bed inputs, and its
# plain <output>.tsv for the authoritative list of cell-type labels (raw
# count column names are "<cell_type>_<sample>_meth"/"..._cov", and a
# cell-type label can itself contain underscores, so cell types are
# matched by known-label prefix against atlas.tsv's own column names
# rather than parsed positionally, which would be ambiguous).
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
# REFINEMENT ON THE ORIGINAL SPEC WORTH FLAGGING EXPLICITLY: the spec's own
# formula (a_kj = mu_kj*(1-phi_k)/phi_k, applied directly) gives every CpG
# for a cell type the SAME total pseudo-count weight (1-phi_k)/phi_k
# regardless of that specific CpG's own pooled coverage -- which does not
# actually produce the "Var[beta_kj] shrinks toward 0 as coverage grows"
# behavior the spec's own narrative promises, since a real per-CpG
# coverage difference never enters a_kj/b_kj under that formula. This
# script instead treats phi_k's implied pseudo-count weight as an
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
p <- add_argument(p, "--raw_counts_file", type = "character", help = "build_atlas.r's <output>_raw_counts.tsv (requires --get_raw_atlas true when building the atlas)")
p <- add_argument(p, "--output_prefix",   type = "character", default = "atlas_beta_binomial", help = "Output file prefix [atlas_beta_binomial]")
p <- add_argument(p, "--min_coverage",    type = "double",  default = 1, help = "Minimum pooled coverage for a CpG's mean/variance to be reported (NA below this), matching build_atlas.r's own --min_coverage convention [1]")
p <- add_argument(p, "--min_df_for_dispersion", type = "integer", default = 30, help = "Minimum sum(total_cov_j - n_j) coverage-based degrees of freedom across multi-replicate CpGs required to estimate a cell type's own phi_k; below this, fall back to the global pooled phi [30]")
p <- add_argument(p, "--min_phi", type = "double", default = 1e-4, help = "Floor on phi_k (rho), keeps a_kj/b_kj finite and Var[beta_kj] > 0 [1e-4]")
p <- add_argument(p, "--max_phi", type = "double", default = 0.999, help = "Ceiling on phi_k (rho) [0.999]")

args <- parse_args(p)

atlas <- fread(args$atlas_file, data.table = TRUE)
cell_types <- setdiff(colnames(atlas), c("chr", "start", "end"))
if (length(cell_types) == 0) stop("--atlas_file has no cell-type columns beyond chr/start/end.")

raw_counts <- fread(args$raw_counts_file, data.table = TRUE)
count_cols <- setdiff(colnames(raw_counts), c("chr", "start", "end"))
meth_cols  <- grep("_meth$", count_cols, value = TRUE)
cov_cols   <- grep("_cov$",  count_cols, value = TRUE)
if (length(meth_cols) == 0 || length(cov_cols) == 0) {
  stop("No *_meth/*_cov columns found in --raw_counts_file. Expected build_atlas.r's <output>_raw_counts.tsv (--get_raw_atlas true).")
}

# Match each meth/cov column pair to its cell type by known-label prefix
# (longest label first, so e.g. "loyfer_Blood-T-CD4-Memory" claims its own
# columns before the shorter "loyfer_Blood-T-CD4" can).
cell_types_by_length <- cell_types[order(-nchar(cell_types))]
match_cell_type <- function(col, suffix) {
  for (ct in cell_types_by_length) {
    prefix <- paste0(ct, "_")
    if (startsWith(col, prefix) && endsWith(col, suffix)) return(ct)
  }
  NA_character_
}
meth_ct <- vapply(meth_cols, match_cell_type, character(1), suffix = "_meth")
cov_ct  <- vapply(cov_cols,  match_cell_type, character(1), suffix = "_cov")

unmatched <- c(meth_cols[is.na(meth_ct)], cov_cols[is.na(cov_ct)])
if (length(unmatched) > 0) {
  message(sprintf("Note: %d raw-count column(s) didn't match any --atlas_file cell type by prefix, dropped: %s",
                   length(unmatched), paste(head(unmatched, 5), collapse = ", ")))
}

pos <- raw_counts[, .(chr, start, end)]

# -----------------------------------------------
# Per cell type: pool meth/cov this script computes independently from the
# raw per-replicate columns (self-contained -- does not assume row order
# agrees with --atlas_file), fit phi_k, then form per-CpG posterior a_kj/b_kj.
# -----------------------------------------------

fit_one_cell_type <- function(ct) {
  ct_meth_cols <- meth_cols[!is.na(meth_ct) & meth_ct == ct]
  ct_cov_cols  <- cov_cols[!is.na(cov_ct) & cov_ct == ct]
  # pair by sample key (column name with the cell-type prefix and _meth/_cov
  # suffix stripped)
  meth_keys <- sub(paste0("^", ct, "_"), "", sub("_meth$", "", ct_meth_cols))
  cov_keys  <- sub(paste0("^", ct, "_"), "", sub("_cov$",  "", ct_cov_cols))
  keys <- intersect(meth_keys, cov_keys)
  if (length(keys) == 0) stop("Cell type '", ct, "' has no matched meth/cov column pairs in --raw_counts_file.")

  meth_mat <- as.matrix(raw_counts[, ct_meth_cols[match(keys, meth_keys)], with = FALSE])
  cov_mat  <- as.matrix(raw_counts[, ct_cov_cols[match(keys, cov_keys)],  with = FALSE])

  total_meth <- rowSums(meth_mat, na.rm = TRUE)
  total_cov  <- rowSums(cov_mat,  na.rm = TRUE)
  mu <- ifelse(total_cov > 0, total_meth / total_cov, NA_real_)

  n_replicates <- length(keys)
  n_covered <- rowSums(!is.na(cov_mat) & cov_mat > 0)

  list(n_replicates = n_replicates, mu = mu, total_meth = total_meth, total_cov = total_cov,
       n_covered = n_covered, meth_mat = meth_mat, cov_mat = cov_mat)
}

fits <- setNames(lapply(cell_types, fit_one_cell_type), cell_types)

# -----------------------------------------------
# Method-of-moments phi_k per cell type
# -----------------------------------------------

# For observation (replicate s, CpG j) with coverage cov_sj, the
# beta-binomial model gives Var(Y_sj) = cov_sj*mu_j*(1-mu_j)*[1+(cov_sj-1)*phi],
# so the Pearson residual cov_sj*(p_sj-mu_j)^2/[mu_j(1-mu_j)] has expectation
# 1+(cov_sj-1)*phi PER OBSERVATION -- the (n-1) inflation term is driven by
# each observation's own COVERAGE (reads within that replicate), not by how
# many replicate samples there are. Pooling S_j = sum_s of that residual
# across a CpG's covered replicates:
#   E[S_j] = n_j - 1 [1 df lost estimating mu_j from the same data] +
#            phi * (total_cov_j - n_j) [overdispersion contribution, scales
#            with total coverage, not replicate count]
# so phi is recovered from the EXCESS of S_j over its pure-binomial
# expectation, divided by the coverage-based term it multiplies:
#   phi_hat = [sum_j(S_j - (n_j-1))] / [sum_j(total_cov_j - n_j)]
# (An earlier version of this estimator divided by (n_bar_replicates - 1)
# instead -- i.e. used the number of replicate SAMPLES as the beta-binomial
# "trial size" -- which is a different, much smaller quantity than each
# observation's own coverage and inflated phi_hat by roughly
# coverage/replicate_count; caught by checking recovered phi against a
# known-truth synthetic simulation before this script was used anywhere.)
estimate_phi <- function(fit) {
  mu <- fit$mu
  meth_mat <- fit$meth_mat
  cov_mat  <- fit$cov_mat
  n_covered <- fit$n_covered
  total_cov <- fit$total_cov

  qualifying <- which(n_covered >= 2 & !is.na(mu) & mu > 0 & mu < 1)
  if (length(qualifying) == 0) return(list(phi = NA_real_, sum_S = NA_real_, sum_meandf = 0, sum_covdf = 0))

  S <- numeric(length(qualifying))
  for (i in seq_along(qualifying)) {
    r <- qualifying[i]
    covered <- which(!is.na(cov_mat[r, ]) & cov_mat[r, ] > 0)
    p_s <- meth_mat[r, covered] / cov_mat[r, covered]
    m <- mu[r]
    S[i] <- sum(cov_mat[r, covered] * (p_s - m)^2) / (m * (1 - m))
  }
  n_j <- n_covered[qualifying]
  cov_j <- total_cov[qualifying]
  sum_S <- sum(S)
  sum_meandf <- sum(n_j - 1)     # df lost estimating each CpG's own mu_j
  sum_covdf  <- sum(cov_j - n_j) # coverage-based term overdispersion acts on

  list(phi = NA_real_, sum_S = sum_S, sum_meandf = sum_meandf, sum_covdf = sum_covdf)
}

phi_stats <- lapply(fits, estimate_phi)

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
  n_replicates = vapply(fits, function(f) f$n_replicates, integer(1)),
  dispersion_covdf = round(vapply(phi_stats, function(s) s$sum_covdf, numeric(1)), 1),
  used_global_fallback = used_global_fallback
)
fwrite(dispersion_dt, paste0(args$output_prefix, "_dispersion.tsv"), sep = "\t")
cat("Per-cell-type dispersion fit:\n"); print(dispersion_dt)

# -----------------------------------------------
# Per-CpG posterior a_kj, b_kj, E[beta_kj], Var[beta_kj], one column per
# cell type, same row order/shape as --atlas_file
# -----------------------------------------------

mean_dt <- copy(pos)
var_dt  <- copy(pos)

for (ct in cell_types) {
  f <- fits[[ct]]
  mu <- f$mu
  m_k <- (1 - phi_k[[ct]]) / phi_k[[ct]]   # prior pseudo-count weight from phi_k
  a0 <- mu * m_k
  b0 <- (1 - mu) * m_k
  a_post <- a0 + f$total_meth
  b_post <- b0 + (f$total_cov - f$total_meth)

  e_beta <- a_post / (a_post + b_post)
  var_beta <- (a_post * b_post) / ((a_post + b_post)^2 * (a_post + b_post + 1))

  below_min_cov <- is.na(f$total_cov) | f$total_cov < args$min_coverage
  e_beta[below_min_cov] <- NA_real_
  var_beta[below_min_cov] <- NA_real_

  mean_dt[[ct]] <- round(e_beta, 4)
  var_dt[[ct]]  <- signif(var_beta, 4)
}

fwrite(mean_dt, paste0(args$output_prefix, "_mean.tsv"), sep = "\t", quote = FALSE)
fwrite(var_dt,  paste0(args$output_prefix, "_var.tsv"),  sep = "\t", quote = FALSE)

cat(sprintf("\nWrote %s_mean.tsv, %s_var.tsv, %s_dispersion.tsv (%d CpGs x %d cell types)\n",
            args$output_prefix, args$output_prefix, args$output_prefix, nrow(pos), length(cell_types)))
