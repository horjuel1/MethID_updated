#!/usr/bin/env Rscript

# Marker selection via Gap Specificity Score (GSS), replacing the ANOVA
# F-statistic, per doc 10's design -- grounded in a 2026 Communications
# Biology paper ("Guidelines on optimizing DNA methylation reference
# panels for cell-type deconvolution") that found effect-size/specificity
# selection beats both plain significance testing and machine learning
# (IDOL/Elastic Net/Random Forest overfit at typical reference-panel
# sample sizes).
#
# Four stages, matching doc 10 exactly:
#
#   1. PRIMARY CRITERION -- GSS, one-vs-rest worst-case gap. For CpG i,
#      candidate cell type k: GSS = max(own_beta - max(everyone else),
#      min(everyone else) - own_beta) -- the guaranteed minimum separation
#      from the closest competitor, in whichever direction (hyper- or
#      hypomethylated) is actually informative. Structurally the one-vs-
#      rest version of what atlas_separability.r already computes
#      pairwise (min_diff = pmin(diff_if_methylated, diff_if_unmethylated));
#      this reuses that same worst-case-gap idea as the primary SELECTION
#      criterion up front, not a post-hoc diagnostic.
#
#   2. SECONDARY CONSISTENCY FILTER -- a count-based (proportions) test of
#      whether a high-GSS gap is actually well-supported by real coverage,
#      not a point estimate from thin data. Needs --coverage_file
#      (bin/build_atlas.r's coverage sidecar: <cell_type>_cov columns) --
#      without it, this stage is skipped and every candidate reports
#      secondary_test_pvalue = NA. Reported, not a hard exclusion, unless
#      --require_secondary_pass is set -- matches doc 10's framing of this
#      as a consistency check on top of GSS, not a competing primary
#      criterion (doc 9's ANOVA-replacement role, demoted here to
#      secondary per doc 10).
#
#   3. CpG-DENSITY-AWARE FINAL PASS -- among GSS-qualifying candidates,
#      prefer genomic clusters of multiple qualifying CpGs within
#      --fragment_window_bp of each other over scattered singletons, since
#      deconvolve_em.r scores a molecule by summing evidence over whatever
#      CpGs one read actually covers -- doc 20 measured directly what this
#      fixes (1.5-2x fragment-length inflation among assigned molecules,
#      since a molecule needs enough overlapping markers to be assigned at
#      all, and short fragments have less genomic room to find several
#      scattered singletons). A preference (tie-break/bonus in ranking),
#      not a hard requirement -- an isolated singleton with a much larger
#      GSS still beats a mediocre cluster.
#
#   4. PAIRWISE GAP-FILLING -- one-vs-rest GSS can leave specific pairs
#      under-resolved even when each cell type individually looks
#      well-specified against the pooled rest. Optional second pass: for
#      any cell type whose nearest-neighbor separability
#      (atlas_separability.r's --nn_table_out) falls below
#      --weak_pair_threshold, compute the PAIRWISE worst-case gap
#      specifically against that neighbor (same formula as
#      atlas_separability.r, not one-vs-rest) and add top scorers even if
#      they didn't clear the generic one-vs-rest bar -- tagged
#      selection_reason="gap_fill" so they're distinguishable from the
#      primary panel.

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(argparser)))

p <- arg_parser("Marker selection via one-vs-rest Gap Specificity Score (GSS), replacing the ANOVA F-statistic, per doc 10")
p <- add_argument(p, "--map_file", type = "character", help = "Atlas map file (TSV/TSV.gz): chr, start, end, <cell_type_beta...> (e.g. bin/build_atlas.r's main output)")
p <- add_argument(p, "--coverage_file", type = "character", default = NA, help = "bin/build_atlas.r's coverage sidecar (<output>_coverage.tsv: <cell_type>_cov, <cell_type>_nsamples columns) -- needed for the secondary consistency filter; omit to skip it")
p <- add_argument(p, "--beta_eps", type = "double", default = 0.01, help = "Clip atlas beta values to [eps, 1-eps], matching deconvolve_em.r [0.01]")
p <- add_argument(p, "--min_gss", type = "double", default = 0.3, help = "Minimum GSS for primary-pass candidacy, per doc 10's numeric guidance [0.3]")
p <- add_argument(p, "--n_markers_per_celltype", type = "integer", default = 50, help = "Target markers per cell type, per doc 10's numeric guidance [50]")
p <- add_argument(p, "--candidate_pool_multiplier", type = "integer", default = 5, help = "How many GSS-ranked candidates per cell type (as a multiple of --n_markers_per_celltype) to carry into the density-aware/secondary-filter stages before final selection -- keeps those slower stages off the full atlas [5]")
p <- add_argument(p, "--min_coverage", type = "integer", default = 10, help = "Minimum pooled coverage (own cell type AND pooled rest) for the secondary consistency test to be computed at all; below this the test is skipped (reported NA) rather than run on too little evidence [10]")
p <- add_argument(p, "--secondary_test_alpha", type = "double", default = 0.05, help = "p-value threshold for passes_secondary_test [0.05]")
p <- add_argument(p, "--require_secondary_pass", type = "flag", help = "Make the secondary consistency test a hard filter (drop candidates that fail it) rather than just reported -- only meaningful with --coverage_file given")
p <- add_argument(p, "--fragment_window_bp", type = "integer", default = 200, help = "Genomic window (bp) for the CpG-density-aware clustering pass, sized to the empirical cfDNA mononucleosome mode (doc 20: ~150-200bp), not an arbitrary spacing [200]")
p <- add_argument(p, "--min_cluster_size", type = "integer", default = 2, help = "Minimum candidates within --fragment_window_bp to count as a preferred cluster [2]")
p <- add_argument(p, "--nn_table", type = "character", default = NA, help = "atlas_separability.r's --nn_table_out CSV -- enables stage 4 (pairwise gap-filling). Omit to skip.")
p <- add_argument(p, "--weak_pair_threshold", type = "double", default = 0.05, help = "Nearest-neighbor separability score (from --nn_table) below which a cell type triggers pairwise gap-filling [0.05]")
p <- add_argument(p, "--n_gap_fill_markers", type = "integer", default = 10, help = "Extra pairwise-gap-filling markers to add per weak pair, per cell type in the pair [10]")
p <- add_argument(p, "--output_file", type = "character", default = "gss_markers.tsv", help = "Output TSV: selected markers, one row per (chrom, start, cell_type)")
p <- add_argument(p, "--summary_file", type = "character", default = NA, help = "Output TSV: per-cell-type selection summary (n selected, fraction clearing GSS>0.5, etc.) [<output_file> with _summary.tsv suffix]")

args <- parse_args(p)
eps <- args$beta_eps

# -----------------------------------------------
# Load atlas
# -----------------------------------------------

atlas <- fread(args$map_file, data.table = TRUE)
cell_types <- setdiff(colnames(atlas), c("chr", "start", "end"))
K <- length(cell_types)
if (K < 2) stop("Need at least 2 cell types in --map_file to compute a one-vs-rest gap.")

complete_idx <- rowSums(is.na(atlas[, ..cell_types])) == 0
atlas <- atlas[complete_idx]
cat(sprintf("Atlas: %d cell types, %d complete-row CpGs\n", K, nrow(atlas)))

beta_mat <- as.matrix(atlas[, ..cell_types])
beta_mat <- pmin(pmax(beta_mat, eps), 1 - eps)
colnames(beta_mat) <- cell_types

cov_mat <- NULL
if (!is.na(args$coverage_file)) {
  cov_dt <- fread(args$coverage_file, data.table = TRUE)
  cov_cols <- paste0(cell_types, "_cov")
  missing_cov <- setdiff(cov_cols, colnames(cov_dt))
  if (length(missing_cov) > 0) {
    stop("--coverage_file is missing expected column(s): ", paste(missing_cov, collapse = ", "))
  }
  cov_dt <- merge(atlas[, .(chr, start, end)], cov_dt, by = c("chr", "start", "end"))
  if (nrow(cov_dt) != nrow(atlas)) stop("--coverage_file doesn't cover every row of --map_file after merging on chr/start/end -- must be build_atlas.r's own sidecar for this exact atlas.")
  cov_mat <- as.matrix(cov_dt[, ..cov_cols])
  colnames(cov_mat) <- cell_types
  cat("Coverage file loaded -- secondary consistency filter enabled.\n")
} else {
  cat("No --coverage_file given -- secondary consistency filter skipped (all candidates get secondary_test_pvalue = NA).\n")
}

# -----------------------------------------------
# Stage 1: one-vs-rest GSS, vectorized (fast row-wise max/min over K-1
# columns per cell type, via do.call(pmax/pmin, ...) rather than apply()).
# -----------------------------------------------

cat("Computing one-vs-rest GSS ...\n")
gss_list <- vector("list", K)
for (k in seq_len(K)) {
  other_cols <- setdiff(seq_len(K), k)
  max_others <- do.call(pmax, as.data.frame(beta_mat[, other_cols, drop = FALSE]))
  min_others <- do.call(pmin, as.data.frame(beta_mat[, other_cols, drop = FALSE]))
  own <- beta_mat[, k]

  gap_hyper <- own - max_others
  gap_hypo  <- min_others - own
  gss <- pmax(gap_hyper, gap_hypo)
  direction <- ifelse(gap_hyper >= gap_hypo, "hyper", "hypo")

  gss_list[[k]] <- data.table(
    chrom = atlas$chr, start = atlas$start, end = atlas$end,
    cell_type = cell_types[k], gss = gss, direction = direction,
    own_beta = own
  )
}
gss_dt <- rbindlist(gss_list)

# -----------------------------------------------
# Primary candidate pool: per cell type, top (n_markers_per_celltype *
# candidate_pool_multiplier) by GSS among those clearing --min_gss --
# keeps the slower per-candidate stages (secondary test, clustering) off
# the full atlas while still giving stage 3's clustering pass real room to
# find genomic neighbors.
# -----------------------------------------------

gss_dt <- gss_dt[gss >= args$min_gss]
if (nrow(gss_dt) == 0) stop("No CpG clears --min_gss=", args$min_gss, " for any cell type -- lower --min_gss or check --map_file's beta values.")

setorder(gss_dt, cell_type, -gss)
pool_size <- args$n_markers_per_celltype * args$candidate_pool_multiplier
candidates <- gss_dt[, .SD[seq_len(min(.N, pool_size))], by = cell_type]
cat(sprintf("Primary candidate pool: %d CpG-cell_type pairs across %d cell type(s) (top %d by GSS each, min_gss=%.2f)\n", nrow(candidates), uniqueN(candidates$cell_type), pool_size, args$min_gss))

# -----------------------------------------------
# Stage 2: secondary consistency filter -- two-proportion test, own cell
# type's pooled (meth, cov) vs. the SUM of everyone else's pooled (meth,
# cov), at each candidate CpG. Only run if --coverage_file was given; only
# meaningful (and only run) where both sides clear --min_coverage.
# -----------------------------------------------

if (!is.null(cov_mat)) {
  cat("Running secondary consistency test on the candidate pool ...\n")
  row_lookup <- match(paste(candidates$chrom, candidates$start), paste(atlas$chr, atlas$start))
  candidates[, secondary_test_pvalue := NA_real_]
  for (i in seq_len(nrow(candidates))) {
    k <- match(candidates$cell_type[i], cell_types)
    r <- row_lookup[i]
    own_cov <- cov_mat[r, k]
    rest_cov <- sum(cov_mat[r, -k])
    if (is.na(own_cov) || is.na(rest_cov) || own_cov < args$min_coverage || rest_cov < args$min_coverage) next
    own_meth <- round(beta_mat[r, k] * own_cov)
    # pooled "rest" methylated count: coverage-weighted mean beta across
    # the other cell types at this row, times pooled rest coverage
    rest_beta_mean <- weighted.mean(beta_mat[r, -k], w = cov_mat[r, -k])
    rest_meth <- round(rest_beta_mean * rest_cov)
    pv <- tryCatch(
      prop.test(c(own_meth, rest_meth), c(own_cov, rest_cov), correct = TRUE)$p.value,
      error = function(e) NA_real_, warning = function(w) suppressWarnings(prop.test(c(own_meth, rest_meth), c(own_cov, rest_cov))$p.value)
    )
    candidates[i, secondary_test_pvalue := pv]
  }
  candidates[, passes_secondary_test := !is.na(secondary_test_pvalue) & secondary_test_pvalue < args$secondary_test_alpha]
  n_tested <- sum(!is.na(candidates$secondary_test_pvalue))
  n_passed <- sum(candidates$passes_secondary_test, na.rm = TRUE)
  cat(sprintf("Secondary test: %d/%d candidates had enough coverage to test, %d passed at alpha=%.2f\n", n_tested, nrow(candidates), n_passed, args$secondary_test_alpha))
  if (isTRUE(args$require_secondary_pass)) {
    n_before <- nrow(candidates)
    candidates <- candidates[is.na(secondary_test_pvalue) == FALSE & passes_secondary_test == TRUE]
    cat(sprintf("--require_secondary_pass: dropped %d/%d candidates that failed or had too little coverage to test\n", n_before - nrow(candidates), n_before))
  }
} else {
  candidates[, `:=`(secondary_test_pvalue = NA_real_, passes_secondary_test = NA)]
}

# -----------------------------------------------
# Stage 3: CpG-density-aware pass. Within each cell type's candidate pool,
# cluster by genomic proximity (--fragment_window_bp), flag membership in
# a qualifying cluster (>= --min_cluster_size candidates within the
# window), and use that as a ranking tie-break/bonus ahead of raw GSS --
# not a hard requirement.
# -----------------------------------------------

cat("Running CpG-density-aware clustering pass ...\n")
candidates[, in_dense_cluster := FALSE]
candidates[, cluster_id := NA_character_]
for (ct in unique(candidates$cell_type)) {
  idx <- which(candidates$cell_type == ct)
  sub <- candidates[idx]
  setorder(sub, chrom, start)
  # Greedy genomic clustering: group consecutive (by position) candidates
  # on the same chromosome that fall within --fragment_window_bp of the
  # cluster's own start.
  sub[, grp := 0L]
  grp_id <- 0L
  cluster_start <- -Inf
  cur_chrom <- ""
  for (i in seq_len(nrow(sub))) {
    if (sub$chrom[i] != cur_chrom || sub$start[i] - cluster_start > args$fragment_window_bp) {
      grp_id <- grp_id + 1L
      cluster_start <- sub$start[i]
      cur_chrom <- sub$chrom[i]
    }
    sub[i, grp := grp_id]
  }
  grp_sizes <- sub[, .N, by = grp]
  qualifying_grps <- grp_sizes[N >= args$min_cluster_size, grp]
  sub[, in_dense_cluster := grp %in% qualifying_grps]
  sub[, cluster_id := paste0(ct, "_cl", grp)]
  candidates[idx[order(sub$chrom, sub$start)], `:=`(in_dense_cluster = sub$in_dense_cluster, cluster_id = sub$cluster_id)]
}
n_clustered <- sum(candidates$in_dense_cluster)
cat(sprintf("%d/%d candidates fall in a qualifying dense cluster (>= %d candidates within %dbp)\n", n_clustered, nrow(candidates), args$min_cluster_size, args$fragment_window_bp))

# Final ranking within each cell type: dense-cluster membership first
# (preference, not requirement), then GSS descending; take top
# --n_markers_per_celltype.
setorder(candidates, cell_type, -in_dense_cluster, -gss)
primary_markers <- candidates[, .SD[seq_len(min(.N, args$n_markers_per_celltype))], by = cell_type]
primary_markers[, selection_reason := "primary_gss"]

# -----------------------------------------------
# Stage 4 (optional): pairwise gap-filling for weak nearest-neighbor pairs.
# -----------------------------------------------

gap_fill_markers <- data.table()
if (!is.na(args$nn_table)) {
  nn <- fread(args$nn_table, data.table = TRUE)
  weak <- nn[separability_score < args$weak_pair_threshold]
  if (nrow(weak) == 0) {
    cat(sprintf("No cell type falls below --weak_pair_threshold=%.3f in --nn_table -- no pairwise gap-filling needed.\n", args$weak_pair_threshold))
  } else {
    cat(sprintf("Pairwise gap-filling for %d weak pair(s) (nearest-neighbor separability < %.3f):\n", nrow(weak), args$weak_pair_threshold))
    log_b <- log(beta_mat); log_1mb <- log(1 - beta_mat)
    done_pairs <- character(0)
    gap_fill_list <- list()
    for (i in seq_len(nrow(weak))) {
      ct_a <- weak$cell_type[i]; ct_b <- weak$nearest_neighbor[i]
      pair_key <- paste(sort(c(ct_a, ct_b)), collapse = "||")
      if (pair_key %in% done_pairs) next
      done_pairs <- c(done_pairs, pair_key)
      ia <- match(ct_a, cell_types); ib <- match(ct_b, cell_types)
      if (is.na(ia) || is.na(ib)) { cat(sprintf("  Skipping %s vs %s -- not both in --map_file\n", ct_a, ct_b)); next }

      diff_if_1 <- abs(log_b[, ia] - log_b[, ib])
      diff_if_0 <- abs(log_1mb[, ia] - log_1mb[, ib])
      pair_min_diff <- pmin(diff_if_1, diff_if_0)
      top_idx <- head(order(-pair_min_diff), args$n_gap_fill_markers)

      cat(sprintf("  %s vs %s: adding top %d pairwise markers (median score %.3f) to both\n", ct_a, ct_b, length(top_idx), median(pair_min_diff[top_idx])))
      for (ct in c(ct_a, ct_b)) {
        k <- match(ct, cell_types)
        gap_fill_list[[length(gap_fill_list) + 1]] <- data.table(
          chrom = atlas$chr[top_idx], start = atlas$start[top_idx], end = atlas$end[top_idx],
          cell_type = ct, gss = pair_min_diff[top_idx], direction = NA_character_, own_beta = beta_mat[top_idx, k],
          secondary_test_pvalue = NA_real_, passes_secondary_test = NA, in_dense_cluster = NA, cluster_id = NA_character_,
          selection_reason = sprintf("gap_fill(%s_vs_%s)", ct_a, ct_b)
        )
      }
    }
    if (length(gap_fill_list) > 0) gap_fill_markers <- rbindlist(gap_fill_list)
  }
} else {
  cat("No --nn_table given -- pairwise gap-filling (stage 4) skipped.\n")
}

# -----------------------------------------------
# Combine, write output
# -----------------------------------------------

final_markers <- rbindlist(list(primary_markers, gap_fill_markers), use.names = TRUE)
# A CpG-cell_type pair added by gap-filling that was already in the
# primary panel stays as "primary_gss" -- gap-filling is meant to ADD
# coverage for weak pairs, not relabel markers already earned on merit.
final_markers[, is_primary := selection_reason == "primary_gss"]
setorder(final_markers, cell_type, -is_primary, -gss)
final_markers <- unique(final_markers, by = c("chrom", "start", "cell_type"))
final_markers[, is_primary := NULL]

fwrite(final_markers, args$output_file, sep = "\t")

summary_file <- if (is.na(args$summary_file)) sub("\\.[^.]+$", "_summary.tsv", args$output_file) else args$summary_file
summary_dt <- final_markers[, .(
  n_selected = .N,
  n_primary = sum(selection_reason == "primary_gss"),
  n_gap_fill = sum(selection_reason != "primary_gss"),
  frac_gss_above_0.5 = round(mean(gss > 0.5, na.rm = TRUE), 3),
  frac_in_dense_cluster = round(mean(in_dense_cluster, na.rm = TRUE), 3),
  median_gss = round(median(gss), 3)
), by = cell_type]
setorder(summary_dt, cell_type)
fwrite(summary_dt, summary_file, sep = "\t")

cat("\nPer-cell-type selection summary:\n")
print(summary_dt)
cat(sprintf("\nWrote %s (%d markers) and %s\n", args$output_file, nrow(final_markers), summary_file))
