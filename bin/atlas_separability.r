#!/usr/bin/env Rscript

# Per-CpG cell-type separability in a methylation atlas, used to sanity-check
# whether the EM's single-molecule confidence claims are plausible for a
# given cell-type pair (see deconvolve_em.r: molecules are scored as a
# per-CpG independent-Bernoulli log-likelihood against each cell type's
# atlas beta).
#
# For every pair of cell types (A, B) and every CpG, we compute the
# worst-case per-call log-likelihood difference a single observed call at
# that site could contribute toward telling A from B:
#   diff_if_methylated   = |log(betaA)   - log(betaB)|
#   diff_if_unmethylated = |log(1-betaA) - log(1-betaB)|
#   min_diff              = min(diff_if_methylated, diff_if_unmethylated)
# min_diff is "how much evidence does this CpG guarantee, regardless of
# which way the call actually goes" -- a site can look discriminating for
# one call direction and useless for the other, so we take the floor, not
# the average. Aggregating min_diff's median across all CpGs gives one
# number per pair: how confusable A and B are, structurally, in this atlas.
# "Nearest neighbor" per cell type is the pair with the smallest such score
# -- the hardest other cell type to rule out with a single CpG.
#
# Also renders a 2D density plot of raw beta values for the most-confusable
# pairs (or --density_pairs, if given), to visualize what "highly
# correlated across CpGs" looks like for a specific pair, plus a full
# pairwise confusability map (all cell types x all cell types, ordered by
# hierarchical clustering on the same score so confusable clusters of
# related cell types group together visually).

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(argparser)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(patchwork)))

p <- arg_parser("Per-CpG cell-type separability diagnostics for a methylation atlas")
p <- add_argument(p, "--map_file",           type = "character", help = "Atlas map file (TSV/TSV.gz): chr, start, end, <cell_type_beta...>")
p <- add_argument(p, "--beta_eps",           type = "double",  default = 0.01, help = "Clip atlas beta values to [eps, 1-eps], matching deconvolve_em.r [0.01]")
p <- add_argument(p, "--marker_file",        type = "character", default = NA, help = "Optional: marker_selection_gss.r's --output_file. Restricts the pairwise distance computation (and, if requested, the lineage map) to the union of GSS-selected marker CpGs across all cell types, per the unified-EM-framework spec's repurposing of GSS as a marker-restricted distance metric for lineage clustering rather than a likelihood filter. Omit to use every complete-row atlas CpG, as before.")
p <- add_argument(p, "--lineage_map_out",    type = "character", default = NA, help = "Optional output TSV (cell_type, lineage_id, lineage_label): cuts the confusability dendrogram already computed below into lineage groups. Requires exactly one of --n_lineages/--lineage_cut_height.")
p <- add_argument(p, "--n_lineages",         type = "integer", default = NA, help = "Cut the confusability dendrogram into exactly this many lineage groups. Mutually exclusive with --lineage_cut_height.")
p <- add_argument(p, "--lineage_cut_height", type = "double",  default = NA, help = "Cut the confusability dendrogram at this distance (median worst-case per-CpG log-lik difference) instead of a fixed cluster count. Mutually exclusive with --n_lineages.")
p <- add_argument(p, "--lineage_dendrogram_out", type = "character", default = NA, help = "Optional output PNG: the confusability dendrogram with the lineage cut marked [<lineage_map_out> with _dendrogram.png suffix, if --lineage_map_out given]")
p <- add_argument(p, "--pair_table_out",     type = "character", default = "atlas_pair_separability.csv", help = "Output CSV: full pairwise separability table")
p <- add_argument(p, "--nn_table_out",       type = "character", default = "atlas_nearest_neighbor.csv",  help = "Output CSV: per-cell-type nearest-neighbor summary")
p <- add_argument(p, "--density_pairs",      type = "character", default = NA, help = "Comma-separated cellTypeA:cellTypeB pairs to plot [default: the 2 most-confusable pairs found]")
p <- add_argument(p, "--n_density_pairs",    type = "integer", default = 2, help = "How many most-confusable pairs to auto-plot if --density_pairs is not given [2]")
p <- add_argument(p, "--density_output_file", type = "character", default = "atlas_density_pairs.png", help = "Output PNG for the 2D density plot(s)")
p <- add_argument(p, "--n_bins",             type = "integer", default = 60, help = "Bins per axis for the 2D density plot [60]")
p <- add_argument(p, "--confusability_output_file", type = "character", default = "atlas_confusability_map.png", help = "Output PNG for the full pairwise confusability heatmap")
p <- add_argument(p, "--correlation_output_file", type = "character", default = "atlas_correlation_map.png", help = "Output PNG for the full pairwise beta-correlation heatmap")

#/dcs11/scharpf/data/horjuela/horjuela/Methylation_Atlas_Dev/publications/Orjuela2026/data/methid_references/brain_liver_pancreas_blood.tsv

args <- parse_args(p)

if (!is.na(args$n_lineages) && !is.na(args$lineage_cut_height)) {
  stop("Pass at most one of --n_lineages / --lineage_cut_height, not both.")
}
if (!is.na(args$lineage_map_out) && is.na(args$n_lineages) && is.na(args$lineage_cut_height)) {
  stop("--lineage_map_out requires exactly one of --n_lineages / --lineage_cut_height.")
}
if ((!is.na(args$n_lineages) || !is.na(args$lineage_cut_height)) && is.na(args$lineage_map_out)) {
  stop("--n_lineages/--lineage_cut_height given without --lineage_map_out -- nowhere to write the result.")
}

atlas <- fread(args$map_file, data.table = TRUE)
cell_types <- setdiff(colnames(atlas), c("chr", "start", "end"))
setnames(atlas, "chr", "chrom")

complete_idx <- rowSums(is.na(atlas[, ..cell_types])) == 0
atlas <- atlas[complete_idx]
cat(sprintf("Atlas: %d cell types, %d complete-row CpGs\n", length(cell_types), nrow(atlas)))

if (!is.na(args$marker_file)) {
  markers <- fread(args$marker_file, data.table = TRUE)
  marker_pos <- unique(markers[, .(chrom, start)])
  before_n <- nrow(atlas)
  atlas <- merge(atlas, marker_pos, by = c("chrom", "start"))
  cat(sprintf("--marker_file given: restricted to the union of %d marker CpG(s) across %d cell type(s) -- %d/%d atlas CpGs retained\n",
              nrow(marker_pos), uniqueN(markers$cell_type), nrow(atlas), before_n))
  if (nrow(atlas) == 0) stop("No atlas CpGs remain after restricting to --marker_file positions -- check --marker_file matches this --map_file's coordinates.")
}

eps <- args$beta_eps
beta_mat <- as.matrix(atlas[, ..cell_types])
beta_mat <- pmin(pmax(beta_mat, eps), 1 - eps)
log_b   <- log(beta_mat)
log_1mb <- log(1 - beta_mat)

# -----------------------------------------------
# Per-CpG worst-case log-likelihood difference, aggregated per pair
# -----------------------------------------------

K <- length(cell_types)
pair_idx <- which(upper.tri(matrix(NA, K, K)), arr.ind = TRUE)

pair_list <- vector("list", nrow(pair_idx))
for (r in seq_len(nrow(pair_idx))) {
  i <- pair_idx[r, 1]; j <- pair_idx[r, 2]
  diff_if_1 <- abs(log_b[, i]   - log_b[, j])
  diff_if_0 <- abs(log_1mb[, i] - log_1mb[, j])
  min_diff  <- pmin(diff_if_1, diff_if_0)
  pair_list[[r]] <- data.table(
    cell_type_a = cell_types[i],
    cell_type_b = cell_types[j],
    median_min_loglik_diff = median(min_diff),
    mean_min_loglik_diff   = mean(min_diff),
    pearson_r = cor(beta_mat[, i], beta_mat[, j])
  )
}
pair_df <- rbindlist(pair_list)
setorder(pair_df, median_min_loglik_diff)

fwrite(pair_df, args$pair_table_out)
cat(sprintf("\nWrote %s (%d pairs)\n", args$pair_table_out, nrow(pair_df)))
cat("\nMost confusable pairs (lowest median worst-case per-CpG log-lik difference):\n")
print(head(pair_df, 10))

# -----------------------------------------------
# Aggregate per cell type: distance to nearest (most confusable) neighbor
# -----------------------------------------------

nn_list <- vector("list", K)
for (k in seq_len(K)) {
  ct <- cell_types[k]
  rows <- pair_df[cell_type_a == ct | cell_type_b == ct]
  best <- rows[which.min(median_min_loglik_diff)]
  neighbor <- if (best$cell_type_a == ct) best$cell_type_b else best$cell_type_a
  nn_list[[k]] <- data.table(
    cell_type = ct,
    nearest_neighbor = neighbor,
    separability_score = best$median_min_loglik_diff,
    pearson_r_with_neighbor = best$pearson_r
  )
}
nn_df <- rbindlist(nn_list)
setorder(nn_df, separability_score)

fwrite(nn_df, args$nn_table_out)
cat(sprintf("\nWrote %s (%d cell types)\n", args$nn_table_out, nrow(nn_df)))
cat("\nPer-cell-type nearest-neighbor separability (lowest = hardest to resolve):\n")
print(nn_df)

# -----------------------------------------------
# Confusability map: full pairwise heatmap, cell types ordered by
# hierarchical clustering on the score so confusable clusters group together
# -----------------------------------------------

full_mat <- matrix(NA_real_, K, K, dimnames = list(cell_types, cell_types))
for (r in seq_len(nrow(pair_df))) {
  a <- pair_df$cell_type_a[r]; b <- pair_df$cell_type_b[r]
  full_mat[a, b] <- pair_df$median_min_loglik_diff[r]
  full_mat[b, a] <- pair_df$median_min_loglik_diff[r]
}

dist_for_clust <- full_mat
diag(dist_for_clust) <- 0
hc <- hclust(as.dist(dist_for_clust), method = "average")
ct_order <- cell_types[hc$order]

# -----------------------------------------------
# Lineage map (layer D prerequisite, unified-EM-framework spec Build Item
# 1a): cut this SAME confusability dendrogram -- not the separate
# correlation-based one below -- into lineage groups, since it directly
# reflects how confusable cell types are under the EM's own per-CpG
# Bernoulli scoring (the quantity layer D's backoff rule actually cares
# about), not just raw beta correlation.
# -----------------------------------------------

if (!is.na(args$lineage_map_out)) {
  lineage_id <- if (!is.na(args$n_lineages)) {
    cutree(hc, k = args$n_lineages)
  } else {
    cutree(hc, h = args$lineage_cut_height)
  }
  lineage_dt <- data.table(
    cell_type = names(lineage_id),
    lineage_id = as.integer(lineage_id)
  )
  lineage_dt[, lineage_label := paste0("lineage_", lineage_id)]
  # A cluster label is only meaningful if it groups >1 cell type; singleton
  # "lineages" are just that cell type on its own -- flagged, not hidden,
  # since a singleton lineage means layer D's same-lineage backoff branch
  # can never trigger for that cell type (there's no sibling to back off to).
  lineage_dt[, lineage_size := .N, by = lineage_id]
  setorder(lineage_dt, lineage_id, cell_type)
  fwrite(lineage_dt, args$lineage_map_out, sep = "\t")
  n_singleton <- sum(lineage_dt$lineage_size == 1)
  cat(sprintf("\nWrote %s: %d cell type(s) into %d lineage group(s) (%d singleton)\n",
              args$lineage_map_out, nrow(lineage_dt), uniqueN(lineage_dt$lineage_id), n_singleton))
  print(lineage_dt[, .(cell_type, lineage_id, lineage_size)])

  dendro_file <- if (!is.na(args$lineage_dendrogram_out)) args$lineage_dendrogram_out else sub("(\\.[^.]+)?$", "_dendrogram.png", args$lineage_map_out, perl = TRUE)
  cut_height_used <- if (!is.na(args$lineage_cut_height)) {
    args$lineage_cut_height
  } else {
    # height that first yields exactly --n_lineages clusters, for the plot's cut line
    sorted_heights <- sort(hc$height)
    k_at_height <- (length(hc$height) + 1) - seq_along(sorted_heights)
    idx <- which(k_at_height == args$n_lineages)
    if (length(idx) > 0) sorted_heights[max(idx)] else NA_real_
  }

  png(dendro_file, width = max(6, 0.35 * K + 3) * 300, height = 5 * 300, res = 300)
  plot(hc, main = "Cell-type lineage clustering (confusability distance)",
       xlab = "", ylab = "Median worst-case per-CpG log-lik difference", sub = "")
  if (!is.na(cut_height_used)) abline(h = cut_height_used, col = "red", lty = 2)
  dev.off()
  cat(sprintf("Wrote %s\n", dendro_file))
}

full_dt <- as.data.table(full_mat, keep.rownames = "cell_type_1")
conf_df <- melt(full_dt, id.vars = "cell_type_1", variable.name = "cell_type_2", value.name = "score")
conf_df[, cell_type_2 := as.character(cell_type_2)]
conf_df[, cell_type_1 := factor(cell_type_1, levels = ct_order)]
conf_df[, cell_type_2 := factor(cell_type_2, levels = rev(ct_order))]

show_labels <- K <= 20
label_layer <- if (show_labels) {
  geom_text(aes(label = ifelse(is.na(score), "", sprintf("%.3f", score))), size = 2.6, color = "white")
} else {
  NULL
}

p_conf <- ggplot(conf_df, aes(x = cell_type_1, y = cell_type_2, fill = score)) +
  geom_tile(color = "white", linewidth = 0.5) +
  label_layer +
  scale_fill_viridis_c(name = "Median worst-case\nlog-lik diff / CpG", trans = "log10", direction = -1, na.value = "grey85") +
  coord_fixed() +
  labs(
    x = NULL, y = NULL,
    title = "Cell-type confusability map",
    subtitle = "Brighter = harder to distinguish with a single CpG (lower worst-case log-lik gap); clustered by similarity"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(args$confusability_output_file, p_conf,
       width = max(6, 0.5 * K + 3), height = max(5, 0.5 * K + 2),
       dpi = 300, bg = "white", limitsize = FALSE)
cat(sprintf("\nWrote %s\n", args$confusability_output_file))

# -----------------------------------------------
# Correlation map: same idea as the confusability map, but the plain
# Pearson correlation (already computed per pair above) rather than the
# log-lik-gap score -- easier to explain to a non-technical audience, and
# the standard input for a correlation-based cell-type hierarchy (doc 18:
# AEGIS built their 34->23->20 tissue hierarchy by collapsing on empirical
# correlation structure; this is the first artifact scoping the same
# approach for smMethID's own EM would need). Clustered on ITS OWN
# structure (1 - r as distance), not reusing the confusability map's
# ordering -- a correlation-based hierarchy should reflect correlation
# structure, even if in practice the two orderings often end up similar.
# -----------------------------------------------

full_mat_r <- matrix(NA_real_, K, K, dimnames = list(cell_types, cell_types))
diag(full_mat_r) <- 1
for (r in seq_len(nrow(pair_df))) {
  a <- pair_df$cell_type_a[r]; b <- pair_df$cell_type_b[r]
  full_mat_r[a, b] <- pair_df$pearson_r[r]
  full_mat_r[b, a] <- pair_df$pearson_r[r]
}

dist_for_clust_r <- 1 - full_mat_r
diag(dist_for_clust_r) <- 0
hc_r <- hclust(as.dist(dist_for_clust_r), method = "average")
ct_order_r <- cell_types[hc_r$order]

corr_dt <- as.data.table(full_mat_r, keep.rownames = "cell_type_1")
corr_df <- melt(corr_dt, id.vars = "cell_type_1", variable.name = "cell_type_2", value.name = "r")
corr_df[, cell_type_2 := as.character(cell_type_2)]
corr_df[, cell_type_1 := factor(cell_type_1, levels = ct_order_r)]
corr_df[, cell_type_2 := factor(cell_type_2, levels = rev(ct_order_r))]

corr_label_layer <- if (show_labels) {
  geom_text(aes(label = sprintf("%.2f", r)), size = 2.6, color = "black")
} else {
  NULL
}

p_corr <- ggplot(corr_df, aes(x = cell_type_1, y = cell_type_2, fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  corr_label_layer +
  scale_fill_gradient2(name = "Pearson r", low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-1, 1)) +
  coord_fixed() +
  labs(
    x = NULL, y = NULL,
    title = "Cell-type beta correlation map",
    subtitle = "Pearson correlation of per-CpG beta across cell types; clustered on correlation structure (1-r)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(args$correlation_output_file, p_corr,
       width = max(6, 0.5 * K + 3), height = max(5, 0.5 * K + 2),
       dpi = 300, bg = "white", limitsize = FALSE)
cat(sprintf("\nWrote %s\n", args$correlation_output_file))

# -----------------------------------------------
# 2D density plots for cell types of interest
# -----------------------------------------------

if (!is.na(args$density_pairs)) {
  pair_specs <- strsplit(args$density_pairs, ",")[[1]]
  plot_pairs <- lapply(pair_specs, function(s) strsplit(trimws(s), ":")[[1]])
} else {
  n_plot <- min(args$n_density_pairs, nrow(pair_df))
  plot_pairs <- lapply(seq_len(n_plot), function(i) c(pair_df$cell_type_a[i], pair_df$cell_type_b[i]))
}

make_density_plot <- function(ct_a, ct_b) {
  stopifnot(ct_a %in% cell_types, ct_b %in% cell_types)
  df <- data.table(a = beta_mat[, ct_a], b = beta_mat[, ct_b])
  r <- cor(df$a, df$b)
  score <- pair_df[(cell_type_a == ct_a & cell_type_b == ct_b) | (cell_type_a == ct_b & cell_type_b == ct_a), median_min_loglik_diff]

  ggplot(df, aes(x = a, y = b)) +
    geom_bin2d(bins = args$n_bins) +
    scale_fill_viridis_c(name = "CpGs", trans = "log10", labels = scales::label_comma()) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 0.4) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      x = ct_a, y = ct_b,
      title = sprintf("r = %.3f, median worst-case log-lik diff/CpG = %.3f", r, score)
    ) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank())
}

plots <- lapply(plot_pairs, function(pr) make_density_plot(pr[1], pr[2]))
combined <- wrap_plots(plots, ncol = min(2, length(plots)))
ggsave(args$density_output_file, combined,
       width = 6 * min(2, length(plots)), height = 6 * ceiling(length(plots) / 2),
       dpi = 300, bg = "white", limitsize = FALSE)
cat(sprintf("\nWrote %s (%d pair%s)\n", args$density_output_file, length(plots), if (length(plots) == 1) "" else "s"))

quit(save = "no")
