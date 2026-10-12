#!/usr/bin/env Rscript

# Pairwise cell-type separability in a methylation atlas, computed on every
# complete-row atlas CpG (the same set deconvolve_em.r scores molecules
# against as a per-CpG independent-Bernoulli log-likelihood).
#
# For each pair of cell types, the separability score is the expected
# evidence one CpG gives for the true cell type over the other: the mean
# per-CpG Bernoulli KL divergence, taking the weaker direction. A molecule
# with n CpGs accumulates about n times this many nats.
#
# Also writes a marker balance table (atlas_marker_balance.tsv) for choosing
# deconvolve_em.r --markers_per_cell_type, and, for each N in
# --marker_sets, the same separability recomputed on only the CpGs that
# deconvolve_em.r --markers_per_cell_type N (gap ranking) would keep: extra
# columns in the pair table and one heatmap per N.

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(argparser)))
suppressWarnings(suppressMessages(library(ggplot2)))

p <- arg_parser("Pairwise cell-type separability and marker balance for a methylation atlas")
p <- add_argument(p, "--map_file",        type = "character", help = "Atlas map file (TSV/TSV.gz): chr, start, end, <cell_type_beta...>")
p <- add_argument(p, "--beta_eps",        type = "double",  default = 0.01, help = "Clip atlas beta values to [eps, 1-eps], matching deconvolve_em.r [0.01]")
p <- add_argument(p, "--marker_balance_out", type = "character", default = "atlas_marker_balance.tsv", help = "Output TSV: per cell type, how many CpGs set it apart from its closest other cell type, and how strong its Nth-best marker is")
p <- add_argument(p, "--pair_table_out",  type = "character", default = "atlas_pair_separability.csv", help = "Output CSV: pairwise separability table")
p <- add_argument(p, "--marker_sets", type = "character", default = "5000,20000", help = "Comma-separated markers-per-cell-type values to also compute separability for (gap ranking, as deconvolve_em.r); '' to skip [5000,20000]")
p <- add_argument(p, "--confusability_output_file", type = "character", default = "atlas_confusability_map.png", help = "Output PNG: pairwise separability heatmap")

args <- parse_args(p)

atlas <- fread(args$map_file, data.table = TRUE)
cell_types <- setdiff(colnames(atlas), c("chr", "start", "end"))
atlas <- atlas[rowSums(is.na(atlas[, ..cell_types])) == 0]
cat(sprintf("Atlas: %d cell types, %d complete-row CpGs\n", length(cell_types), nrow(atlas)))

eps <- args$beta_eps
beta_mat <- pmin(pmax(as.matrix(atlas[, ..cell_types]), eps), 1 - eps)
log_b   <- log(beta_mat)
log_1mb <- log(1 - beta_mat)

# -----------------------------------------------
# Per-pair expected evidence per CpG (min directional KL, nats)
# -----------------------------------------------

K <- length(cell_types)
pair_idx <- which(upper.tri(matrix(0, K, K)), arr.ind = TRUE)
separability <- function(rows) {
  sm <- matrix(0, K, K, dimnames = list(cell_types, cell_types))
  b <- beta_mat[rows, , drop = FALSE]; lb <- log_b[rows, , drop = FALSE]; l1b <- log_1mb[rows, , drop = FALSE]
  for (r in seq_len(nrow(pair_idx))) {
    i <- pair_idx[r, 1]; j <- pair_idx[r, 2]
    kl_ij <- mean(b[, i] * (lb[, i] - lb[, j]) + (1 - b[, i]) * (l1b[, i] - l1b[, j]))
    kl_ji <- mean(b[, j] * (lb[, j] - lb[, i]) + (1 - b[, j]) * (l1b[, j] - l1b[, i]))
    sm[i, j] <- sm[j, i] <- min(kl_ij, kl_ji)
  }
  sm
}
score_mat <- separability(seq_len(nrow(beta_mat)))

pair_df <- data.table(
  cell_type_a = cell_types[pair_idx[, 1]],
  cell_type_b = cell_types[pair_idx[, 2]],
  evidence_per_cpg = score_mat[pair_idx]
)

# Average-linkage clustering, only to order the heatmap so similar cell
# types sit together
hc <- hclust(as.dist(score_mat), method = "average")

# -----------------------------------------------
# Marker balance: per cell type, the gap to its closest other cell type at
# every CpG. Shows whether some cell types have far more distinctive CpGs
# than others (an all-CpG atlas then leans toward them), and how strong the
# Nth-best marker is for deconvolve_em.r --markers_per_cell_type N.
# -----------------------------------------------

marker_sets <- as.integer(trimws(strsplit(args$marker_sets, ",")[[1]]))
marker_sets <- marker_sets[!is.na(marker_sets) & marker_sets > 0]
max_n <- if (length(marker_sets)) max(marker_sets) else 0L
top_idx <- vector("list", K)   # per cell type, CpG rows ranked by gap (as deconvolve_em.r)

balance_rows <- lapply(seq_len(K), function(k) {
  gap <- abs(beta_mat[, k] - beta_mat[, -k, drop = FALSE])
  gap <- do.call(pmin, lapply(seq_len(ncol(gap)), function(j) gap[, j]))
  ord <- order(gap, decreasing = TRUE, method = "radix")
  if (max_n > 0) top_idx[[k]] <<- ord[seq_len(min(max_n, length(ord)))]
  sorted <- gap[ord]
  at <- function(n) if (length(sorted) >= n) round(sorted[n], 3) else NA_real_
  data.table(cell_type = cell_types[k],
             n_gap_ge_0.2 = sum(gap >= 0.2), n_gap_ge_0.3 = sum(gap >= 0.3), n_gap_ge_0.5 = sum(gap >= 0.5),
             gap_at_5000 = at(5000), gap_at_20000 = at(20000), gap_at_50000 = at(50000))
})
balance_dt <- rbindlist(balance_rows)
fwrite(balance_dt, args$marker_balance_out, sep = "\t")
cat(sprintf("\nWrote %s (CpGs that set each cell type apart from its closest other):\n", args$marker_balance_out))
print(balance_dt)

# -----------------------------------------------
# Separability on balanced marker sets: the CpGs deconvolve_em.r
# --markers_per_cell_type N would keep (each cell type's top N by gap, pooled)
# -----------------------------------------------

marker_scores <- list()
for (n in marker_sets) {
  rows <- sort(unique(unlist(lapply(top_idx, function(ix) ix[seq_len(min(n, length(ix)))]))))
  marker_scores[[as.character(n)]] <- separability(rows)
  pair_df[, (paste0("evidence_per_cpg_markers_", n)) := marker_scores[[as.character(n)]][pair_idx]]
  cat(sprintf("Markers per cell type %d: %d CpGs\n", n, length(rows)))
}

setorder(pair_df, evidence_per_cpg)
fwrite(pair_df, args$pair_table_out)
cat(sprintf("\nWrote %s (%d pairs). Hardest pairs to tell apart (all CpGs):\n", args$pair_table_out, nrow(pair_df)))
print(head(pair_df, 10))

# -----------------------------------------------
# Heatmaps, cell types in the same order (clustering on all CpGs) so the
# all-CpG and marker versions can be compared cell by cell
# -----------------------------------------------

ct_order <- cell_types[hc$order]
plot_heatmap <- function(sm, out_file, subtitle) {
  heat_df <- as.data.table(as.table(sm))
  setnames(heat_df, c("cell_type_1", "cell_type_2", "score"))
  heat_df[cell_type_1 == cell_type_2, score := NA_real_]
  heat_df[, cell_type_1 := factor(cell_type_1, levels = ct_order)]
  heat_df[, cell_type_2 := factor(cell_type_2, levels = rev(ct_order))]
  p_conf <- ggplot(heat_df, aes(x = cell_type_1, y = cell_type_2, fill = score)) +
    geom_tile(color = "white", linewidth = 0.5) +
    { if (K <= 20) geom_text(aes(label = ifelse(is.na(score), "", sprintf("%.2f", score))), size = 2.6, color = "white") } +
    scale_fill_viridis_c(name = "Evidence\nper CpG (nats)", trans = "log10", direction = -1, na.value = "grey85") +
    coord_fixed() +
    labs(x = NULL, y = NULL, title = "Cell-type separability", subtitle = subtitle) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  ggsave(out_file, p_conf, width = max(6, 0.5 * K + 3), height = max(5, 0.5 * K + 2),
         dpi = 300, bg = "white", limitsize = FALSE)
  cat(sprintf("Wrote %s\n", out_file))
}

plot_heatmap(score_mat, args$confusability_output_file,
             "All atlas CpGs. Brighter = harder to tell apart; similar cell types grouped")
for (n in names(marker_scores)) {
  plot_heatmap(marker_scores[[n]],
               sub("(\\.[^.]+)?$", paste0("_markers", n, "\\1"), args$confusability_output_file, perl = TRUE),
               sprintf("Balanced markers, %s per cell type. Brighter = harder to tell apart", n))
}
