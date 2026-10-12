#!/usr/bin/env Rscript

# Subset an atlas (chr, start, end, <cell_type betas...>) to the CpGs chosen by
# marker selection, giving a MethID map that main.nf can use directly as --methid_map.

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(argparser)))

p <- arg_parser("Subset an atlas to marker CpGs to make a MethID map")
p <- add_argument(p, "--atlas_file", type = "character", help = "Atlas TSV: chr, start, end, one beta column per cell type")
p <- add_argument(p, "--marker_file", type = "character", help = "Marker TSV with chrom and start columns (e.g. marker_selection_gss.r output)")
p <- add_argument(p, "--output_file", type = "character", default = "methid_map_gss.tsv", help = "Output MethID map")

args <- parse_args(p)

atlas <- fread(args$atlas_file)
markers <- fread(args$marker_file)

# A CpG can be a marker for several cell types; keep each position once
markers <- unique(markers[, .(chr = chrom, start)])

methid_map <- atlas[markers, on = .(chr, start), nomatch = 0]
setorder(methid_map, chr, start)

message("Marker CpGs: ", nrow(markers), "; found in atlas: ", nrow(methid_map))
if (nrow(methid_map) == 0) stop("No marker CpGs matched the atlas -- check that both use the same coordinates")

fwrite(methid_map, args$output_file, sep = "\t", quote = FALSE)

quit(save = "no")
