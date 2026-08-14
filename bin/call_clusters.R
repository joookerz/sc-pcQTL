#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))[1L]])))
source(file.path(script_dir, "scpcqtl_utils.R"))
args <- parse_cli()
stage <- required_arg(args, "stage")
pairs_file <- required_arg(args, "pairs")
outdir <- required_arg(args, "outdir")
celltype <- required_arg(args, "celltype")
max_window <- as.integer(args$max_cluster_genes %||% 50L)
min_window <- as.integer(args$min_cluster_genes %||% 2L)
density_threshold <- as.numeric(args$cluster_density %||% 0.70)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

genes <- fread(file.path(stage, "gene_filtering.tsv"))[chromosome %in% 1:22]
pairs <- fread(pairs_file)
clusters <- list()
assignments <- list()
k <- 0L
for (chr in 1:22) {
  order_columns <- intersect(c("start", "annotation_order"), names(genes))
  chr_info <- copy(genes[chromosome == chr])
  setorderv(chr_info, order_columns)
  gene_names <- chr_info$gene_name
  n <- length(gene_names)
  if (n < min_window) next
  matrix <- matrix(FALSE, n, n, dimnames = list(gene_names, gene_names))
  diag(matrix) <- TRUE
  chr_pairs <- pairs[chromosome == chr & significant == TRUE]
  if (nrow(chr_pairs)) {
    left <- match(chr_pairs$Gene1, gene_names)
    right <- match(chr_pairs$Gene2, gene_names)
    valid <- !is.na(left) & !is.na(right)
    matrix[cbind(left[valid], right[valid])] <- TRUE
    matrix[cbind(right[valid], left[valid])] <- TRUE
  }
  assigned <- setNames(rep(FALSE, n), gene_names)
  chr_index <- 0L
  for (size in seq.int(min(max_window, n), min_window)) {
    for (first in seq_len(n - size + 1L)) {
      window <- gene_names[first:(first + size - 1L)]
      if (any(assigned[window])) next
      density <- mean(matrix[window, window, drop = FALSE][upper.tri(matrix[window, window, drop = FALSE])])
      if (is.finite(density) && density >= density_threshold) {
        k <- k + 1L; chr_index <- chr_index + 1L
        id <- sprintf("SC_chr%d_cluster_%03d", chr, chr_index)
        info <- chr_info[gene_name %in% window]
        clusters[[k]] <- data.table(
          celltype = celltype, cluster_id = id, chromosome = chr,
          cluster_size = size, start_position = min(info$start),
          end_position = max(info$end), cluster_span_bp = max(info$end) - min(info$start),
          edge_density = density, genes = paste(window, collapse = ",")
        )
        assignments[[k]] <- data.table(celltype = celltype, cluster_id = id,
                                       chromosome = chr, gene_name = window)
        assigned[window] <- TRUE
      }
    }
  }
}
cluster_dt <- if (length(clusters)) rbindlist(clusters) else data.table(
  celltype = character(), cluster_id = character(), chromosome = integer(),
  cluster_size = integer(), start_position = integer(), end_position = integer(),
  cluster_span_bp = integer(), edge_density = numeric(), genes = character())
assignment_dt <- if (length(assignments)) rbindlist(assignments) else data.table(
  celltype = character(), cluster_id = character(), chromosome = integer(), gene_name = character())
fwrite(cluster_dt, file.path(outdir, "clusters.tsv"), sep = "\t")
fwrite(assignment_dt, file.path(outdir, "cluster_genes.tsv"), sep = "\t")
fwrite(data.table(celltype = celltype, n_clusters = nrow(cluster_dt),
                  n_cluster_genes = uniqueN(assignment_dt$gene_name)),
       file.path(outdir, "cluster_summary.tsv"), sep = "\t")
