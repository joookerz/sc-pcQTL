#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))[1L]])))
source(file.path(script_dir, "scpcqtl_utils.R"))
args <- parse_cli()
stage <- required_arg(args, "stage")
clusters_file <- required_arg(args, "clusters")
counts_file <- required_arg(args, "counts")
outdir <- required_arg(args, "outdir")
celltype <- required_arg(args, "celltype")
variance_threshold <- as.numeric(args$pca_variance %||% 0.95)
cis_window <- as.integer(args$cis_window %||% 500000L)
covariate_names <- split_csv_arg(args$covariates %||% "age,sex,pc1,pc2,pc3,pc4,pc5,pc6,pf1,pf2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

clusters <- fread(clusters_file)
covars <- readRDS(file.path(stage, "covariates.rds"))
missing_covars <- setdiff(covariate_names, names(covars))
if (length(missing_covars)) stop("Missing PCA/QTL covariates: ", paste(missing_covars, collapse = ", "))

all_cluster_genes <- unique(unlist(strsplit(clusters$genes, ",", fixed = TRUE)))
header <- names(fread(counts_file, nrows = 0L))
available_cluster_genes <- intersect(all_cluster_genes, header)
expression <- if (length(available_cluster_genes)) {
  fread(counts_file, select = available_cluster_genes)
} else {
  data.table()
}
if (nrow(expression) && nrow(expression) != nrow(covars)) {
  stop("Count and covariate row counts differ during cluster PCA")
}

summary_rows <- list()
task_rows <- list()
k <- 0L
for (i in seq_len(nrow(clusters))) {
  cluster <- clusters[i]
  id <- cluster$cluster_id
  chr <- as.integer(cluster$chromosome)
  genes <- strsplit(cluster$genes, ",", fixed = TRUE)[[1L]]
  available <- genes[genes %in% names(expression)]
  cluster_dir <- file.path(outdir, id)
  dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)
  if (length(available) < 2L) {
    summary_rows[[i]] <- data.table(cluster_id = id, status = "failed", reason = "fewer_than_two_genes")
    next
  }
  matrix <- as.matrix(expression[, ..available])
  keep <- complete.cases(matrix)
  matrix <- matrix[keep, , drop = FALSE]
  if (nrow(matrix) < ncol(matrix)) {
    summary_rows[[i]] <- data.table(cluster_id = id, status = "failed", reason = "cells_fewer_than_genes")
    next
  }
  pca <- tryCatch(prcomp(matrix, center = TRUE, scale. = FALSE), error = function(e) e)
  if (inherits(pca, "error")) {
    summary_rows[[i]] <- data.table(cluster_id = id, status = "failed", reason = conditionMessage(pca))
    next
  }
  variance <- pca$sdev^2 / sum(pca$sdev^2)
  cumulative <- cumsum(variance)
  n_retained <- which(cumulative >= variance_threshold)[1L]
  if (!is.finite(n_retained)) n_retained <- length(variance)
  saveRDS(pca, file.path(cluster_dir, "pca.rds"), compress = "gzip", version = 3)
  fwrite(data.table(PC = paste0("PC", seq_along(variance)), eigenvalue = pca$sdev^2,
                    variance_explained = variance, cumulative_variance = cumulative,
                    retained = seq_along(variance) <= n_retained),
         file.path(cluster_dir, "variance_explained.tsv"), sep = "\t")
  fwrite(data.table(gene_name = rownames(pca$rotation), as.data.table(pca$rotation)),
         file.path(cluster_dir, "gene_loadings.tsv"), sep = "\t")
  pheno <- cbind(
    covars[keep, c("individual", covariate_names), with = FALSE],
    as.data.table(pca$x[, seq_len(n_retained), drop = FALSE])
  )
  pheno_file <- file.path(cluster_dir, "phenotypes.tsv")
  fwrite(pheno, pheno_file, sep = "\t")
  region_file <- file.path(cluster_dir, "cis_region.tsv")
  writeLines(paste(chr, max(0L, as.integer(cluster$start_position) - cis_window),
                   as.integer(cluster$end_position) + cis_window, sep = "\t"), region_file)
  summary_rows[[i]] <- data.table(
    cluster_id = id, status = "ok", reason = "", n_genes = length(available),
    n_cells = nrow(matrix), n_donors = uniqueN(pheno$individual),
    n_pcs_total = ncol(pca$x), n_pcs_retained = n_retained,
    retained_variance = cumulative[n_retained]
  )
  for (pc in seq_len(n_retained)) {
    k <- k + 1L
    task_rows[[k]] <- data.table(
      task_id = k, celltype = celltype, cluster_id = id,
      phenotype_id = paste0("PC", pc), chromosome = chr,
      phenotype_file = file.path(id, "phenotypes.tsv"),
      region_file = file.path(id, "cis_region.tsv")
    )
  }
}
summary <- if (length(summary_rows)) rbindlist(summary_rows, fill = TRUE) else data.table(
  cluster_id = character(), status = character(), reason = character())
tasks <- if (length(task_rows)) rbindlist(task_rows) else data.table(
  task_id = integer(), celltype = character(), cluster_id = character(),
  phenotype_id = character(), chromosome = integer(), phenotype_file = character(), region_file = character())
fwrite(summary, file.path(outdir, "pca_summary.tsv"), sep = "\t")
fwrite(tasks, file.path(outdir, "qtl_tasks.tsv"), sep = "\t")
