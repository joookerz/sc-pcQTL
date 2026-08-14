#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))[1L]])))
source(file.path(script_dir, "scpcqtl_utils.R"))
args <- parse_cli()
input_list <- required_arg(args, "input_list")
outdir <- required_arg(args, "outdir")
fdr <- as.numeric(args$qtl_fdr %||% 0.05)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dirs <- fread(input_list, header = FALSE)[[1L]]
dirs <- sort(dirs[dir.exists(dirs)])

variant_output <- file.path(outdir, "all_variant_results.tsv.gz")
phenotype_rows <- list()
acat_rows <- list()
variant_written <- FALSE

for (directory in dirs) {
  metadata_file <- file.path(directory, "metadata.tsv")
  association_file <- file.path(directory, "association.tsv")
  if (!file.exists(metadata_file) || !file.exists(association_file)) next
  metadata <- fread(metadata_file)
  if (nrow(metadata) != 1L) stop("Expected one metadata row in ", metadata_file)
  association <- fread(association_file)
  p_col <- intersect(c("p.value", "pvalue", "Pvalue", "P"), names(association))[1L]
  if (is.na(p_col)) stop("Cannot identify SAIGE-QTL p-value column in ", association_file)
  association[, `:=`(
    celltype = metadata$celltype[[1L]], cluster_id = metadata$cluster_id[[1L]],
    phenotype_id = metadata$phenotype_id[[1L]], chromosome = metadata$chromosome[[1L]],
    qtl_type = "pcQTL", pvalue = as.numeric(get(p_col))
  )]
  association[, qvalue := p.adjust(pvalue, method = "BH")]
  association[, significant := is.finite(qvalue) & qvalue < fdr]

  fwrite(
    association, variant_output, sep = "\t", compress = "gzip",
    append = variant_written, col.names = !variant_written
  )
  variant_written <- TRUE

  n_value <- if ("N" %in% names(association)) suppressWarnings(as.numeric(association$N)) else numeric()
  phenotype_rows[[length(phenotype_rows) + 1L]] <- data.table(
    celltype = metadata$celltype[[1L]], cluster_id = metadata$cluster_id[[1L]],
    phenotype_id = metadata$phenotype_id[[1L]], n_variants = nrow(association),
    n_significant_variants = sum(association$significant),
    min_pvalue = if (any(is.finite(association$pvalue))) min(association$pvalue[is.finite(association$pvalue)]) else NA_real_,
    min_qvalue = if (any(is.finite(association$qvalue))) min(association$qvalue[is.finite(association$qvalue)]) else NA_real_,
    significant = any(association$significant),
    donor_sample_size = if (any(is.finite(n_value))) max(n_value[is.finite(n_value)]) else NA_real_
  )

  acat_file <- file.path(directory, "acat.tsv")
  if (file.exists(acat_file)) {
    acat <- fread(acat_file)
    acat[, `:=`(
      celltype = metadata$celltype[[1L]], cluster_id = metadata$cluster_id[[1L]],
      phenotype_id = metadata$phenotype_id[[1L]]
    )]
    acat_rows[[length(acat_rows) + 1L]] <- acat
  }
}

if (!variant_written) {
  fwrite(data.table(
    celltype = character(), cluster_id = character(), phenotype_id = character(),
    chromosome = integer(), qtl_type = character(), pvalue = numeric(),
    qvalue = numeric(), significant = logical()
  ), variant_output, sep = "\t", compress = "gzip")
}

phenotype <- if (length(phenotype_rows)) rbindlist(phenotype_rows, fill = TRUE) else data.table(
  celltype = character(), cluster_id = character(), phenotype_id = character(),
  n_variants = integer(), n_significant_variants = integer(), min_pvalue = numeric(),
  min_qvalue = numeric(), significant = logical(), donor_sample_size = numeric()
)
acat <- if (length(acat_rows)) rbindlist(acat_rows, fill = TRUE) else data.table()
if (nrow(phenotype)) setorder(phenotype, celltype, cluster_id, phenotype_id)
if (nrow(acat)) setorder(acat, celltype, cluster_id, phenotype_id)
if (nrow(acat)) {
  acat_p_col <- intersect(c("ACAT_p", "acat_p", "pvalue", "Pvalue", "P"), names(acat))[1L]
  if (!is.na(acat_p_col)) {
    acat[, acat_pvalue := as.numeric(get(acat_p_col))]
    acat[, acat_qvalue := p.adjust(acat_pvalue, method = "BH")]
    acat[, acat_significant := is.finite(acat_qvalue) & acat_qvalue < fdr]
  }
}

fwrite(phenotype, file.path(outdir, "phenotype_summary.tsv"), sep = "\t")
fwrite(phenotype[significant == TRUE], file.path(outdir, "significant_pcqtl_phenotypes.tsv"), sep = "\t")
fwrite(acat, file.path(outdir, "region_acat_results.tsv"), sep = "\t")
fwrite(
  phenotype[, .(tested_phenotypes = .N, significant_phenotypes = sum(significant)), by = celltype],
  file.path(outdir, "pcqtl_counts_by_celltype.tsv"), sep = "\t"
)
