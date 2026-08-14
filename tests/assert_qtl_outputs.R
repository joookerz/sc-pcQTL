#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: assert_qtl_outputs.R <outdir>")
outdir <- args[[1L]]
phenotype <- fread(file.path(outdir, "summary", "phenotype_summary.tsv"))
stopifnot(nrow(phenotype) >= 1L, phenotype$significant[[1L]], phenotype$donor_sample_size[[1L]] == 10)
variants <- fread(file.path(outdir, "summary", "all_variant_results.tsv.gz"))
stopifnot(nrow(variants) >= 2L, "pvalue" %in% names(variants), "qvalue" %in% names(variants))
acat <- fread(file.path(outdir, "summary", "region_acat_results.tsv"))
stopifnot(nrow(acat) >= 1L, "acat_qvalue" %in% names(acat), acat$acat_significant[[1L]])
published_null_models <- list.files(
  file.path(outdir, "qtl", "tasks"), pattern = "^saige_null_model",
  recursive = TRUE, full.names = TRUE
)
stopifnot(length(published_null_models) == 0L)
message("Mock SAIGE-QTL integration assertions passed")
