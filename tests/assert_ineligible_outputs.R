#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: assert_ineligible_outputs.R <outdir>")
outdir <- args[[1L]]
qc <- fread(file.path(outdir, "qc", "celltypes", "example", "celltype_qc.tsv"))
stopifnot(!qc$eligible[[1L]], qc$n_filtered_genes[[1L]] == 0L)
stopifnot(!dir.exists(file.path(outdir, "pairs", "example")))
stopifnot(file.exists(file.path(outdir, "pipeline_info", "analysis_parameters.json")))
message("Ineligible-cell-type integration assertions passed")
