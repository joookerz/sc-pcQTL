#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))[1L]])))
source(file.path(script_dir, "scpcqtl_utils.R"))
args <- parse_cli()
input <- required_arg(args, "input")
base_dir <- required_arg(args, "base_dir")
out <- required_arg(args, "out")

samples <- fread(input, colClasses = "character")
required <- c("celltype", "counts")
missing <- setdiff(required, names(samples))
if (length(missing)) stop("Samplesheet is missing columns: ", paste(missing, collapse = ", "))
samples <- samples[, ..required]
if (!nrow(samples)) stop("Samplesheet has no data rows")
if (anyNA(samples) || any(!nzchar(samples$celltype)) || any(!nzchar(samples$counts))) {
  stop("Samplesheet celltype and counts values must be non-empty")
}
invisible(vapply(samples$celltype, safe_celltype, character(1L)))
if (anyDuplicated(samples$celltype)) {
  duplicate <- unique(samples$celltype[duplicated(samples$celltype)])
  stop("Samplesheet contains duplicate cell types: ", paste(duplicate, collapse = ", "))
}
samples[, counts := vapply(counts, function(path) {
  expanded <- path.expand(path)
  candidate <- if (grepl("^/", expanded)) expanded else file.path(base_dir, expanded)
  normalizePath(candidate, mustWork = TRUE)
}, character(1L))]
fwrite(samples, out, sep = ",", quote = TRUE)
