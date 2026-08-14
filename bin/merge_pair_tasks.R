#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))[1L]])))
source(file.path(script_dir, "scpcqtl_utils.R"))
args <- parse_cli()
stage <- required_arg(args, "stage")
input_list <- required_arg(args, "input_list")
outdir <- required_arg(args, "outdir")
test <- args$pair_test %||% "component_union"
scope <- args$pair_scope %||% "fast"
count_family <- args$count_family %||% "poisson"
alpha <- as.numeric(args$pair_alpha %||% 0.05)
save_directional <- as_flag(args$save_directional_tests %||% "false")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

files <- fread(input_list, header = FALSE)[[1L]]
files <- files[file.exists(files)]
directional <- if (length(files)) rbindlist(lapply(files, fread), fill = TRUE) else empty_pair_table()
genes <- fread(file.path(stage, "gene_filtering.tsv"))[keep == TRUE & chromosome %in% 1:22]
counts <- genes[, .N, by = chromosome]
m <- counts[, sum(N * (N - 1) / 2)]

if (nrow(directional)) {
  directional[, `:=`(Gene1 = pmin(response, predictor), Gene2 = pmax(response, predictor))]
  if (test == "component_union") {
    merged <- directional[, {
      ic <- if (any(is.finite(p_count))) which.min(replace(p_count, !is.finite(p_count), Inf)) else NA_integer_
      iz <- if (any(is.finite(p_detection))) which.min(replace(p_detection, !is.finite(p_detection), Inf)) else NA_integer_
      list(
        p_count = if (is.na(ic)) NA_real_ else p_count[ic],
        count_direction = if (is.na(ic)) NA_character_ else paste0(response[ic], "<-", predictor[ic]),
        beta_count = if (is.na(ic)) NA_real_ else beta_count[ic],
        se_count = if (is.na(ic)) NA_real_ else se_count[ic],
        p_detection = if (is.na(iz)) NA_real_ else p_detection[iz],
        detection_direction = if (is.na(iz)) NA_character_ else paste0(response[iz], "<-", predictor[iz]),
        beta_detection = if (is.na(iz)) NA_real_ else beta_detection[iz],
        se_detection = if (is.na(iz)) NA_real_ else se_detection[iz]
      )
    }, by = .(celltype, chromosome, Gene1, Gene2)]
    threshold <- alpha / m
    merged[, significant := (is.finite(p_count) & p_count < threshold) |
                              (is.finite(p_detection) & p_detection < threshold)]
    merged[, threshold := threshold]
  } else {
    merged <- directional[, .SD[which.min(replace(log_p_joint, !is.finite(log_p_joint), Inf))],
                          by = .(celltype, chromosome, Gene1, Gene2)]
    threshold <- alpha / (2 * m)
    merged[, significant := is.finite(log_p_joint) & log_p_joint < log(threshold)]
    merged[, `:=`(best_direction = paste0(response, "<-", predictor), threshold = threshold)]
  }
} else {
  merged <- data.table(celltype = character(), chromosome = integer(), Gene1 = character(), Gene2 = character(), significant = logical())
  threshold <- if (test == "joint_score") alpha / (2 * m) else alpha / m
}

fwrite(merged, file.path(outdir, "all_computed_pairs.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(merged[significant == TRUE], file.path(outdir, "significant_pairs.tsv.gz"), sep = "\t", compress = "gzip")
if (save_directional) fwrite(directional, file.path(outdir, "directional_tests.tsv.gz"), sep = "\t", compress = "gzip")
fwrite(data.table(global_unordered_pairs = m, computed_unordered_pairs = nrow(merged),
                  significant_pairs = sum(merged$significant), threshold = threshold,
                  pair_alpha = alpha, pair_scope = scope, pair_test = test,
                  count_family = count_family,
                  fasthurdle_version = as.character(packageVersion("fasthurdle"))),
       file.path(outdir, "pair_summary.tsv"), sep = "\t")
