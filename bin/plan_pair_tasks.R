#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))[1L]])))
source(file.path(script_dir, "scpcqtl_utils.R"))
args <- parse_cli()
stage <- required_arg(args, "stage")
out <- required_arg(args, "out")
scope <- args$pair_scope %||% "fast"
max_genes <- as.integer(args$max_cluster_genes %||% 50L)
responses_per_task <- as.integer(args$responses_per_task %||% 10L)
if (!scope %in% c("fast", "complete")) stop("pair_scope must be fast or complete")

qc <- fread(file.path(stage, "celltype_qc.tsv"))
genes <- fread(file.path(stage, "gene_filtering.tsv"))[keep == TRUE & chromosome %in% 1:22]
header <- data.table(
  task_id = integer(), chromosome = integer(), response_start = integer(),
  response_end = integer(), block_start = integer(), block_end = integer()
)
summary_file <- sub("\\.tsv$", "_summary.tsv", out)
if (!isTRUE(qc$eligible[[1L]]) || nrow(genes) < 2L) {
  fwrite(header, out, sep = "\t")
  fwrite(data.table(
    pair_scope = scope, global_unordered_pairs = 0, computed_unordered_pairs = 0,
    coverage_fraction = NA_real_, n_tasks = 0L
  ), summary_file, sep = "\t")
  quit(status = 0L)
}

rows <- list()
k <- 0L
for (chr in 1:22) {
  n <- genes[chromosome == chr, .N]
  if (n < 2L) next
  if (scope == "fast") {
    n_blocks <- ceiling(n / max_genes)
    genes_per_block <- ceiling(n / n_blocks)
    for (block in seq_len(n_blocks)) {
      first <- (block - 1L) * genes_per_block + 1L
      last <- min(block * genes_per_block, n)
      if (last <= first) next
      for (response_first in seq.int(first, last, by = responses_per_task)) {
        k <- k + 1L
        rows[[k]] <- data.table(
          task_id = k, chromosome = chr,
          response_start = response_first,
          response_end = min(last, response_first + responses_per_task - 1L),
          block_start = first, block_end = last
        )
      }
    }
  } else {
    for (first in seq.int(1L, n, by = responses_per_task)) {
      k <- k + 1L
      rows[[k]] <- data.table(task_id = k, chromosome = chr,
                              response_start = first,
                              response_end = min(n, first + responses_per_task - 1L),
                              block_start = 1L, block_end = n)
    }
  }
}
tasks <- if (length(rows)) rbindlist(rows) else header
fwrite(tasks, out, sep = "\t")

counts <- genes[, .(n_genes = .N), by = chromosome]
m <- counts[, sum(n_genes * (n_genes - 1) / 2)]
computed <- if (scope == "fast") {
  unique(tasks[, .(chromosome, block_start, block_end)])[
    , sum((block_end - block_start + 1) * (block_end - block_start) / 2)
  ]
} else {
  counts[, sum(vapply(n_genes, function(n) sum(pmin(max_genes - 1L, n - seq_len(n))), numeric(1)))]
}
fwrite(data.table(pair_scope = scope, global_unordered_pairs = m,
                  computed_unordered_pairs = computed,
                  coverage_fraction = computed / m, n_tasks = nrow(tasks)),
       summary_file, sep = "\t")
