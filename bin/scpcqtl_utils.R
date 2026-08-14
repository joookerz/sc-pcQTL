suppressPackageStartupMessages(library(data.table))

parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) stop("Unexpected argument: ", token)
    token <- substring(token, 3L)
    if (grepl("=", token, fixed = TRUE)) {
      bits <- strsplit(token, "=", fixed = TRUE)[[1L]]
      out[[bits[[1L]]]] <- paste(bits[-1L], collapse = "=")
      i <- i + 1L
    } else {
      if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
        out[[token]] <- "true"
        i <- i + 1L
      } else {
        out[[token]] <- args[[i + 1L]]
        i <- i + 2L
      }
    }
  }
  out
}

required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(value)) stop("Missing --", name)
  value
}

as_flag <- function(value) {
  tolower(as.character(value)) %in% c("1", "true", "t", "yes", "y")
}

split_csv_arg <- function(value) {
  values <- trimws(strsplit(value %||% "", ",", fixed = TRUE)[[1L]])
  values[nzchar(values)]
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

normalize_chr <- function(x) {
  value <- sub("^chr", "", as.character(x), ignore.case = TRUE)
  suppressWarnings(as.integer(value))
}

safe_celltype <- function(x) {
  if (!grepl("^[A-Za-z0-9_.-]+$", x)) {
    stop("Cell-type IDs may contain only letters, numbers, '.', '_' and '-': ", x)
  }
  x
}

active_covariates <- function(data, variables) {
  variables[vapply(variables, function(name) {
    values <- data[[name]]
    length(unique(values[!is.na(values)])) > 1L
  }, logical(1L))]
}

load_staged_count_range <- function(stage, chromosome, first_rank, last_rank) {
  manifest_file <- file.path(stage, "count_blocks.tsv")
  if (!file.exists(manifest_file)) stop("Missing staged count-block manifest: ", manifest_file)
  manifest <- fread(manifest_file)
  selected <- manifest[
    chromosome == as.integer(chromosome) & rank_end >= as.integer(first_rank) &
      rank_start <= as.integer(last_rank)
  ][order(rank_start)]
  if (!nrow(selected)) {
    stop("No staged count blocks cover chromosome ", chromosome, " ranks ", first_rank, "-", last_rank)
  }
  matrices <- lapply(selected$file, function(name) readRDS(file.path(stage, name)))
  counts <- do.call(cbind, matrices)
  expected <- seq.int(as.integer(first_rank), as.integer(last_rank))
  genes <- fread(file.path(stage, "gene_filtering.tsv"))[
    keep == TRUE & chromosome == as.integer(chromosome) & chromosome_rank %in% expected
  ][order(chromosome_rank), gene_name]
  missing <- setdiff(genes, colnames(counts))
  if (length(missing)) stop("Staged count blocks are missing genes: ", paste(head(missing, 20L), collapse = ", "))
  counts[, genes, drop = FALSE]
}

load_staged_count_genes <- function(stage, chromosome, genes) {
  gene_table <- fread(file.path(stage, "gene_filtering.tsv"))[
    keep == TRUE & chromosome == as.integer(chromosome) & gene_name %in% genes
  ]
  if (nrow(gene_table) != length(unique(genes))) {
    missing <- setdiff(genes, gene_table$gene_name)
    stop("Requested genes are absent from staged counts: ", paste(missing, collapse = ", "))
  }
  counts <- load_staged_count_range(
    stage, chromosome, min(gene_table$chromosome_rank), max(gene_table$chromosome_rank)
  )
  counts[, genes, drop = FALSE]
}

write_status <- function(path, status, fields = list()) {
  row <- as.data.table(c(list(status = status), fields))
  fwrite(row, path, sep = "\t")
}

canonical_pairs <- function(gene1, gene2) {
  data.table(Gene1 = pmin(gene1, gene2), Gene2 = pmax(gene1, gene2))
}

empty_pair_table <- function() {
  data.table(
    celltype = character(), chromosome = integer(), response = character(),
    predictor = character(), response_rank = integer(), predictor_rank = integer(),
    beta_count = numeric(), se_count = numeric(), p_count = numeric(),
    beta_detection = numeric(), se_detection = numeric(), p_detection = numeric(),
    stat_count = numeric(), stat_detection = numeric(), stat_joint = numeric(),
    log_p_joint = numeric(), neglog10_p_joint = numeric(), p_joint = numeric(),
    test_status = character()
  )
}
