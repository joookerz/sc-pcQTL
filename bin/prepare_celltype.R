#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))[1L]])))
source(file.path(script_dir, "scpcqtl_utils.R"))
args <- parse_cli()

celltype <- safe_celltype(required_arg(args, "celltype"))
counts_file <- required_arg(args, "counts")
annotation_file <- required_arg(args, "annotation")
outdir <- required_arg(args, "outdir")
donor_col <- args$donor_col %||% "individual"
cell_id_col <- args$cell_id_col %||% "barcode"
gene_col <- args$gene_col %||% "gene_name"
chr_col <- args$chromosome_col %||% "chromosome"
start_col <- args$start_col %||% "start"
end_col <- args$end_col %||% "end"
covariates <- split_csv_arg(args$covariates %||% "age,sex,pc1,pc2,pc3,pc4,pc5,pc6,pf1,pf2")
categorical <- split_csv_arg(args$categorical_covariates %||% "sex")
total_col <- args$total_library_col %||% "total_read_counts"
log_col <- args$log_library_col %||% "log_total_read_counts"
min_cells <- as.integer(args$min_cells %||% 10000L)
nonzero_cutoff <- as.numeric(args$min_nonzero_fraction %||% 0.01)
genome_build <- args$genome_build %||% "GRCh38"
count_block_size <- as.integer(args$count_block_size %||% 100L)
if (!is.finite(count_block_size) || count_block_size < 1L) stop("count_block_size must be a positive integer")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
ann <- fread(annotation_file)
required_ann <- c(gene_col, chr_col, start_col, end_col)
missing_ann <- setdiff(required_ann, names(ann))
if (length(missing_ann)) stop("Missing annotation columns: ", paste(missing_ann, collapse = ", "))
setnames(ann, required_ann, c("gene_name", "chromosome", "start", "end"))
ann[, annotation_order := .I]
ann[, chromosome := normalize_chr(chromosome)]
ann[, `:=`(start = as.integer(start), end = as.integer(end))]
ann <- ann[chromosome %in% 1:22 & is.finite(start) & is.finite(end)]
if (anyDuplicated(ann$gene_name)) stop("Gene annotation contains duplicate gene names")

header <- names(fread(counts_file, nrows = 0L))
required_cols <- c(donor_col, cell_id_col, covariates)
missing_cols <- setdiff(required_cols, header)
if (length(missing_cols)) stop("Missing count/covariate columns: ", paste(missing_cols, collapse = ", "))
gene_names <- intersect(ann$gene_name, header)
if (length(gene_names) < 2L) stop("Fewer than two annotated genes occur in the count table")

metadata_cols <- unique(c(donor_col, cell_id_col, covariates, total_col, log_col))
metadata <- fread(counts_file, select = intersect(metadata_cols, header))
if (anyDuplicated(metadata[[cell_id_col]])) stop("Cell IDs are not unique in ", celltype)
if (any(!nzchar(as.character(metadata[[cell_id_col]]))) || any(!nzchar(as.character(metadata[[donor_col]])))) {
  stop("Cell and donor identifiers must be non-empty")
}
varying_within_donor <- covariates[vapply(covariates, function(name) {
  any(metadata[, uniqueN(get(name)), by = donor_col]$V1 > 1L)
}, logical(1L))]
if (length(varying_within_donor)) {
  stop("Donor-level covariates vary within donor: ", paste(varying_within_donor, collapse = ", "))
}
n_cells <- nrow(metadata)
eligible <- n_cells >= min_cells

if (!eligible) {
  genes <- copy(ann)
  genes[, `:=`(
    present_in_counts = gene_name %in% gene_names,
    nonzero_fraction = NA_real_, keep = FALSE, chromosome_rank = NA_integer_
  )]
  setorder(genes, chromosome, start, annotation_order)
  fwrite(genes, file.path(outdir, "gene_filtering.tsv"), sep = "\t")
  fwrite(data.table(
    celltype = celltype, genome_build = genome_build,
    n_cells = n_cells, n_donors = uniqueN(metadata[[donor_col]]),
    n_annotation_genes = nrow(ann), n_annotated_genes = length(gene_names),
    n_filtered_genes = 0L,
    min_cells = min_cells, eligible = FALSE,
    exclusion_reason = sprintf("n_cells_below_%d", min_cells),
    library_size_source = "not_computed_for_ineligible_celltype"
  ), file.path(outdir, "celltype_qc.tsv"), sep = "\t")
  fwrite(data.table(
    chromosome = integer(), block_id = integer(), rank_start = integer(), rank_end = integer(),
    n_genes = integer(), file = character(), first_gene = character(), last_gene = character()
  ), file.path(outdir, "count_blocks.tsv"), sep = "\t")
  writeLines("OK", file.path(outdir, "COMPLETE"))
  quit(status = 0L)
}

has_total <- total_col %in% names(metadata)
has_log <- log_col %in% names(metadata)
total_library <- if (has_total) as.numeric(metadata[[total_col]]) else rep.int(0, n_cells)
log_library <- if (has_log) as.numeric(metadata[[log_col]]) else NULL
if (has_total && (any(!is.finite(total_library)) || any(total_library < 0))) {
  stop("Total library size must be finite and non-negative")
}
if (has_log && any(!is.finite(log_library))) stop("Log library size must be finite")

gene_rows <- list()
block_rows <- list()
gene_row_index <- 0L
block_index <- 0L
for (chr in 1:22) {
  # The manuscript workflow evaluated cluster windows over the complete
  # ordered annotation, while pair tests used only expression-filtered genes.
  chr_info <- copy(ann[chromosome == chr][order(start, annotation_order)])
  if (!nrow(chr_info)) next
  chr_info[, `:=`(
    present_in_counts = gene_name %in% gene_names,
    nonzero_fraction = NA_real_, keep = FALSE,
    chromosome_rank = NA_integer_
  )]
  chr_names <- chr_info[present_in_counts == TRUE, gene_name]
  if (!length(chr_names)) {
    gene_row_index <- gene_row_index + 1L
    gene_rows[[gene_row_index]] <- chr_info
    next
  }
  chr_dt <- fread(counts_file, select = chr_names)
  if (nrow(chr_dt) != n_cells) stop("Row count changed while reading chromosome ", chr)
  chr_mat <- as.matrix(chr_dt)
  storage.mode(chr_mat) <- "double"
  if (any(!is.finite(chr_mat)) || any(chr_mat < 0)) {
    stop("Expression values must be finite and non-negative on chromosome ", chr)
  }
  if (!has_total) total_library <- total_library + rowSums(chr_mat)
  chr_info[match(chr_names, gene_name), nonzero_fraction := colMeans(chr_mat > 0)]
  chr_info[present_in_counts == TRUE, keep := nonzero_fraction >= nonzero_cutoff]
  kept_names <- chr_info[keep == TRUE, gene_name]
  if (length(kept_names)) chr_info[keep == TRUE, chromosome_rank := seq_len(.N)]
  gene_row_index <- gene_row_index + 1L
  gene_rows[[gene_row_index]] <- chr_info

  if (length(kept_names)) {
    for (first in seq.int(1L, length(kept_names), by = count_block_size)) {
      last <- min(length(kept_names), first + count_block_size - 1L)
      names_block <- kept_names[first:last]
      matrix_block <- chr_mat[, match(names_block, chr_names), drop = FALSE]
      colnames(matrix_block) <- names_block
      block_index <- block_index + 1L
      filename <- sprintf("chr%d_block%04d_counts.rds", chr, block_index)
      saveRDS(matrix_block, file.path(outdir, filename), compress = FALSE, version = 3)
      block_rows[[block_index]] <- data.table(
        chromosome = chr, block_id = block_index, rank_start = first, rank_end = last,
        n_genes = length(names_block), file = filename,
        first_gene = names_block[[1L]], last_gene = names_block[[length(names_block)]]
      )
    }
  }
  rm(chr_dt, chr_mat)
  gc(verbose = FALSE)
}

genes <- rbindlist(gene_rows)
setorder(genes, chromosome, start, annotation_order)
kept <- genes[keep == TRUE]
blocks <- if (length(block_rows)) rbindlist(block_rows) else data.table(
  chromosome = integer(), block_id = integer(), rank_start = integer(), rank_end = integer(),
  n_genes = integer(), file = character(), first_gene = character(), last_gene = character()
)
fwrite(genes, file.path(outdir, "gene_filtering.tsv"), sep = "\t")
fwrite(blocks, file.path(outdir, "count_blocks.tsv"), sep = "\t")

if (!has_log) log_library <- log(pmax(total_library, 1))
cov_dt <- metadata[, c(donor_col, cell_id_col, covariates), with = FALSE]
cov_dt[, `:=`(total_read_counts = total_library, log_total_read_counts = log_library)]
setnames(cov_dt, c(donor_col, cell_id_col), c("individual", "cell_id"))
for (name in intersect(categorical, names(cov_dt))) cov_dt[[name]] <- as.factor(cov_dt[[name]])
for (name in setdiff(covariates, categorical)) cov_dt[[name]] <- suppressWarnings(as.numeric(cov_dt[[name]]))
if (any(!complete.cases(cov_dt))) stop("Missing or non-numeric values in IDs, covariates, or library size")
saveRDS(cov_dt, file.path(outdir, "covariates.rds"), compress = "gzip", version = 3)

library_source <- if (has_total && has_log) "input" else if (!has_total && !has_log) {
  "computed_from_annotated_expression_columns"
} else {
  "partially_input_and_partially_computed"
}
fwrite(data.table(
  celltype = celltype, genome_build = genome_build,
  n_cells = n_cells, n_donors = uniqueN(cov_dt$individual),
  n_annotation_genes = nrow(ann), n_annotated_genes = length(gene_names),
  n_filtered_genes = nrow(kept),
  min_cells = min_cells, eligible = TRUE, exclusion_reason = "",
  library_size_source = library_source
), file.path(outdir, "celltype_qc.tsv"), sep = "\t")
writeLines("OK", file.path(outdir, "COMPLETE"))
