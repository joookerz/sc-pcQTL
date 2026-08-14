suppressPackageStartupMessages({
  library(data.table)
  library(testthat)
})

root <- normalizePath(Sys.getenv("SCPCQTL_ROOT", unset = "."))
source(file.path(root, "bin", "scpcqtl_utils.R"))

test_that("pair canonicalization is deterministic", {
  result <- canonical_pairs(c("B", "A"), c("A", "C"))
  expect_equal(result$Gene1, c("A", "A"))
  expect_equal(result$Gene2, c("B", "C"))
})

test_that("chromosomes and active covariates are normalized", {
  expect_equal(normalize_chr(c("chr1", "2", "X")), c(1L, 2L, NA_integer_))
  covars <- data.table(variable = 1:3, constant = 1, group = factor(c("a", "b", "a")))
  expect_equal(active_covariates(covars, c("variable", "constant", "group")), c("variable", "group"))
  expect_error(safe_celltype("not safe"), "Cell-type IDs")
})

test_that("staged count blocks support range and gene reads", {
  stage <- tempfile("scpcqtl-stage-")
  dir.create(stage)
  genes <- data.table(
    gene_name = paste0("G", 1:5), chromosome = 1L, start = seq(100L, 500L, 100L),
    end = seq(150L, 550L, 100L), nonzero_fraction = 1, keep = TRUE,
    chromosome_rank = 1:5
  )
  fwrite(genes, file.path(stage, "gene_filtering.tsv"), sep = "\t")
  left <- matrix(1:9, nrow = 3, dimnames = list(NULL, paste0("G", 1:3)))
  right <- matrix(10:15, nrow = 3, dimnames = list(NULL, paste0("G", 4:5)))
  saveRDS(left, file.path(stage, "block1.rds"), compress = FALSE)
  saveRDS(right, file.path(stage, "block2.rds"), compress = FALSE)
  fwrite(data.table(
    chromosome = 1L, block_id = 1:2, rank_start = c(1L, 4L), rank_end = c(3L, 5L),
    n_genes = c(3L, 2L), file = c("block1.rds", "block2.rds"),
    first_gene = c("G1", "G4"), last_gene = c("G3", "G5")
  ), file.path(stage, "count_blocks.tsv"), sep = "\t")
  expect_equal(colnames(load_staged_count_range(stage, 1, 2, 4)), c("G2", "G3", "G4"))
  expect_equal(colnames(load_staged_count_genes(stage, 1, c("G2", "G4"))), c("G2", "G4"))
})

test_that("fast and complete schedulers have the documented coverage", {
  stage <- tempfile("scpcqtl-plan-")
  dir.create(stage)
  fwrite(data.table(eligible = TRUE), file.path(stage, "celltype_qc.tsv"), sep = "\t")
  fwrite(data.table(
    gene_name = paste0("G", 1:6), chromosome = 1L, start = 1:6, end = 1:6,
    nonzero_fraction = 1, keep = TRUE, chromosome_rank = 1:6
  ), file.path(stage, "gene_filtering.tsv"), sep = "\t")
  script <- file.path(root, "bin", "plan_pair_tasks.R")
  fast <- file.path(stage, "fast.tsv")
  complete <- file.path(stage, "complete.tsv")
  expect_equal(system2("Rscript", c(script, "--stage", stage, "--out", fast,
                                     "--pair_scope", "fast", "--max_cluster_genes", "4")), 0)
  expect_equal(system2("Rscript", c(script, "--stage", stage, "--out", complete,
                                     "--pair_scope", "complete", "--max_cluster_genes", "4",
                                     "--responses_per_task", "10")), 0)
  fast_summary <- fread(file.path(stage, "fast_summary.tsv"))
  complete_summary <- fread(file.path(stage, "complete_summary.tsv"))
  expect_equal(fast_summary$global_unordered_pairs, 15)
  expect_equal(fast_summary$computed_unordered_pairs, 6)
  expect_equal(complete_summary$computed_unordered_pairs, 12)
})

test_that("cluster windows retain filtered-out annotated genes as genomic separators", {
  stage <- tempfile("scpcqtl-cluster-stage-")
  outdir <- tempfile("scpcqtl-cluster-out-")
  dir.create(stage)
  fwrite(data.table(
    gene_name = c("G1", "GAP", "G2"), chromosome = 1L,
    start = c(100L, 200L, 300L), end = c(150L, 250L, 350L),
    annotation_order = 1:3, present_in_counts = c(TRUE, FALSE, TRUE),
    nonzero_fraction = c(1, NA, 1), keep = c(TRUE, FALSE, TRUE),
    chromosome_rank = c(1L, NA_integer_, 2L)
  ), file.path(stage, "gene_filtering.tsv"), sep = "\t")
  pairs <- tempfile(fileext = ".tsv")
  fwrite(data.table(
    chromosome = 1L, Gene1 = "G1", Gene2 = "G2", significant = TRUE
  ), pairs, sep = "\t")
  status <- system2("Rscript", c(
    file.path(root, "bin", "call_clusters.R"),
    "--stage", stage, "--pairs", pairs, "--outdir", outdir,
    "--celltype", "example", "--max_cluster_genes", "2",
    "--min_cluster_genes", "2", "--cluster_density", "0.7"
  ))
  expect_equal(status, 0)
  expect_equal(nrow(fread(file.path(outdir, "clusters.tsv"))), 0L)
})

test_that("library size is reused or deterministically computed", {
  source_counts <- fread(file.path(root, "tests", "fixtures", "counts.tsv"))
  expected_total <- rowSums(as.matrix(source_counts[, .(G1, G2, G3, G4)]))
  input <- tempfile(fileext = ".tsv")
  fwrite(source_counts[, !c("total_read_counts", "log_total_read_counts")], input, sep = "\t")
  outdir <- tempfile("scpcqtl-prepare-")
  script <- file.path(root, "bin", "prepare_celltype.R")
  status <- system2("Rscript", c(
    script, "--celltype", "example", "--counts", input,
    "--annotation", file.path(root, "tests", "fixtures", "genes.tsv"),
    "--outdir", outdir, "--covariates", "sex", "--categorical_covariates", "sex",
    "--min_cells", "10", "--min_nonzero_fraction", "0"
  ))
  expect_equal(status, 0)
  covars <- readRDS(file.path(outdir, "covariates.rds"))
  expect_equal(covars$total_read_counts, expected_total)
  expect_equal(covars$log_total_read_counts, log(pmax(expected_total, 1)))
  qc <- fread(file.path(outdir, "celltype_qc.tsv"))
  expect_equal(qc$library_size_source, "computed_from_annotated_expression_columns")
})

test_that("cell types below the cell-count threshold stop before expression staging", {
  outdir <- tempfile("scpcqtl-ineligible-")
  status <- system2("Rscript", c(
    file.path(root, "bin", "prepare_celltype.R"),
    "--celltype", "small", "--counts", file.path(root, "tests", "fixtures", "counts.tsv"),
    "--annotation", file.path(root, "tests", "fixtures", "genes.tsv"),
    "--outdir", outdir, "--covariates", "sex", "--categorical_covariates", "sex",
    "--min_cells", "21"
  ))
  expect_equal(status, 0)
  qc <- fread(file.path(outdir, "celltype_qc.tsv"))
  expect_false(qc$eligible)
  expect_false(file.exists(file.path(outdir, "covariates.rds")))
  expect_equal(nrow(fread(file.path(outdir, "count_blocks.tsv"))), 0)
})

test_that("samplesheet validation rejects duplicate cell types", {
  sheet <- tempfile(fileext = ".csv")
  fwrite(data.table(
    celltype = c("duplicate", "duplicate"),
    counts = c("counts.tsv", "counts.tsv")
  ), sheet, sep = ",")
  status <- suppressWarnings(system2("Rscript", c(
    file.path(root, "bin", "validate_samplesheet.R"), "--input", sheet,
    "--base_dir", file.path(root, "tests", "fixtures"), "--out", tempfile(fileext = ".csv")
  ), stdout = FALSE, stderr = FALSE))
  expect_false(status == 0)
})

test_that("SAIGE-QTL overrides are resolved while workflow-owned fields are protected", {
  resolver <- file.path(root, "bin", "resolve_saige_params.R")
  defaults <- file.path(root, "assets", "saigeqtl_defaults.tsv")
  user <- tempfile(fileext = ".tsv")
  output <- tempfile(fileext = ".tsv")
  fwrite(data.table(step = "step1", parameter = "maxiter", value = "30"), user, sep = "\t")
  expect_equal(system2("Rscript", c(
    resolver, "--defaults", defaults, "--user", user, "--out", output,
    "--covariates", "age,sex", "--qtl_maf", "0.1"
  )), 0)
  resolved <- fread(output, colClasses = "character")
  expect_equal(resolved[step == "step1" & parameter == "maxiter", value], "30")
  expect_equal(resolved[step == "step1" & parameter == "covarColList", value], "age,sex")
  expect_equal(resolved[step == "step2" & parameter == "minMAF", value], "0.1")

  fwrite(data.table(step = "step1", parameter = "traitType", value = "binary"), user, sep = "\t")
  status <- suppressWarnings(system2("Rscript", c(
    resolver, "--defaults", defaults, "--user", user, "--out", output
  ), stdout = FALSE, stderr = FALSE))
  expect_false(status == 0)
})
