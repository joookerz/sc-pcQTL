#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x

parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) {
      stop("Arguments must be supplied as --name value pairs")
    }
    out[[sub("^--", "", x[[i]])]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}

required <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(value)) stop("Missing --", name)
  value
}

read_table <- function(path, header = TRUE) {
  read.table(path, header = header, stringsAsFactors = FALSE, check.names = FALSE,
             comment.char = "", quote = "")
}

scale_numeric <- function(x) as.numeric(scale(as.numeric(x)))

sample_zero_truncated_poisson <- function(mu) {
  value <- rpois(length(mu), mu)
  zero <- which(value == 0L)
  while (length(zero)) {
    value[zero] <- rpois(length(zero), mu[zero])
    zero <- zero[value[zero] == 0L]
  }
  value
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
fam_file <- required(args, "fam")
eigenvec_file <- required(args, "eigenvec")
dosage_file <- required(args, "dosage")
causal_file <- required(args, "causal")
outdir <- required(args, "outdir")
n_cells <- as.integer(args[["n-cells"]] %||% 12000L)
seed <- as.integer(args$seed %||% 20260814L)
set.seed(seed)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

fam <- read_table(fam_file, header = FALSE)
names(fam) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")
if (nrow(fam) != 503L) stop("Expected 503 unrelated EUR donors; found ", nrow(fam))

pcs <- read_table(eigenvec_file)
names(pcs) <- sub("^#", "", names(pcs))
if (!all(c("FID", "IID") %in% names(pcs))) {
  pcs <- read_table(eigenvec_file, header = FALSE)
  if (ncol(pcs) != 8L) stop("Unexpected PLINK eigenvector format")
  names(pcs) <- c("FID", "IID", paste0("PC", 1:6))
}
pc_cols <- grep("^PC[1-6]$", names(pcs), value = TRUE)
if (length(pc_cols) != 6L) stop("Expected genotype PC1-PC6")

dosage <- read_table(dosage_file)
causal <- read_table(causal_file)
fixed_raw <- intersect(c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"), names(dosage))
dosage_cols <- setdiff(names(dosage), fixed_raw)
if (nrow(causal) != 4L) stop("Expected four causal variants")
for (i in seq_len(nrow(causal))) {
  hit <- dosage_cols[dosage_cols == causal$snp_id[[i]] |
                       startsWith(dosage_cols, paste0(causal$snp_id[[i]], "_"))]
  if (length(hit) != 1L) stop("Could not identify dosage column for ", causal$snp_id[[i]])
  names(dosage)[names(dosage) == hit] <- paste0("G", i)
}

donors <- merge(fam[, c("FID", "IID")], pcs[, c("FID", "IID", pc_cols)], by = c("FID", "IID"))
donors <- merge(donors, dosage[, c("FID", "IID", paste0("G", 1:4))], by = c("FID", "IID"))
donors <- donors[match(fam$IID, donors$IID), , drop = FALSE]
if (nrow(donors) != 503L || anyDuplicated(donors$IID)) stop("Donor metadata merge failed")
for (name in pc_cols) donors[[name]] <- scale_numeric(donors[[name]])
for (name in paste0("G", 1:4)) donors[[name]] <- scale_numeric(donors[[name]])
donors$age <- pmin(80, pmax(18, round(rnorm(nrow(donors), 48, 15))))
donors$sex <- sample(0L:1L, nrow(donors), replace = TRUE)
donors$pf1 <- rnorm(nrow(donors))
donors$pf2 <- 0.25 * donors$pf1 + sqrt(1 - 0.25^2) * rnorm(nrow(donors))
names(donors)[match(paste0("PC", 1:6), names(donors))] <- paste0("pc", 1:6)

minimum_cells <- 10L
if (n_cells < nrow(donors) * minimum_cells) stop("n-cells is too small")
weights <- exp(rnorm(nrow(donors), 0, 0.55))
extra <- as.vector(rmultinom(1L, n_cells - nrow(donors) * minimum_cells, weights / sum(weights)))
cell_counts <- minimum_cells + extra
donor_index <- rep(seq_len(nrow(donors)), cell_counts)
cell_meta <- donors[donor_index, , drop = FALSE]
cell_meta$barcode <- sprintf("sim_immune_cell_%05d", seq_len(nrow(cell_meta)))
cell_meta$individual <- cell_meta$IID

log_library <- pmin(log(15000), pmax(log(1800), rnorm(n_cells, log(5200), 0.35)))
cell_meta$total_read_counts <- as.integer(round(exp(log_library)))
cell_meta$log_total_read_counts <- log_library
library_z <- scale_numeric(log_library)
age_z <- scale_numeric(cell_meta$age)
sex_numeric <- as.numeric(cell_meta$sex)

module_table <- rbind(
  data.frame(module = "C1_count", gene_name = c("SIM_C1_G1", "SIM_C1_G2"),
             position = c(17970000L, 18030000L), loading = c(1.00, 0.90)),
  data.frame(module = "C2_detection", gene_name = c("SIM_C2_G1", "SIM_C2_G2"),
             position = c(24970000L, 25030000L), loading = c(1.00, 0.90)),
  data.frame(module = "C3_joint", gene_name = paste0("SIM_C3_G", 1:4),
             position = c(31925000L, 31975000L, 32025000L, 32075000L),
             loading = c(1.00, 0.90, 0.80, 0.70)),
  data.frame(module = "C4_mixed", gene_name = paste0("SIM_C4_G", 1:5),
             position = c(38900000L, 38950000L, 39000000L, 39050000L, 39100000L),
             loading = c(1.00, 0.80, -0.80, -1.00, 0.60)),
  data.frame(module = "C5_null", gene_name = c("SIM_C5_G1", "SIM_C5_G2"),
             position = c(45970000L, 46030000L), loading = c(1.00, 0.90))
)

background_positions <- c(
  15000000L, 15400000L, 15800000L,
  20000000L, 20400000L, 20800000L, 21200000L, 21600000L, 22000000L,
  27000000L, 27500000L, 28000000L, 28500000L, 29000000L,
  34000000L, 34500000L, 35000000L, 35500000L, 36000000L,
  41000000L, 41600000L, 42200000L, 42800000L,
  48000000L, 48500000L
)
background <- data.frame(
  module = "background",
  gene_name = sprintf("SIM_BG%02d", seq_along(background_positions)),
  position = background_positions,
  loading = 0
)
genes <- rbind(module_table, background)
genes <- genes[order(genes$position), , drop = FALSE]
rownames(genes) <- NULL
if (nrow(genes) != 40L || anyDuplicated(genes$gene_name)) stop("Synthetic gene design is invalid")

gene_annotation <- data.frame(
  gene_name = genes$gene_name,
  chromosome = 22L,
  start = genes$position - 1000L,
  end = genes$position + 1000L
)
write.table(gene_annotation, file.path(outdir, "gene_annotation.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

module_index <- setNames(seq_len(5L), c("C1_count", "C2_detection", "C3_joint", "C4_mixed", "C5_null"))
module_state <- matrix(sample(c(-1, 1), n_cells * 5L, replace = TRUE),
                       nrow = n_cells, ncol = 5L)
donor_detection_effect <- matrix(
  rnorm(nrow(donors) * nrow(genes), sd = 0.55),
  nrow = nrow(donors), ncol = nrow(genes)
)
donor_count_effect <- matrix(
  rnorm(nrow(donors) * nrow(genes), sd = 0.15),
  nrow = nrow(donors), ncol = nrow(genes)
)
module_parameters <- data.frame(
  module = names(module_index),
  detection_loading = c(0.05, 0.00, 4.00, 4.00, 0.00),
  genotype_detection = c(0.04, 0.85, 0.48, 0.55, 0.00),
  genotype_count = c(0.16, 0.00, 0.00, 0.00, 0.00)
)

expression <- matrix(0L, nrow = n_cells, ncol = nrow(genes),
                     dimnames = list(NULL, genes$gene_name))
baseline_detection <- seq(0.28, 0.52, length.out = nrow(genes))
target_positive_mean <- rep(1.0, nrow(genes))
target_positive_mean[genes$module == "C1_count"] <- 2.0
target_positive_mean[genes$module == "C2_detection"] <- 0.03
target_positive_mean[genes$module == "C3_joint"] <- 1.20
target_positive_mean[genes$module == "C4_mixed"] <- 0.03
target_positive_mean[genes$module == "C5_null"] <- 0.03

for (j in seq_len(nrow(genes))) {
  gene <- genes[j, ]
  module_name <- gene$module
  shared_detection <- shared_count <- genotype_detection <- genotype_count <- rep(0, n_cells)
  if (module_name != "background") {
    m <- module_index[[module_name]]
    parameter <- module_parameters[module_parameters$module == module_name, , drop = FALSE]
    loading <- gene$loading
    shared_detection <- loading * parameter$detection_loading *
      module_state[, m]
    if (m <= 4L) {
      genotype_detection <- loading * parameter$genotype_detection * cell_meta[[paste0("G", m)]]
      genotype_count <- loading * parameter$genotype_count * cell_meta[[paste0("G", m)]]
    }
  }

  direct_detection_effect <- rep(0, n_cells)
  if (gene$gene_name == "SIM_C2_G2") {
    direct_detection_effect <- 8.0 * (expression[, "SIM_C2_G1"] > 0) - 4.0
  } else if (gene$gene_name == "SIM_C5_G2") {
    direct_detection_effect <- 8.0 * (expression[, "SIM_C5_G1"] > 0) - 4.0
  }

  eta_detection <- qlogis(baseline_detection[[j]]) +
    0.28 * library_z + 0.04 * age_z + 0.08 * sex_numeric +
    0.08 * cell_meta$pf1 - 0.05 * cell_meta$pf2 +
    donor_detection_effect[donor_index, j] + shared_detection +
    genotype_detection + direct_detection_effect + rnorm(n_cells, sd = 0.05)
  detected <- rbinom(n_cells, 1L, plogis(eta_detection)) == 1L

  direct_count_effect <- rep(0, n_cells)
  if (gene$gene_name == "SIM_C1_G2") {
    direct_count_effect <- 0.04 * expression[, "SIM_C1_G1"]
  } else if (module_name == "C3_joint" && gene$gene_name != "SIM_C3_G1") {
    direct_count_effect <- 0.015 * gene$loading * expression[, "SIM_C3_G1"]
  }
  eta_count <- log_library + log(target_positive_mean[[j]] / 5200) +
    0.03 * age_z + 0.05 * sex_numeric + 0.05 * cell_meta$pf1 +
    donor_count_effect[donor_index, j] + shared_count + genotype_count +
    direct_count_effect
  mu <- pmin(35, pmax(0.02, exp(eta_count)))
  expression[detected, j] <- sample_zero_truncated_poisson(mu[detected])
}

metadata_columns <- c(
  "individual", "barcode", "age", "sex", paste0("pc", 1:6), "pf1", "pf2",
  "total_read_counts", "log_total_read_counts"
)
counts <- cbind(cell_meta[, metadata_columns, drop = FALSE], as.data.frame(expression, check.names = FALSE))
connection <- gzfile(file.path(outdir, "counts.tsv.gz"), open = "wt")
write.table(counts, connection, sep = "\t", quote = FALSE, row.names = FALSE)
close(connection)
write.table(data.frame(celltype = "sim_immune", counts = "counts.tsv.gz"),
            file.path(outdir, "samplesheet.csv"), sep = ",", quote = FALSE, row.names = FALSE)

message(
  "Generated ", length(unique(counts$individual)), " donors, ", nrow(counts),
  " cells, and ", nrow(genes), " genes; gene nonzero fractions range from ",
  format(min(colMeans(expression > 0)), digits = 3), " to ",
  format(max(colMeans(expression > 0)), digits = 3), "."
)
