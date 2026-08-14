#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(data.table); library(fasthurdle)})
script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))[1L]])))
source(file.path(script_dir, "scpcqtl_utils.R"))
args <- parse_cli()
stage <- required_arg(args, "stage")
out <- required_arg(args, "out")
celltype <- required_arg(args, "celltype")
chr <- as.integer(required_arg(args, "chromosome"))
response_start <- as.integer(required_arg(args, "response_start"))
response_end <- as.integer(required_arg(args, "response_end"))
block_start <- as.integer(required_arg(args, "block_start"))
block_end <- as.integer(required_arg(args, "block_end"))
scope <- args$pair_scope %||% "fast"
test <- args$pair_test %||% "component_union"
family <- args$count_family %||% "poisson"
fasthurdle_family <- if (family == "negative_binomial") "negbin" else family
max_genes <- as.integer(args$max_cluster_genes %||% 50L)
covariate_names <- split_csv_arg(args$covariates %||% "age,sex,pc1,pc2,pc3,pc4,pc5,pc6,pf1,pf2")
if (packageVersion("fasthurdle") < package_version("1.1.1")) {
  stop("sc-pcQTL requires fasthurdle >= 1.1.1")
}
if (test == "joint_score" && family != "poisson") stop("joint_score currently requires count_family=poisson")
if (test == "joint_score" && packageVersion("fasthurdle") < package_version("1.2.0")) {
  stop("joint_score requires fasthurdle >= 1.2.0")
}

genes <- fread(file.path(stage, "gene_filtering.tsv"))[keep == TRUE & chromosome == chr]
setorder(genes, chromosome_rank)
gene_names <- genes$gene_name
covars <- readRDS(file.path(stage, "covariates.rds"))
task_first <- if (scope == "fast") block_start else max(1L, response_start - max_genes + 1L)
task_last <- if (scope == "fast") block_end else min(length(gene_names), response_end + max_genes - 1L)
counts <- load_staged_count_range(stage, chr, task_first, task_last)
if (nrow(counts) != nrow(covars)) stop("Staged count and covariate rows differ")

base_covars <- active_covariates(covars, covariate_names)
count_null_formula <- reformulate(base_covars)
zero_null_formula <- reformulate(c("log_total_read_counts", base_covars))
X_null <- model.matrix(count_null_formula, covars)
Z_null <- model.matrix(zero_null_formula, covars)
offsetx <- covars$log_total_read_counts
weights <- rep.int(1, nrow(covars))

wald_one <- function(response_name, predictor_name) {
  y <- as.numeric(counts[, response_name])
  x <- as.numeric(counts[, predictor_name])
  data <- as.data.frame(covars[, c(base_covars, "log_total_read_counts"), with = FALSE])
  data$response <- y
  data$predictor <- x
  count_rhs <- paste(c("predictor", base_covars), collapse = " + ")
  zero_rhs <- paste(c("predictor", "log_total_read_counts", base_covars), collapse = " + ")
  formula <- as.formula(sprintf("response ~ %s + offset(log_total_read_counts) | %s", count_rhs, zero_rhs))
  fit <- tryCatch(fasthurdle(formula, data = data, dist = fasthurdle_family, zero.dist = "binomial"), error = function(e) NULL)
  if (is.null(fit)) return(c(NA, NA, NA, NA, NA, NA))
  co <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  extract <- function(part) {
    table <- co[[part]]
    if (is.null(table) || !"predictor" %in% rownames(table)) return(c(NA, NA, NA))
    c(table["predictor", "Estimate"], table["predictor", "Std. Error"], table["predictor", "Pr(>|z|)"])
  }
  c(extract("count"), extract("zero"))
}

score_response <- function(response_name, predictor_names) {
  y <- as.numeric(counts[, response_name])
  predictors <- counts[, predictor_names, drop = FALSE]
  count_null <- tryCatch(fit_null_count(X_null, y, offsetx = offsetx, weights = weights, dist = "poisson"), error = function(e) NULL)
  zero_null <- tryCatch(fit_null_zero(Z_null, y, weights = weights), error = function(e) NULL)
  count_beta <- count_se <- count_stat <- rep(NA_real_, length(predictor_names))
  zero_beta <- zero_se <- zero_stat <- rep(NA_real_, length(predictor_names))
  if (!is.null(count_null) && isTRUE(count_null$convergence == 0)) {
    count_null <- tryCatch(prepare_score_cache_count(count_null, y, X_null, offsetx, weights), error = function(e) NULL)
    cache <- if (!is.null(count_null)) count_null$score_cache else NULL
    if (!is.null(cache) && isTRUE(cache$valid)) {
      batch <- getFromNamespace("score_test_count_batch_cpp", "fasthurdle")
      res <- tryCatch(batch(
        predictors[cache$Y1 + 1L, , drop = FALSE], cache$Y1, cache$grad_weights,
        cache$v_ee, cache$Y_pos, cache$I_nn_inv, cache$I_nn_beta_inv,
        cache$beta_inv_ok, cache$Xnull_vee_t, cache$X_null_pos, cache$w_pos,
        cache$theta, cache$beta_null, cache$eta_null_pos, cache$mu_pos,
        cache$p0_pos, cache$log_p1_pos, cache$kx_null, FALSE, FALSE, 1e30, NULL
      ), error = function(e) NULL)
      if (!is.null(res)) {
        count_beta <- as.numeric(res$beta); count_se <- as.numeric(res$se)
        count_stat <- as.numeric(res$statistic)
      }
    }
  }
  if (!is.null(zero_null) && isTRUE(zero_null$convergence == 0)) {
    zero_null <- tryCatch(prepare_score_cache_zero(zero_null, y, Z_null, weights), error = function(e) NULL)
    cache <- if (!is.null(zero_null)) zero_null$score_cache else NULL
    if (!is.null(cache) && isTRUE(cache$valid)) {
      for (j in seq_along(predictor_names)) {
        res <- tryCatch(score_test_zero(Z_null, predictors[, j], y, weights = weights,
                                        null_fit_zero = zero_null, spa_cutoff = NULL), error = function(e) NULL)
        if (!is.null(res)) {
          zero_beta[j] <- res$beta[[1L]]; zero_se[j] <- res$se[[1L]]
          zero_stat[j] <- res$statistic[[1L]]
        }
      }
    }
  }
  stat_joint <- count_stat + zero_stat
  log_p_joint <- pchisq(stat_joint, 2, lower.tail = FALSE, log.p = TRUE)
  data.table(beta_count = count_beta, se_count = count_se,
             p_count = pchisq(count_stat, 1, lower.tail = FALSE),
             beta_detection = zero_beta, se_detection = zero_se,
             p_detection = pchisq(zero_stat, 1, lower.tail = FALSE),
             stat_count = count_stat, stat_detection = zero_stat,
             stat_joint = stat_joint,
             log_p_joint = log_p_joint,
             neglog10_p_joint = -log_p_joint / log(10),
             p_joint = exp(log_p_joint))
}

rows <- list()
k <- 0L
for (response_idx in seq.int(response_start, response_end)) {
  predictor_idx <- if (scope == "fast") {
    setdiff(seq.int(block_start, block_end), response_idx)
  } else {
    setdiff(seq.int(max(1L, response_idx - max_genes + 1L),
                    min(length(gene_names), response_idx + max_genes - 1L)), response_idx)
  }
  if (!length(predictor_idx)) next
  base <- data.table(
    celltype = celltype, chromosome = chr, response = gene_names[response_idx],
    predictor = gene_names[predictor_idx], response_rank = response_idx,
    predictor_rank = predictor_idx
  )
  if (test == "component_union") {
    values <- t(vapply(gene_names[predictor_idx], function(name) {
      wald_one(gene_names[response_idx], name)
    }, numeric(6L)))
    stat_count <- (values[, 1L] / values[, 2L])^2
    stat_detection <- (values[, 4L] / values[, 5L])^2
    result <- data.table(beta_count = values[, 1L], se_count = values[, 2L], p_count = values[, 3L],
                         beta_detection = values[, 4L], se_detection = values[, 5L], p_detection = values[, 6L],
                         stat_count = stat_count, stat_detection = stat_detection,
                         stat_joint = stat_count + stat_detection,
                         log_p_joint = NA_real_, neglog10_p_joint = NA_real_, p_joint = NA_real_)
  } else {
    result <- score_response(gene_names[response_idx], gene_names[predictor_idx])
  }
  k <- k + 1L
  rows[[k]] <- cbind(base, result)
}
result <- if (length(rows)) rbindlist(rows, fill = TRUE) else empty_pair_table()
result[, test_status := ifelse(is.finite(p_count) | is.finite(p_detection) | is.finite(p_joint), "ok", "fit_failed")]
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
fwrite(result, out, sep = "\t", compress = "gzip")
