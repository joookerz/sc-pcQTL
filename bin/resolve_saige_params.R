#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))[1L]])))
source(file.path(script_dir, "scpcqtl_utils.R"))
args <- parse_cli()
defaults_file <- required_arg(args, "defaults")
out <- required_arg(args, "out")
user_file <- args$user %||% ""
covariates <- args$covariates %||% "age,sex,pc1,pc2,pc3,pc4,pc5,pc6,pf1,pf2"
qtl_maf <- as.numeric(args$qtl_maf %||% 0.05)
allow_unknown <- as_flag(args$allow_unknown %||% "false")

defaults <- fread(defaults_file, colClasses = "character")
required <- c("step", "parameter", "value")
if (!identical(names(defaults), required)) stop("SAIGE defaults must have columns: ", paste(required, collapse = ","))
resolved <- copy(defaults)

set_value <- function(step, parameter, new_value) {
  hit <- resolved$step == step & resolved$parameter == parameter
  if (any(hit)) resolved[hit, value := as.character(new_value)]
  else resolved <<- rbind(resolved, data.table(step = step, parameter = parameter, value = as.character(new_value)))
}
set_value("step1", "covarColList", covariates)
set_value("step1", "sampleCovarColList", covariates)
set_value("step2", "minMAF", qtl_maf)

protected <- c(
  "phenoFile", "phenoCol", "sampleIDColinphenoFile", "outputPrefix", "plinkFile",
  "bedFile", "bimFile", "famFile", "SAIGEOutputFile", "chrom",
  "GMMATmodelFile", "varianceRatioFile", "rangestoIncludeFile",
  "assocFile", "geneName", "genePval_outputFile", "traitType",
  "covarColList", "sampleCovarColList", "minMAF"
)
known <- unique(c(defaults$parameter, "covarColList", "sampleCovarColList",
  "nThreads", "maxiter", "memoryChunk", "minMAC", "minInfo", "SPAcutoff",
  "IsOutputAFinCaseCtrl", "IsOutputNinCaseCtrl", "IsDropMissingDosages",
  "isCovariateTransform", "markers_per_chunk", "numLinesOutput"))

if (nzchar(user_file)) {
  user <- fread(user_file, colClasses = "character")
  if (!identical(names(user), required)) stop("User SAIGE table must have columns: step,parameter,value")
  if (anyDuplicated(user[, .(step, parameter)])) stop("User SAIGE table contains duplicate step/parameter rows")
  if (any(!nzchar(user$parameter)) || any(!nzchar(user$value))) stop("SAIGE parameter names and values must be non-empty")
  if (any(!user$step %in% c("step1", "step2", "step3"))) stop("SAIGE step must be step1, step2, or step3")
  bad_protected <- intersect(user$parameter, protected)
  if (length(bad_protected)) stop("Workflow-owned SAIGE parameters cannot be overridden: ", paste(bad_protected, collapse = ", "))
  unknown <- setdiff(user$parameter, known)
  if (length(unknown) && !allow_unknown) {
    stop("Unknown SAIGE-QTL parameters: ", paste(unknown, collapse = ", "),
         ". Use --allow_unknown_saige_params to pass parameters for another compatible image.")
  }
  for (i in seq_len(nrow(user))) set_value(user$step[i], user$parameter[i], user$value[i])
}
if (anyDuplicated(resolved[, .(step, parameter)])) stop("Duplicate resolved SAIGE parameters")
setorder(resolved, step, parameter)
fwrite(resolved, out, sep = "\t")
