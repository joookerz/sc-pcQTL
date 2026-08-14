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

evenly_spaced <- function(x, n) {
  x <- unique(x)
  if (length(x) < n) stop("Only ", length(x), " variants are available; need ", n)
  if (length(x) == n) return(x)
  index <- unique(as.integer(round(seq.int(1L, length(x), length.out = n))))
  if (length(index) < n) {
    index <- c(index, setdiff(seq_along(x), index)[seq_len(n - length(index))])
  }
  x[sort(index[seq_len(n)])]
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
mode <- required(args, "mode")
bim <- read_table(required(args, "bim"), header = FALSE)
names(bim) <- c("chromosome", "snp_id", "cm", "position", "allele1", "allele2")
bim$position <- as.integer(bim$position)

if (mode == "background") {
  prune_ids <- readLines(required(args, "prune"), warn = FALSE)
  n <- as.integer(required(args, "n"))
  candidates <- bim[bim$snp_id %in% prune_ids, , drop = FALSE]
  candidates <- candidates[order(candidates$position), , drop = FALSE]
  writeLines(evenly_spaced(candidates$snp_id, n), required(args, "output"))
  quit(status = 0L)
}

if (mode != "signals") stop("--mode must be background or signals")

freq <- read_table(required(args, "freq"))
missing <- read_table(required(args, "missing"))
regions <- read_table(required(args, "regions"))
names(freq) <- sub("^#", "", names(freq))
names(missing) <- sub("^#", "", names(missing))
if (!all(c("SNP", "MAF") %in% names(freq))) stop("Unexpected PLINK .frq columns")
if (!all(c("SNP", "F_MISS") %in% names(missing))) stop("Unexpected PLINK .lmiss columns")

n_per_locus <- as.integer(required(args, "n-per-locus"))
near_n <- min(as.integer(args[["near-n"]] %||% 200L), n_per_locus - 1L)
freq_small <- data.frame(snp_id = freq$SNP, maf = as.numeric(freq$MAF))
missing_small <- data.frame(snp_id = missing$SNP, call_rate = 1 - as.numeric(missing$F_MISS))
variants <- merge(merge(bim, freq_small, by = "snp_id"), missing_small, by = "snp_id")
allele_pair <- paste(pmin(variants$allele1, variants$allele2),
                     pmax(variants$allele1, variants$allele2), sep = "/")
variants$palindromic <- allele_pair %in% c("A/T", "C/G")

selected_rows <- vector("list", nrow(regions))
causal_rows <- vector("list", nrow(regions))
for (i in seq_len(nrow(regions))) {
  region <- regions[i, ]
  pool <- variants[variants$position >= region$start & variants$position <= region$end, , drop = FALSE]
  pool <- pool[order(pool$position), , drop = FALSE]
  if (nrow(pool) < n_per_locus) {
    stop(region$module, " contains only ", nrow(pool), " eligible variants; need ", n_per_locus)
  }
  causal_pool <- pool[
    is.finite(pool$maf) & pool$maf >= 0.20 & pool$maf <= 0.40 &
      is.finite(pool$call_rate) & pool$call_rate >= 0.99 & !pool$palindromic,
    , drop = FALSE
  ]
  if (!nrow(causal_pool)) stop("No eligible causal SNP in ", region$module)
  order_index <- order(
    abs(causal_pool$position - region$center),
    abs(causal_pool$maf - 0.30),
    causal_pool$snp_id
  )
  causal <- causal_pool[order_index[1L], , drop = FALSE]

  nearest <- pool[order(abs(pool$position - causal$position), pool$snp_id), , drop = FALSE]
  near_ids <- unique(c(causal$snp_id, nearest$snp_id))[seq_len(near_n + 1L)]
  remaining <- pool[!pool$snp_id %in% near_ids, , drop = FALSE]
  remaining <- remaining[order(remaining$position), , drop = FALSE]
  spread_ids <- evenly_spaced(remaining$snp_id, n_per_locus - length(near_ids))
  selected_rows[[i]] <- data.frame(module = region$module, snp_id = c(near_ids, spread_ids))
  if (isTRUE(as.integer(region$genetic) == 1L)) {
    causal_rows[[i]] <- data.frame(
      module = region$module,
      snp_id = causal$snp_id,
      chromosome = causal$chromosome,
      position = causal$position,
      allele1 = causal$allele1,
      allele2 = causal$allele2,
      maf = causal$maf,
      call_rate = causal$call_rate
    )
  }
}

selected <- do.call(rbind, selected_rows)
causal <- do.call(rbind, Filter(Negate(is.null), causal_rows))
writeLines(unique(selected$snp_id), required(args, "output"))
write.table(causal, required(args, "causal-output"), sep = "\t", quote = FALSE,
            row.names = FALSE)
