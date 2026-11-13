#!/usr/bin/env Rscript

source("script/validation_20261018/roster.R")
args <- commandArgs(trailingOnly = TRUE)
receipt_dir <- if (length(args)) args[[1L]] else
  "docs/dev-log/validation-20261018/receipts"
jobs <- read.csv(gzfile(file.path(receipt_dir, "jobs.csv.gz")), stringsAsFactors = FALSE)
intervals <- read.csv(gzfile(file.path(receipt_dir,
  "production-posterior_intervals_by_replicate.csv.gz")), stringsAsFactors = FALSE)
posterior_jobs <- jobs[tolower(jobs$adapter) %in% c("posterior", "mi_gls", "posterior_mi"), , drop = FALSE]
if (!nrow(posterior_jobs)) stop("no planned posterior jobs found")
if (!all(c("regime_id", "rep") %in% names(posterior_jobs)) ||
    anyNA(posterior_jobs$regime_id) || anyNA(posterior_jobs$rep))
  stop("posterior jobs lack planned regime and replicate identities")

suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
suppressPackageStartupMessages(library(ape))
dgp_env <- new.env(parent = globalenv())
sys.source("script/mi_gls/regimes.R", envir = dgp_env)
expressions <- parse("script/mi_gls/dgp_v2.R")
if (length(expressions) < 2L) stop("unexpected dgp_v2.R structure")
for (expr in as.list(expressions)[-1L]) eval(expr, envir = dgp_env)

planned <- unique(posterior_jobs[, c("job_id", "regime_id", "rep")])
if (anyDuplicated(planned$job_id)) stop("duplicate posterior job IDs in the plan")
expected <- vector("list", nrow(planned) * 2L)
k <- 0L
for (i in seq_len(nrow(planned))) {
  cell <- dgp_env$simulate_regime_cell(as.integer(planned$regime_id[[i]]),
                                       as.integer(planned$rep[[i]]))
  for (trait in c("x", "y")) {
    k <- k + 1L
    mask <- cell[[paste0("mask_", trait)]]
    expected[[k]] <- data.frame(job_id = planned$job_id[[i]], trait = trait,
      n_expected_truth = sum(mask), stringsAsFactors = FALSE)
  }
}
expected <- do.call(rbind, expected)
audit <- validation_compare_posterior_denominators(intervals, expected)
if (any(audit$status != "match")) {
  print(utils::head(audit[audit$status != "match", , drop = FALSE], 20L))
  stop("posterior interval denominators do not match the independently regenerated masks")
}
tampered <- intervals
first_key <- match(paste(expected$job_id[[1L]], expected$trait[[1L]], sep = "\034"),
                   paste(tampered$job_id, tampered$trait, sep = "\034"))
tampered$n_expected[[first_key]] <- tampered$n_expected[[first_key]] + 1L
negative_control <- validation_compare_posterior_denominators(tampered, expected)
if (!any(negative_control$status == "mismatch"))
  stop("denominator audit failed its deliberately corrupted negative control")

out <- if (length(args) >= 2L) args[[2L]] else
  file.path(receipt_dir, "production-posterior-denominator-audit.csv.gz")
con <- gzfile(out, "wt")
utils::write.csv(audit, con, row.names = FALSE, quote = FALSE)
close(con)
print(table(audit$status))
cat("POSTERIOR_DENOMINATORS_VERIFIED\n")
