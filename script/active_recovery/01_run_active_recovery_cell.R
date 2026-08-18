#!/usr/bin/env Rscript
# Usage: Rscript script/active_recovery/01_run_active_recovery_cell.R n lambda family replicate out_dir [epochs]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L) stop("expected: n lambda family replicate out_dir [epochs]", call. = FALSE)
if (!file.exists("DESCRIPTION")) {
  stop("run this receipt script from the pigauto repository root", call. = FALSE)
}
devtools::load_all(quiet = TRUE)
source("script/active_recovery/00_prepare_active_recovery.R")
n <- as.integer(args[[1L]]); lambda <- as.numeric(args[[2L]])
family <- args[[3L]]; replicate <- as.integer(args[[4L]]); out_dir <- args[[5L]]
epochs <- if (length(args) >= 6L) as.integer(args[[6L]]) else 100L
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, sprintf("active-recovery-n%s-lambda%s-%s-rep%05d.rds", n, lambda, family, replicate))
if (file.exists(out)) stop("receipt already exists; verify SHA/configuration before reuse: ", out, call. = FALSE)
receipt <- tryCatch(
  active_recovery_one(n, lambda, family, replicate, master_seed = 20260818L, epochs = epochs),
  error = function(e) data.frame(n = n, lambda = lambda, family = family,
    replicate = replicate, policy = NA_character_, source_sha = active_recovery_source_sha(),
    status = "error", error = conditionMessage(e), stringsAsFactors = FALSE)
)
saveRDS(receipt, out)
message(out)
