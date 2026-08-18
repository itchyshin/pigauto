#!/usr/bin/env Rscript
# Usage: Rscript script/active_recovery/02_summarise_active_recovery.R results_dir summary.rds
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("expected: results_dir summary.rds", call. = FALSE)
files <- list.files(args[[1L]], pattern = "\\.rds$", full.names = TRUE)
if (!length(files)) stop("no per-replicate receipts found", call. = FALSE)
raw <- do.call(rbind, lapply(files, readRDS))
ok <- raw[raw$status == "ok", , drop = FALSE]
key <- c("n", "lambda", "family", "replicate")
wide <- reshape(ok[, c(key, "policy", "primary")], idvar = key, timevar = "policy", direction = "wide")
paired <- transform(wide,
  active_minus_random = primary.active - primary.random,
  active_minus_uncertainty = primary.active - primary.uncertainty)
summarise <- function(d, contrast) {
  x <- d[[contrast]]
  data.frame(n = d$n[[1L]], lambda = d$lambda[[1L]], family = d$family[[1L]],
    contrast = contrast, replicates = sum(is.finite(x)), mean = mean(x, na.rm = TRUE),
    mcse = stats::sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))),
    ci_low = mean(x, na.rm = TRUE) - 1.96 * stats::sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))),
    ci_high = mean(x, na.rm = TRUE) + 1.96 * stats::sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x))))
}
groups <- split(paired, interaction(paired$n, paired$lambda, paired$family, drop = TRUE))
summary <- do.call(rbind, unlist(lapply(groups, function(d) list(
  summarise(d, "active_minus_random"), summarise(d, "active_minus_uncertainty"))), recursive = FALSE))
failure <- aggregate(I(status != "ok") ~ n + lambda + family + policy, raw, mean)
names(failure)[ncol(failure)] <- "failure_rate"
saveRDS(list(raw = raw, paired = paired, summary = summary, failure = failure), args[[2L]])
print(summary)
print(failure)
