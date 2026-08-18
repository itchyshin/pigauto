#!/usr/bin/env Rscript
# Usage: Rscript 05_summarise_active_recovery_combined.R pilot_dir extension_dir summary.rds
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("expected: pilot_dir extension_dir summary.rds", call. = FALSE)
read_receipts <- function(dir) {
  files <- list.files(file.path(dir, "results"), pattern = "\\.rds$", full.names = TRUE)
  if (!length(files)) stop("no receipts in ", dir, call. = FALSE)
  do.call(rbind, lapply(files, readRDS))
}
raw <- rbind(read_receipts(args[[1L]]), read_receipts(args[[2L]]))
if (nrow(raw) != 32000L) stop("expected 32,000 policy rows", call. = FALSE)
if (!all(raw$status == "ok")) stop("combined summary retains errored receipts", call. = FALSE)
key <- c("n", "lambda", "family", "replicate")
wide <- reshape(raw[, c(key, "policy", "primary")], idvar = key, timevar = "policy", direction = "wide")
if (nrow(wide) != 8000L) stop("expected 8,000 complete replicates", call. = FALSE)
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
saveRDS(list(raw = raw, paired = paired, summary = summary, failure = failure), args[[3L]])
print(summary)
print(failure)
