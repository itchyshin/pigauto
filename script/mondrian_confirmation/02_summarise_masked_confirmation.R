#!/usr/bin/env Rscript
# Usage: Rscript 02_summarise_masked_confirmation.R run_dir summary.rds
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("expected: run_dir summary.rds", call. = FALSE)
x <- lapply(c("split", "mondrian"), function(m) readRDS(file.path(args[[1L]], paste0(m, ".rds"))))
names(x) <- c("split", "mondrian")
if (!all(vapply(x, function(z) identical(z$status, "ok"), logical(1)))) stop("retained method error", call. = FALSE)
tab <- do.call(rbind, Map(function(z, method) transform(z$metrics, method = method), x, names(x)))
saveRDS(list(methods = x, metrics = tab), args[[2L]])
print(tab)
