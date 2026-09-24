# script/rubin_convcheck.R
#
# Why does BACE's own convergence check fail 70 to 90% of pre-run fits while ESS is 750 to 1,600?
# BACE's check (installed build: assess_convergence(method = "summary") -> .assess_summary_convergence) runs on
# the convergence phase, whose fills are deterministic (sample = FALSE). For each trait with missing data it
# takes the per-iteration MEAN of the imputed values and applies three criteria (Geweke is disabled, NA):
#   acf   : |lag-1 autocorrelation of the series| < 0.3
#   trend : cor.test(iteration, series) p > 0.05
#   pct   : mean |percent change| over the last half of the series < 5%   (divides by the series value)
# A trait passes if more than half of its non-NA criteria pass; BACE declares convergence if more than 70% of
# traits pass. This script recomputes each criterion from the summary_stats stored in every pre-run rds and
# reports pass rates per trait and criterion.
#
#   Rscript script/rubin_convcheck.R --dir <prerun dir>

args <- commandArgs(trailingOnly = TRUE)
dir <- args[match("--dir", args) + 1L]
fs <- list.files(dir, "^rubin_.*\\.rds$", recursive = TRUE, full.names = TRUE)
rows <- list()
for (f in fs) {
  x <- readRDS(f); ss <- x$diag$bace$summary_stats
  if (!is.data.frame(ss) || nrow(ss) < 3L) next
  for (v in setdiff(names(ss), "iteration")) {
    s <- as.numeric(ss[[v]]); n_it <- length(s)
    zero_var <- is.na(stats::var(s)) || stats::var(s) < 1e-10
    acf_ok <- if (zero_var) TRUE else { a <- suppressWarnings(stats::cor(s[-n_it], s[-1])); if (is.finite(a)) abs(a) < 0.3 else NA }
    trend_ok <- if (zero_var) TRUE else { p <- tryCatch(suppressWarnings(stats::cor.test(seq_len(n_it), s)$p.value), error = function(e) NA)
                                          if (is.na(p)) NA else p > 0.05 }
    pc <- abs(diff(s) / (s[-n_it] + 1e-10)) * 100
    pc_ok <- mean(pc[ceiling(length(pc) / 2):length(pc)], na.rm = TRUE) < 5
    crit <- c(acf_ok, trend_ok, pc_ok)
    rows[[length(rows) + 1L]] <- data.frame(n = x$n, runs = x$bace[["runs"]], nitt = x$bace[["nitt"]],
      trait = v, n_iter = n_it, mean_abs_level = mean(abs(s)), acf = acf_ok, trend = trend_ok, pct = pc_ok,
      trait_pass = mean(crit, na.rm = TRUE) > 0.5, bace_converged = isTRUE(x$diag$bace$converged))
  }
}
d <- do.call(rbind, rows)
cat(sprintf("%d fits, %d trait series\n\n", length(fs), nrow(d)))
cat("Pass rate per trait and criterion (all fits):\n")
agg <- aggregate(cbind(acf, trend, pct, trait_pass, mean_abs_level) ~ trait, data = d, FUN = mean, na.action = na.pass)
print(agg, digits = 3, row.names = FALSE)
cat("\nTrait pass rate by runs (does more runs help?):\n")
print(stats::xtabs(trait_pass ~ runs + trait, data = aggregate(trait_pass ~ runs + trait, data = d, FUN = mean)), digits = 2)
