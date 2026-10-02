#!/usr/bin/env Rscript
# script/mi_gls/02_summarise.R
#
# Summarise script/mi_gls/01_cell.R output into a bias / empirical-SD /
# SE-ratio / coverage / MCSE table per regime x method x downstream model.
# Usage:
#   Rscript script/mi_gls/02_summarise.R <indir> <out.md>
#
# <indir> holds regime_<id>_rep_<rep>.rds files from 01_cell.R.

args   <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("expected: indir out.md", call. = FALSE)
indir  <- args[[1L]]
out_md <- args[[2L]]

source(file.path("script", "mi_gls", "regimes.R"))

files <- list.files(indir, pattern = "^regime_[0-9]+_rep_[0-9]+\\.rds$",
                    full.names = TRUE)
if (length(files) == 0L) stop("no .rds files found in ", indir, call. = FALSE)

rows <- lapply(files, function(f) {
  r <- readRDS(f)
  d <- r$results
  d$regime_id <- r$regime_id
  d$rep       <- r$rep
  d$true_beta <- r$true_beta
  d
})
all_df <- do.call(rbind, rows)

# One row per regime x method x downstream: bias, empirical SD, mean
# reported SE, SE ratio (mean SE / empirical SD -- < 1 means the method
# UNDERSTATES uncertainty), coverage of the nominal-95% CI, and the
# Monte Carlo SE of that coverage estimate (binomial SE over the number of
# non-NA reps, sqrt(p(1-p)/R)).
summarise_cell <- function(sub) {
  est <- sub$estimate
  ok  <- is.finite(est)
  R   <- sum(ok)
  bias    <- if (R > 0) mean(est[ok]) - sub$true_beta[1L] else NA_real_
  emp_sd  <- if (R > 1) stats::sd(est[ok]) else NA_real_
  mean_se <- if (any(is.finite(sub$se))) mean(sub$se, na.rm = TRUE) else NA_real_
  se_ratio <- if (is.finite(mean_se) && is.finite(emp_sd) && emp_sd > 0) {
    mean_se / emp_sd
  } else NA_real_
  cov_ok <- is.finite(sub$covered)
  Rc     <- sum(cov_ok)
  coverage <- if (Rc > 0) mean(sub$covered[cov_ok]) else NA_real_
  mcse     <- if (Rc > 0 && is.finite(coverage)) {
    sqrt(coverage * (1 - coverage) / Rc)
  } else NA_real_
  data.frame(R = R, bias = bias, emp_sd = emp_sd, mean_se = mean_se,
            se_ratio = se_ratio, coverage = coverage, coverage_mcse = mcse)
}

keys <- unique(all_df[, c("regime_id", "method", "downstream")])
summary_rows <- lapply(seq_len(nrow(keys)), function(i) {
  k <- keys[i, ]
  sub <- all_df[all_df$regime_id == k$regime_id & all_df$method == k$method &
                 all_df$downstream == k$downstream, ]
  cbind(k, summarise_cell(sub))
})
summary_df <- do.call(rbind, summary_rows)
summary_df <- merge(summary_df, regimes, by = "regime_id", all.x = TRUE)
summary_df <- summary_df[order(summary_df$regime_id, summary_df$method,
                               summary_df$downstream), ]

fmt <- function(x, digits = 3) ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))

md <- c(
  "# MI-GLS attenuation campaign summary",
  "",
  sprintf("Generated from %d rep file(s) in `%s`.", length(files), indir),
  "",
  "| regime | lambda | n | mechanism | missing | method | downstream | R | bias | emp_sd | mean_se | se_ratio | coverage | coverage_mcse |",
  "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"
)
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  md <- c(md, sprintf(
    "| %d | %s | %d | %s | %s | %s | %s | %d | %s | %s | %s | %s | %s | %s |",
    r$regime_id, r$lambda, r$n, r$mechanism, r$missing, r$method, r$downstream,
    r$R, fmt(r$bias), fmt(r$emp_sd), fmt(r$mean_se), fmt(r$se_ratio),
    fmt(r$coverage), fmt(r$coverage_mcse)
  ))
}

dir.create(dirname(out_md), recursive = TRUE, showWarnings = FALSE)
writeLines(md, out_md)
cat("wrote", out_md, "(", nrow(summary_df), "rows )\n")
