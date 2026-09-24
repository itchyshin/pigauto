#!/usr/bin/env Rscript
# script/mi_gls/04_acceptance.R
#
# G6 (.unlazy/mi-posterior/GATES.md): simulation acceptance gate for
# multi_impute(draws_method = "posterior"), method = "posterior_full",
# against the "complete" reference, using the summary CSV written by
# script/mi_gls/03_summarise_v2.R.
#
# Rules (every regime x downstream unless noted):
#   1. |paired bias| <= max(0.02, 2.5 * paired_bias_mcse)
#   2. per-cell: coverage >= complete_coverage - 0.05; mean shortfall <= 0.02
#      across all regime x downstream rows
#   3. SE ratio in [0.90, 1.15] under the lambda (phylolm) analysis
#   4. mean SE ratio of posterior_full > posterior_none over the
#      both-missing (missing == "both") regimes
#   5. fit failures <= 2%
#   6. (design review B3) non-converged fits <= 2% of a regime's posterior_full
#      fits; any regime over that threshold FAILS the gate and prints a loud
#      "NONCONVERGED regime <id>: <k>/<n>" line. Rules 1-5 above are already
#      computed from converged reps only (03_summarise_v2.R), so rule 6 is
#      the loud, blocking surface for that filtering.
#
# Usage:
#   Rscript script/mi_gls/04_acceptance.R <summary.csv>
#   Rscript script/mi_gls/04_acceptance.R --selftest

args <- commandArgs(trailingOnly = TRUE)

source(file.path("script", "mi_gls", "regimes.R"))

NONCONVERGED_THRESHOLD <- 0.02

check_nonconverged <- function(df) {
  fails <- character(0)
  pf <- unique(df[df$method == "posterior_full",
                  c("regime_id", "n_fits", "n_converged")])
  pf <- pf[is.finite(pf$n_fits) & pf$n_fits > 0, ]
  for (i in seq_len(nrow(pf))) {
    n <- pf$n_fits[i]; k <- pf$n_converged[i]
    nonconv_frac <- 1 - k / n
    if (is.finite(nonconv_frac) && nonconv_frac > NONCONVERGED_THRESHOLD) {
      cat(sprintf("NONCONVERGED regime %d: %d/%d non-converged (> %.0f%%)\n",
                 pf$regime_id[i], n - k, n, 100 * NONCONVERGED_THRESHOLD))
      fails <- c(fails, sprintf("regime %d: non-converged fraction %.4f > %.2f",
                                pf$regime_id[i], nonconv_frac, NONCONVERGED_THRESHOLD))
    }
  }
  fails
}

run_gate <- function(df) {
  fails <- character(0)

  df <- merge(df, regimes[, c("regime_id", "missing")], by = "regime_id", all.x = TRUE)

  pf <- df[df$method == "posterior_full", ]
  pn <- df[df$method == "posterior_none", ]

  # Rule 6 first (loud, blocking) so NONCONVERGED prints regardless of
  # whether the other rules also fail.
  fails <- c(fails, check_nonconverged(df))

  # Rule 1: |paired bias| <= max(0.02, 2.5 * MCSE).
  bad1 <- pf[is.finite(pf$paired_bias) & is.finite(pf$paired_bias_mcse) &
            abs(pf$paired_bias) > pmax(0.02, 2.5 * pf$paired_bias_mcse), ]
  if (nrow(bad1)) fails <- c(fails, sprintf(
    "bias: regime %d %s |bias|=%.4f > max(0.02, 2.5*MCSE=%.4f)",
    bad1$regime_id, bad1$downstream, abs(bad1$paired_bias), bad1$paired_bias_mcse))

  # Rule 2: per-cell coverage shortfall <= 0.05; mean shortfall <= 0.02.
  shortfall <- pf$complete_coverage - pf$coverage
  bad2 <- pf[is.finite(shortfall) & shortfall > 0.05, ]
  if (nrow(bad2)) fails <- c(fails, sprintf(
    "coverage: regime %d %s shortfall=%.4f > 0.05",
    bad2$regime_id, bad2$downstream, (bad2$complete_coverage - bad2$coverage)))
  mean_shortfall <- mean(shortfall, na.rm = TRUE)
  if (is.finite(mean_shortfall) && mean_shortfall > 0.02) fails <- c(fails, sprintf(
    "coverage: mean shortfall=%.4f > 0.02", mean_shortfall))

  # Rule 3: SE ratio in [0.90, 1.15] under phylolm (lambda) analysis.
  r3 <- pf[pf$downstream == "phylolm", ]
  bad3 <- r3[is.finite(r3$se_ratio) & (r3$se_ratio < 0.90 | r3$se_ratio > 1.15), ]
  if (nrow(bad3)) fails <- c(fails, sprintf(
    "se_ratio: regime %d phylolm se_ratio=%.4f outside [0.90, 1.15]",
    bad3$regime_id, bad3$se_ratio))

  # Rule 4: mean SE ratio of posterior_full > posterior_none, both-missing regimes.
  both_pf <- pf[pf$missing == "both" & is.finite(pf$se_ratio), ]
  both_pn <- pn[pn$missing == "both" & is.finite(pn$se_ratio), ]
  if (nrow(both_pf) && nrow(both_pn)) {
    m_pf <- mean(both_pf$se_ratio); m_pn <- mean(both_pn$se_ratio)
    if (!(m_pf > m_pn)) fails <- c(fails, sprintf(
      "se_ratio ordering: mean posterior_full=%.4f not > posterior_none=%.4f over both-missing regimes",
      m_pf, m_pn))
  } else {
    fails <- c(fails, "se_ratio ordering: insufficient both-missing rows for posterior_full/posterior_none")
  }

  # Rule 5: fit failures <= 2%.
  bad5 <- pf[is.finite(pf$fit_failure_rate) & pf$fit_failure_rate > 0.02, ]
  if (nrow(bad5)) fails <- c(fails, sprintf(
    "fit_failure: regime %d %s failure_rate=%.4f > 0.02",
    bad5$regime_id, bad5$downstream, bad5$fit_failure_rate))

  fails
}

if (length(args) >= 1L && identical(args[[1L]], "--selftest")) {
  base <- data.frame(
    regime_id = rep(c(17, 21), each = 2), method = "posterior_full",
    downstream = rep(c("gls", "phylolm"), 2),
    paired_bias = 0.005, paired_bias_mcse = 0.01,
    coverage = 0.94, complete_coverage = 0.95,
    se_ratio = 1.0, fit_failure_rate = 0.0,
    n_fits = 200L, n_converged = 200L
  )
  base_pn <- base; base_pn$method <- "posterior_none"; base_pn$se_ratio <- 0.8
  base_pn$n_fits <- NA_integer_; base_pn$n_converged <- NA_integer_
  pass_df <- rbind(base, base_pn)

  fail_bias_df <- pass_df
  fail_bias_df$paired_bias[fail_bias_df$method == "posterior_full" &
                           fail_bias_df$regime_id == 17] <- 0.5

  fail_nonconv_df <- pass_df
  fail_nonconv_df$n_converged[fail_nonconv_df$method == "posterior_full" &
                              fail_nonconv_df$regime_id == 17] <- 190L  # 10/200 = 5% > 2%

  pass_fails         <- run_gate(pass_df)
  fail_bias_fails     <- run_gate(fail_bias_df)
  fail_nonconv_fails  <- run_gate(fail_nonconv_df)

  ok <- length(pass_fails) == 0L && length(fail_bias_fails) > 0L &&
    length(fail_nonconv_fails) > 0L &&
    any(grepl("non-converged", fail_nonconv_fails))
  cat("pass-fixture failures:", length(pass_fails), "\n")
  cat("fail-bias-fixture failures:", length(fail_bias_fails), "\n")
  cat("fail-nonconverged-fixture failures:", length(fail_nonconv_fails), "\n")
  if (ok) {
    cat("SELFTEST_OK\n")
  } else {
    cat("SELFTEST FAILED\n")
    if (length(pass_fails)) cat(" pass fixture unexpectedly failed:\n",
                                paste(" -", pass_fails, collapse = "\n"), "\n")
    if (!length(fail_bias_fails)) cat(" fail-bias fixture unexpectedly passed\n")
    if (!length(fail_nonconv_fails) || !any(grepl("non-converged", fail_nonconv_fails))) {
      cat(" fail-nonconverged fixture did not trigger the non-converged rule\n")
    }
  }
  quit(save = "no", status = if (ok) 0L else 1L)
}

if (length(args) < 1L) stop("expected: <summary.csv> or --selftest", call. = FALSE)
summary_csv <- args[[1L]]
df <- utils::read.csv(summary_csv, stringsAsFactors = FALSE)

fails <- run_gate(df)
if (length(fails)) {
  cat("G6 FAILURES:\n")
  for (f in fails) cat(" -", f, "\n")
  quit(save = "no", status = 1L)
} else {
  cat("SIM_ACCEPT_PASS\n")
}
