#!/usr/bin/env Rscript
# script/mi_gls/05_cell_coverage.R
#
# G7 (.unlazy/mi-posterior/GATES.md): per-cell 95% predictive coverage gate
# for multi_impute(draws_method = "posterior"), method = "posterior_full",
# using the cell-coverage CSV written by script/mi_gls/03_summarise_v2.R
# (its 4th argument). Columns expected: regime_id, method, trait,
# mechanism, n, covered_sum, mean_width, n_fits, n_converged,
# median_max_rhat, median_min_ess.
#
# MCAR (exchangeable) masks are GATED to per-regime x trait coverage in
# [0.92, 0.98]. MAR_phylo (clade-biased) masks are reported beside the
# MCAR numbers but never gated (design.md section 5 / design review item
# B3: clade-biased coverage stays descriptive).
#
# Design review item B3: gate only on converged fits (03_summarise_v2.R
# already restricts covered_sum/n to converged reps for posterior_full);
# a regime with > 2% non-converged posterior_full fits FAILS this gate and
# prints a loud "NONCONVERGED regime <id>: <k>/<n>" line, regardless of
# whether its coverage numbers happen to look fine.
#
# Usage:
#   Rscript script/mi_gls/05_cell_coverage.R <cell_coverage.csv>
#   Rscript script/mi_gls/05_cell_coverage.R --selftest

args <- commandArgs(trailingOnly = TRUE)

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
  fails <- c(check_nonconverged(df))
  sub <- df[df$method == "posterior_full", ]

  mcar <- sub[sub$mechanism == "MCAR" & is.finite(sub$n) & sub$n > 0, ]
  mcar$coverage <- mcar$covered_sum / mcar$n
  bad <- mcar[mcar$coverage < 0.92 | mcar$coverage > 0.98, ]
  if (nrow(bad)) fails <- c(fails, sprintf(
    "regime %d trait %s: MCAR coverage=%.4f (n=%d) outside [0.92, 0.98]",
    bad$regime_id, bad$trait, bad$coverage, bad$n))

  attr(fails, "mcar_report") <- mcar
  attr(fails, "mar_phylo_report") <- {
    mar <- sub[sub$mechanism == "MAR_phylo" & is.finite(sub$n) & sub$n > 0, ]
    mar$coverage <- mar$covered_sum / mar$n
    mar
  }
  fails
}

report_side_by_side <- function(fails) {
  mcar <- attr(fails, "mcar_report")
  mar  <- attr(fails, "mar_phylo_report")
  cat("Per-cell coverage by mechanism (posterior_full):\n")
  if (!is.null(mcar) && nrow(mcar)) {
    cat("  MCAR (gated to [0.92, 0.98]):\n")
    for (i in seq_len(nrow(mcar))) cat(sprintf(
      "    regime %d trait %s: coverage=%.4f (n=%d)\n",
      mcar$regime_id[i], mcar$trait[i], mcar$coverage[i], mcar$n[i]))
  }
  if (!is.null(mar) && nrow(mar)) {
    cat("  MAR_phylo (reported, not gated):\n")
    for (i in seq_len(nrow(mar))) cat(sprintf(
      "    regime %d trait %s: coverage=%.4f (n=%d)\n",
      mar$regime_id[i], mar$trait[i], mar$coverage[i], mar$n[i]))
  }
}

if (length(args) >= 1L && identical(args[[1L]], "--selftest")) {
  pass_df <- data.frame(
    regime_id = c(1, 2, 17, 19), mechanism = c("MCAR", "MCAR", "MCAR", "MAR_phylo"),
    method = "posterior_full", trait = c("x", "x", "y", "x"),
    n = c(100, 100, 100, 100), covered_sum = c(95, 94, 96, 70),
    n_fits = 200L, n_converged = 200L
  )
  fail_cov_df <- pass_df
  fail_cov_df$covered_sum[1] <- 60   # coverage 0.60, well outside [0.92, 0.98]

  fail_nonconv_df <- pass_df
  fail_nonconv_df$n_converged[fail_nonconv_df$regime_id == 1] <- 190L  # 5% > 2%

  pass_fails        <- run_gate(pass_df)
  fail_cov_fails     <- run_gate(fail_cov_df)
  fail_nonconv_fails <- run_gate(fail_nonconv_df)

  ok <- length(pass_fails) == 0L && length(fail_cov_fails) > 0L &&
    length(fail_nonconv_fails) > 0L &&
    any(grepl("non-converged", fail_nonconv_fails))
  cat("pass-fixture failures:", length(pass_fails), "\n")
  cat("fail-coverage-fixture failures:", length(fail_cov_fails), "\n")
  cat("fail-nonconverged-fixture failures:", length(fail_nonconv_fails), "\n")
  if (ok) {
    cat("SELFTEST_OK\n")
  } else {
    cat("SELFTEST FAILED\n")
    if (length(pass_fails)) cat(" pass fixture unexpectedly failed\n")
    if (!length(fail_cov_fails)) cat(" fail-coverage fixture unexpectedly passed\n")
    if (!length(fail_nonconv_fails) || !any(grepl("non-converged", fail_nonconv_fails))) {
      cat(" fail-nonconverged fixture did not trigger the non-converged rule\n")
    }
  }
  quit(save = "no", status = if (ok) 0L else 1L)
}

if (length(args) < 1L) stop("expected: <cell_coverage.csv> or --selftest", call. = FALSE)
df <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE)

fails <- run_gate(df)
report_side_by_side(fails)
if (length(fails)) {
  cat("G7 FAILURES:\n")
  for (f in fails) cat(" -", f, "\n")
  quit(save = "no", status = 1L)
} else {
  cat("CELL_COVERAGE_PASS\n")
}
