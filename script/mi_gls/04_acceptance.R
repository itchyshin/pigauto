#!/usr/bin/env Rscript
# script/mi_gls/04_acceptance.R
#
# G6 (.unlazy/mi-posterior/GATES.md): simulation acceptance gate for
# multi_impute(draws_method = "posterior"), method = "posterior_full",
# against the "complete" reference rows, using the summary CSV written by
# script/mi_gls/03_summarise_v2.R. Decisions D1-D6 are recorded in
# docs/dev-log/mi-posterior/design.md section 5c.
#
# Completeness, fail-closed (D3): the expected rows come from
# script/mi_gls/regimes.R, never from the CSV: every regime x downstream
# (gls, phylolm) needs a "complete" row and a posterior_full row, plus a
# posterior_none row where missing == "both" (D6). A missing row, a row
# built for a different rep count (n_expected != MI_N_REPS), a row with
# R <= 0, or a non-finite gated value FAILS; nothing is skipped.
#
# Rules (every expected regime x downstream unless noted):
#   1. |paired bias| <= max(0.02, 2.5 * paired_bias_mcse)
#   2. coverage >= complete coverage - 0.05 per row, and mean shortfall
#      <= 0.02, where shortfall = pmax(complete - coverage, 0) (the
#      positive part, D5), averaged over all posterior_full rows. The
#      coverage truth is rho = 0.7 in regimes 1-16 and the complete-data
#      pseudo-truth in 17-24 (D1, computed in 03).
#   3. SE ratio in [0.90, 1.15] under the lambda (phylolm) analysis. Env
#      MI_SE_RULE selects the reading (D2): "relative" (DEFAULT since CP1,
#      Shinichi 2026-09-24) gates posterior_full ratio / complete-data
#      ratio; "absolute" (the original plan) gates the posterior_full ratio
#      itself. Both numbers are always printed, plus an
#      ANALYSIS_MODEL_SE_RATIO line whenever the complete-data ratio is
#      itself outside the band.
#   4. proper vs improper SE ratio (D4): REPORTED, NOT GATED since CP1
#      (Shinichi 2026-09-24). With Sigma fixed at the posterior mean the
#      plug-in intervals come out about 1.5% wider at these sample sizes
#      (Jensen; S1 measurement), so "proper > improper" need not hold when
#      everything is right. Mean posterior_full vs posterior_none se_ratio
#      is printed within regimes 9-16 and 17-24 per downstream model, plus
#      per-regime pairs. A missing posterior_none row still fails through
#      the completeness check.
#   5. fit failures <= 2% (missing rep files count as failures, D3)
#   6. (design review B3) non-converged posterior_full fits <= 2% of the
#      expected reps (a missing rep file counts as non-converged); a
#      regime over that prints a loud "NONCONVERGED regime <id>: <k>/<n>"
#      line and FAILS. Rules 1-4 use converged reps only (03).
#   Rules 5 and 6 fail only ABOVE 2%: rule 6 compares counts
#   (non-converged > floor(0.02 n)) and rule 5 allows 1e-9 of rounding, so
#   exactly 4/200 passes and 5/200 fails (1 - 196/200 is 0.02000000000000002
#   in floating point).
#
# Env: MI_N_REPS (default 200), MI_REGIMES (default all; a restricted run
# can pass its rules but never prints the G6 token), MI_SE_RULE
# (relative | absolute, default relative).
#
# Usage:
#   Rscript script/mi_gls/04_acceptance.R <summary.csv>
#   Rscript script/mi_gls/04_acceptance.R --selftest

args <- commandArgs(trailingOnly = TRUE)

source(file.path("script", "mi_gls", "regimes.R"))

NONCONVERGED_THRESHOLD <- 0.02
SE_BAND <- c(0.90, 1.15)

fmt4 <- function(x) ifelse(is.finite(x), sprintf("%.4f", x), "NA")

required_rows <- function(regime_ids) {
  rg <- regimes[regimes$regime_id %in% regime_ids, ]
  do.call(rbind, lapply(seq_len(nrow(rg)), function(i) {
    expand.grid(regime_id = rg$regime_id[i],
                method = c("complete", mi_gls_v2_methods(rg$missing[i])),
                downstream = c("gls", "phylolm"),
                KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  }))
}

# Returns the failure messages; attr "report" holds the descriptive lines.
run_gate <- function(df, regime_ids = regimes$regime_id,
                     n_reps = mi_gls_v2_planned_reps, se_rule = "relative") {
  if (!(se_rule %in% c("absolute", "relative"))) {
    stop("MI_SE_RULE must be 'absolute' or 'relative'; got '", se_rule, "'", call. = FALSE)
  }
  fails <- character(0)
  report <- character(0)
  key <- function(d) paste(d$regime_id, d$method, d$downstream)

  need_cols <- c("regime_id", "method", "downstream", "n_expected", "R", "paired_bias",
                 "paired_bias_mcse", "se_ratio", "coverage", "fit_failure_rate",
                 "n_fits", "n_converged")
  miss_cols <- setdiff(need_cols, names(df))
  if (length(miss_cols)) {
    return(structure(sprintf("summary CSV lacks column(s): %s (rebuild it with 03_summarise_v2.R)",
                             paste(miss_cols, collapse = ", ")), report = report))
  }

  # ---- completeness (D3) ----------------------------------------------------
  req <- required_rows(regime_ids)
  df <- df[df$regime_id %in% regime_ids, , drop = FALSE]
  k_df <- key(df)
  for (i in seq_len(nrow(req))) {
    kk <- key(req[i, ])
    hit <- df[k_df == kk, , drop = FALSE]
    lab <- sprintf("regime %d %s %s", req$regime_id[i], req$method[i], req$downstream[i])
    if (nrow(hit) == 0L) { fails <- c(fails, sprintf("missing row: %s", lab)); next }
    if (nrow(hit) > 1L)  { fails <- c(fails, sprintf("duplicate rows (%d): %s", nrow(hit), lab)); next }
    if (!is.finite(hit$n_expected) || hit$n_expected != n_reps) {
      fails <- c(fails, sprintf("%s: n_expected=%s but the gate expects %d reps",
                                lab, hit$n_expected, n_reps))
    }
    if (!is.finite(hit$R) || hit$R <= 0) {
      fails <- c(fails, sprintf("%s: R=%s, no usable reps", lab, hit$R))
    }
  }

  df <- merge(df, regimes[, c("regime_id", "missing")], by = "regime_id", all.x = TRUE)
  cmp <- df[df$method == "complete", c("regime_id", "downstream", "se_ratio", "coverage")]
  names(cmp)[3:4] <- c("c_se_ratio", "c_coverage")
  pf <- merge(df[df$method == "posterior_full", ], cmp, by = c("regime_id", "downstream"), all.x = TRUE)
  pn <- df[df$method == "posterior_none", ]
  pf <- pf[order(pf$regime_id, pf$downstream), ]
  lab_pf <- sprintf("regime %d %s", pf$regime_id, pf$downstream)
  nonfinite <- function(rule, what, bad) {
    if (any(bad)) sprintf("%s: %s non-finite for %s", rule, what, lab_pf[bad]) else character(0)
  }

  # ---- rule 6: non-convergence (loud, blocking) -----------------------------
  nc <- unique(pf[, c("regime_id", "n_fits", "n_converged")])
  for (i in seq_len(nrow(nc))) {
    n <- nc$n_fits[i]; k <- nc$n_converged[i]
    if (!is.finite(n) || !is.finite(k) || n <= 0) {
      fails <- c(fails, sprintf("non-converged: regime %d n_fits/n_converged non-finite", nc$regime_id[i]))
      next
    }
    frac <- 1 - k / n
    if ((n - k) > floor(NONCONVERGED_THRESHOLD * n + 1e-9)) {
      cat(sprintf("NONCONVERGED regime %d: %d/%d non-converged (> %.0f%%)\n",
                  nc$regime_id[i], n - k, n, 100 * NONCONVERGED_THRESHOLD))
      fails <- c(fails, sprintf("regime %d: non-converged fraction %.4f > %.2f",
                                nc$regime_id[i], frac, NONCONVERGED_THRESHOLD))
    }
  }

  # ---- rule 1: paired bias --------------------------------------------------
  ok1 <- is.finite(pf$paired_bias) & is.finite(pf$paired_bias_mcse)
  fails <- c(fails, nonfinite("bias", "paired_bias/paired_bias_mcse", !ok1))
  bad1 <- ok1 & abs(pf$paired_bias) > pmax(0.02, 2.5 * pf$paired_bias_mcse)
  if (any(bad1)) fails <- c(fails, sprintf(
    "bias: %s |bias|=%.4f > max(0.02, 2.5*MCSE=%.4f)",
    lab_pf[bad1], abs(pf$paired_bias[bad1]), pf$paired_bias_mcse[bad1]))

  # ---- rule 2: coverage vs complete (D1, D5) --------------------------------
  shortfall <- pf$c_coverage - pf$coverage
  ok2 <- is.finite(shortfall)
  fails <- c(fails, nonfinite("coverage", "coverage/complete coverage", !ok2))
  bad2 <- ok2 & shortfall > 0.05
  if (any(bad2)) fails <- c(fails, sprintf(
    "coverage: %s coverage=%.4f complete=%.4f shortfall=%.4f > 0.05",
    lab_pf[bad2], pf$coverage[bad2], pf$c_coverage[bad2], shortfall[bad2]))
  if (any(ok2)) {
    mean_shortfall <- mean(pmax(shortfall[ok2], 0))
    report <- c(report, sprintf("MEAN_SHORTFALL %.4f over %d posterior_full rows (positive part; signed mean %.4f)",
                                mean_shortfall, sum(ok2), mean(shortfall[ok2])))
    if (mean_shortfall > 0.02) fails <- c(fails, sprintf(
      "coverage: mean shortfall (positive part) = %.4f > 0.02", mean_shortfall))
  }

  # ---- rule 3: SE ratio under phylolm (D2) ----------------------------------
  rel <- pf$se_ratio / pf$c_se_ratio
  for (i in seq_len(nrow(pf))) {
    gated <- pf$downstream[i] == "phylolm"
    report <- c(report, sprintf("SE_RATIO regime=%d downstream=%s posterior_full=%s complete=%s relative=%s%s",
                                pf$regime_id[i], pf$downstream[i], fmt4(pf$se_ratio[i]),
                                fmt4(pf$c_se_ratio[i]), fmt4(rel[i]),
                                if (gated) sprintf(" (gated, rule=%s)", se_rule) else " (descriptive)"))
    if (is.finite(pf$c_se_ratio[i]) &&
        (pf$c_se_ratio[i] < SE_BAND[1] || pf$c_se_ratio[i] > SE_BAND[2])) {
      report <- c(report, sprintf("ANALYSIS_MODEL_SE_RATIO regime=%d complete=%.4f downstream=%s",
                                  pf$regime_id[i], pf$c_se_ratio[i], pf$downstream[i]))
    }
  }
  r3 <- pf$downstream == "phylolm"
  val3 <- if (se_rule == "absolute") pf$se_ratio else rel
  what3 <- if (se_rule == "absolute") "se_ratio" else "se_ratio/complete se_ratio"
  ok3 <- is.finite(val3)
  fails <- c(fails, nonfinite("se_ratio", what3, r3 & !ok3))
  bad3 <- r3 & ok3 & (val3 < SE_BAND[1] | val3 > SE_BAND[2])
  if (any(bad3)) fails <- c(fails, sprintf(
    "se_ratio: %s %s=%.4f outside [%.2f, %.2f] (rule=%s)",
    lab_pf[bad3], what3, val3[bad3], SE_BAND[1], SE_BAND[2], se_rule))

  # ---- rule 4: proper > improper, per block and downstream (D4) -------------
  both_ids <- regimes$regime_id[regimes$missing == "both" & regimes$regime_id %in% regime_ids]
  blocks <- list("9-16" = both_ids[both_ids <= 16L], "17-24" = both_ids[both_ids > 16L])
  for (bn in names(blocks)) {
    ids <- blocks[[bn]]
    if (!length(ids)) next
    for (ds in c("gls", "phylolm")) {
      full <- pf$se_ratio[match(paste(ids, ds), paste(pf$regime_id, pf$downstream))]
      none <- pn$se_ratio[match(paste(ids, ds), paste(pn$regime_id, pn$downstream))]
      for (j in seq_along(ids)) report <- c(report, sprintf(
        "SE_RATIO_PAIR regime=%d downstream=%s full=%s none=%s", ids[j], ds, fmt4(full[j]), fmt4(none[j])))
      okp <- is.finite(full) & is.finite(none)
      report <- c(report, sprintf(
        "SE_RATIO_ORDER (reported, not gated) regimes %s %s: mean posterior_full=%s posterior_none=%s proper_wider=%s (pairs %d of %d finite)",
        bn, ds, fmt4(mean(full[okp])), fmt4(mean(none[okp])),
        if (any(okp)) mean(full[okp]) > mean(none[okp]) else NA, sum(okp), length(ids)))
    }
  }

  # ---- rule 5: fit failures -------------------------------------------------
  ok5 <- is.finite(pf$fit_failure_rate)
  fails <- c(fails, nonfinite("fit_failure", "fit_failure_rate", !ok5))
  bad5 <- ok5 & pf$fit_failure_rate > 0.02 + 1e-9
  if (any(bad5)) fails <- c(fails, sprintf(
    "fit_failure: %s failure_rate=%.4f > 0.02", lab_pf[bad5], pf$fit_failure_rate[bad5]))

  structure(fails, report = report)
}

if (length(args) >= 1L && identical(args[[1L]], "--selftest")) {
  # A passing summary over the FULL planned grid, then fixtures that must fail.
  make_pass <- function(n_reps = mi_gls_v2_planned_reps) {
    req <- required_rows(regimes$regime_id)
    d <- req
    d$n_expected <- n_reps
    d$R <- n_reps
    d$paired_bias <- ifelse(d$method == "complete", NA_real_, 0.005)
    d$paired_bias_mcse <- ifelse(d$method == "complete", NA_real_, 0.01)
    d$se_ratio <- ifelse(d$method == "posterior_none", 0.8, 1.0)
    d$coverage <- ifelse(d$method == "complete", 0.95, 0.94)
    d$fit_failure_rate <- 0
    d$n_fits <- ifelse(d$method == "complete", NA_integer_, n_reps)
    d$n_converged <- ifelse(d$method == "complete", NA_integer_, n_reps)
    d
  }
  at <- function(d, rid, meth, ds = c("gls", "phylolm")) {
    d$regime_id %in% rid & d$method %in% meth & d$downstream %in% ds
  }
  fx <- list()
  d <- make_pass(); d$paired_bias[at(d, 17, "posterior_full", "gls")] <- 0.5
  fx$bias <- d
  d <- make_pass(); d$n_converged[at(d, 17, "posterior_full")] <- 190L
  fx$nonconverged <- d                                    # 10/200 = 5% > 2%
  d <- make_pass(); fx$missing_regime <- d[d$regime_id != 5, ]
  d <- make_pass(); d$n_converged[at(d, 3, "posterior_full")] <- 150L
  d$fit_failure_rate[at(d, 3, "posterior_full")] <- 0.25
  fx$short_reps <- d                                      # 150/200 rep files present
  fx$wrong_n_expected <- make_pass(n_reps = 100L)         # summary built for 100 reps
  d <- make_pass(); d$se_ratio[at(d, 21, "posterior_full", "phylolm")] <- NA
  fx$nonfinite <- d
  d <- make_pass(); d$se_ratio[at(d, 17:24, "posterior_none", "gls")] <- 1.2
  rule4_reported <- d                                     # improper wider: reported, must NOT fail
  d <- make_pass(); pfr <- which(d$method == "posterior_full")
  d$coverage[pfr] <- rep(c(0.99, 0.905), length.out = length(pfr))
  fx$shortfall_cancel <- d                                # signed mean ~0.0025, positive part 0.0225
  d <- make_pass(); d <- d[!at(d, 12, "posterior_none", "gls"), ]
  fx$missing_pair <- d
  d <- make_pass(); d$n_converged[at(d, 22, "posterior_full")] <- 195L
  d$fit_failure_rate[at(d, 22, "posterior_full")] <- 1 - 195 / 200
  fx$nonconverged_5of200 <- d                             # 2.5% > 2%: both rules fail
  d <- make_pass(); d$se_ratio[at(d, 21, c("complete", "posterior_full"), "phylolm")] <- 0.80
  analysis_model <- d                                     # fails absolute, passes relative
  d <- make_pass(); d$n_converged[at(d, 22, "posterior_full")] <- 196L
  d$fit_failure_rate[at(d, 22, "posterior_full")] <- 1 - 196 / 200
  edge_2pct <- d                                          # exactly 4/200 = 2%: passes

  expect_pattern <- c(bias = "^bias:", nonconverged = "non-converged",
                      missing_regime = "missing row: regime 5 ",
                      short_reps = "non-converged|fit_failure",
                      wrong_n_expected = "n_expected=100",
                      nonfinite = "non-finite",
                      shortfall_cancel = "mean shortfall", missing_pair = "regime 12 posterior_none gls",
                      nonconverged_5of200 = "regime 22: non-converged fraction 0.0250")

  ok <- TRUE
  pass_fails <- run_gate(make_pass())
  cat(sprintf("pass fixture: %d failure(s)\n", length(pass_fails)))
  if (length(pass_fails)) { ok <- FALSE; cat(paste(" -", pass_fails), sep = "\n") }
  pass_rel <- run_gate(make_pass(), se_rule = "relative")
  cat(sprintf("pass fixture (MI_SE_RULE=relative): %d failure(s)\n", length(pass_rel)))
  if (length(pass_rel)) ok <- FALSE
  edge <- run_gate(edge_2pct)
  cat(sprintf("pass fixture, exactly 4/200 non-converged and failed in regime 22: %d failure(s) (want 0)\n",
              length(edge)))
  if (length(edge)) { ok <- FALSE; cat(paste(" -", edge), sep = "\n") }
  for (nm in names(fx)) {
    f <- suppressMessages(capture.output(r <- run_gate(fx[[nm]]), type = "output"))
    hit <- any(grepl(expect_pattern[[nm]], r))
    cat(sprintf("fail fixture %-18s: %d failure(s), expected reason %s: %s\n",
                nm, length(r), shQuote(expect_pattern[[nm]]), if (hit) "yes" else "NO"))
    if (!length(r) || !hit) { ok <- FALSE; cat(paste(" -", r), sep = "\n") }
  }
  r4 <- run_gate(rule4_reported)
  r4_line <- any(grepl("^SE_RATIO_ORDER .*regimes 17-24 gls.*proper_wider=FALSE", attr(r4, "report")))
  cat(sprintf("rule-4 fixture (improper wider in 17-24): %d failure(s) (want 0), reported %s\n",
              length(r4), if (r4_line) "yes" else "NO"))
  if (length(r4) || !r4_line) ok <- FALSE
  am_abs <- run_gate(analysis_model, se_rule = "absolute")
  am_rel <- run_gate(analysis_model, se_rule = "relative")
  am_line <- any(grepl("^ANALYSIS_MODEL_SE_RATIO regime=21 complete=0.8000", attr(am_rel, "report")))
  cat(sprintf("analysis-model fixture: absolute %d failure(s) (want >0), relative %d (want 0), ANALYSIS_MODEL line %s\n",
              length(am_abs), length(am_rel), if (am_line) "printed" else "MISSING"))
  if (!length(am_abs) || length(am_rel) || !am_line) ok <- FALSE
  bad_rule <- tryCatch({ run_gate(make_pass(), se_rule = "loose"); FALSE }, error = function(e) TRUE)
  cat(sprintf("invalid MI_SE_RULE rejected: %s\n", if (bad_rule) "yes" else "NO"))
  if (!bad_rule) ok <- FALSE

  cat(if (ok) "SELFTEST_OK\n" else "SELFTEST FAILED\n")
  quit(save = "no", status = if (ok) 0L else 1L)
}

if (length(args) < 1L) stop("expected: <summary.csv> or --selftest", call. = FALSE)
summary_csv <- args[[1L]]
df <- utils::read.csv(summary_csv, stringsAsFactors = FALSE)
ex <- mi_gls_v2_expected()
se_rule <- Sys.getenv("MI_SE_RULE", "relative")
cat(sprintf("EXPECTED_GRID %s se_rule=%s\n", ex$label, se_rule))
if ("code_sha" %in% names(df)) cat(sprintf("PROVENANCE code_sha=%s\n", paste(unique(df$code_sha), collapse = "|")))

fails <- run_gate(df, regime_ids = ex$regime_ids, n_reps = ex$n_reps, se_rule = se_rule)
rep_lines <- attr(fails, "report")
if (length(rep_lines)) cat(rep_lines, sep = "\n")
if (length(fails)) {
  cat("G6 FAILURES:\n")
  for (f in fails) cat(" -", f, "\n")
  quit(save = "no", status = 1L)
} else if (!ex$full) {
  cat(sprintf("SIM_ACCEPT_SUBSET_OK (%s): rules hold on a partial grid; the G6 token needs all regimes x %d reps\n",
              ex$label, mi_gls_v2_planned_reps))
} else {
  cat("SIM_ACCEPT_PASS\n")
}
