#!/usr/bin/env Rscript
# script/mi_gls/04_acceptance.R
#
# G6 (.unlazy/mi-posterior/GATES.md): simulation acceptance gate for
# multi_impute(draws_method = "posterior"), method = "posterior_full",
# against the "complete" reference rows, using the summary CSV written by
# script/mi_gls/03_summarise_v2.R. Decisions D1-D6 are recorded in
# docs/dev-log/mi-posterior/design.md section 5c.
#
# Gated regimes (CP2 follow-up, Shinichi 2026-09-24; design.md section
# 5e): every rule is evaluated for every regime, but only the gated
# (in-model) regimes 17-40 can FAIL the gate: the Kronecker regimes 17-24
# and the in-model twins 25-40. Regimes 1-16 simulate from the raw
# covariance of a non-ultrametric tree, outside the sampler's model; their
# rule outcomes are printed under a "STRESS TEST (reported, not gated)"
# block. Which regimes are gated comes from regimes.R (regimes$gated),
# never from the CSV.
#
# Completeness, fail-closed (D3), for EVERY regime of the grid, stress
# regimes included (the report must be complete): the expected rows come
# from script/mi_gls/regimes.R, never from the CSV: every regime x
# downstream (gls, phylolm) needs a "complete" row and a posterior_full
# row, plus a posterior_none row where missing == "both" (D6). A missing
# row, a duplicate row, a row built for a different rep count
# (n_expected != MI_N_REPS), a row with fewer rep files present than
# expected (n_present < MI_N_REPS, e.g. 150 of 200) or a row with R <= 0
# FAILS; nothing is skipped. A regime with no rows at all also prints a
# loud MISSING_REGIME line. A non-finite rule value is a rule outcome: it
# fails in a gated regime and is reported in a stress regime, and so are
# non-converged or failed fits when every rep file is present.
#
# Rules (every expected regime x downstream unless noted):
#   1. |paired bias| <= max(0.02, 2.5 * paired_bias_mcse)
#   2. coverage >= complete coverage - 0.05 per row, and mean shortfall
#      <= 0.02, where shortfall = pmax(complete - coverage, 0) (the
#      positive part, D5), averaged over the gated posterior_full rows (the
#      stress rows' mean is reported separately). The coverage truth is
#      rho = 0.7 in regimes 1-16 and 25-40 and the complete-data
#      pseudo-truth in 17-24 (D1, computed in 03).
#   3. SE ratio in [0.90, 1.15] under the lambda (phylolm) analysis. Env
#      MI_SE_RULE selects the reading (D2): "pooled_relative" (DEFAULT since
#      Shinichi's G6 decision, 2026-09-25, option 3) gates the MEAN of
#      posterior_full ratio / complete-data ratio over the gated phylolm
#      rows in the tighter band [0.95, 1.10], and reports each row's
#      relative ratio (rows outside [0.90, 1.15] get an
#      SE_RATIO_ROW_OUTSIDE line, not gated); "relative" (CP1,
#      Shinichi 2026-09-24) gates that relative ratio per row; "absolute"
#      (the original plan) gates the posterior_full ratio itself per row.
#      A non-finite gated row fails under every reading. Both numbers are
#      always printed, plus an ANALYSIS_MODEL_SE_RATIO line whenever the
#      complete-data ratio is itself outside the band.
#   4. proper vs improper SE ratio (D4): REPORTED, NOT GATED since CP1
#      (Shinichi 2026-09-24). On the 30-tip S1 test fixture, fixing Sigma
#      at its posterior mean gave a per-cell predictive variance 1.4-1.6%
#      LARGER than the proper one (intervals about 0.7-0.8% wider; Jensen),
#      so "proper > improper" need not hold when everything is right
#      (design.md 5d.4). At campaign sizes the direction reversed: plug-in
#      per-cell widths were 0.002-0.18% narrower in 48 of 48 rows
#      (cell_coverage.csv) and the plug-in SE ratio was lower in 48 of 48
#      pairs (sim_summary.csv; both-missing regimes 9-24 and 33-40).
#      Mean posterior_full vs posterior_none se_ratio
#      is printed per downstream model within each both-missing block
#      (9-16 stress test, 17-24 Kronecker, 33-40 twins), plus per-regime
#      pairs. A missing posterior_none row still fails through the
#      completeness check.
#   5. fit failures <= 2% (missing rep files count as failures, D3)
#   6. (design review B3) non-converged posterior_full fits <= 2% of the
#      expected reps (a missing rep file counts as non-converged); a
#      regime over that prints a loud "NONCONVERGED regime <id>: <k>/<n>"
#      line and FAILS (a stress regime's line says "stress test, not
#      gated"). Rules 1-4 use converged reps only (03).
#   Rules 5 and 6 fail only ABOVE 2%: rule 6 compares counts
#   (non-converged > floor(0.02 n)) and rule 5 allows 1e-9 of rounding, so
#   exactly 4/200 passes and 5/200 fails (1 - 196/200 is 0.02000000000000002
#   in floating point). Rule 2 likewise fails only ABOVE 0.05 per row and
#   above 0.02 in the mean, with the same 1e-9 allowance: a shortfall of
#   exactly 10/200 passes (0.905 - 0.855 is 0.05000000000000004 in floating
#   point; M2 review 2026-09-24).
#
# Env: MI_N_REPS (default 200), MI_REGIMES (default all; a restricted run
# can pass its rules but never prints the G6 token), MI_SE_RULE
# (pooled_relative | relative | absolute, default pooled_relative).
#
# Usage:
#   Rscript script/mi_gls/04_acceptance.R <summary.csv>
#   Rscript script/mi_gls/04_acceptance.R --selftest

args <- commandArgs(trailingOnly = TRUE)

source(file.path("script", "mi_gls", "regimes.R"))

NONCONVERGED_THRESHOLD <- 0.02
SE_BAND <- c(0.90, 1.15)
POOLED_BAND <- c(0.95, 1.10)   # G6 option 3 (Shinichi, 2026-09-25; results.md option c)

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

# Returns the gate failure messages. Attributes: "report" (descriptive
# lines), "gated_rules" (per-row rule outcomes of the gated regimes) and
# "stress" (the STRESS TEST block: rule outcomes and violations of the
# regimes 1-16, reported, never failing).
run_gate <- function(df, regime_ids = regimes$regime_id,
                     n_reps = mi_gls_v2_planned_reps, se_rule = "pooled_relative") {
  if (!(se_rule %in% c("absolute", "relative", "pooled_relative"))) {
    stop("MI_SE_RULE must be 'pooled_relative', 'relative' or 'absolute'; got '", se_rule, "'",
         call. = FALSE)
  }
  fails <- character(0)
  report <- character(0)
  key <- function(d) paste(d$regime_id, d$method, d$downstream)
  is_gated <- function(rid) mi_gls_v2_is_gated(rid) %in% TRUE

  need_cols <- c("regime_id", "method", "downstream", "n_expected", "n_present", "R", "paired_bias",
                 "paired_bias_mcse", "se_ratio", "coverage", "fit_failure_rate",
                 "n_fits", "n_converged")
  miss_cols <- setdiff(need_cols, names(df))
  if (length(miss_cols)) {
    return(structure(sprintf("summary CSV lacks column(s): %s (rebuild it with 03_summarise_v2.R)",
                             paste(miss_cols, collapse = ", ")),
                     report = report, gated_rules = character(0), stress = character(0)))
  }

  # ---- completeness (D3), every regime including the stress test -------------
  req <- required_rows(regime_ids)
  df <- df[df$regime_id %in% regime_ids, , drop = FALSE]
  for (rid in regime_ids[!(regime_ids %in% df$regime_id)]) {
    cat(sprintf("MISSING_REGIME regime %d (%s): no summary rows\n", rid,
                if (is_gated(rid)) "gated" else "stress test: not gated, but the report must be complete"))
  }
  k_df <- key(df)
  for (i in seq_len(nrow(req))) {
    kk <- key(req[i, ])
    hit <- df[k_df == kk, , drop = FALSE]
    lab <- sprintf("regime %d %s %s", req$regime_id[i], req$method[i], req$downstream[i])
    tag <- if (is_gated(req$regime_id[i])) "" else " (stress-test regime; the report must be complete)"
    if (nrow(hit) == 0L) { fails <- c(fails, sprintf("missing row: %s%s", lab, tag)); next }
    if (nrow(hit) > 1L)  { fails <- c(fails, sprintf("duplicate rows (%d): %s%s", nrow(hit), lab, tag)); next }
    if (!is.finite(hit$n_expected) || hit$n_expected != n_reps) {
      fails <- c(fails, sprintf("%s: n_expected=%s but the gate expects %d reps%s",
                                lab, hit$n_expected, n_reps, tag))
    }
    if (!is.finite(hit$n_present) || hit$n_present < n_reps) {
      fails <- c(fails, sprintf("%s: n_present=%s of %d expected rep files%s",
                                lab, hit$n_present, n_reps, tag))
    }
    if (!is.finite(hit$R) || hit$R <= 0) {
      fails <- c(fails, sprintf("%s: R=%s, no usable reps%s", lab, hit$R, tag))
    }
  }

  df <- merge(df, regimes[, c("regime_id", "missing")], by = "regime_id", all.x = TRUE)
  cmp <- df[df$method == "complete", c("regime_id", "downstream", "se_ratio", "coverage")]
  names(cmp)[3:4] <- c("c_se_ratio", "c_coverage")
  pf <- merge(df[df$method == "posterior_full", ], cmp, by = c("regime_id", "downstream"), all.x = TRUE)
  pn <- df[df$method == "posterior_none", ]
  pf <- pf[order(pf$regime_id, pf$downstream), ]
  pf$gated <- is_gated(pf$regime_id)
  lab_pf <- sprintf("regime %d %s", pf$regime_id, pf$downstream)

  # Rule violations of every regime; split into gate failures (gated
  # regimes) and stress-test lines (regimes 1-16) at the end.
  v_id <- integer(0); v_msg <- character(0)
  flag <- function(bad, msg) {
    bad <- bad %in% TRUE
    if (any(bad)) { v_id <<- c(v_id, pf$regime_id[bad]); v_msg <<- c(v_msg, msg[bad]) }
  }
  nonfinite <- function(rule, what, bad) flag(bad, sprintf("%s: %s non-finite for %s", rule, what, lab_pf))
  outcome <- function(ok, bad) ifelse(!ok, "NONFINITE", ifelse(bad, "FAIL", "PASS"))

  # ---- rule 6: non-convergence (loud; blocking in gated regimes) ------------
  nc <- unique(pf[, c("regime_id", "n_fits", "n_converged")])
  conv_out <- character(0)
  for (i in seq_len(nrow(nc))) {
    rid <- nc$regime_id[i]; n <- nc$n_fits[i]; k <- nc$n_converged[i]
    if (!is.finite(n) || !is.finite(k) || n <= 0) {
      v_id <- c(v_id, rid)
      v_msg <- c(v_msg, sprintf("non-converged: regime %d n_fits/n_converged non-finite", rid))
      conv_out[as.character(rid)] <- "NONFINITE"
      next
    }
    frac <- 1 - k / n
    if ((n - k) > floor(NONCONVERGED_THRESHOLD * n + 1e-9)) {
      cat(sprintf("NONCONVERGED regime %d: %d/%d non-converged (> %.0f%%)%s\n",
                  rid, n - k, n, 100 * NONCONVERGED_THRESHOLD,
                  if (is_gated(rid)) "" else " [stress test, not gated]"))
      v_id <- c(v_id, rid)
      v_msg <- c(v_msg, sprintf("regime %d: non-converged fraction %.4f > %.2f",
                                rid, frac, NONCONVERGED_THRESHOLD))
      conv_out[as.character(rid)] <- sprintf("FAIL(%d/%d)", k, n)
    } else {
      conv_out[as.character(rid)] <- sprintf("PASS(%d/%d)", k, n)
    }
  }

  # ---- rule 1: paired bias --------------------------------------------------
  ok1 <- is.finite(pf$paired_bias) & is.finite(pf$paired_bias_mcse)
  nonfinite("bias", "paired_bias/paired_bias_mcse", !ok1)
  bad1 <- ok1 & abs(pf$paired_bias) > pmax(0.02, 2.5 * pf$paired_bias_mcse)
  flag(bad1, sprintf("bias: %s |bias|=%.4f > max(0.02, 2.5*MCSE=%.4f)",
                     lab_pf, abs(pf$paired_bias), pf$paired_bias_mcse))

  # ---- rule 2: coverage vs complete (D1, D5) --------------------------------
  shortfall <- pf$c_coverage - pf$coverage
  ok2 <- is.finite(shortfall)
  nonfinite("coverage", "coverage/complete coverage", !ok2)
  bad2 <- ok2 & shortfall > 0.05 + 1e-9
  flag(bad2, sprintf("coverage: %s coverage=%.4f complete=%.4f shortfall=%.4f > 0.05",
                     lab_pf, pf$coverage, pf$c_coverage, shortfall))
  mean_sf <- list(); stress_sf_line <- character(0)
  for (grp in c("gated", "stress")) {
    sel <- ok2 & (pf$gated == (grp == "gated"))
    if (!any(sel)) next
    ms <- mean(pmax(shortfall[sel], 0))
    mean_sf[[grp]] <- ms
    line <- sprintf("MEAN_SHORTFALL %s %.4f over %d posterior_full rows (positive part; signed mean %.4f)%s",
                    grp, ms, sum(sel), mean(shortfall[sel]),
                    if (grp == "gated") "" else " [stress test, not gated]")
    if (grp == "gated") report <- c(report, line) else stress_sf_line <- line
  }
  if (!is.null(mean_sf$gated) && mean_sf$gated > 0.02 + 1e-9) fails <- c(fails, sprintf(
    "coverage: mean shortfall (positive part) over the gated regimes = %.4f > 0.02", mean_sf$gated))
  stress_mean_sf <- if (!is.null(mean_sf$stress) && mean_sf$stress > 0.02 + 1e-9) sprintf(
    "coverage: mean shortfall (positive part) over the stress-test regimes = %.4f > 0.02", mean_sf$stress) else character(0)

  # ---- rule 3: SE ratio under phylolm (D2) ----------------------------------
  rel <- pf$se_ratio / pf$c_se_ratio
  for (i in seq_len(nrow(pf))) {
    gated3 <- pf$downstream[i] == "phylolm"
    report <- c(report, sprintf("SE_RATIO regime=%d downstream=%s posterior_full=%s complete=%s relative=%s%s",
                                pf$regime_id[i], pf$downstream[i], fmt4(pf$se_ratio[i]),
                                fmt4(pf$c_se_ratio[i]), fmt4(rel[i]),
                                if (!gated3) " (descriptive)"
                                else if (pf$gated[i]) sprintf(" (gated, rule=%s)", se_rule)
                                else sprintf(" (stress test, rule=%s, not gated)", se_rule)))
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
  nonfinite("se_ratio", what3, r3 & !ok3)
  out3 <- r3 & ok3 & (val3 < SE_BAND[1] | val3 > SE_BAND[2])
  if (se_rule == "pooled_relative") {
    # Per-row relative ratios are reported; the gate is their mean over the
    # gated rows (Shinichi, 2026-09-25). The stress rows' mean is reported.
    bad3 <- rep(FALSE, nrow(pf))
    for (i in which(out3)) report <- c(report, sprintf(
      "SE_RATIO_ROW_OUTSIDE (reported, not gated) %s relative=%.4f outside [%.2f, %.2f]",
      lab_pf[i], val3[i], SE_BAND[1], SE_BAND[2]))
    for (grp in c("gated", "stress")) {
      sel <- r3 & ok3 & (pf$gated == (grp == "gated"))
      if (!any(sel)) next
      pm <- mean(val3[sel])
      report <- c(report, sprintf(
        "POOLED_SE_RATIO %s mean relative=%.4f over %d phylolm rows, band [%.2f, %.2f] (range %.4f-%.4f; %d row(s) outside [%.2f, %.2f])%s",
        grp, pm, sum(sel), POOLED_BAND[1], POOLED_BAND[2], min(val3[sel]), max(val3[sel]),
        sum(out3 & sel), SE_BAND[1], SE_BAND[2],
        if (grp == "gated") "" else " [stress test, not gated]"))
      if (grp == "gated" && (pm < POOLED_BAND[1] || pm > POOLED_BAND[2])) fails <- c(fails, sprintf(
        "se_ratio: pooled mean relative ratio over the gated phylolm rows = %.4f outside [%.2f, %.2f]",
        pm, POOLED_BAND[1], POOLED_BAND[2]))
    }
  } else {
    bad3 <- out3
    flag(bad3, sprintf("se_ratio: %s %s=%.4f outside [%.2f, %.2f] (rule=%s)",
                       lab_pf, what3, val3, SE_BAND[1], SE_BAND[2], se_rule))
  }

  # ---- rule 4: proper > improper, per both-missing block (D4), reported ----
  both <- regimes[regimes$missing == "both", ]
  for (dg in c("tree_raw", "kronecker", "twin")) {
    all_ids <- both$regime_id[both$dgp == dg]
    ids <- all_ids[all_ids %in% regime_ids]
    if (!length(ids)) next
    bn <- sprintf("%d-%d", min(all_ids), max(all_ids))
    for (ds in c("gls", "phylolm")) {
      full <- pf$se_ratio[match(paste(ids, ds), paste(pf$regime_id, pf$downstream))]
      none <- pn$se_ratio[match(paste(ids, ds), paste(pn$regime_id, pn$downstream))]
      for (j in seq_along(ids)) report <- c(report, sprintf(
        "SE_RATIO_PAIR regime=%d downstream=%s full=%s none=%s", ids[j], ds, fmt4(full[j]), fmt4(none[j])))
      okp <- is.finite(full) & is.finite(none)
      report <- c(report, sprintf(
        "SE_RATIO_ORDER (reported, not gated) regimes %s %s [%s]: mean posterior_full=%s posterior_none=%s proper_wider=%s (pairs %d of %d finite)",
        bn, ds, dg, fmt4(mean(full[okp])), fmt4(mean(none[okp])),
        if (any(okp)) mean(full[okp]) > mean(none[okp]) else NA, sum(okp), length(ids)))
    }
  }

  # ---- rule 5: fit failures -------------------------------------------------
  ok5 <- is.finite(pf$fit_failure_rate)
  nonfinite("fit_failure", "fit_failure_rate", !ok5)
  bad5 <- ok5 & pf$fit_failure_rate > 0.02 + 1e-9
  flag(bad5, sprintf("fit_failure: %s failure_rate=%.4f > 0.02", lab_pf, pf$fit_failure_rate))

  # ---- per-row outcomes, gate failures and the stress-test block ------------
  cv <- conv_out[as.character(pf$regime_id)]
  rules <- sprintf("RULES regime=%d downstream=%s bias=%s coverage=%s se_ratio=%s fit_fail=%s converged=%s",
                   pf$regime_id, pf$downstream, outcome(ok1, bad1), outcome(ok2, bad2),
                   ifelse(r3, ifelse(ok3 & !bad3 & out3, "OUTSIDE(pooled)", outcome(ok3, bad3)), "n/a"),
                   outcome(ok5, bad5),
                   ifelse(is.na(cv), "NA", cv))
  v_gated <- is_gated(v_id)
  fails <- c(fails, v_msg[v_gated])
  stress_ids <- regime_ids[!is_gated(regime_ids)]
  stress <- if (length(stress_ids)) {
    n_bad <- length(unique(v_id[!v_gated]))
    c(rules[!pf$gated],
      stress_sf_line,
      if (length(v_msg[!v_gated]) || length(stress_mean_sf)) paste("STRESS_VIOLATION", c(v_msg[!v_gated], stress_mean_sf)),
      sprintf("STRESS_SUMMARY %d rule violation(s) in %d of %d stress-test regime(s); reported, not gated",
              sum(!v_gated) + length(stress_mean_sf), n_bad, length(stress_ids)))
  } else character(0)
  structure(fails, report = report, gated_rules = rules[pf$gated], stress = stress)
}

if (length(args) >= 1L && identical(args[[1L]], "--selftest")) {
  # A passing summary over the FULL planned grid (regimes 1-40), then
  # fixtures that must fail (gated regimes 17-40, or incomplete data in any
  # regime) and stress-test fixtures (regimes 1-16) that must NOT fail.
  make_pass <- function(n_reps = mi_gls_v2_planned_reps) {
    req <- required_rows(regimes$regime_id)
    d <- req
    d$n_expected <- n_reps
    d$n_present <- n_reps
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
  d <- make_pass(); d$paired_bias[at(d, 27, "posterior_full", "phylolm")] <- 0.5
  fx$twin_bias <- d                                       # gated twin (of 3)
  d <- make_pass(); d$coverage[at(d, 35, "posterior_full", "gls")] <- 0.80
  fx$twin_coverage <- d
  d <- make_pass(); d$n_converged[at(d, 17, "posterior_full")] <- 190L
  fx$nonconverged <- d                                    # 10/200 = 5% > 2%
  d <- make_pass(); d$n_converged[at(d, 36, "posterior_full")] <- 190L
  fx$twin_nonconverged <- d
  d <- make_pass(); fx$missing_regime <- d[d$regime_id != 5, ]     # stress regime missing: FAILS
  d <- make_pass(); fx$missing_twin <- d[d$regime_id != 33, ]      # twin regime missing: FAILS
  d <- make_pass(); d$R[at(d, 3, "posterior_full")] <- 0L
  fx$stress_no_usable_reps <- d                           # stress row with R = 0: report incomplete
  d <- make_pass(); d$n_present[d$regime_id == 27] <- 150L; d$R[d$regime_id == 27] <- 150L
  d$n_converged[at(d, 27, "posterior_full")] <- 150L
  d$fit_failure_rate[at(d, 27, "posterior_full")] <- 0.25
  fx$short_reps <- d                                      # 150/200 rep files present (twin 27)
  d <- make_pass(); d$n_present[d$regime_id == 3] <- 150L; d$R[d$regime_id == 3] <- 150L
  d$n_converged[at(d, 3, "posterior_full")] <- 150L
  d$fit_failure_rate[at(d, 3, "posterior_full")] <- 0.25
  fx$stress_short_files <- d                              # stress regime, 150/200 rep files: FAILS
  fx$wrong_n_expected <- make_pass(n_reps = 100L)         # summary built for 100 reps
  d <- make_pass(); d$se_ratio[at(d, 21, "posterior_full", "phylolm")] <- NA
  fx$nonfinite <- d
  d <- make_pass(); d$se_ratio[at(d, 17:24, "posterior_none", "gls")] <- 1.2
  rule4_reported <- d                                     # improper wider: reported, must NOT fail
  d <- make_pass(); pfr <- which(d$method == "posterior_full")
  d$coverage[pfr] <- rep(c(0.99, 0.905), length.out = length(pfr))
  fx$shortfall_cancel <- d                                # signed mean ~0.0025, positive part 0.0225
  d <- make_pass(); d <- d[!at(d, 12, "posterior_none", "gls"), ]
  fx$missing_pair <- d                                    # stress regime row missing: FAILS
  d <- make_pass(); d$n_converged[at(d, 22, "posterior_full")] <- 195L
  d$fit_failure_rate[at(d, 22, "posterior_full")] <- 1 - 195 / 200
  fx$nonconverged_5of200 <- d                             # 2.5% > 2%: both rules fail
  d <- make_pass(); d$se_ratio[at(d, 21, c("complete", "posterior_full"), "phylolm")] <- 0.80
  analysis_model <- d                                     # fails absolute, passes relative
  d <- make_pass(); d$se_ratio[at(d, c(35, 36, 38), "posterior_full", "phylolm")] <- 1.18
  pooled_rows <- d                                        # 3 rows outside: fails relative, passes pooled
  d <- make_pass(); d$se_ratio[at(d, 17:40, "posterior_full", "phylolm")] <- 1.12
  fx$pooled_high <- d                                     # pooled mean 1.12 > 1.10: FAILS (rows in [0.90, 1.15])
  d <- make_pass(); d$n_converged[at(d, 22, "posterior_full")] <- 196L
  d$fit_failure_rate[at(d, 22, "posterior_full")] <- 1 - 196 / 200
  edge_2pct <- d                                          # exactly 4/200 = 2%: passes
  d <- make_pass(); d$coverage[at(d, 27, "complete", "phylolm")] <- 0.905
  d$coverage[at(d, 27, "posterior_full", "phylolm")] <- 0.855
  edge_cov <- d                                           # shortfall exactly 10/200 = 0.05: passes

  # Stress-test fixtures (regimes 1-16): reported, must NOT fail the gate.
  sx <- list()
  d <- make_pass(); d$paired_bias[at(d, 3, "posterior_full", "phylolm")] <- -0.5
  sx$stress_bias <- d
  d <- make_pass(); d$n_converged[at(d, 1, "posterior_full")] <- 190L
  sx$stress_nonconverged <- d
  d <- make_pass(); d$coverage[d$method == "posterior_full" & d$regime_id <= 16] <- 0.80
  sx$stress_coverage <- d                                 # per-row and mean shortfall in 1-16
  d <- make_pass(); d$n_converged[at(d, 3, "posterior_full")] <- 150L
  d$fit_failure_rate[at(d, 3, "posterior_full")] <- 0.25
  sx$stress_failed_fits <- d                              # every rep file present, 25% failed fits
  d <- make_pass(); d$se_ratio[at(d, 2, "posterior_full", "phylolm")] <- NA
  sx$stress_nonfinite <- d

  expect_pattern <- c(bias = "^bias: regime 17 gls",
                      twin_bias = "^bias: regime 27 phylolm",
                      twin_coverage = "^coverage: regime 35 gls coverage=0.8000",
                      nonconverged = "regime 17: non-converged",
                      twin_nonconverged = "regime 36: non-converged",
                      missing_regime = "missing row: regime 5 .*stress-test regime",
                      missing_twin = "missing row: regime 33 ",
                      stress_no_usable_reps = "regime 3 posterior_full gls: R=0, no usable reps",
                      short_reps = "^regime 27 posterior_full gls: n_present=150 of 200 expected rep files$",
                      stress_short_files = "^regime 3 complete gls: n_present=150 of 200 expected rep files \\(stress-test regime",
                      wrong_n_expected = "n_expected=100",
                      nonfinite = "se_ratio: .*non-finite for regime 21 phylolm",
                      pooled_high = "^se_ratio: pooled mean relative ratio .* = 1.1200 outside \\[0.95, 1.10\\]",
                      shortfall_cancel = "mean shortfall .*gated regimes", missing_pair = "regime 12 posterior_none gls",
                      nonconverged_5of200 = "regime 22: non-converged fraction 0.0250")
  stress_pattern <- c(stress_bias = "^STRESS_VIOLATION bias: regime 3 phylolm",
                      stress_nonconverged = "^STRESS_VIOLATION regime 1: non-converged fraction 0.0500",
                      stress_coverage = "^STRESS_VIOLATION coverage: mean shortfall .*stress-test regimes",
                      stress_failed_fits = "^STRESS_VIOLATION fit_failure: regime 3 ",
                      stress_nonfinite = "^STRESS_VIOLATION se_ratio: .*non-finite for regime 2 phylolm")

  ok <- TRUE
  pass_fails <- run_gate(make_pass())
  cat(sprintf("pass fixture: %d failure(s)\n", length(pass_fails)))
  if (length(pass_fails)) { ok <- FALSE; cat(paste(" -", pass_fails), sep = "\n") }
  st0 <- attr(pass_fails, "stress")
  st_ok <- any(grepl("^STRESS_SUMMARY 0 rule violation\\(s\\) in 0 of 16 stress-test regime", st0)) &&
    sum(grepl("^RULES regime=", st0)) == 32L && length(attr(pass_fails, "gated_rules")) == 48L
  cat(sprintf("pass fixture: 32 stress RULES lines + clean STRESS_SUMMARY, 48 gated RULES lines: %s\n",
              if (st_ok) "yes" else "NO"))
  if (!st_ok) ok <- FALSE
  pass_rel <- run_gate(make_pass(), se_rule = "relative")
  cat(sprintf("pass fixture (MI_SE_RULE=relative): %d failure(s)\n", length(pass_rel)))
  if (length(pass_rel)) ok <- FALSE
  edge <- run_gate(edge_2pct)
  cat(sprintf("pass fixture, exactly 4/200 non-converged and failed in regime 22: %d failure(s) (want 0)\n",
              length(edge)))
  if (length(edge)) { ok <- FALSE; cat(paste(" -", edge), sep = "\n") }
  edge_c <- run_gate(edge_cov)
  cat(sprintf("pass fixture, coverage shortfall exactly 10/200 = 0.05 in regime 27: %d failure(s) (want 0)\n",
              length(edge_c)))
  if (length(edge_c)) { ok <- FALSE; cat(paste(" -", edge_c), sep = "\n") }
  for (nm in names(fx)) {
    f <- suppressMessages(capture.output(r <- run_gate(fx[[nm]]), type = "output"))
    hit <- any(grepl(expect_pattern[[nm]], r))
    cat(sprintf("fail fixture %-22s: %d failure(s), expected reason %s: %s\n",
                nm, length(r), shQuote(expect_pattern[[nm]]), if (hit) "yes" else "NO"))
    if (!length(r) || !hit) { ok <- FALSE; cat(paste(" -", r), sep = "\n") }
  }
  for (nm in names(sx)) {
    f <- capture.output(r <- run_gate(sx[[nm]]), type = "output")
    st <- attr(r, "stress")
    hit <- any(grepl(stress_pattern[[nm]], st))
    cat(sprintf("stress fixture %-20s: %d gate failure(s) (want 0), reported in STRESS TEST block %s: %s\n",
                nm, length(r), shQuote(stress_pattern[[nm]]), if (hit) "yes" else "NO"))
    if (length(r) || !hit) { ok <- FALSE; cat(paste(" -", c(r, st)), sep = "\n") }
  }
  f <- capture.output(r <- run_gate(sx$stress_nonconverged), type = "output")
  loud <- any(grepl("^NONCONVERGED regime 1: 10/200 non-converged .*stress test, not gated", f))
  cat(sprintf("stress NONCONVERGED line labelled 'stress test, not gated': %s\n", if (loud) "yes" else "NO"))
  if (!loud) ok <- FALSE
  f <- capture.output(r <- run_gate(fx$missing_regime), type = "output")
  loud <- any(grepl("^MISSING_REGIME regime 5 \\(stress test", f))
  cat(sprintf("missing stress regime prints a loud MISSING_REGIME line: %s\n", if (loud) "yes" else "NO"))
  if (!loud) ok <- FALSE
  r4 <- run_gate(rule4_reported)
  r4_line <- any(grepl("^SE_RATIO_ORDER .*regimes 17-24 gls.*proper_wider=FALSE", attr(r4, "report")))
  r4_twin <- any(grepl("^SE_RATIO_ORDER .*regimes 33-40 phylolm \\[twin\\]", attr(r4, "report")))
  cat(sprintf("rule-4 fixture (improper wider in 17-24): %d failure(s) (want 0), reported %s, twin block 33-40 %s\n",
              length(r4), if (r4_line) "yes" else "NO", if (r4_twin) "yes" else "NO"))
  if (length(r4) || !r4_line || !r4_twin) ok <- FALSE
  am_abs <- run_gate(analysis_model, se_rule = "absolute")
  am_rel <- run_gate(analysis_model, se_rule = "relative")
  am_line <- any(grepl("^ANALYSIS_MODEL_SE_RATIO regime=21 complete=0.8000", attr(am_rel, "report")))
  cat(sprintf("analysis-model fixture: absolute %d failure(s) (want >0), relative %d (want 0), ANALYSIS_MODEL line %s\n",
              length(am_abs), length(am_rel), if (am_line) "printed" else "MISSING"))
  if (!length(am_abs) || length(am_rel) || !am_line) ok <- FALSE
  pr_rel <- run_gate(pooled_rows, se_rule = "relative")
  pr_pool <- run_gate(pooled_rows, se_rule = "pooled_relative")
  pr_lines <- sum(grepl("^SE_RATIO_ROW_OUTSIDE .*regime (35|36|38) phylolm relative=1.1800", attr(pr_pool, "report")))
  pr_mean <- any(grepl("^POOLED_SE_RATIO gated mean relative=1.0225 over 24 phylolm rows", attr(pr_pool, "report")))
  cat(sprintf("pooled-rows fixture: relative %d failure(s) (want 3), pooled_relative %d (want 0), ROW_OUTSIDE lines %d (want 3), POOLED line %s\n",
              length(pr_rel), length(pr_pool), pr_lines, if (pr_mean) "printed" else "MISSING"))
  if (length(pr_rel) != 3L || length(pr_pool) || pr_lines != 3L || !pr_mean) ok <- FALSE
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
se_rule <- Sys.getenv("MI_SE_RULE", "pooled_relative")
cat(sprintf("EXPECTED_GRID %s se_rule=%s\n", ex$label, se_rule))
if ("code_sha" %in% names(df)) cat(sprintf("PROVENANCE code_sha=%s\n", paste(unique(df$code_sha), collapse = "|")))

fails <- run_gate(df, regime_ids = ex$regime_ids, n_reps = ex$n_reps, se_rule = se_rule)
rep_lines <- attr(fails, "report")
if (length(rep_lines)) cat(rep_lines, sep = "\n")
g_ids <- ex$regime_ids[mi_gls_v2_is_gated(ex$regime_ids) %in% TRUE]
cat(sprintf("GATED REGIMES (in-model: Kronecker 17-24, twins 25-40; Shinichi, CP2 follow-up 2026-09-24): %s\n",
            if (length(g_ids)) paste(g_ids, collapse = " ") else "none in this grid"))
if (length(attr(fails, "gated_rules"))) cat(paste(" ", attr(fails, "gated_rules")), sep = "\n")
cat("STRESS TEST (reported, not gated): regimes 1-16, raw tree covariance outside the sampler's model\n")
st <- attr(fails, "stress")
if (length(st)) cat(paste(" ", st), sep = "\n") else cat("  no stress-test regimes in this grid\n")
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
