#!/usr/bin/env Rscript
# script/mi_gls/05_cell_coverage.R
#
# G7 (.unlazy/mi-posterior/GATES.md): per-cell 95% predictive coverage gate
# for multi_impute(draws_method = "posterior"), method = "posterior_full",
# using the cell-coverage CSV written by script/mi_gls/03_summarise_v2.R
# (its 4th argument). Columns expected: regime_id, method, trait,
# mechanism, n, covered_sum, mean_width, n_expected, n_reps_scored, n_fits,
# n_converged, median_max_rhat, median_min_ess.
#
# Gated regimes (CP2 follow-up, Shinichi 2026-09-24; design.md section
# 5e): only the in-model regimes 17-40 (Kronecker 17-24, twins 25-40) can
# FAIL. The gated rows are the MCAR (exchangeable) rows of those regimes:
# 17, 18, 21, 22 and the twins 25, 26, 27, 28, 33, 34, 35, 36, each
# per regime x masked trait, coverage in [0.92, 0.98]. The MCAR rows of
# regimes 1-16 (raw tree covariance, outside the sampler's model) are
# printed under "STRESS TEST (reported, not gated)" with the same rule
# outcomes. Which regimes are gated comes from regimes.R, never the CSV.
#
# Expected MCAR set, fail-closed (design.md section 5c, D3), built from
# regimes.R for every regime of the grid: trait x in x_only regimes, traits
# x and y in both-missing regimes. Completeness fails in ANY regime,
# stress regimes included (the report must be complete): a missing row, a
# duplicate row, a row built for a different rep count
# (n_expected != MI_N_REPS), n <= 0, or a row where fewer reps were scored
# than converged (n_reps_scored < n_converged: a converged rep without cell
# intervals for that trait). Rule outcomes (coverage outside the band, a
# non-finite coverage, non-convergence) fail in gated regimes and are
# reported for stress regimes. MAR_phylo (clade-biased) masks are reported
# but never gated (design.md section 5 / review R5).
#
# Conformal comparator (D7): the conformal per-cell coverage and mean width
# (impute(gnn = FALSE) on the same masked cells) are printed beside
# posterior_full for every regime x trait x mechanism. Descriptive only:
# they never pass or fail the gate.
#
# Design review item B3: gate only on converged fits (03_summarise_v2.R
# restricts covered_sum/n to converged reps); a gated regime (any
# mechanism) with > 2% non-converged posterior_full fits (missing rep
# files count as non-converged) FAILS this gate and prints a loud
# "NONCONVERGED regime <id>: <k>/<n>" line; a stress regime over 2% prints
# the same line marked "stress test, not gated". The 2% is compared on
# counts (non-converged > floor(0.02 n)), so exactly 4/200 passes and
# 5/200 fails.
#
# Env: MI_N_REPS (default 200), MI_REGIMES (default all; a restricted run
# can pass its rules but never prints the G7 token).
#
# Usage:
#   Rscript script/mi_gls/05_cell_coverage.R <cell_coverage.csv>
#   Rscript script/mi_gls/05_cell_coverage.R --selftest

args <- commandArgs(trailingOnly = TRUE)

source(file.path("script", "mi_gls", "regimes.R"))

NONCONVERGED_THRESHOLD <- 0.02
BAND <- c(0.92, 0.98)

fmt4 <- function(x) ifelse(is.finite(x), sprintf("%.4f", x), "NA")
is_gated <- function(rid) mi_gls_v2_is_gated(rid) %in% TRUE

expected_cells <- function(regime_ids, mechanism = c("MCAR", "MAR_phylo")) {
  rg <- regimes[regimes$regime_id %in% regime_ids & regimes$mechanism %in% mechanism, ]
  do.call(rbind, lapply(seq_len(nrow(rg)), function(i) {
    data.frame(regime_id = rg$regime_id[i], trait = mi_gls_v2_traits(rg$missing[i]),
               mechanism = rg$mechanism[i], stringsAsFactors = FALSE)
  }))
}

# Returns list(msgs, ids, out): violation messages with their regime ids
# (the caller decides gated or not) and per-regime outcome labels.
check_nonconverged <- function(sub) {
  msgs <- character(0); ids <- integer(0); out <- character(0)
  pf <- unique(sub[, c("regime_id", "n_fits", "n_converged")])
  for (i in seq_len(nrow(pf))) {
    rid <- pf$regime_id[i]; n <- pf$n_fits[i]; k <- pf$n_converged[i]
    if (!is.finite(n) || !is.finite(k) || n <= 0) {
      msgs <- c(msgs, sprintf("non-converged: regime %d n_fits/n_converged non-finite", rid))
      ids <- c(ids, rid); out[as.character(rid)] <- "NONFINITE"
      next
    }
    nonconv_frac <- 1 - k / n
    # Counts, not 1 - k/n: 1 - 196/200 is 0.02000000000000002 in floating point.
    if ((n - k) > floor(NONCONVERGED_THRESHOLD * n + 1e-9)) {
      cat(sprintf("NONCONVERGED regime %d: %d/%d non-converged (> %.0f%%)%s\n",
                 rid, n - k, n, 100 * NONCONVERGED_THRESHOLD,
                 if (is_gated(rid)) "" else " [stress test, not gated]"))
      msgs <- c(msgs, sprintf("regime %d: non-converged fraction %.4f > %.2f",
                              rid, nonconv_frac, NONCONVERGED_THRESHOLD))
      ids <- c(ids, rid); out[as.character(rid)] <- sprintf("FAIL(%d/%d)", k, n)
    } else {
      out[as.character(rid)] <- sprintf("PASS(%d/%d)", k, n)
    }
  }
  list(msgs = msgs, ids = ids, out = out)
}

# Returns the gate failure messages; attributes "report" (the side-by-side
# lines), "gated_rules" (per-row outcomes of the gated MCAR rows) and
# "stress" (the STRESS TEST block for regimes 1-16, never failing).
run_gate <- function(df, regime_ids = regimes$regime_id, n_reps = mi_gls_v2_planned_reps) {
  need_cols <- c("regime_id", "method", "trait", "mechanism", "n", "covered_sum",
                 "n_expected", "n_reps_scored", "n_fits", "n_converged")
  miss_cols <- setdiff(need_cols, names(df))
  if (length(miss_cols)) {
    return(structure(sprintf("cell-coverage CSV lacks column(s): %s (rebuild it with 03_summarise_v2.R)",
                             paste(miss_cols, collapse = ", ")),
                     report = character(0), gated_rules = character(0), stress = character(0)))
  }
  if (!("mean_width" %in% names(df))) df$mean_width <- NA_real_
  df <- df[df$regime_id %in% regime_ids, , drop = FALSE]
  df$coverage <- ifelse(is.finite(df$n) & df$n > 0, df$covered_sum / df$n, NA_real_)
  sub <- df[df$method == "posterior_full", ]
  fails <- character(0)
  for (rid in regime_ids[!(regime_ids %in% sub$regime_id)]) {
    cat(sprintf("MISSING_REGIME regime %d (%s): no posterior_full rows\n", rid,
                if (is_gated(rid)) "gated" else "stress test: not gated, but the report must be complete"))
  }

  # Rule violations of every regime, split into gate failures (gated
  # regimes) and stress-test lines (regimes 1-16) at the end.
  nc <- check_nonconverged(sub)
  v_id <- nc$ids; v_msg <- nc$msgs

  # ---- expected MCAR set: completeness in every regime, band in gated ------
  exp_mcar <- expected_cells(regime_ids, "MCAR")
  rules <- character(0); rules_gated <- logical(0)
  for (i in seq_len(nrow(exp_mcar))) {
    rid <- exp_mcar$regime_id[i]
    r <- sub[sub$regime_id == rid & sub$trait == exp_mcar$trait[i], ]
    lab <- sprintf("regime %d trait %s", rid, exp_mcar$trait[i])
    tag <- if (is_gated(rid)) "" else " (stress-test regime; the report must be complete)"
    if (nrow(r) == 0L) { fails <- c(fails, sprintf("missing MCAR row: %s%s", lab, tag)); next }
    if (nrow(r) > 1L)  { fails <- c(fails, sprintf("duplicate MCAR rows (%d): %s%s", nrow(r), lab, tag)); next }
    if (!is.finite(r$n_expected) || r$n_expected != n_reps) {
      fails <- c(fails, sprintf("%s: n_expected=%s but the gate expects %d reps%s", lab, r$n_expected, n_reps, tag))
    }
    # Every converged rep must contribute cell intervals for this trait;
    # otherwise reps are dropped silently (D3).
    if (!is.finite(r$n_reps_scored) || !is.finite(r$n_converged) ||
        r$n_reps_scored < r$n_converged) {
      fails <- c(fails, sprintf("%s: n_reps_scored=%s < n_converged=%s (converged reps without cell intervals)%s",
                                lab, r$n_reps_scored, r$n_converged, tag))
    }
    cov_out <- "PASS"
    if (!is.finite(r$n) || r$n <= 0) {
      fails <- c(fails, sprintf("%s: n=%s, no scored cells%s", lab, r$n, tag))
      cov_out <- "NO_CELLS"
    } else if (!is.finite(r$coverage)) {
      v_id <- c(v_id, rid)
      v_msg <- c(v_msg, sprintf("%s: coverage non-finite (covered_sum=%s, n=%s)", lab, r$covered_sum, r$n))
      cov_out <- "NONFINITE"
    } else if (r$coverage < BAND[1] || r$coverage > BAND[2]) {
      v_id <- c(v_id, rid)
      v_msg <- c(v_msg, sprintf("%s: MCAR coverage=%.4f (n=%d) outside [%.2f, %.2f]",
                                lab, r$coverage, as.integer(r$n), BAND[1], BAND[2]))
      cov_out <- "FAIL"
    }
    cv <- nc$out[as.character(rid)]
    rules <- c(rules, sprintf("RULES regime=%d trait=%s MCAR coverage=%s (%s) converged=%s",
                              rid, exp_mcar$trait[i], fmt4(r$coverage), cov_out,
                              if (is.na(cv)) "NA" else cv))
    rules_gated <- c(rules_gated, is_gated(rid))
  }

  # ---- side-by-side report: posterior_full beside conformal (D7) -------------
  report <- "Per-cell coverage and mean width, posterior_full beside conformal (conformal is descriptive):"
  cf <- df[df$method == "conformal", ]
  all_exp <- expected_cells(regime_ids)
  for (mech in c("MCAR", "MAR_phylo")) {
    e <- all_exp[all_exp$mechanism == mech, ]
    if (!nrow(e)) next
    report <- c(report, sprintf("  %s (%s):", mech,
                                if (mech == "MCAR") "posterior_full gated to [0.92, 0.98] in regimes 17-40; regimes 1-16 stress test, not gated"
                                else "reported, not gated"))
    for (i in seq_len(nrow(e))) {
      p <- sub[sub$regime_id == e$regime_id[i] & sub$trait == e$trait[i], ][1, ]
      q <- cf[cf$regime_id == e$regime_id[i] & cf$trait == e$trait[i], ][1, ]
      report <- c(report, sprintf(
        "    regime %d trait %s%s: posterior_full coverage=%s width=%s (n=%s) | conformal coverage=%s width=%s (n=%s)",
        e$regime_id[i], e$trait[i], if (is_gated(e$regime_id[i])) "" else " [stress]",
        fmt4(p$coverage), fmt4(p$mean_width), p$n,
        fmt4(q$coverage), fmt4(q$mean_width), q$n))
    }
  }

  v_gated <- is_gated(v_id)
  fails <- c(fails, v_msg[v_gated])
  stress_ids <- regime_ids[!is_gated(regime_ids)]
  stress <- if (length(stress_ids)) {
    c(rules[!rules_gated],
      if (any(!v_gated)) paste("STRESS_VIOLATION", v_msg[!v_gated]),
      sprintf("STRESS_SUMMARY %d rule violation(s) in %d of %d stress-test regime(s); reported, not gated",
              sum(!v_gated), length(unique(v_id[!v_gated])), length(stress_ids)))
  } else character(0)
  structure(fails, report = report, gated_rules = rules[rules_gated], stress = stress)
}

if (length(args) >= 1L && identical(args[[1L]], "--selftest")) {
  # A passing table over the FULL planned grid, regimes 1-40 (no conformal
  # rows: conformal is descriptive and its absence must not fail), then
  # fixtures that must fail (gated regimes, or incomplete data in any
  # regime) and stress-test fixtures (regimes 1-16) that must NOT fail.
  make_pass <- function(n_reps = mi_gls_v2_planned_reps) {
    e <- expected_cells(regimes$regime_id)
    data.frame(e, method = "posterior_full", n = 1000L,
               covered_sum = ifelse(e$mechanism == "MCAR", 950L, 850L),   # MAR 0.85: not gated
               mean_width = 3.9, n_expected = n_reps, n_reps_scored = n_reps,
               n_fits = n_reps, n_converged = n_reps)
  }
  at <- function(d, rid, tr = c("x", "y")) d$regime_id %in% rid & d$trait %in% tr
  fx <- list()
  d <- make_pass(); d$covered_sum[at(d, 25, "x")] <- 600L
  fx$coverage_band <- d                                        # twin of 1
  d <- make_pass(); d$covered_sum[at(d, 36, "y")] <- 990L
  fx$twin_over_coverage <- d                                   # 0.99 > 0.98
  d <- make_pass(); d$n_converged[d$regime_id == 25] <- 190L
  fx$nonconverged <- d                                         # 10/200 = 5% > 2%
  d <- make_pass(); d$n_converged[d$regime_id == 31] <- 190L
  fx$nonconverged_mar_twin <- d                                # MAR twin: convergence still gated
  d <- make_pass(); fx$missing_regime <- d[d$regime_id != 2, ]   # stress regime missing: FAILS
  d <- make_pass(); fx$missing_twin <- d[d$regime_id != 34, ]    # twin regime missing: FAILS
  d <- make_pass(); fx$missing_trait_y <- d[!at(d, 21, "y"), ]
  d <- make_pass(); d$n_converged[d$regime_id == 22] <- 150L
  fx$short_reps <- d                                           # 150/200 rep files present
  fx$wrong_n_expected <- make_pass(n_reps = 100L)
  d <- make_pass(); d$covered_sum[at(d, 17, "y")] <- NA
  fx$nonfinite <- d
  d <- make_pass(); d$n[at(d, 9, "x")] <- 0L
  fx$zero_n <- d                                               # stress regime without cells: FAILS
  d <- make_pass(); d$n_reps_scored[at(d, 1, "x")] <- 197L
  fx$reps_unscored <- d                                        # 3 converged reps without intervals (stress)
  d <- make_pass(); d$n_reps_scored[at(d, 35, "y")] <- 197L
  fx$reps_unscored_twin <- d
  d <- make_pass(); d$n_reps_scored[at(d, 11, "y")] <- NA
  fx$reps_scored_na <- d
  d <- make_pass(); d$n_converged[d$regime_id == 29] <- 195L; d$n_reps_scored[d$regime_id == 29] <- 195L
  fx$nonconverged_5of200 <- d                                  # 2.5% > 2%
  d <- make_pass(); fx$no_reps_scored_col <- d[, names(d) != "n_reps_scored"]   # CSV from an older 03
  expect_pattern <- c(coverage_band = "regime 25 trait x: MCAR coverage=0.6000",
                      twin_over_coverage = "regime 36 trait y: MCAR coverage=0.9900",
                      nonconverged = "regime 25: non-converged",
                      nonconverged_mar_twin = "regime 31: non-converged",
                      missing_regime = "missing MCAR row: regime 2 trait x (stress-test regime",
                      missing_twin = "missing MCAR row: regime 34 trait x",
                      missing_trait_y = "missing MCAR row: regime 21 trait y",
                      short_reps = "regime 22: non-converged", wrong_n_expected = "n_expected=100",
                      nonfinite = "regime 17 trait y: coverage non-finite", zero_n = "regime 9 trait x: n=0",
                      reps_unscored = "regime 1 trait x: n_reps_scored=197 < n_converged=200",
                      reps_unscored_twin = "regime 35 trait y: n_reps_scored=197 < n_converged=200",
                      reps_scored_na = "regime 11 trait y: n_reps_scored=NA",
                      nonconverged_5of200 = "regime 29: non-converged fraction 0.0250 > 0.02",
                      no_reps_scored_col = "lacks column(s): n_reps_scored")

  # Stress-test fixtures (regimes 1-16): reported, must NOT fail the gate.
  sx <- list()
  d <- make_pass(); d$covered_sum[at(d, 1, "x")] <- 600L
  sx$stress_coverage_band <- d
  d <- make_pass(); d$covered_sum[d$regime_id %in% c(1:4, 9:12) & d$mechanism == "MCAR"] <- 880L
  sx$stress_all_mcar_low <- d                                  # every stress MCAR row at 0.88
  d <- make_pass(); d$n_converged[d$regime_id == 1] <- 190L
  sx$stress_nonconverged <- d
  d <- make_pass(); d$n_converged[d$regime_id == 5] <- 150L; d$n_reps_scored[d$regime_id == 5] <- 150L
  sx$stress_nonconverged_mar <- d
  d <- make_pass(); d$covered_sum[at(d, 12, "y")] <- NA
  sx$stress_nonfinite <- d
  stress_pattern <- c(stress_coverage_band = "^STRESS_VIOLATION regime 1 trait x: MCAR coverage=0.6000",
                      stress_all_mcar_low = "^STRESS_SUMMARY 12 rule violation\\(s\\) in 8 of 16",
                      stress_nonconverged = "^STRESS_VIOLATION regime 1: non-converged fraction 0.0500",
                      stress_nonconverged_mar = "^STRESS_VIOLATION regime 5: non-converged fraction 0.2500",
                      stress_nonfinite = "^STRESS_VIOLATION regime 12 trait y: coverage non-finite")

  ok <- TRUE
  pass_fails <- run_gate(make_pass())
  cat(sprintf("pass fixture (no conformal rows): %d failure(s)\n", length(pass_fails)))
  if (length(pass_fails)) { ok <- FALSE; cat(paste(" -", pass_fails), sep = "\n") }
  n_gr <- length(attr(pass_fails, "gated_rules")); st0 <- attr(pass_fails, "stress")
  gated_set <- sort(unique(as.integer(sub("^RULES regime=([0-9]+) .*", "\\1", attr(pass_fails, "gated_rules")))))
  want_set <- c(17L, 18L, 21L, 22L, 25L, 26L, 27L, 28L, 33L, 34L, 35L, 36L)
  st_ok <- identical(gated_set, want_set) && sum(grepl("^RULES regime=", st0)) == 12L &&
    any(grepl("^STRESS_SUMMARY 0 rule violation", st0))
  cat(sprintf("pass fixture: gated MCAR regimes = {%s} (want {%s}), %d gated rows, 12 stress rows + clean STRESS_SUMMARY: %s\n",
              paste(gated_set, collapse = ","), paste(want_set, collapse = ","), n_gr, if (st_ok) "yes" else "NO"))
  if (!st_ok) ok <- FALSE
  d <- make_pass(); cfr <- d; cfr$method <- "conformal"; cfr$covered_sum <- 500L
  cfr$n_fits <- NA_integer_; cfr$n_converged <- NA_integer_
  with_conf <- run_gate(rbind(d, cfr))
  conf_shown <- any(grepl("conformal coverage=0.5000", attr(with_conf, "report")))
  cat(sprintf("pass fixture + bad conformal rows: %d failure(s) (want 0), conformal printed: %s\n",
              length(with_conf), if (conf_shown) "yes" else "NO"))
  if (length(with_conf) || !conf_shown) ok <- FALSE
  d <- make_pass(); d$n_converged[d$regime_id == 29] <- 196L; d$n_reps_scored[d$regime_id == 29] <- 196L
  edge <- run_gate(d)                                          # exactly 4/200 = 2%: not above 2%
  cat(sprintf("pass fixture, exactly 4/200 non-converged in regime 29: %d failure(s) (want 0)\n", length(edge)))
  if (length(edge)) { ok <- FALSE; cat(paste(" -", edge), sep = "\n") }
  for (nm in names(fx)) {
    invisible(capture.output(r <- run_gate(fx[[nm]])))
    hit <- any(grepl(expect_pattern[[nm]], r, fixed = TRUE))
    cat(sprintf("fail fixture %-22s: %d failure(s), expected reason %s: %s\n",
                nm, length(r), shQuote(expect_pattern[[nm]]), if (hit) "yes" else "NO"))
    if (!length(r) || !hit) { ok <- FALSE; cat(paste(" -", r), sep = "\n") }
  }
  for (nm in names(sx)) {
    invisible(capture.output(r <- run_gate(sx[[nm]])))
    st <- attr(r, "stress")
    hit <- any(grepl(stress_pattern[[nm]], st))
    cat(sprintf("stress fixture %-24s: %d gate failure(s) (want 0), reported in STRESS TEST block %s: %s\n",
                nm, length(r), shQuote(stress_pattern[[nm]]), if (hit) "yes" else "NO"))
    if (length(r) || !hit) { ok <- FALSE; cat(paste(" -", c(r, st)), sep = "\n") }
  }
  f <- capture.output(r <- run_gate(sx$stress_nonconverged))
  loud <- any(grepl("^NONCONVERGED regime 1: 10/200 non-converged .*stress test, not gated", f))
  cat(sprintf("stress NONCONVERGED line labelled 'stress test, not gated': %s\n", if (loud) "yes" else "NO"))
  if (!loud) ok <- FALSE
  f <- capture.output(r <- run_gate(fx$missing_regime))
  loud <- any(grepl("^MISSING_REGIME regime 2 \\(stress test", f))
  cat(sprintf("missing stress regime prints a loud MISSING_REGIME line: %s\n", if (loud) "yes" else "NO"))
  if (!loud) ok <- FALSE
  cat(if (ok) "SELFTEST_OK\n" else "SELFTEST FAILED\n")
  quit(save = "no", status = if (ok) 0L else 1L)
}

if (length(args) < 1L) stop("expected: <cell_coverage.csv> or --selftest", call. = FALSE)
df <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE)
ex <- mi_gls_v2_expected()
cat(sprintf("EXPECTED_GRID %s\n", ex$label))

fails <- run_gate(df, regime_ids = ex$regime_ids, n_reps = ex$n_reps)
cat(attr(fails, "report"), sep = "\n")
g_ids <- intersect(ex$regime_ids[is_gated(ex$regime_ids)], regimes$regime_id[regimes$mechanism == "MCAR"])
cat(sprintf("GATED MCAR REGIMES (in-model: Kronecker 17-24, twins 25-40; Shinichi, CP2 follow-up 2026-09-24): %s\n",
            if (length(g_ids)) paste(g_ids, collapse = " ") else "none in this grid"))
if (length(attr(fails, "gated_rules"))) cat(paste(" ", attr(fails, "gated_rules")), sep = "\n")
cat("STRESS TEST (reported, not gated): regimes 1-16, raw tree covariance outside the sampler's model\n")
st <- attr(fails, "stress")
if (length(st)) cat(paste(" ", st), sep = "\n") else cat("  no stress-test regimes in this grid\n")
if (length(fails)) {
  cat("G7 FAILURES:\n")
  for (f in fails) cat(" -", f, "\n")
  quit(save = "no", status = 1L)
} else if (!ex$full) {
  cat(sprintf("CELL_COVERAGE_SUBSET_OK (%s): rules hold on a partial grid; the G7 token needs all regimes x %d reps\n",
              ex$label, mi_gls_v2_planned_reps))
} else {
  cat("CELL_COVERAGE_PASS\n")
}
