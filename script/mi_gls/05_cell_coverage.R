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
# MCAR (exchangeable) masks are GATED to per-regime x trait coverage in
# [0.92, 0.98]. The expected MCAR set comes from script/mi_gls/regimes.R,
# never from the CSV (design.md section 5c, D3): trait x in x_only
# regimes, traits x and y in both-missing regimes. A missing row, a row
# built for a different rep count (n_expected != MI_N_REPS), n <= 0, a
# non-finite coverage, or a row where fewer reps were scored than
# converged (n_reps_scored < n_converged: a converged rep without cell
# intervals for that trait) FAILS. MAR_phylo (clade-biased) masks are
# reported but never gated (design.md section 5 / review R5).
#
# Conformal comparator (D7): the conformal per-cell coverage and mean width
# (impute(gnn = FALSE) on the same masked cells) are printed beside
# posterior_full for every regime x trait x mechanism. Descriptive only:
# they never pass or fail the gate.
#
# Design review item B3: gate only on converged fits (03_summarise_v2.R
# restricts covered_sum/n to converged reps); a regime with > 2%
# non-converged posterior_full fits (missing rep files count as
# non-converged) FAILS this gate and prints a loud
# "NONCONVERGED regime <id>: <k>/<n>" line. The 2% is compared on counts
# (non-converged > floor(0.02 n)), so exactly 4/200 passes and 5/200 fails.
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

expected_cells <- function(regime_ids, mechanism = c("MCAR", "MAR_phylo")) {
  rg <- regimes[regimes$regime_id %in% regime_ids & regimes$mechanism %in% mechanism, ]
  do.call(rbind, lapply(seq_len(nrow(rg)), function(i) {
    data.frame(regime_id = rg$regime_id[i], trait = mi_gls_v2_traits(rg$missing[i]),
               mechanism = rg$mechanism[i], stringsAsFactors = FALSE)
  }))
}

check_nonconverged <- function(sub) {
  fails <- character(0)
  pf <- unique(sub[, c("regime_id", "n_fits", "n_converged")])
  for (i in seq_len(nrow(pf))) {
    n <- pf$n_fits[i]; k <- pf$n_converged[i]
    if (!is.finite(n) || !is.finite(k) || n <= 0) {
      fails <- c(fails, sprintf("non-converged: regime %d n_fits/n_converged non-finite", pf$regime_id[i]))
      next
    }
    nonconv_frac <- 1 - k / n
    # Counts, not 1 - k/n: 1 - 196/200 is 0.02000000000000002 in floating point.
    if ((n - k) > floor(NONCONVERGED_THRESHOLD * n + 1e-9)) {
      cat(sprintf("NONCONVERGED regime %d: %d/%d non-converged (> %.0f%%)\n",
                 pf$regime_id[i], n - k, n, 100 * NONCONVERGED_THRESHOLD))
      fails <- c(fails, sprintf("regime %d: non-converged fraction %.4f > %.2f",
                                pf$regime_id[i], nonconv_frac, NONCONVERGED_THRESHOLD))
    }
  }
  fails
}

# Returns the failure messages; attr "report" holds the side-by-side lines.
run_gate <- function(df, regime_ids = regimes$regime_id, n_reps = mi_gls_v2_planned_reps) {
  need_cols <- c("regime_id", "method", "trait", "mechanism", "n", "covered_sum",
                 "n_expected", "n_reps_scored", "n_fits", "n_converged")
  miss_cols <- setdiff(need_cols, names(df))
  if (length(miss_cols)) {
    return(structure(sprintf("cell-coverage CSV lacks column(s): %s (rebuild it with 03_summarise_v2.R)",
                             paste(miss_cols, collapse = ", ")), report = character(0)))
  }
  if (!("mean_width" %in% names(df))) df$mean_width <- NA_real_
  df <- df[df$regime_id %in% regime_ids, , drop = FALSE]
  df$coverage <- ifelse(is.finite(df$n) & df$n > 0, df$covered_sum / df$n, NA_real_)
  sub <- df[df$method == "posterior_full", ]
  fails <- check_nonconverged(sub)

  # ---- gated MCAR set, fail-closed (D3) --------------------------------------
  exp_mcar <- expected_cells(regime_ids, "MCAR")
  for (i in seq_len(nrow(exp_mcar))) {
    r <- sub[sub$regime_id == exp_mcar$regime_id[i] & sub$trait == exp_mcar$trait[i], ]
    lab <- sprintf("regime %d trait %s", exp_mcar$regime_id[i], exp_mcar$trait[i])
    if (nrow(r) == 0L) { fails <- c(fails, sprintf("missing MCAR row: %s", lab)); next }
    if (nrow(r) > 1L)  { fails <- c(fails, sprintf("duplicate MCAR rows (%d): %s", nrow(r), lab)); next }
    if (!is.finite(r$n_expected) || r$n_expected != n_reps) {
      fails <- c(fails, sprintf("%s: n_expected=%s but the gate expects %d reps", lab, r$n_expected, n_reps))
    }
    # Every converged rep must contribute cell intervals for this trait;
    # otherwise reps are dropped silently (D3).
    if (!is.finite(r$n_reps_scored) || !is.finite(r$n_converged) ||
        r$n_reps_scored < r$n_converged) {
      fails <- c(fails, sprintf("%s: n_reps_scored=%s < n_converged=%s (converged reps without cell intervals)",
                                lab, r$n_reps_scored, r$n_converged))
    }
    if (!is.finite(r$n) || r$n <= 0) {
      fails <- c(fails, sprintf("%s: n=%s, no scored cells", lab, r$n))
    } else if (!is.finite(r$coverage)) {
      fails <- c(fails, sprintf("%s: coverage non-finite (covered_sum=%s, n=%s)", lab, r$covered_sum, r$n))
    } else if (r$coverage < BAND[1] || r$coverage > BAND[2]) {
      fails <- c(fails, sprintf("%s: MCAR coverage=%.4f (n=%d) outside [%.2f, %.2f]",
                                lab, r$coverage, as.integer(r$n), BAND[1], BAND[2]))
    }
  }

  # ---- side-by-side report: posterior_full beside conformal (D7) -------------
  report <- "Per-cell coverage and mean width, posterior_full beside conformal (conformal is descriptive):"
  cf <- df[df$method == "conformal", ]
  all_exp <- expected_cells(regime_ids)
  for (mech in c("MCAR", "MAR_phylo")) {
    e <- all_exp[all_exp$mechanism == mech, ]
    if (!nrow(e)) next
    report <- c(report, sprintf("  %s (%s):", mech,
                                if (mech == "MCAR") "posterior_full gated to [0.92, 0.98]" else "reported, not gated"))
    for (i in seq_len(nrow(e))) {
      p <- sub[sub$regime_id == e$regime_id[i] & sub$trait == e$trait[i], ][1, ]
      q <- cf[cf$regime_id == e$regime_id[i] & cf$trait == e$trait[i], ][1, ]
      report <- c(report, sprintf(
        "    regime %d trait %s: posterior_full coverage=%s width=%s (n=%s) | conformal coverage=%s width=%s (n=%s)",
        e$regime_id[i], e$trait[i], fmt4(p$coverage), fmt4(p$mean_width), p$n,
        fmt4(q$coverage), fmt4(q$mean_width), q$n))
    }
  }
  structure(fails, report = report)
}

if (length(args) >= 1L && identical(args[[1L]], "--selftest")) {
  # A passing table over the FULL planned grid (no conformal rows: conformal
  # is descriptive and its absence must not fail), then fixtures that must fail.
  make_pass <- function(n_reps = mi_gls_v2_planned_reps) {
    e <- expected_cells(regimes$regime_id)
    data.frame(e, method = "posterior_full", n = 1000L,
               covered_sum = ifelse(e$mechanism == "MCAR", 950L, 850L),   # MAR 0.85: not gated
               mean_width = 3.9, n_expected = n_reps, n_reps_scored = n_reps,
               n_fits = n_reps, n_converged = n_reps)
  }
  at <- function(d, rid, tr = c("x", "y")) d$regime_id %in% rid & d$trait %in% tr
  fx <- list()
  d <- make_pass(); d$covered_sum[at(d, 1, "x")] <- 600L
  fx$coverage_band <- d
  d <- make_pass(); d$n_converged[d$regime_id == 1] <- 190L
  fx$nonconverged <- d                                         # 10/200 = 5% > 2%
  d <- make_pass(); fx$missing_regime <- d[d$regime_id != 2, ]
  d <- make_pass(); fx$missing_trait_y <- d[!at(d, 21, "y"), ]
  d <- make_pass(); d$n_converged[d$regime_id == 22] <- 150L
  fx$short_reps <- d                                           # 150/200 rep files present
  fx$wrong_n_expected <- make_pass(n_reps = 100L)
  d <- make_pass(); d$covered_sum[at(d, 17, "y")] <- NA
  fx$nonfinite <- d
  d <- make_pass(); d$n[at(d, 9, "x")] <- 0L
  fx$zero_n <- d
  d <- make_pass(); d$n_reps_scored[at(d, 1, "x")] <- 197L
  fx$reps_unscored <- d                                        # 3 converged reps without intervals
  d <- make_pass(); d$n_reps_scored[at(d, 11, "y")] <- NA
  fx$reps_scored_na <- d
  d <- make_pass(); d$n_converged[d$regime_id == 5] <- 195L; d$n_reps_scored[d$regime_id == 5] <- 195L
  fx$nonconverged_5of200 <- d                                  # 2.5% > 2%
  d <- make_pass(); fx$no_reps_scored_col <- d[, names(d) != "n_reps_scored"]   # CSV from an older 03
  expect_pattern <- c(coverage_band = "regime 1 trait x: MCAR coverage=0.6000",
                      nonconverged = "non-converged", missing_regime = "missing MCAR row: regime 2 trait x",
                      missing_trait_y = "missing MCAR row: regime 21 trait y",
                      short_reps = "regime 22: non-converged", wrong_n_expected = "n_expected=100",
                      nonfinite = "regime 17 trait y: coverage non-finite", zero_n = "regime 9 trait x: n=0",
                      reps_unscored = "regime 1 trait x: n_reps_scored=197 < n_converged=200",
                      reps_scored_na = "regime 11 trait y: n_reps_scored=NA",
                      nonconverged_5of200 = "regime 5: non-converged fraction 0.0250 > 0.02",
                      no_reps_scored_col = "lacks column(s): n_reps_scored")

  ok <- TRUE
  pass_fails <- run_gate(make_pass())
  cat(sprintf("pass fixture (no conformal rows): %d failure(s)\n", length(pass_fails)))
  if (length(pass_fails)) { ok <- FALSE; cat(paste(" -", pass_fails), sep = "\n") }
  d <- make_pass(); cfr <- d; cfr$method <- "conformal"; cfr$covered_sum <- 500L
  cfr$n_fits <- NA_integer_; cfr$n_converged <- NA_integer_
  with_conf <- run_gate(rbind(d, cfr))
  conf_shown <- any(grepl("conformal coverage=0.5000", attr(with_conf, "report")))
  cat(sprintf("pass fixture + bad conformal rows: %d failure(s) (want 0), conformal printed: %s\n",
              length(with_conf), if (conf_shown) "yes" else "NO"))
  if (length(with_conf) || !conf_shown) ok <- FALSE
  d <- make_pass(); d$n_converged[d$regime_id == 5] <- 196L; d$n_reps_scored[d$regime_id == 5] <- 196L
  edge <- run_gate(d)                                          # exactly 4/200 = 2%: not above 2%
  cat(sprintf("pass fixture, exactly 4/200 non-converged in regime 5: %d failure(s) (want 0)\n", length(edge)))
  if (length(edge)) { ok <- FALSE; cat(paste(" -", edge), sep = "\n") }
  for (nm in names(fx)) {
    invisible(capture.output(r <- run_gate(fx[[nm]])))
    hit <- any(grepl(expect_pattern[[nm]], r, fixed = TRUE))
    cat(sprintf("fail fixture %-16s: %d failure(s), expected reason %s: %s\n",
                nm, length(r), shQuote(expect_pattern[[nm]]), if (hit) "yes" else "NO"))
    if (!length(r) || !hit) { ok <- FALSE; cat(paste(" -", r), sep = "\n") }
  }
  cat(if (ok) "SELFTEST_OK\n" else "SELFTEST FAILED\n")
  quit(save = "no", status = if (ok) 0L else 1L)
}

if (length(args) < 1L) stop("expected: <cell_coverage.csv> or --selftest", call. = FALSE)
df <- utils::read.csv(args[[1L]], stringsAsFactors = FALSE)
ex <- mi_gls_v2_expected()
cat(sprintf("EXPECTED_GRID %s\n", ex$label))

fails <- run_gate(df, regime_ids = ex$regime_ids, n_reps = ex$n_reps)
cat(attr(fails, "report"), sep = "\n")
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
