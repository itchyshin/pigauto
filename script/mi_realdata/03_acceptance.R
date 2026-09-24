#!/usr/bin/env Rscript
# G8 acceptance gate (.unlazy/mi-posterior/GATES.md):
#   "real-data report complete (all planned cells) with model-based vs
#    conformal per-cell coverage and masked downstream slopes"
#
# Fail-closed (orchestrator decision R1, 2026-09-24; tightened in the
# repair round the same day). Prints REALDATA_COMPLETE only if ALL of
# these hold, and otherwise prints one `FAIL ...` line naming every failing
# cell, trait or pair (lib.R::check_acceptance() has the exact rules):
#   1. every planned (dataset, arm, seed) cell has a receipt with status
#      "ok" (a crashed or killed 01_run.R leaves status "running");
#   2. all ok receipts carry ONE non-NA code SHA, recorded from MI_POST_SHA
#      or a clean git HEAD;
#   3. every ok receipt (current schema) has a sampler record showing
#      pigauto defaults (no MI_POST_NITER / MI_POST_BURNIN smoke override);
#   4. every ok receipt records >= 1 kept trait; for every kept trait (and
#      every trait in the coverage table), n_matched == n_masked, where
#      n_matched counts masked cells with a finite model interval, and the
#      model coverage is finite;
#   5. every ok receipt read the split AND the Mondrian conformal results
#      with status "ok", with, for every masked kept trait, a finite
#      coverage and n_interval == n_masked == the model's n_masked;
#   6. every pre-registered pair (pairs.R) has at least one cell whose
#      reference and MI legs both have status "ok", finite slopes and SEs,
#      and an MI pool of all m completions (m_used == m_total), on the
#      pre-registered analysis scale. A pair-cell with status "ok" that
#      misses the finiteness or full-pool condition fails on its own, so a
#      selectively pooled slope can never be reported as if it were whole.
#
# The 5% slope criterion is REPORTED, NOT GATED (orchestrator decision R3,
# 2026-09-24). Reason: the approved plan's Plan-1 real-data check is
# descriptive with a completeness gate (plan table, "Real-data check" row:
# "descriptive, completeness gate"), and a relative difference is unstable
# when the slope is near zero (e.g. the pre-registered Troph ~ Length pair
# is chosen to sit near the detection limit). It is aggregated over ALL
# cells, per pair and overall: count within 5%, max |relative difference|,
# max |absolute difference| and max |difference in reference-SE units|.
# Convergence is reported separately and does not gate either.
#
# Usage:
#   Rscript 03_acceptance.R [outdir]    # default outdir: script/mi_realdata/returned
#   Rscript 03_acceptance.R --selftest  # synthetic fixtures, prints SELFTEST_OK

# Rscript passes a script path containing spaces as `~+~` in --file=
# (e.g. ".../Github~+~Local/..."); undo that before normalizePath().
here <- function() {
  a <- commandArgs(FALSE)
  f <- gsub("~+~", " ", sub("^--file=", "", a[grepl("^--file=", a)]), fixed = TRUE)
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
here_dir <- here()
source(file.path(here_dir, "lib.R"))
source(file.path(here_dir, "pairs.R"))

args <- commandArgs(trailingOnly = TRUE)
selftest <- "--selftest" %in% args
args <- setdiff(args, "--selftest")

show <- function(df) {
  num <- vapply(df, is.double, logical(1))
  df[num] <- lapply(df[num], signif, digits = 4)
  print(df, row.names = FALSE)
}

report <- function(outdir, planned, pairs_list, label) {
  receipts <- collect_receipts(outdir, planned)
  acc <- check_acceptance(receipts, planned, pairs_list)
  cov_tab <- build_coverage_table(receipts)
  slope_tab <- build_slope_table(receipts)
  conv_tab <- build_convergence_table(receipts)

  cat(sprintf("=== mi-posterior real-data acceptance (%s), outdir = %s ===\n", label, outdir))
  cat("\nCell status:\n"); show(acc$cell_status)
  cat("\nPer-trait coverage of masked cells (model-based vs split vs Mondrian):\n")
  if (nrow(cov_tab)) {
    show(cov_tab[c("dataset", "arm", "seed", "trait", "n_masked", "n_matched", "model_coverage",
                   "split_coverage", "mondrian_coverage")])
    cat("\nInterval widths, like for like (mean and median, original scale):\n")
    show(cov_tab[c("dataset", "arm", "seed", "trait", "model_mean_width", "split_mean_width",
                   "mondrian_mean_width", "model_median_width", "split_median_width",
                   "mondrian_median_width")])
  } else {
    cat("  (none yet)\n")
  }
  cat("\nDownstream slope check per cell (5% criterion reported, not gated):\n")
  if (nrow(slope_tab)) {
    show(slope_tab[c("dataset", "arm", "seed", "response", "predictor", "ref_slope", "mi_slope",
                     "diff", "rel_diff", "diff_ref_se", "within_5pct", "ref_status", "mi_status",
                     "m_used", "m_total", "n_nonfinite")])
  } else {
    cat("  (none yet)\n")
  }
  cat("\nPre-registered pairs, aggregated over all cells:\n"); show(acc$pair_status)
  cat(sprintf("SLOPE_5PCT (reported, not gated): %d of %d pair-cells within 5%%\n",
              acc$slope_overall$n_within_5pct, acc$slope_overall$n_pair_cells))

  # Convergence is reported separately from the coverage headline above and
  # never gates REALDATA_COMPLETE (design-review addition, 2026-09-24).
  # Clade-biased ("structured") masks are descriptive only here, same as
  # the coverage/slope tables -- no separate convergence rule for them.
  nonconverged <- conv_tab[!is.na(conv_tab$converged) & !conv_tab$converged, , drop = FALSE]
  cat("\n--- Convergence (reported separately; does not gate REALDATA_COMPLETE) ---\n")
  if (nrow(nonconverged)) {
    for (i in seq_len(nrow(nonconverged))) {
      cat(sprintf("NONCONVERGED %s/%s (max_rhat=%.3f, min_ess=%.0f)\n",
                  nonconverged$dataset[[i]], nonconverged$name[[i]],
                  nonconverged$max_rhat[[i]], nonconverged$min_ess[[i]]))
    }
  } else if (nrow(conv_tab)) {
    cat("All fits with a convergence record report converged.\n")
  } else {
    cat("(no convergence records yet)\n")
  }

  cat("\n--- G8 gate ---\n")
  if (length(acc$failures)) cat(paste0("FAIL ", acc$failures, "\n"), sep = "")
  acc$conv_tab <- conv_tab
  acc$nonconverged <- nonconverged
  acc
}

if (selftest) {
  # Every variant in lib.R::selftest_variants(): "pass" must be complete
  # with no failures (positive control, although one cell is
  # non-converged and one slope is outside 5%); every other variant breaks
  # one G8 condition and must NOT be complete, and must name what failed.
  results <- list()
  for (variant in names(selftest_variants())) {
    spec <- selftest_variants()[[variant]]
    tmp <- tempfile(paste0("mi_realdata_selftest_", variant, "_"))
    dir.create(tmp)
    fx <- build_selftest_fixture(tmp, variant)
    printed <- capture.output(acc <- report(tmp, fx$planned, fx$pairs_list, paste("selftest", variant)))
    unlink(tmp, recursive = TRUE)
    named <- if (all(is.na(spec$must_name))) length(acc$failures) == 0L else
      all(vapply(spec$must_name, function(s) any(startsWith(acc$failures, s)), logical(1)))
    ok <- identical(acc$complete, spec$complete) && named
    cat(sprintf("selftest 03 %-21s expect complete=%-5s got complete=%-5s failures=%d named=%s -> %s\n",
                variant, spec$complete, acc$complete, length(acc$failures), named,
                if (ok) "ok" else "MISMATCH"))
    for (f in acc$failures) cat("    FAIL ", f, "\n", sep = "")
    results[[variant]] <- list(ok = ok, acc = acc, printed = printed)
  }
  pass <- results$pass
  stopifnot(
    all(vapply(results, `[[`, logical(1), "ok")),
    # pass fixture: convergence and the 5% criterion are reported, not gated
    nrow(pass$acc$nonconverged) == 1L,
    any(grepl("^NONCONVERGED synth/synth-mcar-m1 ", pass$printed)),
    pass$acc$pair_status$n_within_5pct == 1L, pass$acc$pair_status$n_rel_finite == 2L,
    isFALSE(pass$acc$pair_status$all_within_5pct),
    any(grepl("^SLOPE_5PCT \\(reported, not gated\\): 1 of 2 pair-cells within 5%", pass$printed)),
    # a failing variant prints its FAIL line in the gate section
    any(grepl("^FAIL cell synth-mcar-m2: status 'missing'", results$missing_cell$printed)),
    # selective pooling is visible in the printed slope table, not only in FAIL
    any(grepl(" ok +ok +2 +20 +18$", results$partial_pool$printed)),
    results$partial_pool$acc$pair_status$n_cells_degraded == 1L,
    results$partial_pool$acc$pair_status$n_rel_finite == 1L   # degraded cell kept out of the 5% count
  )
  cat("\nSELFTEST_OK\n")
  quit(status = 0L)
}

outdir <- if (length(args) >= 1L) args[[1L]] else file.path(here_dir, "returned")
acc <- report(outdir, planned_cells(), mi_realdata_pairs, "real data")
if (acc$complete) {
  cat("\nREALDATA_COMPLETE\n")
} else {
  cat(sprintf("\nNOT complete: %d G8 condition(s) failed (FAIL lines above).\n", length(acc$failures)))
}
