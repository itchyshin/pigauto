#!/usr/bin/env Rscript
# G8 acceptance gate (.unlazy/mi-posterior/GATES.md):
#   "real-data report complete (all planned cells) with model-based vs
#    conformal per-cell coverage and masked downstream slopes"
# Prints REALDATA_COMPLETE iff every planned (dataset, arm, seed) cell has
# a receipt with status "ok" AND every pre-registered pair (pairs.R) has an
# attempted result. Also prints the per-trait coverage table and the slope
# table (with G8's 5% criterion reported per pair, NOT required for
# REALDATA_COMPLETE -- see design.md section 5 / GATES.md G8).
#
# Usage:
#   Rscript 03_acceptance.R [outdir]    # default outdir: script/mi_realdata/returned
#   Rscript 03_acceptance.R --selftest  # synthetic fixture, prints SELFTEST_OK

here <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
here_dir <- here()
source(file.path(here_dir, "lib.R"))
source(file.path(here_dir, "pairs.R"))

args <- commandArgs(trailingOnly = TRUE)
selftest <- "--selftest" %in% args
args <- setdiff(args, "--selftest")

report <- function(outdir, planned, pairs_list, label) {
  receipts <- collect_receipts(outdir, planned)
  acc <- check_acceptance(receipts, planned, pairs_list)
  cov_tab <- build_coverage_table(receipts)
  slope_tab <- build_slope_table(receipts)
  conv_tab <- build_convergence_table(receipts)

  cat(sprintf("=== mi-posterior real-data acceptance (%s), outdir = %s ===\n", label, outdir))
  cat("\nCell status:\n"); print(acc$cell_status, row.names = FALSE)
  cat("\nPer-trait coverage (model-based vs split vs Mondrian):\n")
  if (nrow(cov_tab)) print(cov_tab, row.names = FALSE) else cat("  (none yet)\n")
  cat("\nDownstream slope check (G8's 5% criterion reported, not required):\n")
  if (nrow(slope_tab)) print(slope_tab, row.names = FALSE) else cat("  (none yet)\n")
  cat("\nPre-registered pair status:\n"); print(acc$pair_status, row.names = FALSE)

  if (!acc$complete) {
    if (nrow(acc$missing_cells)) {
      cat("\nMissing/failed planned cells:\n"); print(acc$missing_cells, row.names = FALSE)
    }
    missing_pairs <- acc$pair_status[!acc$pair_status$has_result, , drop = FALSE]
    if (nrow(missing_pairs)) {
      cat("\nPre-registered pairs with no result yet:\n"); print(missing_pairs, row.names = FALSE)
    }
  }

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

  acc$conv_tab <- conv_tab
  acc$nonconverged <- nonconverged
  acc
}

if (selftest) {
  tmp <- tempfile("mi_realdata_selftest_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  fx <- build_selftest_fixture(tmp)
  printed <- capture.output(acc <- report(tmp, fx$planned, fx$pairs_list, "selftest"))
  writeLines(printed)

  # Fixture is deliberately incomplete (seed 2 has no receipt) AND its one
  # present receipt is deliberately non-converged, so both the
  # completeness logic (must stay FALSE) and the NONCONVERGED print path
  # are exercised in one selftest, and convergence must NOT change
  # completeness or the pair/coverage verdicts.
  stopifnot(
    isFALSE(acc$complete),
    nrow(acc$missing_cells) == 1L, identical(acc$missing_cells$status, "missing"),
    nrow(acc$pair_status) == 1L, isTRUE(acc$pair_status$has_result),
    isTRUE(acc$pair_status$has_ok_slope), isTRUE(acc$pair_status$within_5pct),
    nrow(acc$nonconverged) == 1L,
    any(grepl("^NONCONVERGED synth/synth-mcar-m1 ", printed))
  )
  cat("\nSELFTEST_OK\n")
  quit(status = 0L)
}

outdir <- if (length(args) >= 1L) args[[1L]] else file.path(here_dir, "returned")
acc <- report(outdir, planned_cells(), mi_realdata_pairs, "real data")
if (acc$complete) {
  cat("\nREALDATA_COMPLETE\n")
} else {
  cat("\nNOT complete: see missing cells / pairs above.\n")
}
