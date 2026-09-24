#!/usr/bin/env Rscript
# Summarise mi_posterior.rds receipts written by 01_run.R: per-trait
# coverage (model-based vs split-conformal vs Mondrian-conformal), interval
# widths like for like (mean AND median of upper - lower, original trait
# scale, for all three methods), the downstream slope table (reference
# PGLS vs MI-pooled, paired) and its per-pair aggregate over cells, written
# to CSV + a short Markdown report under <outdir>. Descriptive only; the
# G8 gate is 03_acceptance.R.
#
# Usage:
#   Rscript 02_summarise.R [outdir]     # default outdir: script/mi_realdata/returned
#   Rscript 02_summarise.R --selftest   # synthetic fixtures, prints SELFTEST_OK

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

cov_cols <- c("dataset", "arm", "seed", "trait", "n_masked", "n_matched",
              "model_coverage", "split_coverage", "mondrian_coverage")
width_cols <- c("dataset", "arm", "seed", "trait", "model_mean_width", "split_mean_width",
                "mondrian_mean_width", "model_median_width", "split_median_width",
                "mondrian_median_width")
slope_cols <- c("dataset", "arm", "seed", "response", "predictor", "n", "ref_slope", "ref_se",
                "mi_slope", "mi_se", "m_used", "m_total", "n_nonfinite", "diff", "rel_diff",
                "diff_ref_se", "se_ratio", "within_5pct", "ref_status", "mi_status")

show <- function(df) {
  if (!nrow(df)) return("(none)")
  num <- vapply(df, is.double, logical(1))
  df[num] <- lapply(df[num], signif, digits = 4)
  capture.output(print(df, row.names = FALSE))
}

write_report <- function(outdir, receipts, planned, pairs_list, label) {
  cov_tab <- build_coverage_table(receipts)
  slope_tab <- build_slope_table(receipts)
  conv_tab <- build_convergence_table(receipts)
  acc <- check_acceptance(receipts, planned, pairs_list)
  pair_tab <- acc$pair_status
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(cov_tab, file.path(outdir, "coverage_table.csv"), row.names = FALSE)
  utils::write.csv(slope_tab, file.path(outdir, "slope_table.csv"), row.names = FALSE)
  utils::write.csv(pair_tab, file.path(outdir, "pair_summary.csv"), row.names = FALSE)
  utils::write.csv(conv_tab, file.path(outdir, "convergence_table.csv"), row.names = FALSE)

  n_ok <- sum(vapply(receipts, function(r) identical(r$status, "ok"), logical(1)))
  n_nonconverged <- sum(!is.na(conv_tab$converged) & !conv_tab$converged)
  lines <- c(
    sprintf("# mi-posterior real-data summary (%s)", label),
    "",
    sprintf("%d/%d planned cells have an \"ok\" receipt (%d flagged non-converged; see below).",
            n_ok, nrow(planned), n_nonconverged),
    "",
    "## Per-trait coverage of masked cells (model-based vs split vs Mondrian)",
    "",
    show(cov_tab[cov_cols]),
    "",
    "## Interval widths, like for like (mean and median of upper - lower, original scale)",
    "",
    show(cov_tab[width_cols]),
    "",
    "## Downstream slope check, per cell (reference PGLS vs MI-pooled, paired)",
    "",
    "diff = MI - reference; rel_diff = diff / reference; diff_ref_se = diff / reference SE.",
    "m_used of m_total completions entered the Rubin pool; n_nonfinite were non-finite on the",
    "analysis scale. A row with m_used < m_total is a selective pool and fails G8 (03_acceptance.R).",
    "",
    show(slope_tab[slope_cols]),
    "",
    "## Downstream slope check, per pair over all cells (5% criterion reported, not gated)",
    "",
    sprintf("Overall: %d of %d pair-cells with a finite relative difference are within 5%%.",
            acc$slope_overall$n_within_5pct, acc$slope_overall$n_pair_cells),
    "",
    show(pair_tab),
    "",
    "## Convergence (max split R-hat / min bulk ESS over Sigma_P, Sigma_E, lambda)",
    "",
    "Descriptive only; does not gate coverage or REALDATA_COMPLETE (see 03_acceptance.R).",
    "",
    show(conv_tab)
  )
  writeLines(lines, file.path(outdir, "summary.md"))
  list(cov_tab = cov_tab, slope_tab = slope_tab, pair_tab = pair_tab, conv_tab = conv_tab, acc = acc)
}

if (selftest) {
  run_variant <- function(variant) {
    tmp <- tempfile(paste0("mi_realdata_selftest_", variant, "_"))
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    fx <- build_selftest_fixture(tmp, variant)
    receipts <- collect_receipts(tmp, fx$planned)
    out <- write_report(tmp, receipts, fx$planned, fx$pairs_list, paste("selftest", variant))
    out$files <- file.exists(file.path(tmp, c("coverage_table.csv", "slope_table.csv", "pair_summary.csv",
                                              "convergence_table.csv", "summary.md")))
    out
  }

  # Positive control: complete fixture, two cells.
  out <- run_variant("pass")
  wcols <- setdiff(width_cols, c("dataset", "arm", "seed", "trait"))
  stopifnot(
    all(out$files),
    nrow(out$cov_tab) == 4L,                                   # 2 cells x 2 traits
    all(vapply(out$cov_tab[wcols], function(v) all(is.finite(v)), logical(1))),
    nrow(out$slope_tab) == 2L,                                 # 1 pair x 2 cells
    identical(out$slope_tab$within_5pct, c(TRUE, FALSE)),      # +1% and +8%
    all(is.finite(out$slope_tab$diff_ref_se)),
    out$pair_tab$n_within_5pct == 1L, out$pair_tab$n_rel_finite == 2L,
    isFALSE(out$pair_tab$all_within_5pct),                     # aggregated, not last-read
    isTRUE(all.equal(out$pair_tab$max_abs_rel_diff, 0.08)),
    identical(out$slope_tab$m_used, c(20L, 20L)), identical(out$slope_tab$n_nonfinite, c(0L, 0L)),
    nrow(out$conv_tab) == 2L, identical(out$conv_tab$converged, c(FALSE, TRUE))
  )
  cat("selftest 02 pass: 4 coverage rows with 6 finite width columns; 2 slope rows, within_5pct = TRUE, FALSE; pair aggregate 1/2 within 5%\n")

  # Missing conformal: reported as NA with its status, never silently filled.
  out <- run_variant("missing_conformal")
  m2 <- out$cov_tab$seed == 2L
  stopifnot(
    all(out$files), nrow(out$cov_tab) == 4L,
    all(is.na(out$cov_tab$split_coverage[m2])), all(is.na(out$cov_tab$split_mean_width[m2])),
    all(out$cov_tab$split_status[m2] == "error"), all(out$cov_tab$split_status[!m2] == "ok")
  )
  cat("selftest 02 missing_conformal: cell 2 split columns NA with status 'error'\n")

  # Selective pooling (2 of 20 completions): shown in the slope table and
  # summary.md, kept out of the per-pair 5% aggregate, and the gate fails.
  out <- run_variant("partial_pool")
  stopifnot(
    all(out$files),
    identical(out$slope_tab$m_used, c(2L, 20L)), identical(out$slope_tab$n_nonfinite, c(18L, 0L)),
    out$pair_tab$n_cells_degraded == 1L, out$pair_tab$n_cells_ok_slope == 1L,
    out$pair_tab$n_rel_finite == 1L, isFALSE(out$acc$complete)
  )
  cat("selftest 02 partial_pool: m_used = 2, 20 and n_nonfinite = 18, 0 in the slope table; 1 degraded pair-cell kept out of the 5% count; gate not complete\n")
  cat("SELFTEST_OK\n")
  quit(status = 0L)
}

outdir <- if (length(args) >= 1L) args[[1L]] else file.path(here_dir, "returned")
receipts <- collect_receipts(outdir, planned_cells())
out <- write_report(outdir, receipts, planned_cells(), mi_realdata_pairs, "real data")
cat(sprintf("Wrote %s/{coverage_table.csv,slope_table.csv,pair_summary.csv,convergence_table.csv,summary.md}\n",
            outdir))
writeLines(readLines(file.path(outdir, "summary.md")))
