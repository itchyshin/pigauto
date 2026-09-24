#!/usr/bin/env Rscript
# Summarise mi_posterior.rds receipts written by 01_run.R: per-trait
# coverage (model-based vs split-conformal vs Mondrian-conformal) and the
# downstream slope table (reference PGLS vs MI-pooled), written to CSV +
# a short Markdown report under <outdir>.
#
# Usage:
#   Rscript 02_summarise.R [outdir]     # default outdir: script/mi_realdata/returned
#   Rscript 02_summarise.R --selftest   # synthetic fixture, prints SELFTEST_OK

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

write_report <- function(outdir, receipts, planned, pairs_list, label) {
  cov_tab <- build_coverage_table(receipts)
  slope_tab <- build_slope_table(receipts)
  conv_tab <- build_convergence_table(receipts)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(cov_tab, file.path(outdir, "coverage_table.csv"), row.names = FALSE)
  utils::write.csv(slope_tab, file.path(outdir, "slope_table.csv"), row.names = FALSE)
  utils::write.csv(conv_tab, file.path(outdir, "convergence_table.csv"), row.names = FALSE)

  n_ok <- sum(vapply(receipts, function(r) identical(r$status, "ok"), logical(1)))
  n_nonconverged <- sum(!is.na(conv_tab$converged) & !conv_tab$converged)
  lines <- c(
    sprintf("# mi-posterior real-data summary (%s)", label),
    "",
    sprintf("%d/%d planned cells have an \"ok\" receipt (%d flagged non-converged; see below).",
            n_ok, nrow(planned), n_nonconverged),
    "",
    "## Per-trait coverage (model-based vs split vs Mondrian)",
    "",
    capture.output(print(cov_tab, row.names = FALSE)),
    "",
    "## Downstream slope check (reference PGLS vs MI-pooled)",
    "",
    capture.output(print(slope_tab, row.names = FALSE)),
    "",
    "## Convergence (max split R-hat / min bulk ESS over Sigma_P, Sigma_E, lambda)",
    "",
    "Descriptive only -- does not gate coverage or REALDATA_COMPLETE (see 03_acceptance.R).",
    "",
    capture.output(print(conv_tab, row.names = FALSE))
  )
  writeLines(lines, file.path(outdir, "summary.md"))
  list(cov_tab = cov_tab, slope_tab = slope_tab, conv_tab = conv_tab)
}

if (selftest) {
  tmp <- tempfile("mi_realdata_selftest_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  fx <- build_selftest_fixture(tmp)
  receipts <- collect_receipts(tmp, fx$planned)
  out <- write_report(tmp, receipts, fx$planned, fx$pairs_list, "selftest")

  stopifnot(
    file.exists(file.path(tmp, "coverage_table.csv")),
    file.exists(file.path(tmp, "slope_table.csv")),
    file.exists(file.path(tmp, "convergence_table.csv")),
    file.exists(file.path(tmp, "summary.md")),
    nrow(out$cov_tab) == 2L,            # 2 traits from the 1 receipt with data
    nrow(out$slope_tab) == 1L,          # 1 pair from the 1 receipt with data
    isTRUE(out$slope_tab$within_5pct[[1]]),
    nrow(out$conv_tab) == 1L, isFALSE(out$conv_tab$converged[[1]])  # fixture is deliberately non-converged
  )
  cat("SELFTEST_OK\n")
  quit(status = 0L)
}

outdir <- if (length(args) >= 1L) args[[1L]] else file.path(here_dir, "returned")
receipts <- collect_receipts(outdir, planned_cells())
out <- write_report(outdir, receipts, planned_cells(), mi_realdata_pairs, "real data")
cat(sprintf("Wrote %s/{coverage_table.csv,slope_table.csv,summary.md}\n", outdir))
print(out$cov_tab)
print(out$slope_tab)
