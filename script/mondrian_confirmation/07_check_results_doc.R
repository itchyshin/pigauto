#!/usr/bin/env Rscript
# Usage: Rscript 07_check_results_doc.R [returned_dir] [results_md]
# Sanity-checks that docs/dev-log/mondrian-realdata/results.md exists and
# cites every dataset that has a receipt under
# script/mondrian_confirmation/returned/, then regenerates results_table.csv
# (next to results_md) from the same receipts using
# 12_build_results_doc.R's own aggregation and prints ROWS_MATCH only if that
# regeneration is byte-identical to the results_table.csv already on disk --
# i.e. the doc was not hand-edited or built from a stale RESULTS_ROOT. This
# is a consistency check, not a decision-rule implementation -- see
# 08_apply_decision_rule.R for that.
args <- commandArgs(trailingOnly = TRUE)
returned_dir <- if (length(args) >= 1L) args[[1L]] else "script/mondrian_confirmation/returned"
results_md <- if (length(args) >= 2L) args[[2L]] else "docs/dev-log/mondrian-realdata/results.md"

if (!dir.exists(returned_dir)) {
  cat(sprintf("RESULTS_DOC_CHECK=SKIPPED (no such dir: %s)\n", returned_dir))
  quit(status = 0L)
}
# Cell receipts live one level down: returned_dir/<dataset>-<arm>-m<seed>/
# {mask_receipt,mondrian,split}.rds (01_run_masked_confirmation.R's layout).
cell_dirs <- list.dirs(returned_dir, recursive = FALSE, full.names = FALSE)
cell_dirs <- cell_dirs[grepl("^.+-(mcar|structured)-m[0-9]+$", cell_dirs)]
cell_dirs <- cell_dirs[vapply(cell_dirs, function(d) {
  all(file.exists(file.path(returned_dir, d, c("mask_receipt.rds", "mondrian.rds", "split.rds"))))
}, logical(1))]
if (!length(cell_dirs)) {
  cat(sprintf("RESULTS_DOC_CHECK=SKIPPED (no receipts in %s)\n", returned_dir))
  quit(status = 0L)
}
if (!file.exists(results_md)) {
  cat(sprintf("RESULTS_DOC_CHECK=FAIL (missing %s)\n", results_md))
  quit(status = 1L)
}
doc <- paste(readLines(results_md, warn = FALSE), collapse = "\n")
datasets <- unique(sub("-(mcar|structured)-m[0-9]+$", "", cell_dirs))
missing <- datasets[!vapply(datasets, function(d) grepl(d, doc, fixed = TRUE), logical(1))]
if (length(missing)) {
  cat(sprintf("RESULTS_DOC_CHECK=FAIL (results.md does not mention: %s)\n",
              paste(missing, collapse = ", ")))
  quit(status = 1L)
}
cat("RESULTS_DOC_CHECK=OK\n")

# ---------------------------------------------------------------------------
# Regenerate results_table.csv from the same receipts, using
# 12_build_results_doc.R's own build_results_table()/write_results_csv() (a
# single source of truth shared with the script that wrote results.md in the
# first place), and compare byte-for-byte against the CSV already on disk.
# ---------------------------------------------------------------------------
results_csv <- file.path(dirname(results_md), "results_table.csv")
old_source_only <- Sys.getenv("PIGAUTO_12_SOURCE_ONLY", unset = NA)
Sys.setenv(PIGAUTO_12_SOURCE_ONLY = "1")
source("script/mondrian_confirmation/12_build_results_doc.R")
if (is.na(old_source_only)) Sys.unsetenv("PIGAUTO_12_SOURCE_ONLY") else Sys.setenv(PIGAUTO_12_SOURCE_ONLY = old_source_only)

if (!file.exists(results_csv)) {
  cat(sprintf("ROWS_MISMATCH (no such file: %s)\n", results_csv))
  quit(status = 1L)
}
built <- build_results_table(returned_dir)
tmp_csv <- tempfile(fileext = ".csv")
write_results_csv(built$table1, tmp_csv)
rows_match <- identical(readLines(tmp_csv), readLines(results_csv))
unlink(tmp_csv)
if (rows_match) {
  cat("ROWS_MATCH\n")
} else {
  cat("ROWS_MISMATCH\n")
  quit(status = 1L)
}
