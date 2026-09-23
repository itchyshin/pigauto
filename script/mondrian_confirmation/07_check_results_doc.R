#!/usr/bin/env Rscript
# Usage: Rscript 07_check_results_doc.R [returned_dir] [results_md]
# Stub: sanity-checks that docs/dev-log/mondrian-realdata/results.md exists
# and cites every dataset that has a receipt under
# script/mondrian_confirmation/returned/, before 08_apply_decision_rule.R is
# trusted to read either. Not itself a decision-rule implementation -- see
# 08_apply_decision_rule.R for that.
args <- commandArgs(trailingOnly = TRUE)
returned_dir <- if (length(args) >= 1L) args[[1L]] else "script/mondrian_confirmation/returned"
results_md <- if (length(args) >= 2L) args[[2L]] else "docs/dev-log/mondrian-realdata/results.md"

if (!dir.exists(returned_dir)) {
  cat(sprintf("RESULTS_DOC_CHECK=SKIPPED (no such dir: %s)\n", returned_dir))
  quit(status = 0L)
}
files <- list.files(returned_dir, pattern = "\\.rds$", full.names = FALSE)
if (!length(files)) {
  cat(sprintf("RESULTS_DOC_CHECK=SKIPPED (no receipts in %s)\n", returned_dir))
  quit(status = 0L)
}
if (!file.exists(results_md)) {
  cat(sprintf("RESULTS_DOC_CHECK=FAIL (missing %s)\n", results_md))
  quit(status = 1L)
}
doc <- paste(readLines(results_md, warn = FALSE), collapse = "\n")
datasets <- unique(tools::file_path_sans_ext(files))
missing <- datasets[!vapply(datasets, function(d) grepl(d, doc, fixed = TRUE), logical(1))]
if (length(missing)) {
  cat(sprintf("RESULTS_DOC_CHECK=FAIL (results.md does not mention: %s)\n",
              paste(missing, collapse = ", ")))
  quit(status = 1L)
}
cat("RESULTS_DOC_CHECK=OK\n")
