#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork = TRUE)
required <- c(
  "docs/dev-log/gnn-incremental/AFTER-TASK-20261020.md",
  "docs/dev-log/gnn-incremental/PLAN-VS-ACTUAL-20261020.md",
  "docs/dev-log/gnn-incremental/RESULTS-20261020.md",
  "vignettes/validation-results.Rmd",
  "script/gnn_incremental/results/production-20261019-v1-summary/summaries/paired_dataset_metrics.csv",
  "script/gnn_incremental/results/production-20261019-v1-summary/summaries/method_accounting.csv",
  ".unlazy/gnn-incremental/GATES.md"
)
missing <- required[!file.exists(file.path(root, required))]
if (length(missing)) stop("Missing closeout artifacts: ", paste(missing, collapse = ", "))

metrics <- read.csv(file.path(root, required[5L]), stringsAsFactors = FALSE)
accounting <- read.csv(file.path(root, required[6L]), stringsAsFactors = FALSE)
expected_scenarios <- c("nonlinear_primary", "interactive_secondary",
  "nonlinear_low_phylo", "linear_recovery", "null_covariate")
if (length(setdiff(expected_scenarios, unique(metrics$scenario_id))) ||
    length(setdiff(expected_scenarios, unique(accounting$scenario_id)))) {
  stop("Production summaries omit one or more frozen scenarios.")
}
counts <- aggregate(replicate ~ scenario_id + phylo_signal_gate,
  data = unique(metrics[c("scenario_id", "replicate", "phylo_signal_gate")]),
  FUN = function(x) length(unique(x)))
if (nrow(counts) != 10L || any(counts$replicate != 500L)) {
  stop("Expected 500 unique datasets per scenario and gate policy.")
}
if (nrow(accounting) != 60L || any(accounting$requested != 500L) ||
    any(accounting$attempted != 500L) || any(accounting$successful != 500L) ||
    any(accounting$failed != 0L) || any(accounting$nonconverged != 0L)) {
  stop("Method-attempt accounting is incomplete or contains a failure.")
}

report <- readLines(file.path(root, required[3L]), warn = FALSE)
article <- readLines(file.path(root, required[4L]), warn = FALSE)
after <- readLines(file.path(root, required[1L]), warn = FALSE)
plan <- readLines(file.path(root, required[2L]), warn = FALSE)
all_text <- paste(c(report, article, after, plan), collapse = "\n")
for (claim in c("exactly represented", "partly observed", "does not establish")) {
  if (!grepl(claim, all_text, fixed = TRUE)) stop("Required limitation missing: ", claim)
}
headers <- c("## 1. Goal", "## 2. Implemented", "## 3a. Decisions and Rejected Alternatives",
  "## 4. Files Touched", "## 5. Checks Run", "## 6. Tests of the Tests",
  "## 7a. Issue Ledger", "## 8. Consistency Audit", "## 9. What Did Not Go Smoothly",
  "## 10. Known Residuals", "## 11. Team Learning", "## 12. Cross-Product Coverage")
if (length(setdiff(headers, after))) stop("After-task report lacks required sections.")
if (grepl("<fill|<final result", paste(after, collapse = "\n"), ignore.case = TRUE)) {
  stop("After-task report contains an unfinished placeholder.")
}
cat("CLOSEOUT_VALIDATION_PASS\n")
