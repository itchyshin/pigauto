# Oracles for gates G5 and G6 of the lambda-default ledger (rewritten 2026-09-23 after the first
# versions were malformed: G5's EXPECT spanned two output lines, G6 had an R escape error).
# Usage: Rscript docs/dev-log/lambda-default/check_gates_s4_s5.R G5|G6
suppressMessages(devtools::load_all(quiet = TRUE))
which <- commandArgs(trailingOnly = TRUE)[1]
if (identical(which, "G5")) {
  # Tests 1-3 of the per-type guard (the discrete-trait guards) must be byte-identical to origin/main;
  # only test 4 (continuous routing) was amended. Then the whole file must pass.
  head3 <- function(txt) { i <- grep("^test_that\\(", txt); txt[seq_len(i[4] - 1L)] }
  old <- system2("git", c("show", "origin/main:tests/testthat/test-lambda-per-type.R"), stdout = TRUE)
  new <- readLines("tests/testthat/test-lambda-per-type.R")
  same <- identical(head3(old), head3(new))
  df <- as.data.frame(testthat::test_file("tests/testthat/test-lambda-per-type.R", reporter = "silent"))
  cat(sprintf("tests1to3_unchanged=%s pertype_fail=%d pass=%d\n", same, sum(df$failed), sum(df$passed)))
  if (same && sum(df$failed) == 0) cat("G5_OK\n")
} else if (identical(which, "G6")) {
  d <- function(f) { a <- formals(f)$lambda_mode; if (is.character(a)) a[1] else eval(a)[1] }
  vals <- c(fit_baseline = d(fit_baseline), fit_pigauto = d(fit_pigauto), impute = d(impute),
            multi_impute = d(multi_impute))
  print(vals)
  old_fallback <- any(grepl('baseline_arg("lambda_mode", "fixed_1")', readLines("R/multi_impute_trees.R"),
                            fixed = TRUE))
  cat("multi_impute_trees old fallback present:", old_fallback, "\n")
  if (all(vals == "estimate") && !old_fallback) cat("G6_OK\n")
}
