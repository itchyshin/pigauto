# Run only the campaign's public-contract regression file. A skip is not a
# pass: this runner exits non-zero if the file contains any skipped test.
suppressPackageStartupMessages({
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools is required to load the package for contract checks")
  }
  if (!requireNamespace("testthat", quietly = TRUE)) {
    stop("testthat is required to run contract checks")
  }
})

devtools::load_all(".", quiet = TRUE)
test_file <- "tests/testthat/test-validation-contract.R"
result <- testthat::test_file(test_file, reporter = "summary")
expectations <- unlist(lapply(result, `[[`, "results"), recursive = FALSE,
                       use.names = FALSE)
if (!length(expectations)) {
  stop("CONTRACT_FAIL: test file produced no expectations", call. = FALSE)
}
classes <- lapply(expectations, class)
skipped <- sum(vapply(classes, function(x) "expectation_skip" %in% x, logical(1)))
failed <- sum(vapply(classes, function(x) "expectation_failure" %in% x, logical(1)))
errored <- sum(vapply(classes, function(x) "expectation_error" %in% x, logical(1)))
if (skipped > 0L || failed > 0L || errored > 0L) {
  stop(sprintf("CONTRACT_FAIL: %d failure(s), %d error(s), %d skipped test(s)",
               failed, errored, skipped), call. = FALSE)
}
cat(sprintf("CONTRACT_PASS: %d tests, %d expectations, 0 failures/errors/skips\n",
            length(result), length(expectations)))
