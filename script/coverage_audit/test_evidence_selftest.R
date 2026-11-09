options(covr.record_tests = TRUE)
root <- normalizePath(".")
source(file.path(root, "script", "coverage_audit", "test_evidence.R"))
s <- tempfile(fileext = ".R")
t <- tempfile(fileext = ".R")
o <- tempfile(fileext = ".csv")
writeLines(c("tiny <- function(x) {", "  if (x) 1 else 0", "}"), s)
writeLines(c(
  'testthat::test_that("true branch", { testthat::expect_equal(tiny(TRUE), 1) })',
  'testthat::test_that("false branch", { testthat::expect_equal(tiny(FALSE), 0) })'), t)
x <- covr::file_coverage(s, t)
d <- coverage_test_evidence(x, o)
stopifnot(nrow(d) > 0, all(d$named_test),
          any(grepl("true branch", d$test_evidence)),
          any(grepl("false branch", d$test_evidence)))
attr(x, "tests") <- NULL
stopifnot(inherits(tryCatch(coverage_test_evidence(x, o), error = identity), "error"))
unlink(c(s, t, o))
cat("TEST_EVIDENCE_SMOKE_PASS\n")
