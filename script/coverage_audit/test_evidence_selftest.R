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
stopifnot(nrow(d) > 0, !any(d$named_test),
          any(grepl("true branch", d$collective_test_evidence)),
          any(grepl("false branch", d$collective_test_evidence)))
# A single call to the parent must not claim both branches on a shared line.
stopifnot(all(d$single_test_evidence == ""))
writeLines(c("plain <- function(x) {", "  x + 1", "}"), s)
writeLines('testthat::test_that("plain plus", { testthat::expect_equal(plain(1), 2) })', t)
y <- covr::file_coverage(s, t)
positive <- coverage_test_evidence(y, o)
stopifnot(nrow(positive) > 0, all(positive$named_test),
          all(grepl("plain plus", positive$single_test_evidence)))
attr(x, "tests") <- NULL
stopifnot(inherits(tryCatch(coverage_test_evidence(x, o), error = identity), "error"))
unlink(c(s, t, o))
cat("TEST_EVIDENCE_SMOKE_PASS\n")
