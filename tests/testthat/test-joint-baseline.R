test_that("joint_mvn_available() detects Rphylopars correctly", {
  # We expect TRUE in dev environments where Rphylopars is installed
  res <- joint_mvn_available()
  expect_type(res, "logical")
  expect_length(res, 1L)
  expect_equal(res, requireNamespace("Rphylopars", quietly = TRUE))
})
