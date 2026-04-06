test_that("mask_missing returns logical matrix of correct dimensions", {
  X <- matrix(c(1, NA, 3, NA, 5, 6), nrow = 3)
  m <- mask_missing(X)
  expect_true(is.logical(m))
  expect_equal(dim(m), dim(X))
  expect_equal(m[2, 1], FALSE)
  expect_equal(m[1, 1], TRUE)
})

test_that("make_missing_splits val and test indices are disjoint", {
  X      <- matrix(rnorm(200), nrow = 20)
  splits <- make_missing_splits(X, missing_frac = 0.30, seed = 42)
  expect_equal(length(intersect(splits$val_idx, splits$test_idx)), 0)
})

test_that("make_missing_splits produces approximately the right missing fraction", {
  X      <- matrix(rnorm(400), nrow = 20)
  splits <- make_missing_splits(X, missing_frac = 0.25, seed = 1)
  total_missing <- length(splits$val_idx) + length(splits$test_idx)
  actual_frac   <- total_missing / (20 * 20)
  expect_lt(abs(actual_frac - 0.25), 0.02)
})

test_that("make_missing_splits is deterministic given seed", {
  X  <- matrix(rnorm(100), nrow = 10)
  s1 <- make_missing_splits(X, seed = 7)
  s2 <- make_missing_splits(X, seed = 7)
  expect_equal(s1$val_idx,  s2$val_idx)
  expect_equal(s1$test_idx, s2$test_idx)
})

test_that("make_missing_splits mask is TRUE for non-missing cells", {
  X      <- matrix(rnorm(100), nrow = 10)
  splits <- make_missing_splits(X, missing_frac = 0.20, seed = 3)
  expect_false(any(splits$mask[splits$val_idx]))
  expect_false(any(splits$mask[splits$test_idx]))
  n_observed <- sum(splits$mask)
  expect_equal(n_observed, 100 - length(splits$val_idx) - length(splits$test_idx))
})
