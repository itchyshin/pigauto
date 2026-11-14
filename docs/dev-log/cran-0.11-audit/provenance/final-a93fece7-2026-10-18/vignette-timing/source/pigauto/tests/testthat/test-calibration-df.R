# Issue #212: calibration_df() rejects invalid inputs before a silent or
# downstream crash.

test_that("probabilities outside [0, 1] error", {
  expect_error(
    calibration_df(c(1, 0), c(1.5, -0.2), 10),
    "\\[0, 1\\]"
  )
})

test_that("NA pairs warn with the number dropped", {
  expect_warning(
    out <- calibration_df(c(1, 0, 1), c(.5, NA, .7), 2),
    "dropped"
  )
  expect_equal(sum(out$n), 2L)
})

test_that("length mismatch errors", {
  expect_error(
    calibration_df(c(1, 0, 1), c(.5, .2), 2),
    "length"
  )
})

test_that("factor truth errors", {
  expect_error(
    calibration_df(factor(c("a", "b", "a")), c(.5, .2, .7), 2),
    "0/1"
  )
})

test_that("non-binary numeric truth errors", {
  expect_error(
    calibration_df(c(2, 0, 1), c(.5, .2, .7), 2),
    "0/1"
  )
})

test_that("n_bins must be a positive integer", {
  expect_error(
    calibration_df(c(1, 0), c(.5, .2), 0),
    "n_bins"
  )
  expect_error(
    calibration_df(c(1, 0), c(.5, .2), "a"),
    "n_bins"
  )
})

test_that("empty input errors", {
  expect_error(
    calibration_df(numeric(0), numeric(0)),
    "length"
  )
})

test_that("documented example still returns one row per occupied bin", {
  out <- calibration_df(c(0, 0, 1, 1), c(.1, .3, .7, .9), n_bins = 2)
  expect_s3_class(out, "data.frame")
  expect_equal(sum(out$n), 4L)
  expect_true(all(out$obs_freq >= 0 & out$obs_freq <= 1))
})
