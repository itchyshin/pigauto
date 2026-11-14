test_that("CV calibration handles one, capped and uneven validation folds", {
  tm <- list(list(name = "x", type = "continuous", latent_cols = 1L))
  run_cal <- function(n, k, method = "cv_folds", seed = 41L) {
    # Non-validation rows have different truth to guard mask use.
    n_rows <- 2L * n + 1L
    idx <- seq.int(2L, 2L * n, by = 2L)
    truth <- matrix(9, n_rows, 1L)
    truth[idx] <- 0
    mu <- matrix(2, n_rows, 1L)
    delta <- matrix(0, n_rows, 1L)
    mask <- matrix(FALSE, n_rows, 1L)
    mask[idx] <- TRUE
    expect_warning(result <- calibrate_gates(
      tm, mu, delta, truth, mask, gate_grid = c(0, 1), gate_cap = 1,
      gate_method = method, gate_cv_folds = k, min_val_cells = 0, seed = seed),
      "95% split-conformal coverage")
    result
  }
  one_cv <- run_cal(1L, 5L)
  one_split <- run_cal(1L, 5L, method = "single_split")
  expect_equal(one_cv, one_split)
  two_cv <- run_cal(2L, 2L)
  capped <- run_cal(2L, 9L)
  expect_equal(two_cv, capped)
  uneven <- run_cal(7L, 3L)
  expect_equal(uneven, run_cal(7L, 3L))
  for (x in list(one_cv, two_cv, capped, uneven)) {
    expect_equal(unname(x$r_cal_bm), 0)
    expect_equal(unname(x$r_cal_gnn), 1)
    expect_equal(unname(x$r_cal_mean), 0)
  }
})

test_that("CV fold count is a finite integer scalar at least two", {
  tm <- list(list(name = "x", type = "continuous", latent_cols = 1L))
  x <- matrix(0, 2, 1)
  for (bad in list(1, 2.5, NA_real_, Inf, c(2, 3), "3")) {
    expect_error(calibrate_gates(tm, x, x, x, matrix(TRUE, 2, 1),
                                 gate_grid = c(0, 1), gate_cap = 1,
                                 gate_method = "cv_folds", gate_cv_folds = bad),
                 "gate_cv_folds.*integer >= 2")
  }
})
