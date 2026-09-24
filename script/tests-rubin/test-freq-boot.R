# script/tests-rubin/test-freq-boot.R
# Gate G-S2c for mi_freq_A() (parametric bootstrap) and mi_freq_B() (fixed-parameter draws), on the
# same n = 60 cell as G-S2a, M = 20.

source(normalizePath(file.path(testthat::test_path(), "..", "rubin_freq.R")))
source(normalizePath(file.path(testthat::test_path(), "..", "campaign_gnn_off_lib.R")))

testthat::test_that("mi_freq_A: parametric-bootstrap params vary, refits (mostly) succeed", {
  cell <- make_cell("types_mixed", n = 60, seed = 1, lambda = 0.7, rho = 0.5)

  t0 <- Sys.time()
  res <- mi_freq_A(cell, M = 20L)
  wall_a <- as.numeric(Sys.time() - t0, units = "secs")
  cat(sprintf("G-S2c mi_freq_A wall time: %.2f s (M = 20)\n", wall_a))

  testthat::expect_true(res$n_fail <= 2L)
  cat(sprintf("G-S2c mi_freq_A n_fail: %d\n", res$n_fail))

  ok_star <- Filter(Negate(is.null), res$pars_star)
  testthat::expect_gt(length(ok_star), 0L)

  lambdas <- vapply(ok_star, function(x) x$lambda, numeric(1))
  testthat::expect_true(all(is.finite(lambdas)))
  Sigma_p_finite <- vapply(ok_star, function(x) all(is.finite(x$Sigma_p)), logical(1))
  testthat::expect_true(all(Sigma_p_finite))
  Sigma_e_finite <- vapply(ok_star, function(x) is.null(x$Sigma_e) || all(is.finite(x$Sigma_e)), logical(1))
  testthat::expect_true(all(Sigma_e_finite))

  sd_lambda <- stats::sd(lambdas)
  cat(sprintf("G-S2c sd(lambda*): %.4f\n", sd_lambda))
  testthat::expect_gt(sd_lambda, 0)
})

testthat::test_that("mi_freq_B: datasets differ only in missing block cells", {
  cell <- make_cell("types_mixed", n = 60, seed = 1, lambda = 0.7, rho = 0.5)
  block_traits <- default_block_traits(cell)
  df_miss <- cell$df_miss

  t0 <- Sys.time()
  res <- mi_freq_B(cell, M = 20L)
  wall_b <- as.numeric(Sys.time() - t0, units = "secs")
  cat(sprintf("G-S2c mi_freq_B wall time: %.2f s (M = 20)\n", wall_b))

  testthat::expect_length(res$datasets, 20L)

  other_cols <- setdiff(names(df_miss), block_traits)
  other_ok <- vapply(res$datasets, function(d) identical(d[, other_cols], df_miss[, other_cols]), logical(1))
  testthat::expect_true(all(other_ok))

  # observed cells (non-NA in df_miss) are bit-identical in every dataset
  obs_mask <- !is.na(as.matrix(df_miss[, block_traits]))
  obs_identical <- vapply(res$datasets, function(d) {
    dm <- as.matrix(d[, block_traits])
    identical(dm[obs_mask], as.matrix(df_miss[, block_traits])[obs_mask])
  }, logical(1))
  testthat::expect_true(all(obs_identical))

  miss_mask <- is.na(as.matrix(df_miss[, block_traits]))
  first_missing <- as.matrix(res$datasets[[1]][, block_traits])[miss_mask]
  differs <- vapply(res$datasets[-1], function(d) {
    any(as.matrix(d[, block_traits])[miss_mask] != first_missing)
  }, logical(1))
  testthat::expect_true(all(differs))
})
