# script/tests-rubin/test-lib-cellcov.R
# Gate G-S1b for rubin_cell_intervals(). 20,000 replications of M = 20 draws from N(mu, 1) plus one
# fresh draw from the same distribution (the "true" missing value the interval is meant to cover).
# This is the classical iid-normal prediction-interval setup: Xbar +/- t_{M-1} * s * sqrt(1 + 1/M)
# has exact nominal coverage of an independent new draw, so rubin_cell_intervals()'s U = 0 / B-only
# formula (T = (1+1/M) B, df = M-1) should recover ~95% coverage here.
#
# Reps live as the K = 20,000 columns of one M x K draws matrix so the whole simulation is one
# vectorised call to rubin_cell_intervals(), not a 20,000-iteration R loop.

source(normalizePath(file.path(testthat::test_path(), "..", "rubin_lib.R")))

testthat::test_that("rubin_cell_intervals achieves ~95% coverage of a fresh N(mu,1) draw", {
  set.seed(20260924)
  M <- 20L
  K <- 20000L
  mu <- 2.3  # arbitrary, coverage is location-invariant

  draws <- matrix(stats::rnorm(M * K, mean = mu, sd = 1), nrow = M, ncol = K)
  new_draw <- stats::rnorm(K, mean = mu, sd = 1)

  res <- rubin_cell_intervals(draws, conf = 0.95)
  testthat::expect_equal(nrow(res), K)
  testthat::expect_equal(res$df, rep(M - 1L, K))

  covered <- new_draw >= res$lower & new_draw <= res$upper
  coverage <- mean(covered)

  # Percentile interval on the same draws, for comparison only (not asserted tightly; expected
  # to undercover badly at M = 20, historically ~0.87).
  pct <- apply(draws, 2, stats::quantile, probs = c(0.025, 0.975))
  covered_pct <- new_draw >= pct[1, ] & new_draw <= pct[2, ]
  coverage_pct <- mean(covered_pct)

  message(sprintf(
    "G-S1b: rubin_cell_intervals coverage = %.4f (target [0.945, 0.955]); percentile-interval coverage = %.4f (informational, ~0.87 expected)",
    coverage, coverage_pct))

  testthat::expect_gte(coverage, 0.945)
  testthat::expect_lte(coverage, 0.955)
})

testthat::test_that("rubin_cell_intervals flags B = 0 cells with a zero-width interval, not an error", {
  draws <- matrix(rep(1.5, 20), nrow = 20, ncol = 1)  # all M draws identical
  res <- rubin_cell_intervals(draws)
  testthat::expect_true(res$zero_width)
  testthat::expect_equal(res$B, 0)
  testthat::expect_equal(res$lower, res$upper)
  testthat::expect_equal(res$lower, 1.5)
})
