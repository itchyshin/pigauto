test_that("Bayesian lambda weights remain normalized for real likelihood scales", {
  set.seed(603)
  tree <- ape::rtree(12)
  correlation <- stats::cov2cor(ape::vcv.phylo(tree))
  grid <- c(.005, .2, .5, .8, .995)
  for (values in list(sin(seq_len(12)), rep(1, 12),
                     rep(c(1, -1), 6) * 1e-300,
                     rep(c(1, -1), 6) * 1e308)) {
    result <- bayes_lambda_for_col(values, correlation, lambda_grid = grid)
    expect_identical(result$lambda_grid, grid)
    expect_true(all(is.finite(result$weights)))
    expect_true(all(result$weights >= 0 & result$weights <= 1))
    expect_equal(sum(result$weights), 1)
    expect_true(any(result$weights > 0))
    expect_equal(result$lambda_post_mean, sum(grid * result$weights))
    expect_true(result$lambda_post_mean >= min(grid) &&
                  result$lambda_post_mean <= max(grid))
    expect_true(is.finite(result$lambda_post_entropy))
    expect_true(result$lambda_post_entropy >= 0 &&
                  result$lambda_post_entropy <= log(length(grid)) + 1e-12)
  }
})

test_that("Bayesian lambda keeps genuine sparse and empty-grid fallbacks", {
  expected <- list(lambda_grid = 1, weights = 1,
    lambda_post_mean = 1, lambda_post_entropy = 0)
  expect_identical(bayes_lambda_for_col(seq_len(9), diag(9)), expected)
  expect_identical(bayes_lambda_for_col(seq_len(12), diag(12),
                                      lambda_grid = numeric()), expected)
})
