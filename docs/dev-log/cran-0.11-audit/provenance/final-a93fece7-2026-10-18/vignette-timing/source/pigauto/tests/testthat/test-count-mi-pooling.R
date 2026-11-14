# Issue #40 regression tests: count / proportion / zi_count MI pooling
# must not be pulled to the moon by a single outlier dropout draw.
#
# The underlying bug: when n_imputations > 1, pigauto's `pool_imputations`
# used to compute rowMeans across the M decoded draws. For count traits
# the decoding is expm1(latent * sd + mean), so a dropout-noisy latent
# ~+5 SD gets amplified to thousands, contaminating the whole cell's
# pooled estimate.  Fix: median-pool (the default) or opt back into
# mean-pool via pool_method = "mean".

test_that("pool_imputations median is robust to expm1-amplified count outliers", {
  # Synthesize 5 decode results for 3 species x 1 count trait.
  # Four of them agree on values near 2; the 5th has a dropout-amplified
  # outlier at 1000. Mean would be ~202 (moon); median should be near 2.
  tm_count <- list(type = "count", name = "n")
  draws <- list(
    list(imputed = data.frame(n = c(2L, 3L, 2L)),      probabilities = list()),
    list(imputed = data.frame(n = c(2L, 3L, 2L)),      probabilities = list()),
    list(imputed = data.frame(n = c(3L, 2L, 2L)),      probabilities = list()),
    list(imputed = data.frame(n = c(2L, 3L, 3L)),      probabilities = list()),
    # outlier draw on species 1 only (as if dropout made latent huge)
    list(imputed = data.frame(n = c(1000L, 3L, 2L)),   probabilities = list())
  )
  # Fake latent_runs: N(0, 1) 3x1 each, not used for count pooling
  latent_runs <- lapply(seq_along(draws), function(i) matrix(stats::rnorm(3), 3, 1))

  out_med  <- pigauto:::pool_imputations(draws, latent_runs,
                                          list(tm_count),
                                          pool_method = "median")
  out_mean <- pigauto:::pool_imputations(draws, latent_runs,
                                          list(tm_count),
                                          pool_method = "mean")

  # Median is unaffected by the outlier on species 1.
  expect_equal(out_med$imputed$n[1], 2L)
  # Mean is pulled to ~202 by the outlier.
  expect_gt(out_mean$imputed$n[1], 100L)

  # Species 2 and 3 have no outliers and should agree closely across methods.
  expect_lte(abs(as.integer(out_med$imputed$n[2]) -
                  as.integer(out_mean$imputed$n[2])), 1L)
  expect_lte(abs(as.integer(out_med$imputed$n[3]) -
                  as.integer(out_mean$imputed$n[3])), 1L)
})

test_that("pool_imputations median-robustness applies to proportion traits", {
  # plogis-decoded proportion: a dropout-noisy latent ~+5 SD gives
  # plogis(5) ~ 0.993. One outlier can swing the mean toward 1.
  tm_prop <- list(type = "proportion", name = "p")
  draws <- list(
    list(imputed = data.frame(p = c(0.50, 0.40, 0.55)), probabilities = list()),
    list(imputed = data.frame(p = c(0.48, 0.42, 0.53)), probabilities = list()),
    list(imputed = data.frame(p = c(0.52, 0.38, 0.56)), probabilities = list()),
    list(imputed = data.frame(p = c(0.49, 0.41, 0.54)), probabilities = list()),
    # outlier pushes species 1's decoded value to 0.99
    list(imputed = data.frame(p = c(0.99, 0.40, 0.55)), probabilities = list())
  )
  latent_runs <- lapply(seq_along(draws), function(i) matrix(stats::rnorm(3), 3, 1))

  out_med  <- pigauto:::pool_imputations(draws, latent_runs,
                                          list(tm_prop),
                                          pool_method = "median")
  out_mean <- pigauto:::pool_imputations(draws, latent_runs,
                                          list(tm_prop),
                                          pool_method = "mean")

  # Median stays near 0.50 (robust to the outlier)
  expect_lt(abs(out_med$imputed$p[1] - 0.50), 0.05)
  # Mean gets pulled up by the 0.99 outlier
  expect_gt(out_mean$imputed$p[1], 0.55)
})

test_that("pool_imputations continuous traits pool by mean (unchanged)", {
  tm_cont <- list(type = "continuous", name = "x")
  draws <- list(
    list(imputed = data.frame(x = c(1.0, 2.0)), probabilities = list()),
    list(imputed = data.frame(x = c(1.2, 2.2)), probabilities = list()),
    list(imputed = data.frame(x = c(0.8, 1.8)), probabilities = list())
  )
  latent_runs <- lapply(seq_along(draws), function(i) matrix(stats::rnorm(2), 2, 1))
  # default pool_method = "median" still means-pools for continuous
  out <- pigauto:::pool_imputations(draws, latent_runs, list(tm_cont))
  expect_equal(out$imputed$x, c(1.0, 2.0), tolerance = 1e-9)
})

test_that("pool_method = 'mean' restores pre-v0.9.2 behaviour byte-for-byte", {
  # For count traits with no outliers, mean and median agree to within 1.
  tm_count <- list(type = "count", name = "n")
  set.seed(42L)
  draws <- lapply(seq_len(5), function(i) {
    list(imputed = data.frame(n = as.integer(stats::rpois(4, 3))),
         probabilities = list())
  })
  latent_runs <- lapply(seq_along(draws), function(i) matrix(stats::rnorm(4), 4, 1))

  out_mean <- pigauto:::pool_imputations(draws, latent_runs,
                                          list(tm_count),
                                          pool_method = "mean")
  # Mean path compiles & returns sane integers
  expect_type(out_mean$imputed$n, "integer")
  expect_true(all(out_mean$imputed$n >= 0L))
})
