# script/tests-rubin/test-lib-pool.R
# Gate G-S1a for rubin_pool(). Part (a): a fixed hand-worked toy checked to 1e-12. Part (b): an
# oracle check against mice::pool() on 20 bootstrap-perturbed lm fits.

source(normalizePath(file.path(testthat::test_path(), "..", "rubin_lib.R")))

testthat::test_that("rubin_pool matches a hand-worked M = 5 toy to 1e-12", {
  # q = 1:5, u = c(0.1, 0.2, 0.3, 0.4, 0.5); M = 5.
  #   qbar = mean(q)           = 3
  #   W    = mean(u)           = 0.3
  #   B    = var(q)            = sum((q - 3)^2) / 4 = (4+1+0+1+4)/4 = 2.5
  #   T    = W + (1 + 1/5) B   = 0.3 + 1.2 * 2.5     = 3.3
  #   riv  = (1 + 1/5) B / W   = 1.2 * 2.5 / 0.3      = 10
  #   lambda = (1+1/5)B / T    = 3 / 3.3              = 10/11 = 0.909090909090909...
  #   df (df_com = Inf)  = (M-1) / lambda^2 = 4 / (10/11)^2 = 4 * 121/100 = 4.84
  #   fmi  = (riv + 2/(df+3)) / (riv+1) = (10 + 2/7.84) / 11 = 0.9322820037105752
  #   se   = sqrt(T) = sqrt(3.3) = 1.8165902124584949
  #   crit = qt(0.975, 4.84)  = 2.5963527612534527
  #   lower/upper = qbar -/+ crit * se
  q <- 1:5
  u <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  r <- rubin_pool(q, u, df_com = Inf, conf = 0.95)

  testthat::expect_equal(r$estimate, 3, tolerance = 1e-12)
  testthat::expect_equal(r$W, 0.3, tolerance = 1e-12)
  testthat::expect_equal(r$B, 2.5, tolerance = 1e-12)
  testthat::expect_equal(r$T, 3.3, tolerance = 1e-12)
  testthat::expect_equal(r$riv, 10, tolerance = 1e-12)
  testthat::expect_equal(r$lambda, 10 / 11, tolerance = 1e-12)
  testthat::expect_equal(r$df, 4.84, tolerance = 1e-12)
  testthat::expect_equal(r$fmi, 0.93228200371057524, tolerance = 1e-12)
  testthat::expect_equal(r$se, 1.8165902124584949, tolerance = 1e-12)
  testthat::expect_equal(r$lower, -1.71650901418261, tolerance = 1e-10)
  testthat::expect_equal(r$upper, 7.71650901418261, tolerance = 1e-10)
})

testthat::test_that("rubin_pool reproduces mice::pool() on 20 bootstrap lm fits", {
  testthat::skip_if_not_installed("mice")
  set.seed(20260924)
  n <- 60
  x <- stats::rnorm(n)
  y <- 1.5 + 2 * x + stats::rnorm(n)
  M <- 20
  fits <- lapply(seq_len(M), function(i) {
    idx <- sample(n, n, replace = TRUE)
    stats::lm(y[idx] ~ x[idx])
  })

  p <- mice::pool(mice::as.mira(fits))$pooled

  for (j in seq_len(nrow(p))) {
    q <- vapply(fits, function(f) stats::coef(f)[j], numeric(1))
    u <- vapply(fits, function(f) stats::vcov(f)[j, j], numeric(1))
    df_com <- stats::df.residual(fits[[1]])
    r <- rubin_pool(q, u, df_com = df_com)

    testthat::expect_equal(r$estimate, p$estimate[j], tolerance = 1e-8,
                            label = paste("estimate, term", p$term[j]))
    testthat::expect_equal(r$T, p$t[j], tolerance = 1e-8,
                            label = paste("T, term", p$term[j]))
    testthat::expect_equal(r$df, p$df[j], tolerance = 1e-8,
                            label = paste("df, term", p$term[j]))
  }
})
