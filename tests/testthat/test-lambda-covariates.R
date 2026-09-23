# Tests for S3 of the joint-lambda-default lane (docs/dev-log/2026-09-22-joint-lambda-alignment.md):
# Pagel's lambda threaded through the covariate-aware BM path
# (bm_impute_col_with_cov(), build_pagel_nll_cache(y, R, X = )).

# ---- lambda = 1 bit-identical ----------------------------------------------

test_that("[lambda-cov] lambda = 1 is bit-identical to the pre-lambda bm_impute_col_with_cov()", {
  # Fixture and reference values captured from bm_impute_col_with_cov() BEFORE
  # the lambda argument was added (no lambda arg existed at all, i.e. today's
  # implicit lambda = 1). See docs/dev-log/lambda-default/S3-covariates-report.md
  # for the generating script.
  set.seed(20260922L)
  n <- 80L
  tree <- ape::rcoal(n)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]

  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  X <- cbind(1, x1, x2)
  beta_true <- c(0, 0.8, -0.5)
  phylo_e <- as.numeric(t(chol(R)) %*% rnorm(n))
  y <- as.numeric(X %*% beta_true) + phylo_e

  set.seed(123L)
  na_idx <- sample(n, size = round(0.3 * n))
  y[na_idx] <- NA
  miss_idx <- which(is.na(y))

  out <- pigauto:::bm_impute_col_with_cov(y, X, R, nugget = 1e-6,
                                            ridge = 0.0, lrt_threshold = 0.02,
                                            lambda = 1.0)

  ref_mu6 <- c(0.452235020019112, -1.58303316206551, -0.201014777091942,
               -1.3826846565077, 0.852319507648334, 1.10907550612489)
  ref_se6 <- c(0.122583429624254, 0.155041063942026, 0.419359208216456,
               0.175459395566474, 0.199312958058435, 0.0601000761508676)
  ref_beta <- c(0.117578685386471, 0.797173056236766, -0.497680387053763)

  expect_equal(out$mu[miss_idx][1:6], ref_mu6, tolerance = 1e-8)
  expect_equal(out$se[miss_idx][1:6], ref_se6, tolerance = 1e-8)
  expect_equal(out$beta_hat, ref_beta, tolerance = 1e-8)
  expect_true(out$used_cov)
  # No new field leaks in at lambda = 1.
  expect_false("lambda_hat" %in% names(out))
})

# ---- X = NULL cache is untouched --------------------------------------------

test_that("[lambda-cov] build_pagel_nll_cache(y, R) with X = NULL is unchanged", {
  # Reference value captured from the pre-edit closure (no X argument existed).
  set.seed(555L)
  n2 <- 30L
  tree2 <- ape::rcoal(n2)
  R2 <- stats::cov2cor(ape::vcv.phylo(tree2))
  R2 <- R2[tree2$tip.label, tree2$tip.label]
  y2 <- as.numeric(t(chol(R2)) %*% rnorm(n2))
  y2[c(2L, 7L, 15L, 22L)] <- NA

  cache <- pigauto:::build_pagel_nll_cache(y2, R2)
  expect_equal(cache$nll(0.4), -16.2152757282217, tolerance = 1e-10)

  # X = NULL explicitly is the same as omitting it.
  cache_explicit <- pigauto:::build_pagel_nll_cache(y2, R2, X = NULL)
  expect_equal(cache_explicit$nll(0.4), cache$nll(0.4), tolerance = 1e-12)
})

# ---- lambda recovery with covariates ----------------------------------------

# Default n = 1200: at n = 300 a single lambda_hat draw is noisy (MAE ~ 0.14,
# see the S3 report) -- fine for the bias-gated recovery test below (which
# passes n = 300 explicitly and averages over many replicates), but not for
# the single-draw "covariate effect survives" test further down, which keeps
# this default.
.sim_lambda_cov <- function(seed, n = 1200L, lambda_true = 0.3, cov_beta = c(0.8, -0.5)) {
  set.seed(seed)
  tree <- ape::rcoal(n)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  X <- cbind(1, x1, x2)
  beta_true <- c(0, cov_beta)
  phylo_e <- as.numeric(t(chol(R)) %*% rnorm(n))
  iid_e <- stats::rnorm(n)
  y <- as.numeric(X %*% beta_true) +
    sqrt(lambda_true) * phylo_e + sqrt(1 - lambda_true) * iid_e
  na_idx <- sample(n, size = round(0.2 * n))
  y[na_idx] <- NA
  list(y = y, X = X, R = R)
}

test_that("[lambda-cov] recovers lambda with covariates: mean bias at n = 300", {
  skip_on_cran()
  # n = 300 (per the lane brief), gated on the BIAS of mean(lambda_hat), not
  # on the per-draw MAE. At n = 300 a single lambda_hat draw is noisy (MAE
  # ~ 0.14, matching the S3 report's n = 300 measurement) -- expected, and
  # NOT gated below, only reported via message().
  #
  # 20 seeds is NOT enough replicates to evaluate the bias gate reliably:
  # with seeds 1:20, mean(lambda_hat) - 0.3 = -0.1046, which fails a < 0.05
  # gate, but the SE of that mean (~0.03, from sd(lambda_hat) ~ 0.13) is
  # itself close to the 0.05 threshold, so the pass/fail outcome for any
  # fixed 20-seed window is essentially a coin flip (checked: seeds 21:40
  # give bias -0.028, seeds 61:80 give -0.140). Extending to more replicates
  # narrows the estimate of the true bias rather than changing it -- the
  # cumulative mean over seeds 1:n_seeds stabilises around bias ~ -0.035 to
  # -0.04 by n_seeds ~ 400-600 (checked up to 800). 600 replicates at n = 300
  # runs in ~6s and gives bias ~ -0.039, a safe margin under the 0.05 gate.
  lambda_true <- 0.3
  n_seeds <- 600L
  lambda_hats <- vapply(seq_len(n_seeds), function(seed) {
    dat <- .sim_lambda_cov(seed, n = 300L, lambda_true = lambda_true)
    out <- pigauto:::bm_impute_col_with_cov(dat$y, dat$X, dat$R,
                                              lrt_threshold = 0,
                                              lambda = "estimate")
    out$lambda_hat
  }, numeric(1L))

  bias <- mean(lambda_hats) - lambda_true
  mae  <- mean(abs(lambda_hats - lambda_true))
  message(sprintf(
    "[lambda-cov] n = 300, %d seeds: mean(lambda_hat) = %.4f, bias = %.4f, MAE = %.4f (MAE not gated)",
    n_seeds, mean(lambda_hats), bias, mae))

  expect_lt(abs(bias), 0.05)
  # Seed 1's own draw, in a symmetric band around the true value.
  expect_gte(lambda_hats[1L], 0.15)
  expect_lte(lambda_hats[1L], 0.45)
})

test_that("[lambda-cov] covariate effect survives lambda = \"estimate\"", {
  dat <- .sim_lambda_cov(77L, lambda_true = 0.3)
  out <- pigauto:::bm_impute_col_with_cov(dat$y, dat$X, dat$R,
                                            lrt_threshold = 0,
                                            lambda = "estimate")
  # beta_hat[2:3] correspond to x1, x2 with true coefficients 0.8, -0.5.
  expect_true(out$beta_hat[2L] > 0)
  expect_true(out$beta_hat[3L] < 0)
  expect_lt(abs(out$beta_hat[2L] - 0.8), 0.3)
  expect_lt(abs(out$beta_hat[3L] - (-0.5)), 0.3)
})

# ---- numeric lambda interpolates --------------------------------------------

test_that("[lambda-cov] lambda = 0 matches OLS-with-covariates on missing cells", {
  set.seed(88L)
  n <- 60L
  tree <- ape::rcoal(n)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  X <- cbind(1, x1, x2)
  beta_true <- c(1, 0.8, -0.5)
  y <- as.numeric(X %*% beta_true) + stats::rnorm(n) * 0.5
  na_idx <- sample(n, 15L)
  y[na_idx] <- NA

  out0 <- pigauto:::bm_impute_col_with_cov(y, X, R, lambda = 0.0, lrt_threshold = 0)

  obs <- which(!is.na(y))
  fit_ols <- lm.fit(X[obs, , drop = FALSE], y[obs])
  pred_ols <- as.numeric(X %*% fit_ols$coefficients)

  # At lambda = 0, R(lambda) = I: the GLS-with-covariates BLUP collapses to
  # plain OLS with no phylogenetic shrinkage on the residual.
  expect_equal(out0$mu[na_idx], pred_ols[na_idx], tolerance = 1e-6)
})

test_that("[lambda-cov] lambda = 0.5 RMSE <= lambda = 1 RMSE on a lambda = 0.5 DGP", {
  set.seed(1L)
  n <- 150L
  tree <- ape::rcoal(n)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  X <- cbind(1, x1, x2)
  beta_true <- c(0, 0.8, -0.5)
  phylo_e <- as.numeric(t(chol(R)) %*% rnorm(n))
  iid_e <- stats::rnorm(n)
  lambda_true <- 0.5
  y_full <- as.numeric(X %*% beta_true) +
    sqrt(lambda_true) * phylo_e + sqrt(1 - lambda_true) * iid_e
  na_idx <- sample(n, round(0.25 * n))
  y <- y_full
  y[na_idx] <- NA

  out_half <- pigauto:::bm_impute_col_with_cov(y, X, R, lambda = 0.5, lrt_threshold = 0)
  out_one  <- pigauto:::bm_impute_col_with_cov(y, X, R, lambda = 1.0, lrt_threshold = 0)
  rmse_half <- sqrt(mean((out_half$mu[na_idx] - y_full[na_idx])^2))
  rmse_one  <- sqrt(mean((out_one$mu[na_idx] - y_full[na_idx])^2))

  expect_lte(rmse_half, rmse_one)
})

# ---- transform_tree_pagel comment correction --------------------------------

test_that("[lambda-cov] transform_tree_pagel is exact on the correlation scale for a non-ultrametric tree", {
  set.seed(99L)
  tree <- ape::rtree(30L)
  expect_false(ape::is.ultrametric(tree))

  lambda <- 0.4
  out <- pigauto:::transform_tree_pagel(tree, lambda)
  R0 <- stats::cov2cor(ape::vcv.phylo(tree))
  R_out <- stats::cov2cor(ape::vcv.phylo(out))
  tips <- tree$tip.label
  R_out <- R_out[tips, tips]
  expected <- lambda * R0 + (1 - lambda) * diag(nrow(R0))

  expect_equal(R_out, expected, tolerance = 1e-10, ignore_attr = TRUE)
})
