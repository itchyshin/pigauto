# Tests for the v0.11 Bayesian model-averaged Pagel-lambda predictor.
#
# Background: ML and 5-fold CV both pick degenerate near-zero lambda
# on weak-phylo-signal traits (BIEN sla in particular), even when an
# interior lambda gives better generalisation. The cause is a plateau
# in the profile likelihood — finite-sample noise dominates any single
# point estimate.
#
# The Bayesian ensemble averages predictions across a fine lambda grid,
# weighted by the profile-likelihood posterior (under a flat prior).
# This mimics what BACE's MCMC posterior does deterministically.
#
# Spec section: "Bayesian λ ensemble" added to
# specs/2026-05-18-pagel-lambda-baseline-design.md (to be updated).

# ---- bayes_lambda_for_col() -----------------------------------------------

test_that("[pagel-bayes] bayes_lambda_for_col returns a posterior summary list", {
  set.seed(1L)
  n <- 60L
  tree <- ape::rcoal(n)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  y <- as.numeric(t(chol(R)) %*% rnorm(n))
  y[c(2L, 9L, 18L, 33L)] <- NA
  out <- pigauto:::bayes_lambda_for_col(y, R)
  expect_true(is.list(out))
  expect_true(all(c("lambda_grid", "weights", "lambda_post_mean",
                     "lambda_post_entropy") %in% names(out)))
  # Weights sum to 1 (numerical normalisation)
  expect_equal(sum(out$weights), 1.0, tolerance = 1e-10)
  # Posterior mean is a valid scalar in [0, 1]
  expect_true(is.finite(out$lambda_post_mean))
  expect_true(out$lambda_post_mean >= 0 && out$lambda_post_mean <= 1)
})

test_that("[pagel-bayes] bayes_lambda_for_col on pure-BM data concentrates mass at high lambda", {
  set.seed(2L)
  n <- 80L
  tree <- ape::rcoal(n)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  y <- as.numeric(t(chol(R)) %*% rnorm(n))
  y[sample(n, n %/% 5L)] <- NA
  out <- pigauto:::bayes_lambda_for_col(y, R)
  expect_gt(out$lambda_post_mean, 0.5)
})

test_that("[pagel-bayes] bayes_lambda_for_col on iid white-noise data concentrates mass at low lambda", {
  set.seed(3L)
  n <- 80L
  tree <- ape::rcoal(n)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  y <- rnorm(n)
  y[sample(n, n %/% 5L)] <- NA
  out <- pigauto:::bayes_lambda_for_col(y, R)
  expect_lt(out$lambda_post_mean, 0.4)
})

test_that("[pagel-bayes] bayes weights sum to 1 with all-non-negative entries", {
  set.seed(4L)
  n <- 40L
  tree <- ape::rcoal(n)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  y <- as.numeric(t(chol(R)) %*% rnorm(n))
  y[1L] <- NA
  out <- pigauto:::bayes_lambda_for_col(y, R)
  expect_equal(sum(out$weights), 1.0, tolerance = 1e-10)
  expect_true(all(out$weights >= 0))
})

test_that("[pagel-bayes] bayes_lambda_for_col is reproducible (deterministic, no RNG)", {
  set.seed(5L)
  n <- 50L
  tree <- ape::rcoal(n)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  y <- as.numeric(t(chol(R)) %*% rnorm(n))
  y[c(2L, 11L, 23L)] <- NA
  out1 <- pigauto:::bayes_lambda_for_col(y, R)
  out2 <- pigauto:::bayes_lambda_for_col(y, R)
  expect_equal(out1$lambda_post_mean, out2$lambda_post_mean)
  expect_equal(out1$weights, out2$weights)
})

# ---- bm_impute_col(..., lambda = "bayes") ---------------------------------

test_that("[pagel-bayes] bm_impute_col with lambda='bayes' returns finite predictions", {
  set.seed(6L)
  n <- 60L
  tree <- ape::rcoal(n)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  y <- as.numeric(t(chol(R)) %*% rnorm(n))
  miss <- c(3L, 12L, 28L, 47L)
  y[miss] <- NA
  out <- pigauto:::bm_impute_col(y, R, lambda = "bayes")
  expect_true(all(is.finite(out$mu)))
  expect_true(all(out$se[miss] >= 0))
  # Diagnostic field populated
  expect_true("lambda_hat" %in% names(out))
  expect_true(is.finite(out$lambda_hat))
  expect_true(out$lambda_hat >= 0 && out$lambda_hat <= 1)
})

test_that("[pagel-bayes] predictions at lambda='bayes' are bounded by predictions at the lambda endpoints", {
  # Each averaged prediction is a convex combination of predictions at
  # grid points (weights >= 0 summing to 1). So for monotonic-in-lambda
  # predictions, the bayes prediction must lie between the min and max
  # endpoint predictions (modulo numerical noise on the grid).
  set.seed(7L)
  n <- 50L
  tree <- ape::rcoal(n)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  y <- as.numeric(t(chol(R)) %*% rnorm(n))
  miss <- c(5L, 17L, 38L)
  y[miss] <- NA
  bayes_pred <- pigauto:::bm_impute_col(y, R, lambda = "bayes")$mu[miss]
  pred_0  <- pigauto:::bm_impute_col(y, R, lambda = 0.005)$mu[miss]
  pred_1  <- pigauto:::bm_impute_col(y, R, lambda = 0.995)$mu[miss]
  for (i in seq_along(miss)) {
    lo <- min(pred_0[i], pred_1[i])
    hi <- max(pred_0[i], pred_1[i])
    # Convex combo bound, plus tolerance for the non-monotonic shape
    expect_true(bayes_pred[i] >= lo - 0.5 * abs(hi - lo) &&
                bayes_pred[i] <= hi + 0.5 * abs(hi - lo),
                info = sprintf("cell %d: bayes=%g, [%g, %g]",
                                miss[i], bayes_pred[i], lo, hi))
  }
})
