# Tests for the v0.10 Pagel's lambda baseline (spec: 2026-05-18-pagel-lambda-baseline-design.md).

# ---- transform_tree_pagel() ----------------------------------------------

test_that("[pagel] transform_tree_pagel returns a phylo with same tip set", {
  set.seed(1L)
  tree <- ape::rcoal(20L)
  out <- transform_tree_pagel(tree, lambda = 0.5)
  expect_s3_class(out, "phylo")
  expect_setequal(out$tip.label, tree$tip.label)
  expect_equal(ape::Ntip(out), ape::Ntip(tree))
  expect_equal(ape::Nnode(out), ape::Nnode(tree))
})

test_that("[pagel] transform_tree_pagel(lambda = 1) is identity on internal edges", {
  set.seed(2L)
  tree <- ape::rcoal(15L)
  out <- transform_tree_pagel(tree, lambda = 1.0)
  # All edge lengths preserved
  expect_equal(out$edge.length, tree$edge.length, tolerance = 1e-12)
})

test_that("[pagel] transform_tree_pagel(lambda = 0) produces a star tree (zero off-diagonal correlation)", {
  set.seed(3L)
  tree <- ape::rcoal(10L)
  out <- transform_tree_pagel(tree, lambda = 0.0)
  n_tip <- ape::Ntip(tree)
  is_terminal <- tree$edge[, 2L] <= n_tip
  # Internal edges zeroed
  expect_true(all(out$edge.length[!is_terminal] == 0))
  # Correlation matrix is identity (modulo numerical noise)
  R <- stats::cov2cor(ape::vcv.phylo(out))
  expect_equal(R, diag(n_tip), tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("[pagel] transform_tree_pagel(lambda = 0.5) halves off-diagonal correlations", {
  set.seed(4L)
  tree <- ape::rcoal(12L)
  out <- transform_tree_pagel(tree, lambda = 0.5)
  R0 <- stats::cov2cor(ape::vcv.phylo(tree))
  R_lam <- stats::cov2cor(ape::vcv.phylo(out))
  # Diagonal stays at 1
  expect_equal(diag(R_lam), rep(1, ncol(R_lam)), tolerance = 1e-10,
               ignore_attr = TRUE)
  # Off-diagonal halved (modulo tip ordering)
  tips <- tree$tip.label
  R_lam <- R_lam[tips, tips]
  off_R0 <- R0[upper.tri(R0)]
  off_lam <- R_lam[upper.tri(R_lam)]
  expect_equal(off_lam, 0.5 * off_R0, tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("[pagel] transform_tree_pagel rejects out-of-range lambda", {
  tree <- ape::rcoal(8L)
  expect_error(transform_tree_pagel(tree, lambda = -0.1))
  expect_error(transform_tree_pagel(tree, lambda = 1.5))
  expect_error(transform_tree_pagel(tree, lambda = NA_real_))
})

test_that("[pagel] transform_tree_pagel induces R(lambda) = lambda*R + (1-lambda)*I on vcv", {
  set.seed(5L)
  tree <- ape::rcoal(10L)
  R0 <- stats::cov2cor(ape::vcv.phylo(tree))
  lambda <- 0.4
  R_lambda_via_formula <- lambda * R0 + (1 - lambda) * diag(nrow(R0))
  tree_t <- transform_tree_pagel(tree, lambda = lambda)
  R_lambda_via_tree <- stats::cov2cor(ape::vcv.phylo(tree_t))
  # Tip ordering may differ; align on tip.label.
  tips <- tree$tip.label
  R_lambda_via_tree <- R_lambda_via_tree[tips, tips]
  expect_equal(R_lambda_via_tree, R_lambda_via_formula, tolerance = 1e-10)
})

# ---- bm_impute_col(..., lambda = ) ---------------------------------------

test_that("[pagel] bm_impute_col with lambda=1 is identical to no-lambda call (back-compat)", {
  set.seed(10L)
  tree <- ape::rcoal(20L)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  y <- as.numeric(t(chol(R)) %*% rnorm(20L))
  y[c(3L, 7L, 15L)] <- NA
  ref <- pigauto:::bm_impute_col(y, R)
  out <- pigauto:::bm_impute_col(y, R, lambda = 1.0)
  expect_equal(out$mu, ref$mu, tolerance = 1e-12)
  expect_equal(out$se, ref$se, tolerance = 1e-12)
})

test_that("[pagel] bm_impute_col with lambda=0 predicts the GLS grand mean on missing cells", {
  set.seed(11L)
  tree <- ape::rcoal(25L)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  y <- as.numeric(t(chol(R)) %*% rnorm(25L)) + 10  # offset
  miss_idx <- c(2L, 9L, 18L, 22L)
  y[miss_idx] <- NA
  out <- pigauto:::bm_impute_col(y, R, lambda = 0.0)
  # Under lambda = 0, R = I, GLS mean = arithmetic mean of observed cells.
  mu_obs <- mean(y, na.rm = TRUE)
  expect_equal(out$mu[miss_idx], rep(mu_obs, length(miss_idx)),
               tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("[pagel] bm_impute_col rejects out-of-range lambda", {
  tree <- ape::rcoal(10L)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  y <- rnorm(10L); y[1L] <- NA
  expect_error(pigauto:::bm_impute_col(y, R, lambda = -0.1))
  expect_error(pigauto:::bm_impute_col(y, R, lambda = 1.5))
})

test_that("[pagel] bm_impute_col is finite + smooth across lambda in [0, 1]", {
  # The conditional mean is non-linear in lambda (mu_hat itself depends
  # on lambda via the GLS estimator), so we don't require monotone
  # interpolation between endpoints. But predictions should be finite
  # and vary continuously across the lambda range.
  set.seed(12L)
  tree <- ape::rcoal(20L)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  y <- as.numeric(t(chol(R)) %*% rnorm(20L))
  miss <- c(4L, 11L)
  y[miss] <- NA
  lams <- seq(0, 1, length.out = 11)
  mus <- sapply(lams, function(l) pigauto:::bm_impute_col(y, R, lambda = l)$mu[miss])
  expect_true(all(is.finite(mus)))
  # Continuity: adjacent lambda predictions should not jump wildly.
  diffs <- apply(mus, 1L, function(row) max(abs(diff(row))))
  obs_sd <- stats::sd(y, na.rm = TRUE)
  expect_true(all(diffs < obs_sd),
              info = sprintf("max adj-diff = %s, sd(y_obs) = %.3f",
                              paste(round(diffs, 3), collapse = ", "), obs_sd))
})

# ---- bm_impute_col(..., lambda = "estimate") -----------------------------

test_that('[pagel] bm_impute_col(lambda = "estimate") returns a value in [0, 1]', {
  set.seed(20L)
  tree <- ape::rcoal(40L)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  y <- as.numeric(t(chol(R)) %*% rnorm(40L))
  y[c(2L, 8L, 25L)] <- NA
  out <- pigauto:::bm_impute_col(y, R, lambda = "estimate")
  expect_true("lambda_hat" %in% names(out))
  expect_true(is.finite(out$lambda_hat))
  expect_true(out$lambda_hat >= 0 && out$lambda_hat <= 1)
})

test_that('[pagel] bm_impute_col(lambda = "estimate") recovers high lambda on strong-BM data', {
  # Simulate a trait under pure BM (lambda = 1). ML should recover lambda_hat near 1.
  set.seed(21L)
  tree <- ape::rcoal(60L)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  L <- t(chol(R))
  y <- as.numeric(L %*% rnorm(60L))
  y[seq(1L, 60L, by = 5L)] <- NA  # 12 missing
  out <- pigauto:::bm_impute_col(y, R, lambda = "estimate")
  expect_gt(out$lambda_hat, 0.7)
})

test_that('[pagel] bm_impute_col(lambda = "estimate") recovers low lambda on iid (white-noise) data', {
  # Simulate iid white-noise trait (no phylo structure). ML should
  # recover lambda_hat near 0.
  set.seed(22L)
  tree <- ape::rcoal(60L)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  y <- rnorm(60L)
  y[seq(1L, 60L, by = 5L)] <- NA
  out <- pigauto:::bm_impute_col(y, R, lambda = "estimate")
  expect_lt(out$lambda_hat, 0.3)
})

# ---- fit_baseline + dispatcher threading ----------------------------------

test_that("[pagel] fit_baseline(lambda_mode='fixed_1') is back-compatible (no behaviour change)", {
  set.seed(30L)
  tree <- ape::rcoal(30L)
  df <- data.frame(
    x1 = stats::rnorm(30L),
    x2 = stats::rnorm(30L),
    row.names = tree$tip.label
  )
  df$x1[c(2L, 5L, 17L)] <- NA
  df$x2[c(3L, 11L)] <- NA
  pd <- preprocess_traits(df, tree)
  bl_default <- fit_baseline(pd, tree)
  bl_fixed   <- fit_baseline(pd, tree, lambda_mode = "fixed_1")
  expect_equal(bl_fixed$mu, bl_default$mu, tolerance = 1e-12)
  expect_equal(bl_fixed$se, bl_default$se, tolerance = 1e-12)
})

test_that("[pagel] fit_baseline(lambda_mode='estimate') runs end-to-end + produces finite predictions", {
  set.seed(31L)
  tree <- ape::rcoal(40L)
  df <- data.frame(
    x1 = stats::rnorm(40L),    # iid: should pull lambda toward 0
    x2 = stats::rnorm(40L),
    row.names = tree$tip.label
  )
  df$x1[c(2L, 5L, 17L, 28L)] <- NA
  df$x2[c(3L, 11L, 35L)] <- NA
  pd <- preprocess_traits(df, tree)
  bl_est <- fit_baseline(pd, tree, lambda_mode = "estimate")
  expect_true(all(is.finite(bl_est$mu)))
  expect_true(all(is.finite(bl_est$se)))
  expect_true(all(bl_est$se >= 0))
})

test_that("[pagel] fit_baseline lambda_mode='estimate' produces different predictions vs 'fixed_1' on iid data", {
  set.seed(32L)
  tree <- ape::rcoal(50L)
  df <- data.frame(
    x1 = stats::rnorm(50L),       # iid (no phylo signal)
    row.names = tree$tip.label
  )
  miss <- c(2L, 7L, 22L, 41L)
  df$x1[miss] <- NA
  pd <- preprocess_traits(df, tree)
  bl_fixed <- fit_baseline(pd, tree, lambda_mode = "fixed_1")
  bl_est   <- fit_baseline(pd, tree, lambda_mode = "estimate")
  # For iid data, the lambda-estimate path should shrink toward grand mean
  # (different from lambda=1 BM kriging predictions). Predictions on missing
  # cells should differ between modes.
  diff <- bl_fixed$mu[miss, 1L] - bl_est$mu[miss, 1L]
  expect_true(max(abs(diff)) > 1e-4,
              info = sprintf("max abs diff = %.6f (should be > 1e-4)",
                              max(abs(diff))))
})

test_that("[pagel] impute(lambda_mode=...) threads to baseline and model config", {
  skip_if_no_libtorch()
  set.seed(33L)
  tree <- ape::rcoal(24L)
  df <- data.frame(
    x = stats::rnorm(24L),
    row.names = tree$tip.label
  )
  df$x[c(3L, 9L, 17L)] <- NA

  res <- suppressWarnings(impute(
    df, tree,
    lambda_mode = "estimate",
    epochs = 1L,
    missing_frac = 0.2,
    verbose = FALSE,
    seed = 33L
  ))

  expect_equal(res$fit$model_config$lambda_mode, "estimate")
  expect_true(all(is.finite(unlist(res$prediction$imputed, use.names = FALSE))))
})
