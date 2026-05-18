# Tests for the v0.11 eigendecomposition-cached Pagel-λ NLL.
# Spec: specs/2026-05-18-pagel-lambda-eigendecomp-speedup-design.md.

# Dense-Cholesky reference implementation (matches v0.10 ml_lambda_for_col's
# inner nll), used to verify the cached version produces identical results.
.nll_dense_reference <- function(y, R, lambda, nugget = 1e-6) {
  obs <- which(!is.na(y))
  n_o <- length(obs)
  R_oo <- R[obs, obs, drop = FALSE]
  y_o <- y[obs]
  ones <- rep(1, n_o)
  R_l <- lambda * R_oo
  diag(R_l) <- lambda * diag(R_oo) + (1 - lambda) + nugget
  L <- chol(R_l)
  chol_solve <- function(b) backsolve(L, forwardsolve(t(L), b))
  a <- chol_solve(ones)
  b <- chol_solve(y_o)
  mu_hat <- sum(b) / sum(a)
  e <- y_o - mu_hat
  e_solve <- chol_solve(e)
  sigma2 <- as.numeric(crossprod(e, e_solve)) / max(n_o - 1L, 1L)
  log_det <- 2 * sum(log(diag(L)))
  0.5 * ((n_o - 1L) * log(sigma2) + log_det)
}

# ---- correctness ----------------------------------------------------------

test_that("[pagel-perf] build_pagel_nll_cache returns a closure that matches dense Cholesky NLL", {
  set.seed(1L)
  tree <- ape::rcoal(40L)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  y <- as.numeric(t(chol(R)) %*% rnorm(40L))
  y[c(2L, 7L, 15L, 28L)] <- NA

  cache <- pigauto:::build_pagel_nll_cache(y, R)
  expect_true(is.function(cache$nll))
  for (lam in c(0.001, 0.05, 0.2, 0.5, 0.8, 0.99)) {
    ref <- .nll_dense_reference(y, R, lambda = lam)
    out <- cache$nll(lam)
    expect_equal(out, ref, tolerance = 1e-8,
                 info = sprintf("lambda = %.3f: cached %.8f vs dense %.8f",
                                 lam, out, ref))
  }
})

test_that("[pagel-perf] cached NLL matches dense across 30 random datasets", {
  for (seed in seq_len(30L)) {
    set.seed(seed)
    n <- sample(20L:80L, 1L)
    tree <- ape::rcoal(n)
    R <- stats::cov2cor(ape::vcv.phylo(tree))
    R <- R[tree$tip.label, tree$tip.label]
    y <- as.numeric(t(chol(R)) %*% rnorm(n))
    n_miss <- max(2L, n %/% 5L)
    y[sample(n, n_miss)] <- NA
    cache <- pigauto:::build_pagel_nll_cache(y, R)
    for (lam in c(0.01, 0.3, 0.7, 0.99)) {
      ref <- .nll_dense_reference(y, R, lambda = lam)
      out <- cache$nll(lam)
      expect_equal(out, ref, tolerance = 1e-7,
                   info = sprintf("seed=%d lam=%.3f", seed, lam))
    }
  }
})

test_that("[pagel-perf] cache stores n_o and works on all-observed input", {
  set.seed(2L)
  tree <- ape::rcoal(15L)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  y <- as.numeric(t(chol(R)) %*% rnorm(15L))
  cache <- pigauto:::build_pagel_nll_cache(y, R)
  expect_equal(cache$n_o, 15L)
})

# ---- refactored ml_lambda_for_col still produces the same lambda_hat -------

test_that("[pagel-perf] eigendecomp-based ml_lambda_for_col reproduces v0.10 lambda_hat", {
  # Strong-BM trait: lambda_hat should land in upper half regardless of impl.
  set.seed(3L)
  tree <- ape::rcoal(60L)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  L <- t(chol(R))
  y <- as.numeric(L %*% rnorm(60L))
  y[seq(1L, 60L, by = 5L)] <- NA
  lam_hat <- pigauto:::ml_lambda_for_col(y, R)
  expect_gt(lam_hat, 0.5)
  expect_lte(lam_hat, 1.0)

  # Iid white-noise: lambda_hat should land near 0.
  set.seed(4L)
  y_iid <- rnorm(60L)
  y_iid[seq(1L, 60L, by = 5L)] <- NA
  lam_iid <- pigauto:::ml_lambda_for_col(y_iid, R)
  expect_lt(lam_iid, 0.3)
})

# ---- speed regression ------------------------------------------------------

test_that("[pagel-perf] 11-point grid evaluation at n=600 completes in under 3 sec", {
  skip_on_cran()
  set.seed(5L)
  tree <- ape::rcoal(600L)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  y <- as.numeric(t(chol(R)) %*% rnorm(600L))
  y[seq(1L, 600L, by = 4L)] <- NA  # 25% missing
  t0 <- Sys.time()
  cache <- pigauto:::build_pagel_nll_cache(y, R)
  for (lam in seq(0.005, 0.995, length.out = 11L)) cache$nll(lam)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("\n  [pagel-perf] eigendecomp + 11 NLL evals at n=600: %.2f s\n",
              elapsed))
  expect_lt(elapsed, 3.0)
})
