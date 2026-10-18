# Tests for Hadfield-Nakagawa sparse S^{-1} construction.

test_that("build_henderson_S_inv reproduces dense A^{-1} exactly on small trees", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("ape")
  set.seed(2026L)

  for (n in c(10L, 50L)) {
    tree <- ape::rcoal(n)
    H <- build_henderson_S_inv(tree)

    # Henderson's Q is on the raw VCV (depth-weighted) scale.
    A_raw  <- ape::vcv(tree)
    A_raw  <- A_raw[tree$tip.label, tree$tip.label]
    A_inv_dense <- solve(A_raw)

    b <- diag(n)
    A_inv_sparse <- henderson_R_inv_apply(b, H)

    expect_equal(dim(A_inv_sparse), dim(A_inv_dense))
    expect_lt(max(abs(A_inv_sparse - A_inv_dense)) / max(abs(A_inv_dense)),
              1e-6)  # relative tol; tiny at n=10/50
  }
})

test_that("henderson_log_det_R matches log|A| on small trees", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("ape")
  set.seed(2026L)
  tree <- ape::rcoal(20L)
  H <- build_henderson_S_inv(tree)
  A_raw <- ape::vcv(tree)[tree$tip.label, tree$tip.label]
  ld_dense  <- as.numeric(determinant(A_raw, logarithm = TRUE)$modulus)
  ld_sparse <- henderson_log_det_R(H)
  expect_equal(ld_sparse, ld_dense, tolerance = 1e-8)
})

test_that("build_henderson_S_inv returns sparse Q with O(n) nonzeros", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("ape")
  set.seed(2026L)
  for (n in c(50L, 500L)) {
    tree <- ape::rcoal(n)
    H <- build_henderson_S_inv(tree)
    # Each non-root edge contributes at most 4 nonzeros to Q (2 diag
    # accumulations + 2 off-diagonal entries). The diagonal also picks
    # up children contributions but those merge with existing entries.
    nnz <- Matrix::nnzero(H$Q)
    # Loose upper bound: 6 * (n - 1) for binary trees per the paper.
    expect_lt(nnz, 6L * H$N)
  }
})

test_that("build_henderson_S_inv tip ordering matches tree$tip.label", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("ape")
  set.seed(2026L)
  tree <- ape::rcoal(15L)
  H <- build_henderson_S_inv(tree)
  expect_equal(H$tip_labels, tree$tip.label)
  expect_equal(length(H$tip_idx), H$n_tips)
  # tip_sqrt_d is stored to enable correlation-scale Schur (cor_scale=TRUE).
  expect_equal(length(H$tip_sqrt_d), H$n_tips)
  expect_true(all(H$tip_sqrt_d > 0))
  A_raw <- ape::vcv(tree)[tree$tip.label, tree$tip.label]
  expect_equal(H$tip_sqrt_d, sqrt(diag(A_raw)), tolerance = 1e-10)
})

test_that("henderson_R_inv_apply(cor_scale = TRUE) matches solve(R) on a non-ultrametric tree", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("ape")
  set.seed(2026L)
  # ape::rtree() draws non-ultrametric trees (variable tip depths), which
  # is the regime where A-scale vs R-scale matters most. Verifies the
  # cor_scale rescale matches dense solve(R) to ~1e-6 relative.
  for (n in c(20L, 60L)) {
    tree <- ape::rtree(n)
    H <- build_henderson_S_inv(tree)
    R   <- cov2cor(ape::vcv(tree))[tree$tip.label, tree$tip.label]
    R_inv_dense <- solve(R)
    # b = identity gives R_inv as the result column-by-column.
    R_inv_hen <- henderson_R_inv_apply(diag(n), H, cor_scale = TRUE)
    expect_equal(dim(R_inv_hen), dim(R_inv_dense))
    rel_err <- max(abs(R_inv_hen - R_inv_dense)) / max(abs(R_inv_dense))
    expect_lt(rel_err, 1e-6)
  }
})

test_that("henderson_R_inv_apply(cor_scale = FALSE) still returns A^{-1}", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("ape")
  set.seed(2026L)
  tree <- ape::rtree(30L)
  H <- build_henderson_S_inv(tree)
  A <- ape::vcv(tree)[tree$tip.label, tree$tip.label]
  A_inv_dense <- solve(A)
  A_inv_hen <- henderson_R_inv_apply(diag(nrow(A)), H, cor_scale = FALSE)
  rel_err <- max(abs(A_inv_hen - A_inv_dense)) / max(abs(A_inv_dense))
  expect_lt(rel_err, 1e-6)
})

test_that("henderson_bm_predict(cor_scale = TRUE) approximates bm_impute_col(., R)", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("ape")
  set.seed(2026L)
  tree <- ape::rcoal(40L)
  R <- cov2cor(ape::vcv(tree))[tree$tip.label, tree$tip.label]
  y <- as.numeric(t(chol(R)) %*% rnorm(40L))
  names(y) <- tree$tip.label
  miss <- sample(40L, 12L)
  y_obs <- y; y_obs[miss] <- NA

  # Dense reference: bm_impute_col operates on R (correlation scale).
  dense <- bm_impute_col(y_obs, R, nugget = 1e-8)

  H <- build_henderson_S_inv(tree)
  hen <- henderson_bm_predict(y_obs, H, cor_scale = TRUE)

  # Means at the missing tips should agree to ~1% relative. The two paths
  # share the same conditional-MVN math but differ in numerical realisation:
  # dense bm_impute_col forms (R_oo + nugget I)^{-1}, sparse henderson uses
  # the Cholesky factorisation of Q_no_no with a tip-edge-only ridge. Both
  # converge to the same posterior mean as n -> infty, but at small n the
  # discrepancy is dominated by the nugget/ridge gap.
  rel_err <- max(abs(hen$mu[miss] - dense$mu[miss])) /
             max(abs(dense$mu[miss]))
  expect_lt(rel_err, 1e-2)
})
