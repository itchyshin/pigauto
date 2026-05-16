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
})
