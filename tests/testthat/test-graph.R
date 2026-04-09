test_that("build_phylo_graph returns list with adj and coords", {
  set.seed(1)
  tree <- ape::rtree(20)
  g <- build_phylo_graph(tree, k_eigen = 4L)
  expect_type(g, "list")
  expect_named(g, c("adj", "coords", "n", "sigma", "D"))
})

test_that("build_phylo_graph adj is square and n matches ntips", {
  set.seed(2)
  tree <- ape::rtree(15)
  g <- build_phylo_graph(tree, k_eigen = 4L)
  expect_equal(nrow(g$adj), 15L)
  expect_equal(ncol(g$adj), 15L)
  expect_equal(g$n, 15L)
})

test_that("build_phylo_graph coords has k_eigen columns", {
  set.seed(3)
  tree <- ape::rtree(25)
  g <- build_phylo_graph(tree, k_eigen = 6L)
  expect_equal(ncol(g$coords), 6L)
  expect_equal(nrow(g$coords), 25L)
})

test_that("build_phylo_graph adj is symmetric and entries are in (0, 1)", {
  set.seed(4)
  tree <- ape::rtree(20)
  g    <- build_phylo_graph(tree)
  # Symmetric normalisation: entries should be positive and < 1
  expect_true(all(g$adj > 0))
  expect_lte(max(g$adj), 1.0 + 1e-8)
  # Symmetric: A[i,j] == A[j,i]
  expect_equal(g$adj, t(g$adj), tolerance = 1e-8)
})

test_that("build_phylo_graph caching returns identical results", {
  set.seed(5)
  tree  <- ape::rtree(20)
  cache <- tempfile(fileext = ".rds")
  on.exit(unlink(cache))
  g1 <- build_phylo_graph(tree, k_eigen = 4L, cache_path = cache)
  g2 <- build_phylo_graph(tree, k_eigen = 4L, cache_path = cache)
  expect_equal(g1$adj,    g2$adj)
  expect_equal(g1$coords, g2$coords)
})

test_that("build_phylo_graph returns the cophenetic distance matrix in $D", {
  set.seed(6)
  tree <- ape::rtree(30)
  g    <- build_phylo_graph(tree, k_eigen = 4L)

  expect_true("D" %in% names(g))
  expect_true(is.matrix(g$D))
  expect_equal(dim(g$D), c(30L, 30L))

  # Should match a fresh cophenetic() call bit-for-bit.
  D_ref <- ape::cophenetic.phylo(tree)
  expect_equal(unname(g$D), unname(D_ref))
})

test_that("sparse (RSpectra) and dense eigensolvers agree up to sign flips", {
  skip_if_not_installed("RSpectra")

  set.seed(7)
  # Use ape::rtree (asymmetric branch lengths): the resulting graph
  # Laplacian has well-separated eigenvalues, so individual eigenvectors
  # are uniquely defined up to sign. Ultrametric trees (ape::rcoal) give
  # a highly symmetric Laplacian with large degenerate eigenvalue
  # clusters, in which case the sparse and dense solvers can pick
  # different bases for the degenerate subspace -- that is numerically
  # correct behaviour but makes column-by-column comparison ill-posed,
  # so we avoid it in this test.
  #
  # > 500 tips ensures the sparse path is triggered by default.
  tree <- ape::rtree(600)
  D    <- ape::cophenetic.phylo(tree)

  k <- 8L

  # Dense reference (force dense by raising the threshold above n).
  dense_vecs  <- pigauto:::spectral_features_from_D(
    D, k = k, sigma_mult = 0.5, dense_threshold = 10000L
  )
  # Sparse path (Lanczos via RSpectra).
  sparse_vecs <- pigauto:::spectral_features_from_D(
    D, k = k, sigma_mult = 0.5, dense_threshold = 0L
  )

  expect_equal(dim(dense_vecs), dim(sparse_vecs))

  # Column-by-column agreement up to sign (well-defined only when the
  # selected eigenvalues are non-degenerate; rtree gives that).
  for (j in seq_len(k)) {
    s <- sign(sum(dense_vecs[, j] * sparse_vecs[, j]))
    if (s == 0) s <- 1
    expect_equal(s * sparse_vecs[, j], dense_vecs[, j], tolerance = 1e-4,
                 info = paste("eigenvector column", j))
  }

  # Subspace-invariant check: the projector V V^T is basis-independent,
  # so dense and sparse must agree even if the basis rotates inside a
  # degenerate eigenvalue cluster. This is the test that actually
  # matters for downstream GNN inputs, since the model only sees the
  # k-dim subspace the vectors span.
  P_dense  <- dense_vecs  %*% t(dense_vecs)
  P_sparse <- sparse_vecs %*% t(sparse_vecs)
  expect_equal(P_sparse, P_dense, tolerance = 1e-4)
})

test_that("fit_baseline reuses graph$D when supplied instead of recomputing", {
  # Regression test for Fix B: passing `graph` to fit_baseline should give
  # numerically identical mu and se as the old code path (no graph).
  set.seed(8)
  n    <- 50L
  tree <- ape::rcoal(n)
  sp   <- tree$tip.label
  traits <- data.frame(
    row.names = sp,
    x = stats::rnorm(n),
    y = stats::rnorm(n)
  )
  traits[sample.int(n, 5), 1] <- NA
  traits[sample.int(n, 5), 2] <- NA

  pd    <- preprocess_traits(traits, tree, log_transform = FALSE)
  graph <- build_phylo_graph(tree, k_eigen = 4L)

  bl_no_graph   <- fit_baseline(pd, tree)
  bl_with_graph <- fit_baseline(pd, tree, graph = graph)

  expect_equal(bl_no_graph$mu, bl_with_graph$mu, tolerance = 1e-10)
  expect_equal(bl_no_graph$se, bl_with_graph$se, tolerance = 1e-10)
})
