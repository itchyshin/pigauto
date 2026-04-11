# tests/testthat/test-bm-internal.R
# Unit tests for internal BM baseline (R/bm_internal.R)

test_that("phylo_cor_matrix returns correct structure", {
  set.seed(1)
  tree <- ape::rtree(20)
  R <- phylo_cor_matrix(tree)

  expect_true(is.matrix(R))
  expect_equal(nrow(R), 20L)
  expect_equal(ncol(R), 20L)
  # Diagonal = 1 (correlation matrix)
  expect_equal(unname(diag(R)), rep(1, 20), tolerance = 1e-12)
  # Symmetric
  expect_equal(R, t(R))
  # Positive definite (all eigenvalues > 0)
  evals <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(evals > 0))
  # Dimnames match tip labels
  expect_equal(rownames(R), tree$tip.label)
  expect_equal(colnames(R), tree$tip.label)
})

test_that("bm_impute_col recovers observed values exactly", {
  set.seed(2)
  tree <- ape::rtree(30)
  R <- phylo_cor_matrix(tree)

  # Simulate some data
  y <- rnorm(30)
  # Mask 10 species
  miss <- sample(30, 10)
  y_masked <- y
  y_masked[miss] <- NA

  res <- bm_impute_col(y_masked, R)

  # Observed species keep exact values
  obs <- setdiff(seq_len(30), miss)
  expect_equal(res$mu[obs], y_masked[obs])
  # Observed species have se = 0
 expect_equal(res$se[obs], rep(0, length(obs)))
  # Missing species have finite mu and positive se
  expect_true(all(is.finite(res$mu[miss])))
  expect_true(all(res$se[miss] > 0))
})

test_that("bm_impute_col returns correct structure", {
  set.seed(3)
  tree <- ape::rtree(20)
  R <- phylo_cor_matrix(tree)

  y <- rnorm(20)
  y[c(1, 5, 10)] <- NA

  res <- bm_impute_col(y, R)

  expect_true(is.list(res))
  expect_named(res, c("mu", "se"))
  expect_length(res$mu, 20L)
  expect_length(res$se, 20L)
  expect_true(all(is.finite(res$mu)))
  expect_true(all(is.finite(res$se)))
  expect_true(all(res$se >= 0))
})

test_that("bm_impute_col handles all-observed case", {
  set.seed(4)
  tree <- ape::rtree(15)
  R <- phylo_cor_matrix(tree)

  y <- rnorm(15)  # no NAs
  res <- bm_impute_col(y, R)

  expect_equal(res$mu, y)
  expect_equal(res$se, rep(0, 15))
})

test_that("bm_impute_col falls back for < 5 observations", {
  set.seed(5)
  tree <- ape::rtree(20)
  R <- phylo_cor_matrix(tree)

  y <- rep(NA_real_, 20)
  y[c(1, 2, 3)] <- c(1.0, 2.0, 3.0)  # only 3 observed

  res <- bm_impute_col(y, R)

  # Observed keep values
  expect_equal(res$mu[1:3], c(1.0, 2.0, 3.0))
  # Missing get global mean
  expect_equal(res$mu[4:20], rep(mean(c(1, 2, 3)), 17))
  # SE for missing = global sd
  expect_equal(res$se[4:20], rep(sd(c(1, 2, 3)), 17))
  # SE for observed = 0
  expect_equal(res$se[1:3], rep(0, 3))
})

test_that("bm_impute_col handles near-singular R_oo via nugget back-off", {
  # Star tree: all tips at equal distance from root
  n <- 15
  tree <- ape::stree(n, type = "star")
  tree$edge.length <- rep(1, nrow(tree$edge))
  R <- phylo_cor_matrix(tree)

  y <- rnorm(n)
  y[c(1, 2, 3)] <- NA

  # Should not error, thanks to nugget regularisation
  res <- expect_no_error(bm_impute_col(y, R))
  expect_true(all(is.finite(res$mu)))
  expect_true(all(is.finite(res$se)))
})

test_that("bm_impute_col handles zero observations", {
  set.seed(7)
  tree <- ape::rtree(10)
  R <- phylo_cor_matrix(tree)

  y <- rep(NA_real_, 10)  # all missing
  res <- bm_impute_col(y, R)

  # Fallback: mu = 0, se = 1
  expect_equal(res$mu, rep(0, 10))
  expect_equal(res$se, rep(1, 10))
})

test_that("fit_baseline with internal BM produces finite mu and se", {
  skip_if_not_installed("torch")
  skip_if(!torch::torch_is_installed(), "torch backend not installed")

  set.seed(8)
  tree <- ape::rtree(30)
  # Simple continuous data
  n <- 30
  df <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    row.names = tree$tip.label
  )
  df$x1[sample(n, 5)] <- NA
  df$x2[sample(n, 5)] <- NA

  pd <- preprocess_traits(df, tree)
  graph <- build_phylo_graph(tree, k_eigen = 4L)
  bl <- fit_baseline(pd, tree, graph = graph)

  expect_true(all(is.finite(bl$mu)))
  expect_true(all(is.finite(bl$se)))
  expect_equal(nrow(bl$mu), n)
  expect_equal(ncol(bl$mu), ncol(pd$X_scaled))
})

test_that("fit_baseline with internal BM and splits produces finite results", {
  skip_if_not_installed("torch")
  skip_if(!torch::torch_is_installed(), "torch backend not installed")

  set.seed(9)
  tree <- ape::rtree(30)
  n <- 30
  df <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    row.names = tree$tip.label
  )

  pd <- preprocess_traits(df, tree)
  spl <- make_missing_splits(pd$X_scaled, missing_frac = 0.25, seed = 1,
                             trait_map = pd$trait_map)
  graph <- build_phylo_graph(tree, k_eigen = 4L)
  bl <- fit_baseline(pd, tree, splits = spl, graph = graph)

  expect_true(all(is.finite(bl$mu)))
  expect_true(all(is.finite(bl$se)))
  # SE should be > 0 for at least some species (those missing after masking)
  expect_true(any(bl$se > 0))
})

test_that("build_phylo_graph caches R_phy alongside D", {
  set.seed(10)
  tree <- ape::rtree(20)
  g <- build_phylo_graph(tree, k_eigen = 4L)

  expect_true(!is.null(g$R_phy))
  expect_equal(nrow(g$R_phy), 20L)
  expect_equal(ncol(g$R_phy), 20L)
  expect_equal(unname(diag(g$R_phy)), rep(1, 20), tolerance = 1e-12)
  expect_equal(rownames(g$R_phy), tree$tip.label)
})
