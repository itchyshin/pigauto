test_that("build_phylo_graph returns list with adj and coords", {
  set.seed(1)
  tree <- ape::rtree(20)
  g <- build_phylo_graph(tree, k_eigen = 4L)
  expect_type(g, "list")
  expect_named(g, c("adj", "coords", "n", "sigma"))
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
