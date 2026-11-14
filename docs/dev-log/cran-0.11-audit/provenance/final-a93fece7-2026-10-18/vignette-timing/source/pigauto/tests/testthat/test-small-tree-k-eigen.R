# Issue #208: default k_eigen = "auto" completes on 3-, 4-, and 5-tip trees.

small_tree_impute <- function(n) {
  set.seed(208L + n)
  tree <- ape::rtree(n)
  df <- data.frame(
    y = as.numeric(ape::rTraitCont(tree)),
    row.names = tree$tip.label
  )
  df$y[1L] <- NA_real_
  impute(df, tree, epochs = 5L, verbose = FALSE, seed = 208L)
}

test_that("impute() with default k_eigen completes on 3-, 4-, and 5-tip trees", {
  for (n in c(3L, 4L, 5L)) {
    res <- small_tree_impute(n)
    expect_s3_class(res, "pigauto_result")
    expect_equal(nrow(res$completed), n)
  }
})

test_that("auto k_eigen clamps so k + 1 <= n on 3- and 4-tip trees", {
  g3 <- build_phylo_graph(ape::rtree(3))
  g4 <- build_phylo_graph(ape::rtree(4))
  expect_lte(ncol(g3$coords) + 1L, 3L)
  expect_lte(ncol(g4$coords) + 1L, 4L)
})
