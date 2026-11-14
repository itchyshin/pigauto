test_that("phylogenetic masking selects the sole eligible internal node", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  x <- matrix(1:3, ncol = 1, dimnames = list(c("a", "b", "c"), "x"))
  for (seed in c(1L, 7L, 19L)) {
    set.seed(seed)
    p <- generate_missing_probs(x, 3L, 1L, "MAR_phylo",
      list(n_clades = 1L, p_clade = .7, p_base = .1), tree, NULL)
    expect_equal(as.numeric(p), c(.7, .7, .1), info = paste("seed", seed))
  }
  star <- ape::read.tree(text = "(a:1,b:1,c:1);")
  set.seed(7)
  p <- generate_missing_probs(x, 3L, 1L, "MAR_phylo",
    list(n_clades = 2L, p_clade = .7, p_base = .1), star, NULL)
  expect_equal(as.numeric(p), rep(.1, 3))
})
