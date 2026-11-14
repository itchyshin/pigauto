test_that("read_tree returns every tree from a multi-tree Newick file", {
  path <- tempfile(fileext = ".tre")
  on.exit(unlink(path), add = TRUE)

  writeLines(c("(a:1,b:1,c:1);", "(a:2,b:2,c:2);"), path)
  trees <- read_tree(path)

  expect_s3_class(trees, "multiPhylo")
  expect_length(trees, 2L)
  expect_true(all(vapply(trees, inherits, logical(1), "phylo")))
})
