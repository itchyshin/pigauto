test_that("fit_ovr_categorical_fits returns K probability vectors per tip", {
  skip_if_not_installed("Rphylopars")
  set.seed(100)
  tree <- ape::rtree(30)
  df <- data.frame(
    x = rnorm(30),
    z = factor(sample(c("P", "Q", "R"), 30, TRUE), levels = c("P", "Q", "R")),
    row.names = tree$tip.label
  )
  df$z[c(5, 15, 25)] <- NA
  pd <- preprocess_traits(df, tree)

  probs <- fit_ovr_categorical_fits(pd, tree, trait_name = "z", splits = NULL)

  expect_equal(dim(probs), c(30, 3))
  expect_equal(colnames(probs), c("z=P", "z=Q", "z=R"))
  expect_true(all(is.finite(probs) | is.na(probs)))
  finite_vals <- probs[!is.na(probs) & is.finite(probs)]
  expect_true(all(finite_vals >= 0 & finite_vals <= 1))
})
