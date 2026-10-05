# evaluate(fit) with data = NULL uses the training matrix stored on the fit
# (issue #211). gnn = FALSE keeps this off the torch path.

test_that("evaluate(fit) with data = NULL matches evaluate(fit, data = pd)", {
  set.seed(211L)
  tree <- ape::rtree(24L)
  df <- data.frame(
    row.names = tree$tip.label,
    x = ape::rTraitCont(tree, sigma = 0.4),
    y = ape::rTraitCont(tree, sigma = 0.4)
  )
  pd <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.25, seed = 211L,
                                trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, tree, splits = splits, gnn = FALSE, verbose = FALSE,
                     seed = 211L)

  expect_true(is.matrix(fit$X_scaled))
  eval_stored <- evaluate(fit)
  eval_data <- evaluate(fit, data = pd)

  expect_s3_class(eval_stored, "data.frame")
  expect_true(all(c("pigauto", "baseline") %in% eval_stored$method))
  expect_equal(eval_stored, eval_data)
})

test_that("evaluate(fit) still requires data when the fit has no X_scaled", {
  set.seed(212L)
  tree <- ape::rtree(20L)
  df <- data.frame(
    row.names = tree$tip.label,
    x = ape::rTraitCont(tree, sigma = 0.4)
  )
  pd <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.25, seed = 212L,
                                trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, tree, splits = splits, gnn = FALSE, verbose = FALSE,
                     seed = 212L)
  fit$X_scaled <- NULL
  expect_error(evaluate(fit), "data.*required")
})
