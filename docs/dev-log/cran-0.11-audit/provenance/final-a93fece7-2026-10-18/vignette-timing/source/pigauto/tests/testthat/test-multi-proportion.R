test_that("preprocess_traits handles multi_proportion groups (CLR encoding)", {
  set.seed(1)
  tree <- ape::rcoal(30)
  df   <- simulate_multi_proportion_traits(tree, K = 5L, signal = 0.8, seed = 1)

  pd <- preprocess_traits(df, tree,
                          multi_proportion_groups = list(comp = names(df)))

  # Trait map has the group entry and NO per-column entries for members
  expect_true("comp" %in% names(pd$trait_map))
  expect_false(any(names(df) %in% names(pd$trait_map)))

  tm <- pd$trait_map$comp
  expect_equal(tm$type, "multi_proportion")
  expect_equal(tm$n_latent, 5L)
  expect_equal(tm$input_cols, names(df))
  expect_equal(tm$levels, names(df))
  expect_length(tm$mean, 5L)
  expect_length(tm$sd, 5L)

  # Total latent width == K
  expect_equal(ncol(pd$X_scaled), 5L)

  # Row sums of the un-z-scored CLR reconstruction should be ~0
  clr_back <- pd$X_scaled
  for (k in seq_len(tm$n_latent)) {
    clr_back[, k] <- clr_back[, k] * tm$sd[k] + tm$mean[k]
  }
  expect_true(all(abs(rowSums(clr_back)) < 1e-8))
})


test_that("multi_proportion group name cannot collide with a column", {
  set.seed(2)
  tree <- ape::rcoal(10)
  df   <- simulate_multi_proportion_traits(tree, K = 3L)

  expect_error(
    preprocess_traits(df, tree,
                      multi_proportion_groups = list(cat1 = names(df))),
    "collide with existing column"
  )
})


test_that("multi_proportion groups must have >= 2 columns", {
  set.seed(3)
  tree <- ape::rcoal(10)
  df   <- simulate_multi_proportion_traits(tree, K = 3L)

  expect_error(
    preprocess_traits(df, tree,
                      multi_proportion_groups = list(comp = "cat1")),
    "needs >= 2 columns"
  )
})


test_that("fit_baseline produces BM predictions for K CLR columns", {
  set.seed(4)
  tree <- ape::rcoal(30)
  df   <- simulate_multi_proportion_traits(tree, K = 4L, signal = 0.9, seed = 4)
  pd   <- preprocess_traits(df, tree,
                            multi_proportion_groups = list(comp = names(df)))
  splits   <- make_missing_splits(pd$X_scaled, missing_frac = 0.2,
                                  val_frac = 0.25, seed = 1)
  graph    <- build_phylo_graph(tree, k_eigen = 4L)
  baseline <- fit_baseline(pd, tree, splits = splits, graph = graph)

  expect_equal(dim(baseline$mu), c(30L, 4L))
  expect_false(any(is.na(baseline$mu)))
})


test_that("impute() end-to-end on simulated multi_proportion data", {
  skip_on_cran()
  set.seed(5)
  tree <- ape::rcoal(40)
  df   <- simulate_multi_proportion_traits(tree, K = 5L, signal = 0.7, seed = 5)

  # Introduce missingness: 20% of rows have their whole composition dropped
  # (compositional data is either complete for a row or missing — you can't
  # observe 3 of 5 components)
  miss_rows <- sample(nrow(df), round(0.2 * nrow(df)))
  df_miss <- df
  df_miss[miss_rows, ] <- NA

  result <- impute(df_miss, tree,
                   multi_proportion_groups = list(comp = names(df)),
                   epochs = 50L, verbose = FALSE, missing_frac = 0)

  # Output shape
  expect_equal(nrow(result$prediction$imputed), nrow(df))
  # Each component has its own column in $imputed
  expect_true(all(names(df) %in% names(result$prediction$imputed)))
  # Probabilities list has the group matrix
  expect_true("comp" %in% names(result$prediction$probabilities))
  pm <- result$prediction$probabilities$comp
  expect_equal(dim(pm), c(40L, 5L))
  expect_equal(colnames(pm), names(df))
  # Row sums ~ 1
  expect_true(all(abs(rowSums(pm) - 1) < 1e-6))
})
