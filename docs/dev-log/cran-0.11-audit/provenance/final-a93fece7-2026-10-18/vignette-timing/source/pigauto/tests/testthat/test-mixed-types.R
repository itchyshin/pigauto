# Tests for mixed-type trait support

make_mixed_data <- function(n = 30, seed = 100) {
  set.seed(seed)
  tree <- ape::rtree(n)
  sp   <- tree$tip.label

  df <- data.frame(
    row.names = sp,
    mass    = abs(stats::rnorm(n)) + 0.5,              # continuous
    clutch  = as.integer(stats::rpois(n, 3) + 1L),     # count
    migr    = factor(sample(c("no", "yes"), n, replace = TRUE)),  # binary
    diet    = factor(sample(c("herb", "carn", "omni"), n,
                            replace = TRUE)),           # categorical
    threat  = ordered(sample(c("LC", "VU", "EN"), n,
                             replace = TRUE),
                      levels = c("LC", "VU", "EN"))     # ordinal
  )

  list(tree = tree, df = df)
}

test_that("preprocess_traits detects all 5 trait types", {
  td <- make_mixed_data()
  pd <- preprocess_traits(td$df, td$tree)

  expect_equal(pd$trait_map$mass$type,   "continuous")
  expect_equal(pd$trait_map$clutch$type, "count")
  expect_equal(pd$trait_map$migr$type,   "binary")
  expect_equal(pd$trait_map$diet$type,   "categorical")
  expect_equal(pd$trait_map$threat$type, "ordinal")

  # Latent dimensions: 1 + 1 + 1 + 3 + 1 = 7
  expect_equal(pd$p_latent, 7L)
  expect_equal(ncol(pd$X_scaled), 7L)
})

test_that("make_missing_splits with trait_map works at trait level", {
  td <- make_mixed_data()
  pd <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, missing_frac = 0.25,
                             seed = 1, trait_map = pd$trait_map)

  expect_true(length(spl$val_idx) > 0)
  expect_true(length(spl$test_idx) > 0)
  expect_equal(length(intersect(spl$val_idx, spl$test_idx)), 0)
  expect_true(!is.null(spl$val_idx_trait))
})

test_that("fit_baseline handles all trait types", {
  td <- make_mixed_data()
  pd <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 1, trait_map = pd$trait_map)
  bl <- fit_baseline(pd, td$tree, splits = spl)

  expect_true(is.matrix(bl$mu))
  expect_true(is.matrix(bl$se))
  expect_equal(dim(bl$mu), dim(pd$X_scaled))
  expect_true(all(is.finite(bl$mu)))
})

test_that("evaluate_imputation returns type-specific metrics", {
  td <- make_mixed_data()
  pd <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 1, trait_map = pd$trait_map)
  bl <- fit_baseline(pd, td$tree, splits = spl)

  ev <- evaluate_imputation(bl$mu, pd$X_scaled, spl, trait_map = pd$trait_map)
  expect_s3_class(ev, "data.frame")
  expect_true("type" %in% names(ev))

  # Check we get metrics for multiple types
  types_found <- unique(ev$type)
  expect_true(length(types_found) > 1L)
})

test_that("encode and decode round-trip preserves continuous values", {
  set.seed(42)
  tree <- ape::rtree(10)
  df <- data.frame(row.names = tree$tip.label,
                   x = abs(rnorm(10)) + 1)
  pd <- preprocess_traits(df, tree, log_transform = TRUE)

  # Manually decode
  tm <- pd$trait_map$x
  decoded <- pd$X_scaled[, tm$latent_cols] * tm$sd + tm$mean
  if (tm$log_transform) decoded <- exp(decoded)

  # Should be close to original
  orig <- df$x[match(tree$tip.label, rownames(df))]
  ok <- !is.na(orig) & !is.na(decoded)
  expect_equal(unname(decoded[ok]), unname(orig[ok]), tolerance = 1e-10)
})
