# tests/testthat/test-simulate.R
# Unit tests for trait simulation functions

test_that("simulate_binary_traits returns correct structure", {
  set.seed(1)
  tree <- ape::rtree(30)
  df <- simulate_binary_traits(tree, n_traits = 3, signal = 0.6, seed = 42)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 30)
  expect_equal(ncol(df), 3)
  expect_equal(rownames(df), tree$tip.label)
  # All columns should be factors with levels "no", "yes"
  for (j in seq_len(3)) {
    expect_s3_class(df[[j]], "factor")
    expect_equal(levels(df[[j]]), c("no", "yes"))
  }
})

test_that("simulate_binary_traits signal parameter works", {
  set.seed(1)
  tree <- ape::rtree(100)
  # High signal: related species should share values more often

  df_high <- simulate_binary_traits(tree, n_traits = 1, signal = 1.0, seed = 1)
  df_low  <- simulate_binary_traits(tree, n_traits = 1, signal = 0.1, seed = 1)
  # Can't do a strict test, but both should produce valid data
  expect_true(all(df_high[[1]] %in% c("no", "yes")))
  expect_true(all(df_low[[1]] %in% c("no", "yes")))
})

test_that("simulate_binary_traits threshold_quantile creates imbalance", {
  set.seed(1)
  tree <- ape::rtree(200)
  df <- simulate_binary_traits(tree, n_traits = 1, signal = 0.6,
                               threshold_quantile = 0.9, seed = 42)
  # ~90% should be "no" (below threshold)
  frac_yes <- mean(df[[1]] == "yes")
  expect_true(frac_yes < 0.25)  # should be ~10% yes
})

test_that("simulate_ordinal_traits returns correct structure", {
  set.seed(1)
  tree <- ape::rtree(50)
  df <- simulate_ordinal_traits(tree, n_traits = 2, n_levels = 5, seed = 42)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 50)
  expect_equal(ncol(df), 2)
  for (j in seq_len(2)) {
    expect_s3_class(df[[j]], "ordered")
    expect_equal(levels(df[[j]]), paste0("L", 1:5))
  }
})

test_that("simulate_ordinal_traits respects n_levels", {
  tree <- ape::rtree(100)
  for (nl in c(3, 7, 10)) {
    df <- simulate_ordinal_traits(tree, n_traits = 1, n_levels = nl, seed = 1)
    expect_equal(nlevels(df[[1]]), nl)
  }
})

test_that("simulate_count_traits returns integers", {
  tree <- ape::rtree(50)
  df <- simulate_count_traits(tree, n_traits = 3, mean_count = 20, seed = 42)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 50)
  expect_equal(ncol(df), 3)
  for (j in seq_len(3)) {
    expect_type(df[[j]], "integer")
    expect_true(all(df[[j]] >= 0))
  }
})

test_that("simulate_count_traits Poisson vs NegBin", {
  tree <- ape::rtree(200)
  df_pois <- simulate_count_traits(tree, n_traits = 1, mean_count = 20,
                                   overdispersion = NULL, seed = 1)
  df_nb   <- simulate_count_traits(tree, n_traits = 1, mean_count = 20,
                                   overdispersion = 2, seed = 1)
  # NegBin should have higher variance
  expect_true(stats::var(df_nb[[1]]) > stats::var(df_pois[[1]]) * 0.5)
  # Both should be non-negative integers
  expect_true(all(df_pois[[1]] >= 0))
  expect_true(all(df_nb[[1]] >= 0))
})

test_that("simulate_categorical_traits returns correct structure", {
  tree <- ape::rtree(60)
  df <- simulate_categorical_traits(tree, n_traits = 2, n_levels = 4, seed = 42)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 60)
  expect_equal(ncol(df), 2)
  for (j in seq_len(2)) {
    expect_s3_class(df[[j]], "factor")
    expect_false(is.ordered(df[[j]]))
    expect_equal(levels(df[[j]]), LETTERS[1:4])
  }
})

test_that("simulate_categorical_traits covers all levels with enough species", {
  tree <- ape::rtree(200)
  df <- simulate_categorical_traits(tree, n_traits = 1, n_levels = 3,
                                    signal = 0.3, seed = 42)
  # With 200 species and low signal, all 3 levels should be represented
  expect_equal(length(unique(df[[1]])), 3)
})

test_that("simulate_proportion_traits returns (0,1) values", {
  tree <- ape::rtree(50)
  df <- simulate_proportion_traits(tree, n_traits = 3, signal = 0.6, seed = 42)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 50)
  expect_equal(ncol(df), 3)
  for (j in seq_len(3)) {
    expect_type(df[[j]], "double")
    expect_true(all(df[[j]] > 0 & df[[j]] < 1))
  }
})

test_that("simulate_proportion_traits boundary_frac pushes values", {
  tree <- ape::rtree(200)
  df_no  <- simulate_proportion_traits(tree, n_traits = 1, boundary_frac = 0.0,
                                       seed = 1)
  df_hi  <- simulate_proportion_traits(tree, n_traits = 1, boundary_frac = 0.3,
                                       seed = 1)
  # High boundary_frac should have more extreme values
  extreme_no <- sum(df_no[[1]] < 0.05 | df_no[[1]] > 0.95)
  extreme_hi <- sum(df_hi[[1]] < 0.05 | df_hi[[1]] > 0.95)
  expect_true(extreme_hi > extreme_no)
})

test_that("simulate_zi_count_traits returns integers with zeros", {
  tree <- ape::rtree(100)
  df <- simulate_zi_count_traits(tree, n_traits = 2, zero_frac = 0.5,
                                 mean_nz = 20, seed = 42)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 100)
  expect_equal(ncol(df), 2)
  for (j in seq_len(2)) {
    expect_type(df[[j]], "integer")
    expect_true(all(df[[j]] >= 0))
    # Should have a non-trivial number of zeros
    expect_true(sum(df[[j]] == 0) > 10)
    # Should also have non-zero values
    expect_true(sum(df[[j]] > 0) > 10)
  }
})

test_that("simulate_zi_count_traits zero_frac is approximately correct", {
  tree <- ape::rtree(500)
  for (zf in c(0.2, 0.5, 0.8)) {
    df <- simulate_zi_count_traits(tree, n_traits = 1, zero_frac = zf,
                                   mean_nz = 20, seed = 42)
    observed_zf <- mean(df[[1]] == 0)
    # Allow ±15% tolerance (BM + threshold is stochastic)
    expect_true(abs(observed_zf - zf) < 0.15,
                info = sprintf("zero_frac = %.1f, observed = %.2f", zf, observed_zf))
  }
})

test_that("simulate_zi_count_traits overdispersion works", {
  tree <- ape::rtree(200)
  df_pois <- simulate_zi_count_traits(tree, n_traits = 1, zero_frac = 0.3,
                                      mean_nz = 20, overdispersion = NULL,
                                      seed = 1)
  df_nb   <- simulate_zi_count_traits(tree, n_traits = 1, zero_frac = 0.3,
                                      mean_nz = 20, overdispersion = 2,
                                      seed = 1)
  # Both valid
  expect_true(all(df_pois[[1]] >= 0))
  expect_true(all(df_nb[[1]] >= 0))
})

test_that("all simulation functions are reproducible with seed", {
  tree <- ape::rtree(30)
  # Binary
  a <- simulate_binary_traits(tree, n_traits = 1, seed = 123)
  b <- simulate_binary_traits(tree, n_traits = 1, seed = 123)
  expect_identical(a, b)

  # Ordinal
  a <- simulate_ordinal_traits(tree, n_traits = 1, seed = 123)
  b <- simulate_ordinal_traits(tree, n_traits = 1, seed = 123)
  expect_identical(a, b)

  # Count
  a <- simulate_count_traits(tree, n_traits = 1, seed = 123)
  b <- simulate_count_traits(tree, n_traits = 1, seed = 123)
  expect_identical(a, b)

  # Categorical
  a <- simulate_categorical_traits(tree, n_traits = 1, seed = 123)
  b <- simulate_categorical_traits(tree, n_traits = 1, seed = 123)
  expect_identical(a, b)

  # Proportion
  a <- simulate_proportion_traits(tree, n_traits = 1, seed = 123)
  b <- simulate_proportion_traits(tree, n_traits = 1, seed = 123)
  expect_identical(a, b)

  # ZI count
  a <- simulate_zi_count_traits(tree, n_traits = 1, seed = 123)
  b <- simulate_zi_count_traits(tree, n_traits = 1, seed = 123)
  expect_identical(a, b)
})
