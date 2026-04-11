test_that("mask_missing returns logical matrix of correct dimensions", {
  X <- matrix(c(1, NA, 3, NA, 5, 6), nrow = 3)
  m <- mask_missing(X)
  expect_true(is.logical(m))
  expect_equal(dim(m), dim(X))
  expect_equal(m[2, 1], FALSE)
  expect_equal(m[1, 1], TRUE)
})

test_that("make_missing_splits val and test indices are disjoint", {
  X      <- matrix(rnorm(200), nrow = 20)
  splits <- make_missing_splits(X, missing_frac = 0.30, seed = 42)
  expect_equal(length(intersect(splits$val_idx, splits$test_idx)), 0)
})

test_that("make_missing_splits produces approximately the right missing fraction", {
  X      <- matrix(rnorm(400), nrow = 20)
  splits <- make_missing_splits(X, missing_frac = 0.25, seed = 1)
  total_missing <- length(splits$val_idx) + length(splits$test_idx)
  actual_frac   <- total_missing / (20 * 20)
  expect_lt(abs(actual_frac - 0.25), 0.02)
})

test_that("make_missing_splits is deterministic given seed", {
  X  <- matrix(rnorm(100), nrow = 10)
  s1 <- make_missing_splits(X, seed = 7)
  s2 <- make_missing_splits(X, seed = 7)
  expect_equal(s1$val_idx,  s2$val_idx)
  expect_equal(s1$test_idx, s2$test_idx)
})

test_that("make_missing_splits mask is TRUE for non-missing cells", {
  X      <- matrix(rnorm(100), nrow = 10)
  splits <- make_missing_splits(X, missing_frac = 0.20, seed = 3)
  expect_false(any(splits$mask[splits$val_idx]))
  expect_false(any(splits$mask[splits$test_idx]))
  n_observed <- sum(splits$mask)
  expect_equal(n_observed, 100 - length(splits$val_idx) - length(splits$test_idx))
})

# ---- MAR / MNAR mechanism tests --------------------------------------------

test_that("MCAR mechanism is backward-compatible (default)", {
  X  <- matrix(rnorm(200), nrow = 20)
  s1 <- make_missing_splits(X, seed = 42)
  s2 <- make_missing_splits(X, seed = 42, mechanism = "MCAR")
  expect_equal(s1$val_idx, s2$val_idx)
  expect_equal(s1$test_idx, s2$test_idx)
  expect_equal(s2$mechanism, "MCAR")
})

test_that("MAR_trait produces non-uniform missingness correlated with driver", {
  set.seed(123)
  n <- 200; p <- 5
  X <- matrix(rnorm(n * p), nrow = n)
  # Make column 1 the driver — species with high col 1 should be more often missing
  X[, 1] <- seq(-3, 3, length.out = n)  # deterministic gradient

  splits <- make_missing_splits(X, missing_frac = 0.30, seed = 1,
                                mechanism = "MAR_trait",
                                mechanism_args = list(driver_col = 1, beta = 3.0))

  # Val + test indices are disjoint

  expect_equal(length(intersect(splits$val_idx, splits$test_idx)), 0)

  # Mechanism is recorded
  expect_equal(splits$mechanism, "MAR_trait")

  # Check non-uniform: rows with high driver values should have more missing cells
  all_miss <- c(splits$val_idx, splits$test_idx)
  miss_rows <- ((all_miss - 1L) %% n) + 1L
  # Top-quartile rows should have more missing cells than bottom-quartile
  top_rows <- which(X[, 1] > quantile(X[, 1], 0.75))
  bot_rows <- which(X[, 1] < quantile(X[, 1], 0.25))
  n_top <- sum(miss_rows %in% top_rows)
  n_bot <- sum(miss_rows %in% bot_rows)
  expect_gt(n_top, n_bot)  # non-uniform: top rows have more missing
})

test_that("MAR_phylo requires a tree", {
  X <- matrix(rnorm(100), nrow = 10)
  expect_error(
    make_missing_splits(X, mechanism = "MAR_phylo"),
    "tree.*required"
  )
})

test_that("MAR_phylo produces clade-structured missingness", {
  skip_if_not_installed("ape")
  set.seed(99)
  tree <- ape::rtree(100)
  n <- 100; p <- 4
  X <- matrix(rnorm(n * p), nrow = n)
  rownames(X) <- tree$tip.label

  splits <- make_missing_splits(X, missing_frac = 0.30, seed = 1,
                                mechanism = "MAR_phylo",
                                mechanism_args = list(n_clades = 2,
                                                      p_clade = 0.8,
                                                      p_base = 0.05),
                                tree = tree)

  expect_equal(splits$mechanism, "MAR_phylo")
  expect_equal(length(intersect(splits$val_idx, splits$test_idx)), 0)

  # Missingness should be concentrated — some rows have many missing cells,
  # others have few. Under MCAR with 0.30, variance across rows is lower.
  all_miss <- c(splits$val_idx, splits$test_idx)
  miss_rows <- ((all_miss - 1L) %% n) + 1L
  row_counts <- tabulate(miss_rows, nbins = n)
  # High variance = clade structure (some rows have many missing, others few)
  expect_gt(sd(row_counts), 0.5)
})

test_that("MNAR produces value-dependent missingness", {
  set.seed(77)
  n <- 500; p <- 6
  X <- matrix(rnorm(n * p), nrow = n)
  # Put extreme values in a clear top group
  X[1:50, ] <- 4.0   # very high extreme
  X[451:500, ] <- -4.0  # very low extreme

  splits <- make_missing_splits(X, missing_frac = 0.25, seed = 1,
                                mechanism = "MNAR",
                                mechanism_args = list(beta = 3.0))

  expect_equal(splits$mechanism, "MNAR")
  expect_equal(length(intersect(splits$val_idx, splits$test_idx)), 0)

  # Extreme rows should have more missing cells than central rows
  all_miss <- c(splits$val_idx, splits$test_idx)
  miss_rows <- ((all_miss - 1L) %% n) + 1L
  row_counts <- tabulate(miss_rows, nbins = n)
  extreme_rows <- c(1:50, 451:500)
  central_rows <- 200:300
  mean_extreme <- mean(row_counts[extreme_rows])
  mean_central <- mean(row_counts[central_rows])
  expect_gt(mean_extreme, mean_central)
})

test_that("MAR_trait works with trait_map", {
  skip_if_not_installed("ape")
  set.seed(42)
  n <- 50; p <- 3
  X <- matrix(rnorm(n * p), nrow = n)
  X[, 1] <- seq(-2, 2, length.out = n)  # gradient driver

  tm <- list(
    list(name = "t1", type = "continuous", latent_cols = 1L, n_latent = 1L),
    list(name = "t2", type = "continuous", latent_cols = 2L, n_latent = 1L),
    list(name = "t3", type = "continuous", latent_cols = 3L, n_latent = 1L)
  )

  splits <- make_missing_splits(X, missing_frac = 0.25, seed = 1,
                                trait_map = tm,
                                mechanism = "MAR_trait",
                                mechanism_args = list(driver_col = 1, beta = 2.0))

  expect_equal(splits$mechanism, "MAR_trait")
  expect_true(length(splits$val_idx_trait) > 0)
  expect_true(length(splits$test_idx_trait) > 0)
  expect_equal(length(intersect(splits$val_idx, splits$test_idx)), 0)
})

test_that("MNAR works with trait_map", {
  set.seed(55)
  n <- 50; p <- 3
  X <- matrix(rnorm(n * p), nrow = n)

  tm <- list(
    list(name = "t1", type = "continuous", latent_cols = 1L, n_latent = 1L),
    list(name = "t2", type = "continuous", latent_cols = 2L, n_latent = 1L),
    list(name = "t3", type = "continuous", latent_cols = 3L, n_latent = 1L)
  )

  splits <- make_missing_splits(X, missing_frac = 0.25, seed = 1,
                                trait_map = tm,
                                mechanism = "MNAR",
                                mechanism_args = list(beta = 2.0))

  expect_equal(splits$mechanism, "MNAR")
  expect_equal(length(intersect(splits$val_idx, splits$test_idx)), 0)
})
