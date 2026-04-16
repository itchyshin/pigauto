test_that("aggregate_to_species produces n_species rows with correct per-type semantics", {
  skip_if_not_installed("Rphylopars")
  set.seed(1)
  tree <- ape::rtree(6)
  df <- data.frame(
    species = rep(tree$tip.label, each = 3),
    x = rnorm(18),
    y = factor(sample(c("A", "B"), 18, TRUE), levels = c("A", "B")),
    z = factor(sample(c("P", "Q", "R"), 18, TRUE), levels = c("P", "Q", "R"))
  )
  df$x[c(5, 6)] <- NA
  pd <- preprocess_traits(df, tree, species_col = "species")
  expect_true(isTRUE(pd$multi_obs))
  expect_equal(pd$n_species, 6L)
  expect_equal(pd$n_obs, 18L)

  agg <- aggregate_to_species(pd)
  expect_equal(nrow(agg$X_species), 6L)
  expect_equal(ncol(agg$X_species), ncol(pd$X_scaled))

  x_col <- which(colnames(pd$X_scaled) == "x")
  for (s in seq_len(6)) {
    obs_rows <- which(pd$obs_to_species == s)
    expected <- mean(pd$X_scaled[obs_rows, x_col], na.rm = TRUE)
    if (is.nan(expected)) expected <- NA_real_
    expect_equal(agg$X_species[s, x_col], expected, tolerance = 1e-9)
  }

  y_col <- which(colnames(pd$X_scaled) == "y")
  for (s in seq_len(6)) {
    obs_rows <- which(pd$obs_to_species == s)
    prop <- mean(pd$X_scaled[obs_rows, y_col], na.rm = TRUE)
    expected <- if (is.nan(prop)) NA_real_ else as.numeric(prop >= 0.5)
    expect_equal(agg$X_species[s, y_col], expected, tolerance = 1e-9)
  }

  cat_col_names <- grep("^z=", colnames(pd$X_scaled), value = TRUE)
  cat_col_idx   <- which(colnames(pd$X_scaled) %in% cat_col_names)
  for (s in seq_len(6)) {
    obs_rows <- which(pd$obs_to_species == s)
    props <- colMeans(pd$X_scaled[obs_rows, cat_col_idx, drop = FALSE], na.rm = TRUE)
    if (all(is.nan(props))) {
      expect_true(all(is.na(agg$X_species[s, cat_col_idx])))
    } else {
      expected <- rep(0, length(cat_col_idx))
      expected[which.max(props)] <- 1
      expect_equal(unname(agg$X_species[s, cat_col_idx]), expected, tolerance = 1e-9)
    }
  }
})

test_that("build_liability_matrix accepts multi-obs data", {
  skip_if_not_installed("Rphylopars")
  set.seed(2)
  tree <- ape::rtree(10)
  df <- data.frame(
    species = rep(tree$tip.label, each = 3),
    x = rnorm(30),
    y = factor(sample(c("A", "B"), 30, TRUE), levels = c("A", "B"))
  )
  pd <- preprocess_traits(df, tree, species_col = "species")
  out <- build_liability_matrix(pd, splits = NULL)
  expect_equal(nrow(out$X_liab), 10L)
  expect_true(all(is.finite(out$X_liab)))
})

test_that("fit_ovr_categorical_fits runs on multi-obs data", {
  skip_if_not_installed("Rphylopars")
  set.seed(3)
  tree <- ape::rtree(20)
  df <- data.frame(
    species = rep(tree$tip.label, each = 3),
    x = rnorm(60),
    z = factor(sample(c("P", "Q", "R"), 60, TRUE), levels = c("P", "Q", "R"))
  )
  pd <- preprocess_traits(df, tree, species_col = "species")
  probs <- fit_ovr_categorical_fits(pd, tree, trait_name = "z", splits = NULL)
  expect_equal(dim(probs), c(20L, 3L))
  finite_vals <- probs[!is.na(probs) & is.finite(probs)]
  expect_true(all(finite_vals >= 0 & finite_vals <= 1))
})

test_that("fit_joint_mvn_baseline runs on multi-obs data", {
  skip_if_not_installed("Rphylopars")
  set.seed(4)
  tree <- ape::rtree(20)
  df <- data.frame(
    species = rep(tree$tip.label, each = 3),
    x1 = rnorm(60),
    x2 = rnorm(60)
  )
  pd <- preprocess_traits(df, tree, species_col = "species")
  res <- fit_joint_mvn_baseline(pd, tree, splits = NULL)
  expect_equal(nrow(res$mu), 20L)
  expect_equal(ncol(res$mu), 2L)
  expect_true(all(is.finite(res$mu)))
  expect_true(all(res$se >= 0))
})

test_that("fit_baseline uses Level-C paths on multi-obs data", {
  skip_if_not_installed("Rphylopars")
  set.seed(5)
  tree <- ape::rtree(30)
  df <- data.frame(
    species = rep(tree$tip.label, each = 3),
    x = rnorm(90),
    y = factor(sample(c("A", "B"), 90, TRUE), levels = c("A", "B")),
    z = factor(sample(c("P", "Q", "R", "S"), 90, TRUE),
               levels = c("P", "Q", "R", "S"))
  )
  pd <- preprocess_traits(df, tree, species_col = "species")
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.2,
                                 seed = 5, trait_map = pd$trait_map)
  bl <- fit_baseline(pd, tree, splits = splits)
  expect_equal(nrow(bl$mu), 30L)
  expect_true(all(is.finite(bl$mu[, "x"])))
  expect_true(all(bl$mu[, "y"] >= qlogis(0.01) - 1e-6))
  expect_true(all(bl$mu[, "y"] <= qlogis(0.99) + 1e-6))
  cat_cols <- grep("^z=", colnames(pd$X_scaled), value = TRUE)
  row_sums <- rowSums(exp(bl$mu[, cat_cols]))
  expect_true(all(abs(row_sums - 1) < 1e-6))
})
