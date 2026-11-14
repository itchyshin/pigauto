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

test_that("aggregate_to_species with soft_aggregate preserves proportions", {
  skip_if_not_installed("Rphylopars")
  set.seed(10)
  tree <- ape::rtree(6)
  df <- data.frame(
    species = rep(tree$tip.label, each = 4),
    y = factor(c(rep("B", 3), "A",    # sp1: 3/4 = 0.75 class-B
                 rep("A", 4),           # sp2: 0/4 = 0.0 class-B
                 rep("B", 4),           # sp3: 4/4 = 1.0 class-B
                 "A","A","B","B",       # sp4: 2/4 = 0.5 class-B
                 rep("A", 4),           # sp5: 0/4
                 rep("B", 4)),          # sp6: 4/4
               levels = c("A","B"))
  )
  pd <- preprocess_traits(df, tree, species_col = "species")

  agg_hard <- aggregate_to_species(pd)
  agg_soft <- aggregate_to_species(pd, soft_aggregate = TRUE)

  y_col <- which(colnames(pd$X_scaled) == "y")

  # Hard: sp1 (0.75) -> 1, sp4 (0.5) -> 1 (>= threshold)
  expect_equal(agg_hard$X_species[1, y_col], 1, tolerance = 1e-9)
  expect_equal(agg_hard$X_species[4, y_col], 1, tolerance = 1e-9)

  # Soft: sp1 -> 0.75, sp4 -> 0.5 (proportions preserved)
  expect_equal(agg_soft$X_species[1, y_col], 0.75, tolerance = 1e-9)
  expect_equal(agg_soft$X_species[4, y_col], 0.5, tolerance = 1e-9)
  # Pure cases: same in both modes
  expect_equal(agg_soft$X_species[2, y_col], 0, tolerance = 1e-9)
  expect_equal(agg_soft$X_species[3, y_col], 1, tolerance = 1e-9)

  # is_proportion_col flag
  expect_true(agg_soft$is_proportion_col[y_col])
  expect_false(agg_hard$is_proportion_col[y_col])
})

test_that("build_liability_matrix with soft_aggregate uses soft E-step on proportions", {
  skip_if_not_installed("Rphylopars")
  set.seed(11)
  tree <- ape::rtree(8)
  df <- data.frame(
    species = rep(tree$tip.label, each = 4),
    x = rnorm(32),
    y = factor(c("A","A","B","B",   # sp1: p=0.5 (max ambiguity)
                 "B","B","B","B",   # sp2: p=1.0 (pure)
                 sample(c("A","B"), 24, replace = TRUE)),
               levels = c("A","B"))
  )
  pd <- preprocess_traits(df, tree, species_col = "species")

  out_hard <- build_liability_matrix(pd, splits = NULL, soft_aggregate = FALSE)
  out_soft <- build_liability_matrix(pd, splits = NULL, soft_aggregate = TRUE)

  y_col <- which(colnames(out_soft$X_liab) == "y")
  # sp1: p=0.5 with N(0,1) prior => soft posterior mean near 0
  expect_lt(abs(out_soft$X_liab[1, y_col]), 1e-6)
  # sp2: p=1.0 => same as hard (boundary case)
  hard_val <- estep_liability_binary(y = 1L, mu_prior = 0, sd_prior = 1)$mean
  expect_equal(out_soft$X_liab[2, y_col], hard_val, tolerance = 1e-6)
  # Hard mode: sp1 gets thresholded to 1, so estep_liability_binary(1)
  expect_equal(out_hard$X_liab[1, y_col], hard_val, tolerance = 1e-6)
  # Hard and soft differ on sp1
  expect_gt(abs(out_hard$X_liab[1, y_col] - out_soft$X_liab[1, y_col]), 0.5)
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

test_that("fit_joint_threshold_baseline propagates soft_aggregate", {
  skip_if_not_installed("Rphylopars")
  set.seed(20)
  tree <- ape::rtree(20)
  df <- data.frame(
    species = rep(tree$tip.label, each = 3),
    x = rnorm(60),
    y = factor(sample(c("A","B"), 60, TRUE), levels = c("A","B"))
  )
  pd <- preprocess_traits(df, tree, species_col = "species")
  res_hard <- fit_joint_threshold_baseline(pd, tree, splits = NULL, soft_aggregate = FALSE)
  res_soft <- fit_joint_threshold_baseline(pd, tree, splits = NULL, soft_aggregate = TRUE)
  expect_equal(dim(res_hard$mu_liab), dim(res_soft$mu_liab))
  y_col <- which(colnames(res_hard$mu_liab) == "y")
  diff <- max(abs(res_hard$mu_liab[, y_col] - res_soft$mu_liab[, y_col]))
  expect_gt(diff, 1e-4)
})

test_that("fit_ovr_categorical_fits propagates soft_aggregate", {
  skip_if_not_installed("Rphylopars")
  set.seed(21)
  tree <- ape::rtree(20)
  df <- data.frame(
    species = rep(tree$tip.label, each = 3),
    x = rnorm(60),
    z = factor(sample(c("P","Q","R"), 60, TRUE), levels = c("P","Q","R"))
  )
  pd <- preprocess_traits(df, tree, species_col = "species")
  probs_hard <- fit_ovr_categorical_fits(pd, tree, trait_name = "z",
                                           splits = NULL, soft_aggregate = FALSE)
  probs_soft <- fit_ovr_categorical_fits(pd, tree, trait_name = "z",
                                           splits = NULL, soft_aggregate = TRUE)
  expect_equal(dim(probs_hard), dim(probs_soft))
  diff <- max(abs(probs_hard - probs_soft), na.rm = TRUE)
  expect_gt(diff, 1e-4)
})

test_that("fit_baseline respects multi_obs_aggregation = 'soft'", {
  skip_if_not_installed("Rphylopars")
  set.seed(22)
  tree <- ape::rtree(25)
  df <- data.frame(
    species = rep(tree$tip.label, each = 3),
    x = rnorm(75),
    y = factor(sample(c("A","B"), 75, TRUE), levels = c("A","B"))
  )
  pd <- preprocess_traits(df, tree, species_col = "species")
  bl_hard <- fit_baseline(pd, tree, splits = NULL, multi_obs_aggregation = "hard")
  bl_soft <- fit_baseline(pd, tree, splits = NULL, multi_obs_aggregation = "soft")
  expect_equal(dim(bl_hard$mu), dim(bl_soft$mu))
  y_col <- which(colnames(pd$X_scaled) == "y")
  diff <- max(abs(bl_hard$mu[, y_col] - bl_soft$mu[, y_col]))
  expect_gt(diff, 1e-4)
  expect_true(all(is.finite(bl_hard$mu[, y_col])))
  expect_true(all(is.finite(bl_soft$mu[, y_col])))
})

test_that("multi_obs_aggregation default 'hard' is identical to explicit 'hard'", {
  skip_if_not_installed("Rphylopars")
  set.seed(23)
  tree <- ape::rtree(20)
  df <- data.frame(
    species = rep(tree$tip.label, each = 3),
    x = rnorm(60),
    y = factor(sample(c("A","B"), 60, TRUE), levels = c("A","B")),
    z = factor(sample(c("P","Q","R","S"), 60, TRUE), levels = c("P","Q","R","S"))
  )
  pd <- preprocess_traits(df, tree, species_col = "species")
  bl_default <- fit_baseline(pd, tree, splits = NULL)
  bl_hard    <- fit_baseline(pd, tree, splits = NULL, multi_obs_aggregation = "hard")
  expect_equal(bl_default$mu, bl_hard$mu, tolerance = 1e-12)
  expect_equal(bl_default$se, bl_hard$se, tolerance = 1e-12)
})
