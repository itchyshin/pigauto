# Small synthetic dataset helpers
make_test_data <- function(n = 40, p = 2, seed = 42) {
  set.seed(seed)
  tree <- ape::rtree(n)
  sp   <- tree$tip.label
  df   <- data.frame(
    row.names = sp,
    tr1 = abs(stats::rnorm(n)) + 0.5,
    tr2 = abs(stats::rnorm(n)) + 0.5
  )
  list(tree = tree, df = df)
}

test_that("fit_pigauto returns a pigauto_fit object", {
  td  <- make_test_data()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 1, trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     epochs = 20L, eval_every = 10L, patience = 5L,
                     verbose = FALSE, seed = 1)
  expect_s3_class(fit, "pigauto_fit")
  expect_true("model_state" %in% names(fit))
  expect_true("model_config" %in% names(fit))
  expect_true("trait_map" %in% names(fit))
  expect_true("latent_names" %in% names(fit))
  expect_true("val_rmse" %in% names(fit))
  expect_true("test_rmse" %in% names(fit))

  # GPU memory fix (2026-04-21): model_state must be on CPU so the
  # returned fit object doesn't hold GPU tensor refs that pin ~40 GB
  # of training-time activations and break predict() at large n.
  state_devices <- vapply(fit$model_state,
                           function(t) as.character(t$device$type),
                           character(1))
  expect_true(all(state_devices == "cpu"),
              info = paste("model_state tensors must be on CPU; got:",
                            paste(unique(state_devices), collapse = ", ")))
})

test_that("predict.pigauto_fit returns pigauto_pred with correct dimensions", {
  td  <- make_test_data()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 1, trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     epochs = 20L, eval_every = 10L, patience = 5L,
                     verbose = FALSE, seed = 1)
  pred <- predict(fit, return_se = TRUE)

  expect_true(is.data.frame(pred$imputed))
  expect_equal(nrow(pred$imputed), 40L)
  expect_equal(ncol(pred$imputed), 2L)
  expect_true(is.matrix(pred$se))
  expect_equal(dim(pred$se), c(40L, 2L))
})

test_that("predicted values are finite", {
  td  <- make_test_data()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 2, trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     epochs = 20L, eval_every = 10L, patience = 5L,
                     verbose = FALSE, seed = 2)
  pred <- predict(fit, return_se = FALSE)
  imp  <- pred$imputed
  expect_true(all(sapply(imp, function(col) all(is.finite(col)))))
})

test_that("evaluate_imputation returns data.frame with expected columns", {
  td  <- make_test_data()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 3, trait_map = pd$trait_map)
  bl  <- fit_baseline(pd, td$tree, splits = spl)
  ev  <- evaluate_imputation(bl$mu, pd$X_scaled, spl, trait_map = pd$trait_map)
  expect_s3_class(ev, "data.frame")
  expect_true(all(c("split", "trait", "type", "n", "rmse") %in% names(ev)))
})

test_that("predict with n_imputations > 1 returns multiple datasets", {
  td  <- make_test_data()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 4, trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     epochs = 20L, eval_every = 10L, patience = 5L,
                     verbose = FALSE, seed = 4)
  pred <- predict(fit, n_imputations = 3L)

  expect_s3_class(pred, "pigauto_pred")
  expect_equal(pred$n_imputations, 3L)
  expect_true(is.list(pred$imputed_datasets))
  expect_equal(length(pred$imputed_datasets), 3L)
  expect_true(is.matrix(pred$se))
})


# ---- Multi-observation per species tests ------------------------------------

test_that("fit_pigauto and predict work with multi-obs data", {
  set.seed(50)
  tree <- ape::rtree(20)

  # 3 observations per species (60 rows, 20 species)
  df <- data.frame(
    species = rep(tree$tip.label, each = 3),
    tr1 = abs(stats::rnorm(60)) + 0.5,
    tr2 = abs(stats::rnorm(60)) + 0.5
  )

  pd  <- preprocess_traits(df, tree, species_col = "species")
  expect_true(pd$multi_obs)
  expect_equal(pd$n_obs, 60L)
  expect_equal(pd$n_species, 20L)

  spl <- make_missing_splits(pd$X_scaled, seed = 50, trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, tree, splits = spl,
                     epochs = 20L, eval_every = 10L, patience = 5L,
                     verbose = FALSE, seed = 50)

  expect_s3_class(fit, "pigauto_fit")
  expect_true(fit$multi_obs)

  pred <- predict(fit, return_se = TRUE)
  expect_s3_class(pred, "pigauto_pred")
  expect_equal(nrow(pred$imputed), 60L)  # obs-level output
  expect_true(pred$multi_obs)
  expect_equal(length(pred$obs_species), 60L)
  expect_true(is.matrix(pred$se))
  expect_equal(nrow(pred$se), 60L)
  expect_true(all(is.finite(pred$imputed$tr1)))
})


# ---- impute() convenience wrapper tests ------------------------------------

test_that("impute() works end-to-end with rownames (single-obs)", {
  td <- make_test_data(n = 20, seed = 60)
  result <- impute(td$df, td$tree, epochs = 20L, verbose = FALSE, seed = 60)
  expect_s3_class(result, "pigauto_result")
  expect_s3_class(result$prediction, "pigauto_pred")
  expect_s3_class(result$fit, "pigauto_fit")
  expect_true(nrow(result$prediction$imputed) == 20L)
})

test_that("impute() works with species_col (multi-obs)", {
  set.seed(61)
  tree <- ape::rtree(15)
  df <- data.frame(
    sp = rep(tree$tip.label, each = 2),
    x = abs(stats::rnorm(30)) + 0.5
  )
  result <- impute(df, tree, species_col = "sp",
                   epochs = 20L, verbose = FALSE, seed = 61)
  expect_s3_class(result, "pigauto_result")
  expect_equal(nrow(result$prediction$imputed), 30L)  # obs-level
  expect_true(result$prediction$multi_obs)
})


# ---- Tests for new features: attention, calibration, conformal ---------------

test_that("fit_pigauto with attention produces valid fit", {
  td  <- make_test_data()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 1, trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     epochs = 30L, eval_every = 10L, patience = 5L,
                     use_attention = TRUE, verbose = FALSE, seed = 1)
  expect_s3_class(fit, "pigauto_fit")
  expect_true(fit$model_config$use_attention)
  pred <- predict(fit)
  expect_s3_class(pred, "pigauto_pred")
  expect_equal(nrow(pred$imputed), 40L)
})

test_that("fit_pigauto stores calibrated_gates and conformal_scores", {
  td  <- make_test_data()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 1, trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     epochs = 30L, eval_every = 10L, patience = 5L,
                     verbose = FALSE, seed = 1)
  expect_true(!is.null(fit$calibrated_gates))
  expect_length(fit$calibrated_gates, ncol(pd$X_scaled))
  expect_true(all(fit$calibrated_gates >= 0))
  expect_true(!is.null(fit$conformal_scores))
  expect_true(all(fit$conformal_scores[!is.na(fit$conformal_scores)] > 0))
})

test_that("predict with conformal scores returns intervals", {
  td  <- make_test_data()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 1, trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     epochs = 30L, eval_every = 10L, patience = 5L,
                     verbose = FALSE, seed = 1)
  pred <- predict(fit, return_se = TRUE)
  expect_true(!is.null(pred$conformal_lower))
  expect_true(!is.null(pred$conformal_upper))
  expect_equal(dim(pred$conformal_lower), dim(pred$conformal_upper))
  # Lower should be <= upper for all non-NA values
  valid <- !is.na(pred$conformal_lower) & !is.na(pred$conformal_upper)
  expect_true(all(pred$conformal_lower[valid] <= pred$conformal_upper[valid]))
})

test_that("phylo label propagation gives species-specific discrete baselines", {
  # Create mixed-type data
  data(avonet300, tree300, package = "pigauto")
  traits <- avonet300
  rownames(traits) <- traits$Species_Key
  traits$Species_Key <- NULL
  pd <- preprocess_traits(traits, tree300)
  spl <- make_missing_splits(pd$X_scaled, seed = 1, trait_map = pd$trait_map)
  bl <- fit_baseline(pd, tree300, spl)

  # Check that binary/categorical baselines vary across species
  # (old baseline was constant across species)
  for (tm in pd$trait_map) {
    if (tm$type %in% c("binary", "categorical")) {
      lc <- tm$latent_cols[1]
      vals <- bl$mu[, lc]
      # Should have variation (not all the same value)
      expect_true(sd(vals) > 0,
                  info = paste("Trait", tm$name, "should have species-specific baseline"))
    }
  }
})

test_that("attention + calibration + conformal work with mixed types", {
  data(avonet300, tree300, package = "pigauto")
  traits <- avonet300[1:50, ]  # subset for speed
  rownames(traits) <- traits$Species_Key
  traits$Species_Key <- NULL
  subtree <- ape::keep.tip(tree300, rownames(traits))

  pd  <- preprocess_traits(traits, subtree)
  spl <- make_missing_splits(pd$X_scaled, seed = 1, trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, subtree, splits = spl,
                     epochs = 30L, eval_every = 10L, patience = 5L,
                     use_attention = TRUE, verbose = FALSE, seed = 1)

  expect_true(fit$model_config$use_attention)
  expect_true(!is.null(fit$calibrated_gates))
  expect_true(!is.null(fit$conformal_scores))

  pred <- predict(fit, return_se = TRUE)
  expect_s3_class(pred, "pigauto_pred")
  expect_true(!is.null(pred$conformal_lower))
  expect_true(!is.null(pred$conformal_upper))
})
