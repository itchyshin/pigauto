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
