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
  spl <- make_missing_splits(pd$X_scaled, seed = 1)
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     epochs = 20L, eval_every = 10L, patience = 5L,
                     verbose = FALSE, seed = 1)
  expect_s3_class(fit, "pigauto_fit")
  expect_named(fit, c("model_state", "model_config", "graph", "baseline",
                       "norm", "species_names", "trait_names", "splits",
                       "history", "val_rmse", "test_rmse"))
})

test_that("predict.pigauto_fit returns matrix of correct dimensions", {
  td  <- make_test_data()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 1)
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     epochs = 20L, eval_every = 10L, patience = 5L,
                     verbose = FALSE, seed = 1)
  pred <- predict(fit, return_se = FALSE)
  expect_true(is.matrix(pred))
  expect_equal(nrow(pred), 40L)
  expect_equal(ncol(pred), 2L)
})

test_that("predicted values are finite", {
  td  <- make_test_data()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 2)
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     epochs = 20L, eval_every = 10L, patience = 5L,
                     verbose = FALSE, seed = 2)
  pred <- predict(fit, return_se = FALSE)
  expect_true(all(is.finite(pred)))
})

test_that("evaluate_imputation returns data.frame with expected columns", {
  td  <- make_test_data()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 3)
  bl  <- fit_baseline(pd, td$tree, splits = spl)
  ev  <- evaluate_imputation(bl$mu, pd$X_scaled, spl)
  expect_s3_class(ev, "data.frame")
  expect_true(all(c("split", "trait", "n", "rmse", "pearson_r") %in% names(ev)))
})
