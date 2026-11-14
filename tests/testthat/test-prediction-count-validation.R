make_prediction_count_fit <- function(gnn = TRUE) {
  withr::local_seed(241)
  tree <- ape::rtree(12)
  traits <- data.frame(x = sin(seq_len(12) / 3), row.names = tree$tip.label)
  traits$x[c(2, 5)] <- NA_real_
  data <- preprocess_traits(traits, tree, trait_types = c(x = "continuous"),
                            log_transform = FALSE)
  fit_pigauto(data, tree, splits = NULL, gnn = gnn, epochs = 1L,
              eval_every = 1L, hidden_dim = 8L, n_gnn_layers = 1L,
              dropout = 0, edge_dropout = 0, use_transformer_blocks = FALSE,
              refine_steps = 1L, predict_method = "per_column",
              verbose = FALSE, seed = 41)
}

test_that("shared positive integer validation handles integer limits", {
  expect_identical(.check_positive_integer(1, "count"), 1L)
  expect_identical(.check_positive_integer(2L, "count"), 2L)
  expect_identical(.check_positive_integer(.Machine$integer.max, "count"),
                   .Machine$integer.max)
  for (bad in list(0, -1, 1.5, NA_real_, Inf, c(1, 2), "2", integer(),
                   as.double(.Machine$integer.max) + 1, 1e100)) {
    expect_error(.check_positive_integer(bad, "count"),
                 "'count' must be a positive integer", fixed = TRUE)
  }
})

test_that("real GNN and GNN-off valid predictions remain finite", {
  skip_if_no_libtorch()
  skip_if_not_installed("ape")
  # Select a real CPU device at the external device-selection boundary.
  local_mocked_bindings(get_device = function() torch::torch_device("cpu"),
                        .package = "pigauto")
  fit <- make_prediction_count_fit()
  plain <- make_prediction_count_fit(FALSE)
  original <- fit$X_scaled
  for (n in c(1L, 2L)) {
    pred <- predict(fit, n_imputations = n)
    expect_s3_class(pred, "pigauto_pred")
    expect_identical(dim(pred$imputed_latent), c(12L, 1L))
    expect_true(all(is.finite(pred$imputed_latent)))
  }
  two_steps <- fit
  two_steps$model_config$refine_steps <- 2L
  expect_true(all(is.finite(predict(two_steps)$imputed_latent)))
  expect_identical(fit$X_scaled, original)
  plain_pred <- predict(plain)
  expect_s3_class(plain_pred, "pigauto_pred")
  expect_true(all(is.finite(plain_pred$imputed_latent)))
  for (bad in list(0L, NULL)) {
    changed <- plain
    changed$model_config$refine_steps <- bad
    expect_equal(predict(changed)$imputed_latent, plain_pred$imputed_latent)
  }
})

test_that("prediction rejects invalid imputation counts in both arms", {
  skip_if_no_libtorch()
  skip_if_not_installed("ape")
  local_mocked_bindings(get_device = function() torch::torch_device("cpu"),
                        .package = "pigauto")
  for (fit in list(make_prediction_count_fit(), make_prediction_count_fit(FALSE))) {
    for (bad in list(0, -1, 1.5, NA_real_, Inf, c(1, 2), "2", integer())) {
      expect_error(predict(fit, n_imputations = bad),
                   "n_imputations.*positive integer")
    }
  }
})

test_that("prediction rejects invalid stored GNN refinement steps", {
  skip_if_no_libtorch()
  skip_if_not_installed("ape")
  local_mocked_bindings(get_device = function() torch::torch_device("cpu"),
                        .package = "pigauto")
  fit <- make_prediction_count_fit()
  # Zero reproduces the configuration in an actual saved legacy GNN fit.
  for (bad in list(0, NULL, -1, 1.5, NA_real_, Inf, c(1, 2), "2", integer())) {
    changed <- fit
    changed$model_config$refine_steps <- bad
    expect_error(predict(changed), "refine_steps.*positive integer")
  }
})

test_that("invalid imputation counts fail before device selection", {
  skip_if_no_libtorch()
  skip_if_not_installed("ape")
  local_mocked_bindings(get_device = function() torch::torch_device("cpu"),
                        .package = "pigauto")
  fit <- make_prediction_count_fit()
  with_mocked_bindings(
    expect_error(predict(fit, n_imputations = 0L), "n_imputations.*positive integer"),
    get_device = function() stop("DEVICE_SELECTED"), .package = "pigauto"
  )
})

test_that("invalid stored refinement fails before device selection", {
  skip_if_no_libtorch()
  skip_if_not_installed("ape")
  local_mocked_bindings(get_device = function() torch::torch_device("cpu"),
                        .package = "pigauto")
  fit <- make_prediction_count_fit()
  fit$model_config$refine_steps <- 0L
  with_mocked_bindings(
    expect_error(predict(fit), "refine_steps.*positive integer"),
    get_device = function() stop("DEVICE_SELECTED"), .package = "pigauto"
  )
})
