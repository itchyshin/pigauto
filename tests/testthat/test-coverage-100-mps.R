test_that("real MPS prediction retains scalar fixed-effect bias", {
  skip_if_not(isTRUE(torch::torch_is_installed()), "real MPS runtime unavailable")
  skip_if_not("backends_mps_is_available" %in% getNamespaceExports("torch"),
              "real MPS device unavailable")
  skip_if_not(isTRUE(torch::backends_mps_is_available()), "real MPS device unavailable")
  skip_if_not(identical(get_device()$type, "mps"), "real MPS device unavailable")
  old_threads <- torch::torch_get_num_threads()
  on.exit(torch::torch_set_num_threads(old_threads), add = TRUE)
  torch::torch_set_num_threads(1L)
  expect_identical(get_device()$type, "mps")

  withr::local_seed(319)
  tree <- ape::rtree(12)
  i <- seq_len(12)
  traits <- data.frame(x = sin(i / 3) + i / 10, row.names = tree$tip.label)
  traits$x[c(2, 5)] <- NA_real_
  covs <- data.frame(c1 = cos(i * .4), c2 = sin(i * .7))
  data <- preprocess_traits(traits, tree, trait_types = c(x = "continuous"),
                            log_transform = FALSE, covariates = covs)
  fit <- fit_pigauto(data, tree, splits = NULL, gnn = TRUE, epochs = 1L,
    eval_every = 1L, hidden_dim = 8L, n_gnn_layers = 1L, dropout = 0,
    edge_dropout = 0, use_transformer_blocks = FALSE, refine_steps = 1L,
    predict_method = "per_column", verbose = FALSE, seed = 41)
  expect_equal(dim(fit$covariates), c(12L, 2L))
  expect_equal(dim(as.array(fit$model_state[["cov_linear.weight"]])), c(1L, 2L))
  expect_equal(length(as.array(fit$model_state[["cov_linear.bias"]])), 1L)
  # Known model parameters isolate the fixed-effect term while keeping the
  # real fitted model, actual MPS prediction path and preprocessed covariates.
  fit$model_state[["cov_linear.weight"]] <- torch::torch_tensor(
    matrix(c(.4, -.25), 1, 2), dtype = torch::torch_float())
  fit$model_state[["cov_linear.bias"]] <- torch::torch_tensor(
    .7, dtype = torch::torch_float())$reshape(c(1L))
  fit$calibrated_gates <- 0
  fit$r_cal_bm <- 1; fit$r_cal_gnn <- 0; fit$r_cal_mean <- 0
  fit$mean_baseline_per_col <- 0
  first <- predict(fit, return_se = FALSE)
  expected <- as.numeric(fit$baseline$mu[, 1]) +
    as.numeric(as.matrix(fit$covariates) %*% c(.4, -.25)) + .7
  missing <- which(is.na(traits$x))
  expect_identical(dim(first$imputed_latent), c(12L, 1L))
  expect_true(all(is.finite(first$imputed_latent)))
  expect_equal(unname(first$imputed_latent[, 1]), expected, tolerance = 2e-5)
  expect_equal(unname(first$imputed_latent[missing, 1]), expected[missing], tolerance = 2e-5)
  completed <- build_completed(traits, first$imputed, species_col = NULL,
                               input_row_order = NULL)$completed
  expect_identical(completed$x[!is.na(traits$x)], traits$x[!is.na(traits$x)])
  expect_true(all(is.finite(completed$x)))

  fit$model_state[["cov_linear.bias"]] <- torch::torch_tensor(
    1.1, dtype = torch::torch_float())$reshape(c(1L))
  second <- predict(fit, return_se = FALSE)
  expect_true(all(is.finite(second$imputed_latent)))
  expect_equal(unname(second$imputed_latent[missing, 1] - first$imputed_latent[missing, 1]),
               rep(.4, length(missing)), tolerance = 2e-5)
})
