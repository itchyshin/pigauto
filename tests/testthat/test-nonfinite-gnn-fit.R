test_that("GNN fitting retains finite closed-gate and restored predictions", {
  skip_if_no_libtorch()
  testthat::local_mocked_bindings(cuda_is_available = function() FALSE,
    backends_mps_is_available = function() FALSE, .package = "torch")
  expect_identical(get_device()$type, "cpu")
  set.seed(41)
  tree <- ape::rtree(12)
  traits <- data.frame(x = seq_len(12) / 10, row.names = tree$tip.label)
  data <- preprocess_traits(traits, tree)
  graph <- build_phylo_graph(tree)
  run_fit <- function(lr, splits = NULL, epochs = 3L) fit_pigauto(data, tree, graph = graph,
    splits = splits, gnn = TRUE, epochs = epochs, eval_every = 1L, patience = 2L,
    lr = lr, warmup_epochs = 1L, hidden_dim = 8L, n_gnn_layers = 1L,
    use_transformer_blocks = FALSE, use_attention = TRUE,
    corruption_ramp = 0L, edge_dropout = 0, corruption_rate = 0.8,
    dropout = 0, verbose = FALSE, seed = 41, predict_method = "per_column",
    refine_steps = 1L)
  closed <- run_fit(1e10)
  expect_true(all(is.finite(predict(closed, return_se = FALSE)$imputed_latent)))
  good <- run_fit(0.003)
  expect_s3_class(good, "pigauto_fit")
  expect_true(all(is.finite(predict(good, return_se = FALSE)$imputed_latent)))
  untrained <- run_fit(0.003, epochs = 0L)
  expect_equal(nrow(untrained$history), 0L)
  expect_true(all(is.finite(predict(untrained, return_se = FALSE)$imputed_latent)))
  splits <- make_missing_splits(data$X_scaled, missing_frac = 0.25,
                                val_frac = 0.5, seed = 41)
  expect_warning(restored <- run_fit(1e4, splits), "Small validation set")
  expect_true(is.finite(restored$val_rmse))
  expect_true(all(is.finite(predict(restored, return_se = FALSE)$imputed_latent)))
})

test_that("a closed GNN gate keeps finite predictions across refinement", {
  skip_if_no_libtorch()
  testthat::local_mocked_bindings(cuda_is_available = function() FALSE,
    backends_mps_is_available = function() FALSE, .package = "torch")
  set.seed(41)
  tree <- ape::rtree(12)
  traits <- data.frame(x = seq_len(12) / 10, row.names = tree$tip.label)
  traits$x[2] <- NA_real_
  data <- preprocess_traits(traits, tree, covariates = data.frame(z = sin(seq_len(12))))
  graph <- build_phylo_graph(tree)
  fit <- fit_pigauto(data, tree, graph = graph, splits = NULL, gnn = TRUE,
    epochs = 3L, eval_every = 1L, patience = 2L, lr = 1e5, warmup_epochs = 1L,
    hidden_dim = 8L, n_gnn_layers = 1L, use_transformer_blocks = FALSE,
    use_attention = TRUE, corruption_ramp = 0L, edge_dropout = 0,
    corruption_rate = 0.8, dropout = 0, verbose = FALSE, seed = 41,
    predict_method = "per_column", refine_steps = 8L)
  expect_null(fit$calibrated_gates)
  expect_identical(as.numeric(torch::torch_sigmoid(fit$model_state$res_raw)$item()), 0)
  pred <- predict(fit, return_se = FALSE)$imputed_latent
  expect_true(all(is.finite(pred)))
  fixed <- as.matrix(fit$covariates) %*% t(as.matrix(as.array(fit$model_state[["cov_linear.weight"]])))
  fixed <- sweep(fixed, 2L, as.numeric(as.array(fit$model_state[["cov_linear.bias"]])), "+")
  expect_equal(unname(pred), unname(fit$baseline$mu + fixed), tolerance = 1e-5)
  positive <- fit
  positive$model_state$res_raw <- torch::torch_zeros_like(fit$model_state$res_raw)
  expect_error(predict(positive, return_se = FALSE), "non-finite.*learning rate")
  cfg <- positive$model_config
  model <- ResidualPhyloDAE(input_dim = cfg$input_dim, hidden_dim = cfg$hidden_dim,
    coord_dim = cfg$k_eigen, cov_dim = cfg$cov_dim, per_column_rs = cfg$per_column_rs,
    n_gnn_layers = cfg$n_gnn_layers, gate_cap = cfg$gate_cap,
    use_attention = cfg$use_attention, n_user_cov = cfg$n_user_cov,
    dropout = cfg$dropout, use_transformer_blocks = cfg$use_transformer_blocks,
    n_heads = cfg$n_heads, ffn_mult = cfg$ffn_mult)
  model$load_state_dict(positive$model_state)
  model$eval()
  X <- data$X_scaled
  observed <- !is.na(X)
  X[!observed] <- fit$baseline$mu[!observed]
  tensor <- function(x) torch::torch_tensor(x, dtype = torch::torch_float())
  expect_error(pigauto_refine_forward(model, tensor(X), tensor(graph$coords),
    tensor(graph$adj), NULL, tensor(graph$D_sq), tensor(fit$baseline$mu),
    tensor(fit$covariates), torch::torch_tensor(observed, dtype = torch::torch_bool()),
    8L, torch::torch_device("cpu")), "non-finite.*learning rate")
  splits <- make_missing_splits(data$X_scaled, missing_frac = 0.25,
                                val_frac = 0.5, seed = 41)
  expect_warning(validated <- fit_pigauto(data, tree, graph = graph, splits = splits,
    gnn = TRUE, epochs = 3L, eval_every = 1L, patience = 2L, lr = 1e5,
    warmup_epochs = 1L, hidden_dim = 8L, n_gnn_layers = 1L,
    use_transformer_blocks = FALSE, use_attention = TRUE, corruption_ramp = 0L,
    edge_dropout = 0, corruption_rate = 0.8, dropout = 0, verbose = FALSE,
    seed = 41, predict_method = "per_column", refine_steps = 8L),
    "Small validation set")
  expect_true(all(is.finite(predict(validated, return_se = FALSE)$imputed_latent)))
  calibrated <- positive
  calibrated$calibrated_gates <- 0
  calibrated$r_cal_gnn <- 0
  calibrated$r_cal_bm <- 1
  calibrated$r_cal_mean <- 0
  calibrated$mean_baseline_per_col <- 0
  expect_equal(unname(predict(calibrated, return_se = FALSE)$imputed_latent),
               unname(fit$baseline$mu + fixed), tolerance = 1e-5)
  calibrated$calibrated_gates <- calibrated$r_cal_gnn <- 0.4
  calibrated$r_cal_bm <- 0.6
  expect_error(predict(calibrated, return_se = FALSE), "non-finite.*learning rate")
})

test_that("calibration treats an invalid unused delta as zero contribution", {
  truth <- matrix(rep(c(0, 1), 10), ncol = 1)
  mu <- truth
  bad <- matrix(rep(c(Inf, NaN), 10), ncol = 1)
  continuous <- list(type = "continuous", latent_cols = 1L, name = "x")
  binary <- list(type = "binary", latent_cols = 1L, name = "x")
  expect_equal(compute_corner_loss(0, 1:20, continuous, mu, bad, truth, FALSE), 0)
  expect_identical(compute_corner_loss(0.5, 1:20, continuous, mu, bad, truth, FALSE), Inf)
  expect_identical(compute_corner_loss(0.5, 1:20, binary, mu, bad, truth, FALSE), Inf)
  expect_equal(compute_corner_loss(c(1, 0, 0), 1:20, continuous, mu, bad,
                                   truth, TRUE, 0), 0)
  expect_identical(compute_corner_loss(1e-12, 1:20, continuous, mu, bad,
                                      truth, FALSE), Inf)
  expect_identical(compute_corner_loss(0, 1:20, binary, mu, bad, truth, FALSE,
                                      fixed_cal = bad), Inf)
  cal <- calibrate_gates(list(continuous), mu, bad, truth, matrix(TRUE, 20, 1),
    gate_grid = c(0, 0.5), gate_cap = 0.8, gate_method = "single_split", seed = 41)
  expect_equal(unname(cal$r_cal_gnn), 0)
  score <- compute_conformal_scores(list(continuous), 0, mu, bad, truth,
                                     matrix(TRUE, 20, 1))
  expect_equal(unname(score), 0)
})
