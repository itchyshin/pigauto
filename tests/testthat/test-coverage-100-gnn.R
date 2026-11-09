test_that("mixed losses agree with hand calculations for every trait", {
  skip_if_no_libtorch()
  types <- c("continuous", "count", "ordinal", "proportion", "binary",
             "categorical", "zi_count", "multi_proportion")
  cols <- list(1L, 2L, 3L, 4L, 5L, 6:8, 9:10, 11:12)
  tm <- Map(function(type, lc) list(type = type, latent_cols = lc,
                                    n_latent = length(lc)), types, cols)
  truth <- matrix(0, 3, 12)
  truth[, 1:4] <- matrix(rep(1:4, each = 3), 3)
  truth[, 5] <- c(0, 1, 0)
  truth[, 6:8] <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  truth[, 9] <- c(1, 0, 1)
  truth[, 10] <- c(2, NA, 2)
  truth[, 11:12] <- rbind(c(1, -1), c(-1, 1), c(1, -1))
  pred <- torch::torch_zeros(c(3L, 12L), requires_grad = TRUE)
  y <- torch::torch_tensor(truth, dtype = torch::torch_float())
  mask <- matrix(FALSE, 3, 12); mask[1:2, ] <- TRUE
  m <- torch::torch_tensor(mask, dtype = torch::torch_bool())
  expected <- c(1, 4, 9, 16, log(2), log(3), (log(2) + 4) / 2, 1)
  for (i in seq_along(tm)) {
    expect_equal(as.numeric(compute_mixed_loss(pred, y, m, tm[i])$item()),
                 expected[i], tolerance = 1e-6, info = types[i])
    expect_equal(composite_val_loss(pred, y, m, tm[i]), expected[i],
                 tolerance = 1e-6, info = types[i])
  }
  loss <- compute_mixed_loss(pred, y, m, tm)
  expect_equal(as.numeric(loss$item()), mean(expected), tolerance = 1e-6)
  expect_equal(composite_val_loss(pred, y, m, tm), mean(expected), tolerance = 1e-6)
  loss$backward()
  expect_true(all(is.finite(as.array(pred$grad))))
  empty <- torch::torch_zeros(c(3L, 12L), dtype = torch::torch_bool())
  expect_equal(as.numeric(compute_mixed_loss(pred, y, empty, tm)$item()), 0)
  expect_identical(composite_val_loss(pred, y, empty, tm), Inf)
})

test_that("learning rate warmup and cosine endpoints are explicit", {
  expect_equal(cosine_lr(1, 2, 6, 0.1, 0.01), 0.05)
  expect_equal(cosine_lr(2, 2, 6, 0.1, 0.01), 0.1)
  expect_equal(cosine_lr(4, 2, 6, 0.1, 0.01), 0.055)
  expect_equal(cosine_lr(6, 2, 6, 0.1, 0.01), 0.01)
  expect_equal(cosine_lr(3, 2, 2, 0.1), 0.1)
})

test_that("calibration corner losses and empty validation preserve the baseline", {
  mu <- matrix(c(0, 2, 4), ncol = 1)
  delta <- matrix(c(2, 4, 6), ncol = 1)
  truth <- matrix(c(1, 3, 5), ncol = 1)
  tm <- list(type = "continuous", latent_cols = 1L, name = "x")
  expect_equal(compute_corner_loss(0.5, 1:3, tm, mu, delta, truth, FALSE), 0)
  expect_equal(compute_corner_loss(c(0, 0, 1), 1:3, tm, mu, delta, truth,
                                   TRUE, 3), 8 / 3)
  expect_identical(compute_corner_loss(0, integer(), tm, mu, delta, truth, FALSE), Inf)
  expect_error(compute_corner_loss(0, 1:3, tm, mu, delta, truth, FALSE,
                                   fixed_cal = matrix(0, 2, 1)), "same dimensions")
  cal <- calibrate_gates(list(tm), mu, delta, truth, matrix(FALSE, 3, 1),
                         gate_grid = c(0, 0.5), gate_cap = 0.8)
  expect_equal(unname(cal$r_cal_bm), 1)
  expect_equal(unname(cal$r_cal_gnn), 0)
  expect_equal(unname(cal$r_cal_mean), 0)
})

test_that("calibration split methods keep the better baseline corner", {
  tm <- list(type = "continuous", latent_cols = 1L, name = "x")
  truth <- matrix(0, 20, 1)
  mu <- truth
  delta <- matrix(2, 20, 1)
  mask <- matrix(TRUE, 20, 1)
  for (method in c("single_split", "median_splits", "cv_folds")) {
    cal <- calibrate_gates(list(tm), mu, delta, truth, mask,
      gate_grid = c(0, 0.5, 1), gate_cap = 1, gate_method = method,
      gate_splits_B = 3L, gate_cv_folds = 3L, seed = 42)
    expect_equal(unname(cal$r_cal_bm), 1, info = method)
    expect_equal(unname(cal$r_cal_gnn), 0, info = method)
  }
  expect_error(calibrate_gates(list(tm), mu, delta, truth, mask,
    gate_grid = c(0, 1), gate_cap = 1, gate_method = "cv_folds",
    gate_cv_folds = 1), "integer >= 2")
  expect_error(simplex_grid(0.3), "evenly divide")
  expect_equal(nrow(simplex_grid(0.5)), 6L)
  expect_equal(mean_baseline_scalar(c(0, 1), c(TRUE, TRUE), "binary"), 0)
  expect_equal(mean_baseline_scalar(c(NA, 1), c(TRUE, FALSE), "continuous"), 0)
  expect_error(mean_baseline_scalar(1, TRUE, "unknown"), "unsupported")
})

test_that("calibration discrete corners use classification and count scale", {
  mu <- matrix(c(-1, 1, 3, 0, 0, 3, 0, 0), 2, 4)
  truth <- matrix(c(0, 1, 1, 0, 0, 1, 0, 0), 2, 4)
  binary <- list(type = "binary", latent_cols = 1L)
  categorical <- list(type = "categorical", latent_cols = 2:4)
  expect_equal(compute_corner_loss(0, 1:2, binary, mu, mu, truth, FALSE), 0)
  expect_equal(compute_corner_loss(0, 1:2, categorical, mu, mu, truth, FALSE), 0)
  z <- list(type = "zi_count", latent_cols = 1:2, mean = 0, sd = 1)
  zm <- cbind(c(0, 0), c(log(3), log(3)))
  zy <- cbind(c(0, 1), c(NA, log(3)))
  expect_equal(compute_corner_loss(c(0, 0, 1), 1:2, z, zm, zm, zy,
                                   TRUE, c(0, log(3))), 1)
})

test_that("tiny GNN fit returns finite predictions and bounded gates", {
  skip_if_no_libtorch()
  testthat::local_mocked_bindings(cuda_is_available = function() FALSE,
    backends_mps_is_available = function() FALSE, .package = "torch")
  expect_identical(get_device()$type, "cpu")
  set.seed(17)
  tree <- ape::rtree(20)
  traits <- data.frame(x = seq_len(20) / 10, row.names = tree$tip.label)
  traits$x[c(3, 7)] <- NA_real_
  data <- preprocess_traits(traits, tree)
  splits <- make_missing_splits(data$X_scaled, missing_frac = 0.25,
                                seed = 3, trait_map = data$trait_map)
  expect_warning(fit <- fit_pigauto(data, tree, splits = splits, gnn = TRUE,
    hidden_dim = 8L, n_gnn_layers = 1L, epochs = 2L, eval_every = 1L,
    patience = 2L, dropout = 0, refine_steps = 1L, verbose = FALSE, seed = 4),
    "Small validation set")
  expect_s3_class(fit, "pigauto_fit")
  expect_true(isTRUE(fit$model_config$gnn))
  expect_gt(length(fit$model_state), 0L)
  expect_gt(nrow(fit$history), 0L)
  expect_true(all(is.finite(fit$history$loss_rec)))
  expect_true(all(fit$r_cal_gnn >= 0 & fit$r_cal_gnn <= 1))
  expect_equal(unname(fit$r_cal_bm + fit$r_cal_gnn + fit$r_cal_mean), 1,
               tolerance = 1e-6)
  pred <- predict(fit, return_se = FALSE)
  expect_true(all(is.finite(pred$imputed_latent)))
  expect_equal(dim(pred$imputed_latent), c(20L, 1L))
})

test_that("attachment explains how to install a missing torch dependency", {
  attach <- .onAttach
  env <- new.env(parent = environment(attach))
  env$requireNamespace <- function(...) FALSE
  environment(attach) <- env
  expect_message(attach("unused", "pigauto"), "install.packages\\('torch'\\)")
})

test_that("GPU checkpoint reports formatted results only when enabled", {
  testthat::local_mocked_bindings(gpu_mem_str = function() "1.00/2.00 GiB",
                                 .package = "pigauto")
  withr::local_envvar(PIGAUTO_DEBUG_GPU_MEM = "1")
  expect_output(gpu_mem_checkpoint("fit"), "\\[GPU @ fit\\] 1.00/2.00 GiB")
  withr::local_envvar(PIGAUTO_DEBUG_GPU_MEM = "0")
  expect_output(gpu_mem_checkpoint("fit"), NA)
})

test_that("device selection prioritizes available backends (mocked routing)", {
  testthat::local_mocked_bindings(
    cuda_is_available = function() TRUE,
    backends_mps_is_available = function() TRUE,
    torch_device = function(type) type, .package = "torch")
  expect_identical(get_device(), "cuda")
})

test_that("device selection falls back to MPS or CPU (mocked routing)", {
  testthat::local_mocked_bindings(
    cuda_is_available = function() FALSE,
    backends_mps_is_available = function() TRUE,
    torch_device = function(type) type, .package = "torch")
  expect_identical(get_device(), "mps")
  testthat::local_mocked_bindings(backends_mps_is_available = function() FALSE,
                                 .package = "torch")
  expect_identical(get_device(), "cpu")
  expect_identical(gpu_mem_str(), "")
})

test_that("GPU memory formatter reads the structured statistics (mocked API)", {
  testthat::local_mocked_bindings(
    cuda_is_available = function() TRUE,
    cuda_memory_stats = function() list(allocated_bytes.all.current = 1024^3,
                                        reserved_bytes.all.current = 2 * 1024^3),
    .package = "torch")
  expect_identical(gpu_mem_str(), "1.00/2.00 GiB (alloc/reserved, via cuda_memory_stats)")
})

test_that("GPU memory formatter handles alternate keys and subprocess fallback", {
  testthat::local_mocked_bindings(
    cuda_is_available = function() TRUE,
    cuda_memory_stats = function() list(allocated_bytes.foo.current = 2 * 1024^3,
                                        reserved_bytes.foo.current = 3 * 1024^3),
    .package = "torch")
  expect_match(gpu_mem_str(), "2.00/3.00 GiB \\(alloc/reserved, matched keys:")
  testthat::local_mocked_bindings(cuda_memory_stats = function() stop("unavailable"),
                                 .package = "torch")
  formatter <- gpu_mem_str
  env <- new.env(parent = environment(formatter))
  env$system <- function(...) "1024, 4096"
  environment(formatter) <- env
  expect_identical(formatter(), "1.00/4.00 GiB (used/total, via nvidia-smi)")
  env$system <- function(...) "bad"
  expect_identical(formatter(), "(memory query unavailable; CUDA alive)")
})

test_that("masked tensor RMSE uses only selected cells", {
  skip_if_no_libtorch()
  pred <- torch::torch_tensor(c(1, 3, 100))
  truth <- torch::torch_tensor(c(0, 1, 0))
  mask <- torch::torch_tensor(c(TRUE, TRUE, FALSE))
  expect_equal(as.numeric(rmse_torch(pred, truth, mask)$item()), sqrt(2.5),
               tolerance = 1e-6)
})

test_that("small CV calibration falls back safely and reports its limit", {
  tm <- list(type = "continuous", latent_cols = 1L, name = "x")
  mu <- matrix(0, 1, 1); delta <- matrix(2, 1, 1)
  truth <- matrix(0, 1, 1); mask <- matrix(TRUE, 1, 1)
  expect_warning(cal <- calibrate_gates(list(tm), mu, delta, truth, mask,
    gate_grid = c(0, 1), gate_cap = 1, gate_method = "cv_folds",
    gate_cv_folds = 5L, seed = 4), "Small validation set")
  expect_equal(unname(cal$r_cal_bm), 1)
  expect_equal(unname(cal$r_cal_gnn), 0)
})

test_that("three-way calibration retains a perfect mean corner", {
  tm <- list(type = "continuous", latent_cols = 1L, name = "x")
  mu <- matrix(3, 20, 1); delta <- matrix(4, 20, 1)
  truth <- matrix(1, 20, 1); mask <- matrix(TRUE, 20, 1)
  cal <- calibrate_gates(list(tm), mu, delta, truth, mask,
    gate_grid = c(0, 1), gate_cap = 1, safety_floor = TRUE,
    mean_baseline_per_col = 1, simplex_step = 0.5, seed = 11)
  expect_equal(unname(cal$r_cal_bm), 0)
  expect_equal(unname(cal$r_cal_gnn), 0)
  expect_equal(unname(cal$r_cal_mean), 1)
})

test_that("conformal scoring checks dimensions before blending", {
  tm <- list(list(name = "x", type = "continuous", latent_cols = 1L))
  x <- matrix(0, 3, 1); mask <- matrix(TRUE, 3, 1)
  score <- function(...) compute_conformal_scores(tm, 0, x, x, x, mask, ...)
  expect_error(compute_conformal_scores(tm, 0, x, matrix(0, 2, 1), x, mask),
               "delta_cal.*same dimensions")
  expect_error(compute_conformal_scores(tm, 0, x, x, matrix(0, 2, 1), mask),
               "X_truth_r.*same dimensions")
  expect_error(compute_conformal_scores(tm, 0, x, x, x, matrix(TRUE, 2, 1)),
               "val_mask_mat.*same dimensions")
  expect_error(score(fixed_cal = matrix(0, 2, 1)), "fixed_cal.*same dimensions")
  expect_error(compute_conformal_scores(tm, c(0, 1), x, x, x, mask),
               "weights and means must have length")
  expect_error(score(method = "mondrian"), "D_sq.*obs_mask_mat")
  missing <- matrix(NA_real_, 3, 1)
  expect_true(is.na(compute_conformal_scores(tm, 0, x, x, missing, mask)[[1]]))
})

test_that("Mondrian fallback returns the global score without donors", {
  tm <- list(list(name = "x", type = "continuous", latent_cols = 1L))
  mu <- matrix(0, 3, 1); truth <- matrix(c(1, 2, 3), 3, 1)
  expect_message(scores <- compute_conformal_scores(tm, 0, mu, mu, truth,
    matrix(TRUE, 3, 1), method = "mondrian", D_sq = matrix(0, 3, 3),
    obs_mask_mat = matrix(FALSE, 3, 1), verbose = TRUE), "fallback")
  info <- attr(scores, "mondrian")$x
  expect_equal(as.numeric(scores), 3)
  expect_true(info$fallback)
  expect_equal(info$near_score, 3)
  expect_equal(info$far_score, 3)
  expect_true(is.na(info$threshold))
  expect_length(mondrian_locality(matrix(0, 2, 2), integer(), integer()), 0L)
})

test_that("phylogenetic gate overrides choose only supported corners", {
  bm <- c(0.7, 0.5); gnn <- c(0.3, 0.3); mean <- c(0, 0.2)
  unchanged <- apply_phylo_signal_gate_override(bm, gnn, mean, integer(), TRUE, c(1, 1))
  expect_equal(unchanged$r_cal_bm, bm)
  expect_false(unchanged$used_mean_corner)
  mean_out <- apply_phylo_signal_gate_override(bm, gnn, mean, 2L, TRUE, c(1, 1))
  expect_equal(c(mean_out$r_cal_bm[2], mean_out$r_cal_gnn[2], mean_out$r_cal_mean[2]), c(0, 0, 1))
  expect_true(mean_out$used_mean_corner)
  expect_warning(bm_out <- apply_phylo_signal_gate_override(bm, gnn, mean, 2L,
    FALSE, NULL), "falling back to BM")
  expect_equal(c(bm_out$r_cal_bm[2], bm_out$r_cal_gnn[2], bm_out$r_cal_mean[2]), c(1, 0, 0))
})

test_that("validation mask splitting handles empty and sparse columns", {
  mask <- cbind(rep(FALSE, 6), c(TRUE, TRUE, rep(FALSE, 4)), rep(TRUE, 6))
  out <- split_val_cal_conf(mask, TRUE, min_val_cells = 2L, seed = 9L)
  expect_false(any(out$val_mask_cal[, 1]))
  expect_identical(out$val_mask_cal[, 2], mask[, 2])
  expect_identical(out$val_mask_conf[, 2], mask[, 2])
  expect_equal(sum(out$val_mask_cal[, 3]), 3)
  expect_equal(sum(out$val_mask_conf[, 3]), 3)
  expect_false(any(out$val_mask_cal[, 3] & out$val_mask_conf[, 3]))
  expect_identical(out$val_mask_cal | out$val_mask_conf, mask)
})

test_that("fit rejects unsupported combinations before training", {
  expect_error(fit_pigauto(list(), NULL, joint_refine_iter = -1), "non-negative integer")
  expect_error(fit_pigauto(list(), NULL), "pigauto_data")
  set.seed(5)
  tree <- ape::rtree(5)
  df <- data.frame(x = 1:5, row.names = tree$tip.label)
  data <- preprocess_traits(df, tree)
  graph <- list(coords = matrix(0, 5, 2))
  baseline <- list(mu = matrix(0, 5, 1))
  expect_error(fit_pigauto(data, tree, graph = graph, baseline = baseline,
                           gnn = NA, verbose = FALSE), "gnn.*TRUE or FALSE")
  expect_error(fit_pigauto(data, tree, graph = graph, baseline = baseline,
                           gnn = FALSE, conformal_method = "mondrian",
                           verbose = FALSE), "not supported when gnn = FALSE")
  data$multi_obs <- TRUE
  expect_error(fit_pigauto(data, tree, graph = graph, baseline = baseline,
                           gnn = TRUE, conformal_method = "mondrian",
                           verbose = FALSE), "not supported for multi-observation")
})

test_that("safety-floor means respect type-specific training masks", {
  x <- cbind(gate = c(0, 1, 1), magnitude = c(NA, 2, 4),
             comp1 = c(-1, 1, 2), comp2 = c(1, -1, -2))
  data <- list(X_scaled = x, trait_map = list(
    list(type = "zi_count", latent_cols = 1:2),
    list(type = "multi_proportion", latent_cols = 3:4)))
  mask <- matrix(FALSE, 3, 4); mask[3, ] <- TRUE
  means <- compute_safety_floor_means(data, mask, splits = NULL, safety_floor = TRUE)
  expect_equal(unname(means), c(0, 2, 0, 0))
  expect_null(compute_safety_floor_means(data, mask, splits = NULL, safety_floor = FALSE))
})

test_that("zero-inflated loss handles separate gate and finite magnitude masks", {
  skip_if_no_libtorch()
  tm <- list(list(type = "zi_count", latent_cols = 1:2))
  pred <- torch::torch_zeros(c(2L, 2L))
  truth <- torch::torch_tensor(rbind(c(1, 2), c(0, NA_real_)))
  gate_only <- torch::torch_tensor(rbind(c(TRUE, FALSE), c(TRUE, FALSE)))
  mag_only <- torch::torch_tensor(rbind(c(FALSE, TRUE), c(FALSE, TRUE)))
  for (metric in list(function(m) as.numeric(compute_mixed_loss(pred, truth, m, tm)$item()),
                      function(m) composite_val_loss(pred, truth, m, tm))) {
    expect_equal(metric(gate_only), log(2) / 2, tolerance = 1e-6)
    expect_equal(metric(mag_only), 2, tolerance = 1e-6)
  }
})

test_that("unsupported compositional mean corner falls back safely", {
  tm <- list(type = "multi_proportion", latent_cols = 1:2, name = "colour")
  mu <- matrix(0, 20, 2); delta <- matrix(1, 20, 2)
  truth <- matrix(0, 20, 2); mask <- matrix(TRUE, 20, 2)
  expect_identical(compute_corner_loss(c(0, 0, 1), 1:20, tm, mu, delta,
    truth, TRUE, c(0, 0)), Inf)
  expect_error(calibrate_gates(list(tm), mu, delta, truth, mask,
    gate_grid = c(0, 1), gate_cap = 1, safety_floor = TRUE),
    "requires mean_baseline_per_col")
  cal <- calibrate_gates(list(tm), mu, delta, truth, mask,
    gate_grid = c(0, 1), gate_cap = 1, safety_floor = TRUE,
    mean_baseline_per_col = c(0, 0), simplex_step = 0.5, seed = 4)
  expect_equal(unname(cal$r_cal_bm), c(1, 1))
  expect_equal(unname(cal$r_cal_gnn), c(0, 0))
  expect_equal(unname(cal$r_cal_mean), c(0, 0))
})

test_that("unsplit CPU training retains its legacy graph and uncalibrated state", {
  skip_if_no_libtorch()
  testthat::local_mocked_bindings(cuda_is_available = function() FALSE,
    backends_mps_is_available = function() FALSE, .package = "torch")
  set.seed(41)
  tree <- ape::rtree(12)
  traits <- data.frame(x = seq_len(12) / 10, row.names = tree$tip.label)
  data <- preprocess_traits(traits, tree)
  graph <- build_phylo_graph(tree)
  graph$D_sq <- NULL
  messages <- capture.output(fit <- fit_pigauto(data, tree, graph = graph,
    gnn = TRUE, epochs = 2L, eval_every = 1L, patience = 2L,
    hidden_dim = 8L, n_gnn_layers = 1L, use_transformer_blocks = FALSE,
    use_attention = TRUE, corruption_ramp = 0L, edge_dropout = 0,
    corruption_rate = 0.8, dropout = 0, verbose = TRUE, seed = 41),
    type = "message")
  expect_true(isTRUE(fit$model_config$gnn))
  expect_gt(length(fit$model_state), 0L)
  expect_gt(nrow(fit$history), 0L)
  expect_true(all(is.finite(fit$history$loss_rec)))
  expect_true(all(is.na(fit$history$val_loss)))
  expect_null(fit$r_cal_gnn)
  expect_null(fit$conformal_scores)
  expect_null(fit$graph$D_sq)
  expect_true(any(grepl("Using device: cpu", messages, fixed = TRUE)))
  expect_true(any(grepl("epoch", messages, fixed = TRUE)))
})

test_that("no-split baseline fit uses the baseline corner", {
  set.seed(32)
  tree <- ape::rtree(8)
  data <- preprocess_traits(data.frame(x = seq_len(8) / 10,
                                      row.names = tree$tip.label), tree)
  fit <- fit_pigauto(data, tree, gnn = FALSE, verbose = FALSE)
  expect_false(fit$model_config$gnn)
  expect_equal(unname(fit$r_cal_bm), 1)
  expect_equal(unname(fit$r_cal_gnn), 0)
  expect_equal(unname(fit$r_cal_mean), 0)
  expect_null(fit$conformal_scores)
})

test_that("zero-inflated trait trains both latent columns on CPU", {
  skip_if_no_libtorch()
  testthat::local_mocked_bindings(cuda_is_available = function() FALSE,
    backends_mps_is_available = function() FALSE, .package = "torch")
  expect_identical(get_device()$type, "cpu")
  set.seed(52)
  tree <- ape::rtree(16)
  traits <- data.frame(parasites = rep(c(0L, 1L, 2L, 4L), 4),
                       row.names = tree$tip.label)
  data <- preprocess_traits(traits, tree, trait_types = c(parasites = "zi_count"))
  fit <- fit_pigauto(data, tree, gnn = TRUE, hidden_dim = 8L,
    n_gnn_layers = 1L, epochs = 2L, eval_every = 1L, patience = 2L,
    corruption_rate = 0.8, corruption_ramp = 0L, dropout = 0,
    use_transformer_blocks = FALSE, verbose = FALSE, seed = 52)
  expect_identical(data$trait_map[[1]]$type, "zi_count")
  expect_length(data$trait_map[[1]]$latent_cols, 2L)
  expect_true(isTRUE(fit$model_config$gnn))
  expect_gt(nrow(fit$history), 0L)
  expect_true(all(is.finite(fit$history$loss_rec)))
  expect_equal(fit$model_config$input_dim, 2L)
  expect_true(all(is.finite(predict(fit, return_se = FALSE)$imputed_latent)))
})

test_that("discrete split verification keeps an accurate baseline", {
  tm <- list(type = "binary", latent_cols = 1L, name = "binary")
  truth <- matrix(rep(c(0, 1), 10), ncol = 1)
  mu <- matrix(ifelse(truth == 1, 2, -2), ncol = 1)
  delta <- -mu
  mask <- matrix(TRUE, 20, 1)
  for (method in c("single_split", "median_splits", "cv_folds")) {
    cal <- calibrate_gates(list(tm), mu, delta, truth, mask,
      gate_grid = c(0, 0.5, 1), gate_cap = 1, gate_method = method,
      gate_splits_B = 3L, gate_cv_folds = 3L, seed = 23)
    expect_equal(unname(cal$r_cal_bm), 1, info = method)
    expect_equal(unname(cal$r_cal_gnn), 0, info = method)
  }
})

test_that("unsupported corner types return infinite loss and unnamed gates are labelled", {
  x <- matrix(rep(c(0, 1), 10), ncol = 1)
  unsupported <- list(type = "unsupported", latent_cols = 1L, name = "x")
  expect_identical(compute_corner_loss(0, 1:20, unsupported, x, x, x, FALSE), Inf)
  binary <- list(type = "binary", latent_cols = 1L, name = "x")
  messages <- capture.output(cal <- calibrate_gates(list(binary), x, x, x,
    matrix(TRUE, 20, 1), gate_grid = c(0, 1), gate_cap = 1,
    gate_method = "single_split", seed = NULL, verbose = TRUE), type = "message")
  expect_true(any(grepl("col1=", messages, fixed = TRUE)))
  expect_equal(unname(cal$r_cal_bm + cal$r_cal_gnn), 1)
})

test_that("CPU GNN prediction pools reproducible stochastic draws", {
  skip_if_no_libtorch()
  testthat::local_mocked_bindings(cuda_is_available = function() FALSE,
    backends_mps_is_available = function() FALSE, .package = "torch")
  expect_identical(get_device()$type, "cpu")
  set.seed(41)
  tree <- ape::rtree(12)
  traits <- data.frame(x = seq_len(12) / 10, row.names = tree$tip.label)
  traits$x[c(2, 5)] <- NA_real_
  d <- preprocess_traits(traits, tree, trait_types = c(x = "continuous"),
                         log_transform = FALSE,
                         covariates = data.frame(z = cos(seq_len(12))))
  fit <- fit_pigauto(d, tree, gnn = TRUE, epochs = 2L, eval_every = 1L,
    hidden_dim = 8L, n_gnn_layers = 1L, corruption_rate = .8, corruption_ramp = 0L,
    dropout = .2, edge_dropout = 0, use_transformer_blocks = FALSE,
    use_attention = TRUE, verbose = FALSE, seed = 41, predict_method = "per_column",
    refine_steps = 2L)
  first <- predict(fit, n_imputations = 3L, return_se = TRUE)
  torch::torch_manual_seed(1234)
  second <- predict(fit, n_imputations = 3L, return_se = TRUE)
  expect_identical(first$n_imputations, 3L)
  expect_length(first$imputed_datasets, 3L)
  expect_equal(dim(first$imputed_latent), c(12L, 1L))
  expect_true(all(is.finite(first$imputed_latent)))
  expect_equal(first$imputed_latent, second$imputed_latent)
  expect_equal(first$imputed_datasets, second$imputed_datasets)
  draws <- vapply(first$imputed_datasets, function(x) x$x, numeric(12L))
  expect_equal(first$imputed$x, rowMeans(draws), tolerance = 1e-8)
  expect_true(all(is.finite(draws)))
  expect_gt(stats::sd(draws[2, ]), 0)
  expect_gt(first$se_latent[2, 1], 0)
  malformed <- fit
  malformed$covariates <- fit$covariates[-1, , drop = FALSE]
  expect_error(predict(malformed, return_se = FALSE),
               "Stored covariates have 11 row")
})

test_that("legacy trait-map-free training uses finite continuous losses", {
  skip_if_no_libtorch()
  testthat::local_mocked_bindings(cuda_is_available = function() FALSE,
    backends_mps_is_available = function() FALSE, .package = "torch")
  set.seed(42)
  tree <- ape::rtree(20)
  traits <- data.frame(x = seq_len(20) / 10, row.names = tree$tip.label)
  data <- preprocess_traits(traits, tree, trait_types = c(x = "continuous"),
                            log_transform = FALSE)
  splits <- make_missing_splits(data$X_scaled, missing_frac = .25,
                                val_frac = .5, seed = 41)
  graph <- build_phylo_graph(tree)
  baseline <- fit_baseline(data, tree, splits = splits, predict_method = "per_column")
  # Mimic the older continuous-only pigauto_data contract, retaining a real
  # baseline computed from the same observed data and held-out mask.
  data$trait_map <- NULL
  messages <- capture.output(fit <- fit_pigauto(data, tree, graph = graph,
    baseline = baseline, splits = splits, gnn = TRUE, epochs = 3L, eval_every = 1L,
    patience = 1L, warmup_epochs = 10000000L,
    hidden_dim = 8L, n_gnn_layers = 1L, corruption_rate = .8, corruption_ramp = 0L,
    dropout = 0, edge_dropout = 0, use_transformer_blocks = FALSE,
    use_attention = TRUE, verbose = TRUE, seed = 41, predict_method = "per_column",
    refine_steps = 1L), type = "message")
  expect_null(fit$trait_map)
  expect_identical(nrow(fit$history), 2L)
  expect_true(all(is.finite(fit$history$loss_rec)))
  expect_true(all(is.finite(fit$history$val_loss)))
  expect_true(is.finite(fit$val_rmse) && fit$val_rmse >= 0)
  expect_true(is.finite(fit$test_rmse) && fit$test_rmse >= 0)
  expect_true(any(grepl(" | val ", messages, fixed = TRUE)))
  expect_true(any(grepl("Early stopping at epoch 2", messages, fixed = TRUE)))
  pred <- predict(fit, return_se = FALSE)
  expect_equal(dim(pred$imputed), c(20L, 1L))
  expect_true(all(is.finite(pred$imputed)))
})

test_that("phylogenetic signal gate selects the available safe corner", {
  skip_if_no_libtorch()
  skip_if_not_installed("phytools")
  testthat::local_mocked_bindings(cuda_is_available = function() FALSE,
    backends_mps_is_available = function() FALSE, .package = "torch")
  set.seed(71)
  tree <- ape::rtree(30)
  traits <- data.frame(x = sin(seq_len(30) * 1.7), row.names = tree$tip.label)
  data <- preprocess_traits(traits, tree, trait_types = c(x = "continuous"),
                            log_transform = FALSE)
  splits <- make_missing_splits(data$X_scaled, missing_frac = .2,
                                val_frac = .5, seed = 41)
  graph <- build_phylo_graph(tree)
  for (floor in c(TRUE, FALSE)) {
    warnings <- character()
    fit <- withCallingHandlers(fit_pigauto(data, tree, graph = graph, splits = splits,
      gnn = TRUE, epochs = 2L, eval_every = 1L, hidden_dim = 8L, n_gnn_layers = 1L,
      corruption_rate = .8, corruption_ramp = 0L, dropout = 0, edge_dropout = 0,
      use_transformer_blocks = FALSE, use_attention = TRUE, verbose = FALSE,
      seed = 41, predict_method = "per_column", refine_steps = 1L,
      safety_floor = floor, phylo_signal_gate = TRUE, phylo_signal_threshold = 1),
      warning = function(w) {
        msg <- conditionMessage(w)
        if (grepl("Small validation set", msg, fixed = TRUE) ||
            grepl("phylo_signal_gate triggered", msg, fixed = TRUE)) {
          warnings <<- c(warnings, msg)
          invokeRestart("muffleWarning")
        }
      })
    expect_true(is.finite(fit$phylo_signal_per_trait[[1]]))
    expect_lt(fit$phylo_signal_per_trait[[1]], 1)
    expect_true(fit$phylo_gate_triggered[[1]])
    expect_true(any(grepl("Small validation set", warnings, fixed = TRUE)))
    expect_identical(any(grepl("phylo_signal_gate triggered", warnings, fixed = TRUE)), !floor)
    expect_equal(unname(c(fit$r_cal_bm, fit$r_cal_gnn, fit$r_cal_mean)),
                 if (floor) c(0, 0, 1) else c(1, 0, 0))
    latent <- predict(fit, return_se = FALSE)$imputed_latent
    expected <- if (floor) matrix(fit$mean_baseline_per_col[[1]], 30, 1) else fit$baseline$mu
    expect_equal(unname(latent), unname(expected), tolerance = 1e-6)
  }
})

test_that("GPU debug checkpoint stays silent on CPU", {
  testthat::local_mocked_bindings(cuda_is_available = function() FALSE,
                                  .package = "torch")
  withr::local_envvar(c(PIGAUTO_DEBUG_GPU_MEM = "1"))
  expect_identical(gpu_mem_str(), "")
  expect_silent(result <- withVisible(gpu_mem_checkpoint("CPU checkpoint")))
  expect_false(result$visible)
  expect_null(result$value)
})

test_that("GNN compositional safety floor retains a valid simplex", {
  skip_if_no_libtorch()
  testthat::local_mocked_bindings(cuda_is_available = function() FALSE,
    backends_mps_is_available = function() FALSE, .package = "torch")
  set.seed(62)
  tree <- ape::rtree(20)
  p <- seq(.15, .85, length.out = 20)
  traits <- data.frame(p = p, q = 1 - p, row.names = tree$tip.label)
  traits[c(3, 7), ] <- NA_real_
  data <- preprocess_traits(traits, tree, multi_proportion_groups = list(comp = c("p", "q")))
  splits <- make_missing_splits(data$X_scaled, missing_frac = .25, val_frac = .5,
                                seed = 41, trait_map = data$trait_map)
  expect_warning(fit <- fit_pigauto(data, tree, splits = splits, gnn = TRUE,
    epochs = 2L, eval_every = 1L, hidden_dim = 8L, n_gnn_layers = 1L,
    corruption_rate = .8, corruption_ramp = 0L, dropout = 0, edge_dropout = 0,
    use_transformer_blocks = FALSE, use_attention = TRUE, verbose = FALSE,
    seed = 41, predict_method = "per_column", refine_steps = 2L,
    safety_floor = TRUE), "Small validation set")
  expect_identical(fit$trait_map[[1]]$type, "multi_proportion")
  expect_equal(unname(fit$r_cal_bm), c(1, 1))
  expect_equal(unname(fit$r_cal_gnn), c(0, 0))
  expect_equal(unname(fit$r_cal_mean), c(0, 0))
  pred <- predict(fit, return_se = FALSE)$imputed
  expect_identical(names(pred), c("p", "q"))
  expect_equal(dim(pred), c(20L, 2L))
  expect_true(all(is.finite(as.matrix(pred))) && all(as.matrix(pred) >= 0))
  expect_equal(unname(rowSums(pred)), rep(1, 20), tolerance = 1e-8)
})

test_that("three distinct split winners use the safe median fallback", {
  y <- c(-1.66710909640267, 1.08876938532769, -.415842429402306,
         -.36465927826914, 1.85592940783941, .541907004893235,
         1.23523793516361, -2.45768633900871, 1.35019298381003,
         .1148686074746, -1.18289480967858, -.110135472269955)
  bm <- c(1.39086793809816, -1.15463082528809, .653747855091367,
          .476129787419971, 1.15577388971215, .791322284128638,
          -1.54766181366128, -.582397211304587, -1.47146001603008,
          -.501436242579701, -1.10450956093243, .0777423526095165)
  gnn <- c(-.750401151450222, -.523861568989111, .258239919840661,
           .328866790235863, .44891244654269, 1.70632216263598,
           -1.36811575877369, -1.86155299553292, -.146869572432759,
           .766836998051743, -.0036644962697545, -.36653236855842)
  corners <- list(bm, gnn, rep(0, length(y)))
  # Independent calculations of the three half-A winners and half-B checks.
  winners <- vapply(seq_len(3L), function(k) {
    set.seed(100L + 17L + k - 1L)
    perm <- sample(length(y))
    a <- perm[1:6]; b <- perm[7:12]
    loss <- function(idx) vapply(corners, function(pred) mean((pred[idx] - y[idx])^2), numeric(1))
    la <- loss(a); lb <- loss(b)
    best <- which.min(la)
    gain <- (lb[1] - lb[best]) / max(lb[1], 1e-12)
    if (!is.finite(gain) || gain < .02) {
      best <- if (lb[3] <= lb[1]) 3L else 1L
    }
    best
  }, integer(1))
  expect_identical(winners, c(3L, 1L, 2L))
  weights <- diag(3)[winners, , drop = FALSE]
  expect_equal(apply(weights, 2L, stats::median), c(0, 0, 0))
  expect_lte(mean(y^2), mean((bm - y)^2))
  tm <- list(list(name = "x", type = "continuous", latent_cols = 1L))
  expect_warning(cal <- calibrate_gates(tm, matrix(bm), matrix(gnn), matrix(y),
    matrix(TRUE, length(y), 1), gate_grid = c(0, 1), gate_cap = 1,
    gate_method = "median_splits", gate_splits_B = 3L, seed = 100L,
    safety_floor = TRUE, mean_baseline_per_col = 0, simplex_step = 1), "95%")
  expect_equal(unname(c(cal$r_cal_bm, cal$r_cal_gnn, cal$r_cal_mean)), c(0, 0, 1))
})
