test_that('evaluation mode masks GNN input context', {
  skip_if_no_libtorch()
  torch::torch_set_num_threads(1L)
  set.seed(31)
  tree <- ape::rtree(30)
  pd <- preprocess_traits(data.frame(x = rnorm(30), y = rnorm(30),
                                     row.names = tree$tip.label), tree,
                           log_transform = FALSE)
  spl <- make_missing_splits(pd$X_scaled, seed = 3, trait_map = pd$trait_map)
  fit <- suppressWarnings(fit_pigauto(pd, tree, splits = spl, gnn = TRUE,
      epochs = 2L, eval_every = 1L, patience = 2L, verbose = FALSE))
  # Force a nonzero GNN contribution so a closed gate cannot hide leakage.
  p <- ncol(pd$X_scaled)
  fit$calibrated_gates <- rep(1, p)
  fit$r_cal_bm <- rep(0, p)
  fit$r_cal_gnn <- rep(1, p)
  fit$r_cal_mean <- rep(0, p)
  fit$mean_baseline_per_col <- rep(0, p)
  held <- c(spl$val_idx, spl$test_idx)
  pred <- predict(fit, evaluation = TRUE)
  expect_equal(pred$imputed_latent,
               predict(fit, .mask_observed_idx = held)$imputed_latent)
  altered <- fit
  altered$X_scaled[held] <- altered$X_scaled[held] + 100
  expect_equal(pred$imputed_latent,
               predict(altered, evaluation = TRUE)$imputed_latent)
  expect_error(evaluate_imputation(predict(fit), pd$X_scaled, spl), 'evaluation')
})
