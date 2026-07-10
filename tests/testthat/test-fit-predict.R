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

test_that("predict.pigauto_fit uses observed latent cells as DAE context", {
  td  <- make_test_data(n = 16, p = 1, seed = 20260511)
  td$df$tr1[c(3, 9, 12)] <- NA_real_
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, missing_frac = 0.1,
                             seed = 20260511, trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     epochs = 2L, eval_every = 1L, patience = 2L,
                     refine_steps = 1L, dropout = 0,
                     use_transformer_blocks = FALSE,
                     verbose = FALSE, seed = 20260511)

  p <- as.integer(fit$model_config$input_dim)
  fit$calibrated_gates <- rep(1, p)
  fit$r_cal_bm         <- rep(0, p)
  fit$r_cal_gnn        <- rep(1, p)
  fit$r_cal_mean       <- rep(0, p)
  fit$mean_baseline_per_col <- rep(0, p)

  pred_ref <- predict(fit, return_se = FALSE)$imputed_latent

  fit_shifted <- fit
  observed <- !is.na(fit_shifted$X_scaled)
  fit_shifted$X_scaled[observed] <- fit_shifted$X_scaled[observed] + 25
  pred_shifted <- predict(fit_shifted, return_se = FALSE)$imputed_latent

  expect_false(isTRUE(all.equal(pred_ref, pred_shifted, tolerance = 1e-8)),
               info = "prediction should depend on stored observed cells, not only the baseline")
})

test_that("predict.pigauto_fit rejects stale X_scaled dimensions", {
  td  <- make_test_data(n = 12, p = 1, seed = 20260512)
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, missing_frac = 0.1,
                             seed = 20260512, trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     epochs = 1L, eval_every = 1L, patience = 1L,
                     refine_steps = 1L, use_transformer_blocks = FALSE,
                     verbose = FALSE, seed = 20260512)

  fit$X_scaled <- fit$X_scaled[-1, , drop = FALSE]
  expect_error(
    predict(fit, return_se = FALSE),
    regexp = "Stored `X_scaled`"
  )
})

test_that("predict.pigauto_fit can mask held-out cells from DAE context", {
  td  <- make_test_data(n = 16, p = 1, seed = 20260513)
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, missing_frac = 0.1,
                             seed = 20260513, trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     epochs = 2L, eval_every = 1L, patience = 2L,
                     refine_steps = 1L, dropout = 0,
                     use_transformer_blocks = FALSE,
                     verbose = FALSE, seed = 20260513)

  p <- as.integer(fit$model_config$input_dim)
  fit$calibrated_gates <- rep(1, p)
  fit$r_cal_bm         <- rep(0, p)
  fit$r_cal_gnn        <- rep(1, p)
  fit$r_cal_mean       <- rep(0, p)
  fit$mean_baseline_per_col <- rep(0, p)

  hold_idx <- c(spl$val_idx, spl$test_idx)
  pred_ref <- predict(fit, return_se = FALSE,
                      .mask_observed_idx = hold_idx)$imputed_latent

  fit_shifted <- fit
  fit_shifted$X_scaled[hold_idx] <- fit_shifted$X_scaled[hold_idx] + 25
  pred_shifted <- predict(fit_shifted, return_se = FALSE,
                          .mask_observed_idx = hold_idx)$imputed_latent

  expect_equal(pred_shifted, pred_ref, tolerance = 1e-6)
})

test_that("predict.pigauto_fit adds cov_linear fixed effects outside the blend", {
  skip_if_no_libtorch()

  n <- 4L
  sp <- paste0("sp", seq_len(n))
  covariates <- matrix(seq_len(n) - 1, ncol = 1,
                       dimnames = list(sp, "temp"))

  model <- ResidualPhyloDAE(
    input_dim = 1L, hidden_dim = 4L, coord_dim = 1L, cov_dim = 3L,
    per_column_rs = TRUE, n_gnn_layers = 1L, gate_cap = 0.8,
    use_attention = FALSE, n_user_cov = 1L, dropout = 0,
    use_transformer_blocks = FALSE, n_heads = 1L, ffn_mult = 1L
  )
  torch::with_no_grad({
    model$cov_linear$weight$copy_(
      torch::torch_tensor(matrix(2, nrow = 1L, ncol = 1L))
    )
    model$cov_linear$bias$copy_(torch::torch_tensor(1))
  })

  fit <- structure(
    list(
      model_state = lapply(model$state_dict(), function(t) t$detach()$cpu()),
      model_config = list(
        input_dim = 1L, hidden_dim = 4L, k_eigen = 1L, cov_dim = 3L,
        per_column_rs = TRUE, n_gnn_layers = 1L, gate_cap = 0.8,
        use_attention = FALSE, n_user_cov = 1L, dropout = 0,
        refine_steps = 1L, use_transformer_blocks = FALSE,
        n_heads = 1L, ffn_mult = 1L
      ),
      graph = list(coords = matrix(0, n, 1), adj = diag(n), D_sq = NULL),
      baseline = list(
        mu = matrix(0, n, 1, dimnames = list(sp, "y")),
        se = matrix(0, n, 1, dimnames = list(sp, "y"))
      ),
      species_names = sp,
      obs_species = sp,
      obs_to_species = NULL,
      X_scaled = matrix(NA_real_, n, 1, dimnames = list(sp, "y")),
      n_species = n,
      n_obs = n,
      multi_obs = FALSE,
      trait_names = "y",
      latent_names = "y",
      trait_map = list(list(
        name = "y", type = "continuous", latent_cols = 1L,
        mean = 0, sd = 1, log_transform = FALSE
      )),
      calibrated_gates = c(y = 0),
      r_cal_bm = c(y = 0),
      r_cal_gnn = c(y = 0),
      r_cal_mean = c(y = 0),
      mean_baseline_per_col = c(y = 0),
      conformal_scores = NULL,
      covariates = covariates,
      cov_means = 0,
      cov_sds = 1,
      cov_names = "temp"
    ),
    class = "pigauto_fit"
  )

  pred <- predict(fit, return_se = FALSE)
  expected <- as.numeric(1 + 2 * covariates[, 1])

  expect_equal(as.numeric(pred$imputed_latent[, "y"]), expected,
               tolerance = 1e-6)
  expect_equal(pred$imputed$y, expected, tolerance = 1e-6)
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

test_that("impute() multi-obs aligns completed rows when input is shuffled", {
  # Regression test for the row-alignment bug found 2026-04-26.
  # When the input data.frame's species column is NOT in tree-tip order,
  # preprocess_traits internally reorders rows to tree-tip order, but the
  # returned `result$completed` must still align with the user's input
  # row order.  Specifically: each non-NA row must keep its original value,
  # and each NA row's imputed value must be consistent with the species
  # of THAT original row (not some other row).
  set.seed(20260426)
  tree <- ape::rtree(20)
  # Two obs per species, species column SHUFFLED (not tree-tip sorted).
  sp_shuffled <- sample(rep(tree$tip.label, each = 2L))
  # Strong species effect: each species has its own true mean.
  sp_truth <- setNames(rnorm(20, mean = 10, sd = 5), tree$tip.label)
  trait_full <- sp_truth[sp_shuffled] + rnorm(length(sp_shuffled), sd = 0.05)
  df <- data.frame(species = sp_shuffled, trait = trait_full,
                    stringsAsFactors = FALSE)

  # Mask exactly one of each species's two observations (cell-level MCAR).
  mask_idx <- vapply(unique(df$species),
                      function(sp) sample(which(df$species == sp), 1L),
                      integer(1))
  df_obs <- df
  df_obs$trait[mask_idx] <- NA

  result <- impute(df_obs, tree, species_col = "species",
                    epochs = 60L, verbose = FALSE, seed = 1L,
                    missing_frac = 0.0)

  # 1. Observed rows must retain their original values exactly.
  observed_idx <- setdiff(seq_len(nrow(df)), mask_idx)
  expect_equal(result$completed$trait[observed_idx],
               df$trait[observed_idx],
               tolerance = 1e-8)

  # 2. Imputed rows must be close to the OTHER observed obs of the same
  #    species (within ~0.5 of the true species mean).  This is the bug:
  #    if predictions are mis-aligned, masked rows get predictions from a
  #    different species and miss by 5-15 units instead of <1.
  imp_pred <- result$completed$trait[mask_idx]
  imp_truth <- sp_truth[as.character(df$species[mask_idx])]
  expect_true(all(abs(imp_pred - imp_truth) < 1.0),
              info = sprintf(
                "Imputed values must be near species truth. Worst diff: %.3f",
                max(abs(imp_pred - imp_truth))))

  # 3. Pearson correlation between predicted and species-truth must be high.
  expect_gt(cor(imp_pred, imp_truth), 0.95)
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

test_that("predict returns conformal intervals for proportion traits", {
  set.seed(112)
  n <- 50L
  tree <- ape::rtree(n)
  df <- data.frame(
    prop = stats::runif(n, min = 0.05, max = 0.95),
    row.names = tree$tip.label
  )
  df$prop[sample.int(n, 10L)] <- NA_real_

  res <- impute(df, tree,
                trait_types = c(prop = "proportion"),
                epochs = 20L, eval_every = 10L, patience = 5L,
                verbose = FALSE, seed = 112L, missing_frac = 0.30)
  pred <- predict(res$fit, return_se = TRUE)

  expect_true("prop" %in% names(res$fit$conformal_scores))
  expect_true(is.finite(res$fit$conformal_scores["prop"]))
  expect_true("prop" %in% colnames(pred$conformal_lower))
  expect_true("prop" %in% colnames(pred$conformal_upper))

  lo <- pred$conformal_lower[, "prop"]
  hi <- pred$conformal_upper[, "prop"]
  ok <- is.finite(lo) & is.finite(hi)
  expect_true(any(ok))
  expect_true(all(lo[ok] >= 0 & lo[ok] <= 1))
  expect_true(all(hi[ok] >= 0 & hi[ok] <= 1))
  expect_true(all(lo[ok] <= hi[ok]))
})

test_that("cross-validation conformal coverage handles proportion intervals", {
  trait_map <- list(list(
    name = "prop", type = "proportion", latent_cols = 1L,
    mean = 0, sd = 1
  ))
  pred <- list(
    conformal_lower = matrix(c(0.2, 0.2), ncol = 1,
                             dimnames = list(NULL, "prop")),
    conformal_upper = matrix(c(0.8, 0.8), ncol = 1,
                             dimnames = list(NULL, "prop"))
  )
  X_scaled <- matrix(stats::qlogis(c(0.5, 0.9)), ncol = 1)

  cov <- pigauto:::compute_conformal_coverage_fold(
    pred, X_scaled, test_idx = c(1L, 2L), trait_map, fold = 1L, rep = 1L
  )

  expect_equal(cov$type, "proportion")
  expect_equal(cov$coverage, 0.5)
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

test_that("predict.pigauto_fit catches wrong-length calibrated_gates override", {
  # Regression test for Opus C.4 (2026-04-28).
  # Previously: setting calibrated_gates to a vector of the wrong length
  # produced a cryptic torch shape error inside the encoder linear layer
  # ("input and weight.T shapes cannot be multiplied (n x p-1 and p x h)")
  # because (1, 0)-shaped weights broadcast with (n, p) baselines to (n, 0)
  # and corrupted the next refine step's enc1 input width.
  data(avonet300, tree300, package = "pigauto")
  traits <- avonet300[1:40, c("Species_Key", "Mass")]
  rownames(traits) <- traits$Species_Key
  traits$Species_Key <- NULL
  subtree <- ape::keep.tip(tree300, rownames(traits))

  res <- pigauto::impute(traits, subtree, epochs = 20L, verbose = FALSE,
                         seed = 1L, missing_frac = 0.0)
  fit <- res$fit
  p   <- as.integer(fit$model_config$input_dim)

  # 1. Wrong-length numeric override should error with a clear message.
  fit_bad <- fit
  fit_bad$calibrated_gates <- rep(0, p + 1L)
  fit_bad$r_cal_gnn        <- rep(0, p + 1L)
  fit_bad$r_cal_bm         <- rep(1, p + 1L)
  fit_bad$r_cal_mean       <- rep(0, p + 1L)
  expect_error(predict(fit_bad, return_se = FALSE),
               "calibrated_gates has length")

  # 2. Zero-length override (e.g. rep(0, length(NULL))) is treated as "no
  #    calibration" and predict falls through to the learned-gate path.
  fit_empty <- fit
  fit_empty$calibrated_gates <- numeric(0)
  fit_empty$r_cal_gnn        <- numeric(0)
  fit_empty$r_cal_bm         <- numeric(0)
  fit_empty$r_cal_mean       <- numeric(0)
  pred_empty <- predict(fit_empty, return_se = FALSE)
  expect_s3_class(pred_empty, "pigauto_pred")

  # 3. Legacy/scrubbed sub-slots that are zero-length should fall back to
  #    the valid calibrated_gates vector, not create zero-width torch tensors.
  fit_partial_empty <- fit
  fit_partial_empty$calibrated_gates      <- rep(0, p)
  fit_partial_empty$r_cal_gnn             <- numeric(0)
  fit_partial_empty$r_cal_bm              <- numeric(0)
  fit_partial_empty$r_cal_mean            <- numeric(0)
  fit_partial_empty$mean_baseline_per_col <- numeric(0)
  pred_partial_empty <- predict(fit_partial_empty, return_se = FALSE)
  expect_s3_class(pred_partial_empty, "pigauto_pred")

  # 4. Correct-length zero gate (full BM, no GNN) succeeds.
  fit_zero <- fit
  fit_zero$calibrated_gates <- rep(0, p)
  fit_zero$r_cal_gnn        <- rep(0, p)
  fit_zero$r_cal_bm         <- rep(1, p)
  fit_zero$r_cal_mean       <- rep(0, p)
  pred_zero <- predict(fit_zero, return_se = FALSE)
  expect_s3_class(pred_zero, "pigauto_pred")
})

test_that("conformal_split_val=TRUE changes conformal scores (no double-dipping)", {
  # Regression test for Opus C.3 (2026-04-28).
  # When `conformal_split_val = TRUE`, the calibration half and the
  # conformal half must not share cells — otherwise post-selection
  # bias makes the conformal quantile a biased estimator of the true
  # residual quantile.  Default is FALSE (legacy single-set path) to
  # preserve bench-grade RMSE on small-val datasets; this test exercises
  # the opt-in path and verifies that toggling the flag produces a
  # different conformal-score vector (the split path uses different
  # cells and therefore a different empirical quantile).  Use a dataset
  # large enough that each trait has at least
  # `2 * min_val_cells = 40` cells so the per-column split actually
  # fires (smaller traits silently fall back to the single-set path).
  # n=200 + 50% missing + 40% val per masked-trait gives ~40 val cells
  # per column, just above the default `2 * min_val_cells = 40` split
  # threshold so the per-column split actually fires here.
  set.seed(90)
  tree <- ape::rtree(200)
  df   <- data.frame(
    row.names = tree$tip.label,
    tr1 = abs(stats::rnorm(200)) + 0.5,
    tr2 = abs(stats::rnorm(200)) + 0.5
  )
  pd  <- preprocess_traits(df, tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 90, trait_map = pd$trait_map,
                              missing_frac = 0.50, val_frac = 0.40)

  fit_split <- fit_pigauto(pd, tree, splits = spl,
                            epochs = 30L, eval_every = 10L, patience = 5L,
                            verbose = FALSE, seed = 90,
                            conformal_split_val = TRUE)
  fit_full  <- fit_pigauto(pd, tree, splits = spl,
                            epochs = 30L, eval_every = 10L, patience = 5L,
                            verbose = FALSE, seed = 90,
                            conformal_split_val = FALSE)

  expect_true(!is.null(fit_split$conformal_scores))
  expect_true(!is.null(fit_full$conformal_scores))
  # When val is generous enough, both should be finite and the split
  # version should differ from the post-selected estimate (proof that
  # the toggle actually changes the conformal-score computation).
  ok <- is.finite(fit_split$conformal_scores) &
        is.finite(fit_full$conformal_scores)
  if (any(ok)) {
    expect_false(isTRUE(all.equal(fit_split$conformal_scores[ok],
                                    fit_full$conformal_scores[ok])))
  }
})
