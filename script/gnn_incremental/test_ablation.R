#!/usr/bin/env Rscript

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required for this local acceptance check.", call. = FALSE)
}
devtools::load_all(".", quiet = TRUE)
source("script/gnn_incremental/dgp.R")
source("script/gnn_incremental/ablation.R")

case <- make_gnn_case(
  n_species = 80L, obs_per = 3L, n_covs = 5L,
  variance_multiplier = 0.3, f_type = "nonlinear",
  beta_strength = 1, seed = 782031L,
  sp_miss_frac = 0.1, within_miss = 0.05
)
set.seed(782032L)
input_order <- sample(seq_len(nrow(case$df_observed)))
traits <- case$df_observed[input_order, c("species", "y"), drop = FALSE]
covariates <- case$df_observed[input_order,
                               c("species", case$predictor_names), drop = FALSE]

result <- impute(
  traits = traits, tree = case$tree, species_col = "species",
  covariates = covariates, gnn = TRUE, phylo_signal_gate = FALSE,
  missing_frac = 0.5, epochs = 3L, eval_every = 1L, patience = 2L,
  hidden_dim = 16L, trait_embed_dim = 8L, k_eigen = 4L,
  predict_method = "per_column", seed = 782033L, verbose = FALSE
)
fit <- result$fit
ablated <- ablate_gnn_in_fit(fit)

stopifnot(identical(ablated$model_state, fit$model_state))
stopifnot(identical(ablated$model_config, fit$model_config))
stopifnot(identical(ablated$covariates, fit$covariates))
stopifnot(identical(ablated$r_cal_mean, fit$r_cal_mean))
stopifnot(identical(ablated$mean_baseline_per_col,
                    fit$mean_baseline_per_col))
stopifnot(all(ablated$r_cal_gnn == 0))
stopifnot(isTRUE(all.equal(ablated$r_cal_bm,
                            fit$r_cal_bm + fit$r_cal_gnn)))
stopifnot(isTRUE(all.equal(
  ablated$r_cal_bm + ablated$r_cal_gnn + ablated$r_cal_mean,
  fit$r_cal_bm + fit$r_cal_gnn + fit$r_cal_mean
)))

prediction_off <- predict(ablated, return_se = FALSE)
prediction_full <- result$prediction
order_internal_to_input <- result$data$input_row_order
stopifnot(length(order_internal_to_input) == nrow(traits))
stopifnot(!anyNA(order_internal_to_input))
expected <- prediction_full$imputed$y
observed <- prediction_off$imputed$y
stopifnot(length(expected) == length(observed))
stopifnot(all(is.finite(expected)))
stopifnot(all(is.finite(observed)))
input_idx <- match(seq_len(nrow(traits)), order_internal_to_input)
stopifnot(!anyNA(input_idx))
outer_mask <- case$mask[input_order]
stopifnot(all(is.finite(expected[input_idx[outer_mask]])))
stopifnot(all(is.finite(observed[input_idx[outer_mask]])))
stopifnot(isTRUE(all.equal(
  expected[input_idx[outer_mask]],
  result$completed$y[outer_mask], tolerance = 1e-10
)))
truth <- case$df_complete$y[input_order]
stopifnot(is.finite(mean((truth[outer_mask] -
  expected[input_idx[outer_mask]])^2)))
stopifnot(is.finite(mean((truth[outer_mask] -
  observed[input_idx[outer_mask]])^2)))
stopifnot(all(result$completed$y[!outer_mask] == traits$y[!outer_mask],
               na.rm = TRUE))

# Exercise the enabled policy as well. The fitted signal must equal an
# explicit calculation with this fit's validation and test cells masked.
gate_on <- impute(
  traits = traits, tree = case$tree, species_col = "species",
  covariates = covariates, gnn = TRUE, phylo_signal_gate = TRUE,
  missing_frac = 0.5, epochs = 3L, eval_every = 1L, patience = 2L,
  hidden_dim = 16L, trait_embed_dim = 8L, k_eigen = 4L,
  predict_method = "per_column", seed = 782033L, verbose = FALSE
)
expected_signal <- compute_phylo_signal_per_trait(
  data = gate_on$data, tree = gate_on$tree, method = "lambda", min_tips = 20L,
  splits = gate_on$splits
)
stopifnot(isTRUE(all.equal(gate_on$fit$phylo_signal_per_trait,
                            expected_signal)))
stopifnot(nrow(gate_on$fit$history) > 0L)
stopifnot(all(is.finite(gate_on$fit$phylo_signal_per_trait)))
stopifnot(any(gate_on$fit$phylo_gate_triggered))
triggered <- which(gate_on$fit$phylo_gate_triggered)
stopifnot(all(gate_on$fit$r_cal_bm[triggered] == 1))
stopifnot(all(gate_on$fit$r_cal_gnn[triggered] == 0))
stopifnot(all(gate_on$fit$r_cal_mean[triggered] == 0))
delivered_gate_on <- predict(gate_on$fit, return_se = FALSE)
stopifnot(all(is.finite(delivered_gate_on$imputed$y)))
fixed_weight <- as.matrix(as.array(
  gate_on$fit$model_state[["cov_linear.weight"]]$detach()$cpu()
))
fixed_bias <- as.numeric(as.array(
  gate_on$fit$model_state[["cov_linear.bias"]]$detach()$cpu()
))
fixed_effect <- gate_on$fit$covariates %*% t(fixed_weight)
fixed_effect <- sweep(fixed_effect, 2L, fixed_bias, "+")
bm_corner <- gate_on$fit$baseline$mu[gate_on$fit$obs_to_species, , drop = FALSE] +
  fixed_effect
gate_latent <- which(gate_on$fit$phylo_gate_triggered)
stopifnot(length(gate_latent) == 1L)
outer_internal <- input_idx[outer_mask]
stopifnot(isTRUE(all.equal(
  delivered_gate_on$imputed_latent[outer_internal, gate_latent],
  bm_corner[outer_internal, gate_latent], tolerance = 1e-5
)))
cat("ABLATION_VALIDATION_PASS\n")
