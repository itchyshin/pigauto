#' Impute missing traits using a fitted pigauto model
#'
#' Runs a single forward pass through the fitted model and returns imputed
#' trait values back-transformed to the original scale.  Supports all
#' trait types (continuous, binary, categorical, ordinal, count, proportion)
#' and MC dropout for multiple imputation (when \code{n_imputations > 1}). The fitted model is a gated ensemble
#' of a phylogenetic baseline and a graph neural network correction;
#' prediction is the per-trait blend
#' \code{(1 - r_cal) * baseline + r_cal * delta_GNN}.
#'
#' @details
#' When \code{n_imputations > 1}, the model runs in train mode (dropout
#' active) to generate M stochastic forward passes.  Point estimates are
#' computed as the mean (continuous, count) or mode (binary, categorical,
#' ordinal) across passes.  Between-imputation variance provides proper
#' uncertainty estimates following Rubin's rules.  The M complete datasets
#' are returned in \code{imputed_datasets} for downstream pooling.
#'
#' **Decoding per type:**
#' \describe{
#'   \item{continuous}{reverse z-score, then \code{exp()} if log-transformed}
#'   \item{binary}{\code{sigmoid(latent)} to probability, round to 0/1}
#'   \item{count}{reverse z-score of log1p, \code{expm1()}, round, clip >= 0}
#'   \item{ordinal}{reverse z-score, round to nearest valid integer level}
#'   \item{categorical}{\code{softmax()} over K latent columns, argmax}
#' }
#'
#' @param object object of class \code{"pigauto_fit"}.
#' @param newdata \code{NULL} (use the training data) or a
#'   \code{"pigauto_data"} object for new species.
#' @param return_se logical. Compute standard errors? (default \code{TRUE}).
#' @param n_imputations integer. Number of MC dropout forward passes (default
#'   \code{1L}).  Set to e.g. 10 or 20 for proper multiple imputation with
#'   between-imputation variance.
#' @param ... ignored.
#' @return A list of class \code{"pigauto_pred"} with:
#'   \describe{
#'     \item{imputed}{data.frame of imputed values in original scale with
#'       proper R types (numeric, integer, factor, ordered).}
#'     \item{imputed_latent}{Numeric matrix (n x p_latent) of predictions in
#'       latent scale.}
#'     \item{se}{Numeric matrix (n x n_original_traits) of per-cell
#'       uncertainty.  Continuous/count: SE in original scale.
#'       Binary: Bernoulli SD \code{sqrt(p*(1-p))}.
#'       Categorical: entropy of probability distribution.
#'       Ordinal: SE in integer scale.  \code{NULL} if
#'       \code{return_se = FALSE}.}
#'     \item{probabilities}{Named list.  Binary traits: numeric probability
#'       vector.  Categorical traits: n x K probability matrix.  Other
#'       types: not present.}
#'     \item{imputed_datasets}{List of M data.frames when
#'       \code{n_imputations > 1}; \code{NULL} otherwise.}
#'     \item{trait_map}{Trait map from the fitted model.}
#'     \item{species_names}{Character vector.}
#'     \item{trait_names}{Character vector.}
#'     \item{n_imputations}{Integer, number of imputations performed.}
#'   }
#' @examples
#' \dontrun{
#' pred <- predict(fit, return_se = TRUE)
#' pred$imputed        # data.frame, original scale
#' pred$se             # matrix, uncertainty
#' pred$probabilities  # list of prob vectors/matrices
#'
#' # Multiple imputation (MC dropout)
#' pred10 <- predict(fit, n_imputations = 10)
#' pred10$imputed_datasets  # 10 complete data.frames
#' }
#' @importFrom torch torch_tensor torch_float with_no_grad torch_zeros
#' @importFrom torch torch_cat
#' @export
predict.pigauto_fit <- function(object, newdata = NULL, return_se = TRUE,
                                n_imputations = 1L, ...) {
  cfg       <- object$model_config
  device    <- get_device()
  trait_map <- object$trait_map
  has_trait_map <- !is.null(trait_map)

  # Reconstruct model (backward compat with old saves that lack new params)
  per_col      <- isTRUE(cfg$per_column_rs)
  n_gnn_layers <- cfg$n_gnn_layers %||% 1L
  gate_cap     <- cfg$gate_cap %||% 0.5
  use_attention <- cfg$use_attention %||% FALSE
  n_user_cov   <- cfg$n_user_cov %||% 0L
  model <- ResidualPhyloDAE(
    input_dim     = as.integer(cfg$input_dim),
    hidden_dim    = as.integer(cfg$hidden_dim),
    coord_dim     = as.integer(cfg$k_eigen),
    cov_dim       = as.integer(cfg$cov_dim),
    per_column_rs = per_col,
    n_gnn_layers  = as.integer(n_gnn_layers),
    gate_cap      = gate_cap,
    use_attention = use_attention,
    n_user_cov    = as.integer(n_user_cov)
  )
  model$to(device = device)
  model$load_state_dict(object$model_state)

  # Calibrated gates override learned gates
  calibrated_gates <- object$calibrated_gates  # NULL if not available
  use_calibrated   <- !is.null(calibrated_gates)

  # Prepare data
  multi_obs <- isTRUE(object$multi_obs)
  if (is.null(newdata)) {
    MU_species <- object$baseline$mu    # n_species x p
    coords     <- object$graph$coords   # n_species x k
    adj        <- object$graph$adj      # n_species x n_species
    obs_to_sp  <- object$obs_to_species # NULL or integer vector
  } else {
    stop("newdata support not yet implemented.")
  }

  # For multi-obs, expand baseline to observation level
  if (multi_obs) {
    n_obs     <- object$n_obs
    n_species <- object$n_species
    MU <- MU_species[obs_to_sp, , drop = FALSE]
    rownames(MU) <- NULL
  } else {
    n_obs     <- nrow(MU_species)
    n_species <- n_obs
    MU        <- MU_species
    obs_to_sp <- NULL
  }

  n <- n_obs   # rows in output
  p <- ncol(MU)
  n_imp <- as.integer(n_imputations)

  t_X_fill <- torch::torch_tensor(MU, dtype = torch::torch_float(),
                                  device = device)
  t_MU     <- torch::torch_tensor(MU, dtype = torch::torch_float(),
                                  device = device)
  t_coords <- torch::torch_tensor(coords, dtype = torch::torch_float(),
                                  device = device)
  t_adj    <- torch::torch_tensor(adj,    dtype = torch::torch_float(),
                                  device = device)

  # Observation-to-species mapping tensor
  if (multi_obs) {
    t_obs_to_sp <- torch::torch_tensor(
      as.integer(obs_to_sp), dtype = torch::torch_long(), device = device
    )
  } else {
    t_obs_to_sp <- NULL
  }

  # ---- Covariates (environmental conditioners) ------------------------------
  has_covariates <- !is.null(object$covariates)
  t_covariates   <- NULL
  if (has_covariates) {
    t_covariates <- torch::torch_tensor(
      object$covariates, dtype = torch::torch_float(), device = device
    )
  }

  # ---- Inference (single or MC dropout) ------------------------------------
  latent_runs <- vector("list", n_imp)

  # Pre-create calibrated gates tensor (once, outside the loop)
  if (use_calibrated) {
    t_cal_gates <- torch::torch_tensor(
      calibrated_gates, dtype = torch::torch_float(), device = device
    )$unsqueeze(1L)
  }

  # BM SE tensor for MC dropout BM-draw injection (latent / z-score scale).
  # BM_SE = 0 for observed cells → observed values never perturbed.
  # BM_SE > 0 for originally-missing cells → t_BM_draw ~ N(BM_mu, BM_se)
  # per imputation, held fixed across refine steps so each m draws ONE
  # consistent BM posterior sample.  The blend then uses t_BM_draw instead
  # of t_MU in the (1 - gate) term:
  #   pred = (1-r)*t_BM_draw + r*delta_dropout
  # → when gate=0: pred = t_BM_draw  ← proper BM posterior draw, non-zero variance
  # → when gate>0: both BM draws and GNN dropout contribute variance
  if (n_imp > 1L) {
    bse_mat <- object$baseline$se            # n_species x p_latent
    if (multi_obs) bse_mat <- bse_mat[obs_to_sp, , drop = FALSE]
    t_BM_SE <- torch::torch_tensor(
      bse_mat, dtype = torch::torch_float(), device = device
    )
  }

  # Default: deterministic baseline for n_imp == 1 path (same as t_MU).
  # Overwritten each iteration when n_imp > 1.
  t_BM_draw <- t_MU

  for (m in seq_len(n_imp)) {
    if (n_imp == 1L) {
      model$eval()
      X_iter <- t_X_fill$clone()
      # t_BM_draw stays as t_MU  → blend is identical to original single-pass
    } else {
      model$train()   # activates dropout on GNN hidden layers
      # Draw one BM posterior sample for this imputation.
      # BM_SE = 0 for observed cells → X_iter unchanged at observed positions.
      noise     <- torch::torch_randn(c(n, p), dtype = torch::torch_float(),
                                     device = device)
      t_BM_draw <- t_MU + noise * t_BM_SE   # fixed BM draw for this m
      X_iter    <- t_BM_draw$clone()          # GNN input starts from BM draw
    }
    torch::with_no_grad({
      mask_ind0 <- torch::torch_zeros(c(n, 1L), device = device)
      for (step in seq_len(cfg$refine_steps)) {
        cov_parts0 <- list(t_MU, mask_ind0)
        if (has_covariates) cov_parts0[[length(cov_parts0) + 1L]] <- t_covariates
        covs0 <- torch::torch_cat(cov_parts0, dim = 2L)
        out   <- model(X_iter, t_coords, covs0, t_adj, t_obs_to_sp)
        # Use t_BM_draw (the BM posterior sample) in the baseline term so that
        # between-imputation variance is non-zero even when the gate is 0.
        if (use_calibrated) {
          pred <- (1 - t_cal_gates) * t_BM_draw + t_cal_gates * out$delta
        } else {
          pred <- (1 - out$rs) * t_BM_draw + out$rs * out$delta
        }
        X_iter <- pred
      }
    })
    latent_runs[[m]] <- as.matrix(X_iter$cpu())
  }

  # Residual scale (per-column vector or legacy scalar)
  rs_val <- as.numeric(out$rs$cpu()$squeeze())

  # ---- Legacy path: no trait_map (old pigauto_fit objects) -----------------
  if (!has_trait_map) {
    return(decode_continuous_legacy(latent_runs, object, rs_val,
                                   return_se, n_imp))
  }

  # ---- Mixed-type decoding -------------------------------------------------
  row_labels <- if (multi_obs) object$obs_species else object$species_names
  decode_results <- lapply(latent_runs, function(lat) {
    decode_from_latent(lat, trait_map, row_labels)
  })

  if (n_imp == 1L) {
    imputed     <- decode_results[[1]]$imputed
    probs       <- decode_results[[1]]$probabilities
    latent_pred <- latent_runs[[1]]
  } else {
    pool <- pool_imputations(decode_results, latent_runs, trait_map)
    imputed     <- pool$imputed
    probs       <- pool$probabilities
    latent_pred <- pool$latent_mean
  }

  rownames(latent_pred) <- row_labels
  lat_names <- latent_names_from_map(trait_map)
  if (length(lat_names) == ncol(latent_pred)) {
    colnames(latent_pred) <- lat_names
  }

  # ---- SE ------------------------------------------------------------------
  se_mat <- NULL
  se_latent_mat <- NULL
  if (return_se) {
    # Expand species-level baseline SE to obs level if multi_obs
    bse <- object$baseline$se
    if (multi_obs) {
      bse <- bse[obs_to_sp, , drop = FALSE]
      rownames(bse) <- NULL
    }
    # Row names for output
    row_names <- if (multi_obs) object$obs_species else object$species_names

    if (n_imp > 1L) {
      se_mat <- compute_mc_se(decode_results, trait_map,
                              bse, latent_runs[[1]],
                              row_names)
    } else {
      se_mat <- compute_single_se(latent_runs[[1]], probs, trait_map,
                                  bse, row_names)
    }

    se_latent_mat <- compute_latent_se(latent_runs, trait_map,
                                       bse, row_names)
  }

  # ---- Conformal prediction intervals -----------------------------------
  conformal_lower <- NULL
  conformal_upper <- NULL
  conformal_scores_out <- object$conformal_scores

  if (!is.null(conformal_scores_out) && has_trait_map) {
    n_traits <- length(trait_map)
    conformal_lower <- matrix(NA_real_, nrow = n, ncol = n_traits)
    conformal_upper <- matrix(NA_real_, nrow = n, ncol = n_traits)
    colnames(conformal_lower) <- vapply(trait_map, "[[", character(1), "name")
    colnames(conformal_upper) <- vapply(trait_map, "[[", character(1), "name")
    rownames(conformal_lower) <- row_labels
    rownames(conformal_upper) <- row_labels

    for (tm in trait_map) {
      nm <- tm$name
      lc <- tm$latent_cols

      if (!(tm$type %in% c("continuous", "count", "ordinal"))) next
      if (is.na(conformal_scores_out[nm])) next

      q <- conformal_scores_out[nm]

      if (tm$type == "continuous") {
        # Convert conformal score from latent to original scale
        pred_latent <- latent_pred[, lc[1]]
        pred_low  <- (pred_latent - q) * tm$sd + tm$mean
        pred_high <- (pred_latent + q) * tm$sd + tm$mean
        if (isTRUE(tm$log_transform)) {
          conformal_lower[, nm] <- exp(pred_low)
          conformal_upper[, nm] <- exp(pred_high)
        } else {
          conformal_lower[, nm] <- pred_low
          conformal_upper[, nm] <- pred_high
        }

      } else if (tm$type == "count") {
        pred_latent <- latent_pred[, lc[1]]
        pred_low  <- (pred_latent - q) * tm$sd + tm$mean
        pred_high <- (pred_latent + q) * tm$sd + tm$mean
        conformal_lower[, nm] <- pmax(expm1(pred_low), 0)
        conformal_upper[, nm] <- expm1(pred_high)

      } else if (tm$type == "ordinal") {
        pred_latent <- latent_pred[, lc[1]]
        pred_low  <- (pred_latent - q) * tm$sd + tm$mean
        pred_high <- (pred_latent + q) * tm$sd + tm$mean
        K <- length(tm$levels)
        conformal_lower[, nm] <- pmax(round(pred_low), 0)
        conformal_upper[, nm] <- pmin(round(pred_high), K - 1L)
      }
    }
  }

  # ---- Imputed datasets for multiple imputation ----------------------------
  imputed_datasets <- if (n_imp > 1L) {
    lapply(decode_results, "[[", "imputed")
  }

  structure(
    list(
      imputed          = imputed,
      imputed_latent   = latent_pred,
      se               = se_mat,
      se_latent        = se_latent_mat,
      probabilities    = probs,
      imputed_datasets = imputed_datasets,
      conformal_lower  = conformal_lower,
      conformal_upper  = conformal_upper,
      conformal_scores = conformal_scores_out,
      calibrated_gates = calibrated_gates,
      trait_map        = trait_map,
      species_names    = object$species_names,
      obs_species      = if (multi_obs) object$obs_species else NULL,
      obs_to_species   = obs_to_sp,
      multi_obs        = multi_obs,
      trait_names      = object$trait_names,
      n_imputations    = n_imp
    ),
    class = "pigauto_pred"
  )
}


# ---- Internal: legacy all-continuous decoding ---------------------------------

decode_continuous_legacy <- function(latent_runs, object, rs_val,
                                    return_se, n_imp) {
  if (n_imp == 1L) {
    pred_scaled <- latent_runs[[1]]
  } else {
    pred_scaled <- Reduce("+", latent_runs) / n_imp
  }

  rownames(pred_scaled) <- object$species_names
  colnames(pred_scaled) <- object$trait_names

  sds    <- object$norm$sds
  means  <- object$norm$means
  is_log <- isTRUE(object$norm$log_transform)

  # Back-transform predictions
  imputed <- sweep(pred_scaled, 2, sds, "*")
  imputed <- sweep(imputed, 2, means, "+")
  if (is_log) imputed <- exp(imputed)

  se_out <- NULL
  if (return_se) {
    if (n_imp > 1L) {
      all_bt <- lapply(latent_runs, function(lr) {
        bt <- sweep(lr, 2, sds, "*")
        bt <- sweep(bt, 2, means, "+")
        if (is_log) bt <- exp(bt)
        bt
      })
      se_mc <- matrix(NA_real_, nrow(imputed), ncol(imputed),
                      dimnames = dimnames(imputed))
      for (j in seq_len(ncol(imputed))) {
        vals <- sapply(all_bt, function(m) m[, j])
        se_mc[, j] <- apply(vals, 1, stats::sd)
      }
      # Add baseline SE in quadrature
      bm_se <- object$baseline$se
      se_bm_orig <- sweep(bm_se, 2, sds, "*")
      if (is_log) {
        pred_log <- sweep(pred_scaled, 2, sds, "*")
        pred_log <- sweep(pred_log, 2, means, "+")
        se_bm_orig <- exp(pred_log) * se_bm_orig
      }
      se_out <- sqrt(se_bm_orig^2 + se_mc^2)
    } else {
      bm_se   <- object$baseline$se
      se_orig <- sweep(bm_se, 2, sds, "*")
      if (is_log) {
        pred_log <- sweep(pred_scaled, 2, sds, "*")
        pred_log <- sweep(pred_log, 2, means, "+")
        se_orig  <- exp(pred_log) * se_orig  # delta method
      }
      rownames(se_orig) <- object$species_names
      colnames(se_orig) <- object$trait_names
      se_out <- se_orig
    }
  }

  list(imputed = imputed, se = se_out)
}


# ---- Internal: decode latent matrix to original-scale data.frame ---------------

decode_from_latent <- function(latent_mat, trait_map, species_names) {
  n <- nrow(latent_mat)
  # Use make.unique for multi-obs (species_names may have duplicates)
  imputed <- data.frame(row.names = make.unique(species_names, sep = "."))
  probs <- list()

  for (tm in trait_map) {
    nm <- tm$name
    lc <- tm$latent_cols

    if (tm$type == "continuous") {
      vals <- latent_mat[, lc] * tm$sd + tm$mean
      if (tm$log_transform) vals <- exp(vals)
      imputed[[nm]] <- vals

    } else if (tm$type == "count") {
      vals <- latent_mat[, lc] * tm$sd + tm$mean
      vals <- expm1(vals)
      vals <- pmax(round(vals), 0)
      imputed[[nm]] <- as.integer(vals)

    } else if (tm$type == "ordinal") {
      vals <- latent_mat[, lc] * tm$sd + tm$mean
      K <- length(tm$levels)
      int_vals <- as.integer(pmin(pmax(round(vals), 0), K - 1L))
      imputed[[nm]] <- factor(tm$levels[int_vals + 1L],
                              levels = tm$levels, ordered = TRUE)

    } else if (tm$type == "binary") {
      prob <- expit(latent_mat[, lc])
      probs[[nm]] <- prob
      pred_class <- ifelse(prob >= 0.5, tm$levels[2], tm$levels[1])
      imputed[[nm]] <- factor(pred_class, levels = tm$levels)

    } else if (tm$type == "categorical") {
      logits <- latent_mat[, lc, drop = FALSE]
      prob_mat <- softmax_rows(logits)
      colnames(prob_mat) <- tm$levels
      rownames(prob_mat) <- species_names
      probs[[nm]] <- prob_mat
      pred_idx <- apply(prob_mat, 1, which.max)
      imputed[[nm]] <- factor(tm$levels[pred_idx], levels = tm$levels)

    } else if (tm$type == "proportion") {
      # Inverse logit: latent -> logit scale -> (0,1)
      vals_logit <- latent_mat[, lc] * tm$sd + tm$mean
      vals <- stats::plogis(vals_logit)
      imputed[[nm]] <- vals

    } else if (tm$type == "zi_count") {
      # Expected value: E[X] = P(non-zero) * E[count | non-zero]
      p_nz <- expit(latent_mat[, lc[1]])
      count_logscale <- latent_mat[, lc[2]] * tm$sd + tm$mean
      count_hat <- expm1(count_logscale)
      count_hat <- pmax(count_hat, 0)
      ev <- p_nz * count_hat
      imputed[[nm]] <- as.integer(pmax(round(ev), 0L))
      probs[[nm]] <- p_nz  # probability of non-zero
    }
  }

  list(imputed = imputed, probabilities = probs)
}


# ---- Internal: pool multiple imputations ----------------------------------------

pool_imputations <- function(decode_results, latent_runs, trait_map) {
  M <- length(decode_results)
  latent_mean <- Reduce("+", latent_runs) / M
  sp_names <- rownames(decode_results[[1]]$imputed)
  imputed <- data.frame(row.names = sp_names)
  probs <- list()

  for (tm in trait_map) {
    nm <- tm$name

    if (tm$type == "continuous") {
      vals <- rowMeans(sapply(decode_results, function(dr) dr$imputed[[nm]]))
      imputed[[nm]] <- vals

    } else if (tm$type == "count") {
      vals <- rowMeans(sapply(decode_results, function(dr) {
        as.numeric(dr$imputed[[nm]])
      }))
      imputed[[nm]] <- as.integer(pmax(round(vals), 0L))

    } else if (tm$type == "ordinal") {
      K <- length(tm$levels)
      int_vals <- rowMeans(sapply(decode_results, function(dr) {
        as.integer(dr$imputed[[nm]]) - 1L
      }))
      int_vals <- as.integer(pmin(pmax(round(int_vals), 0), K - 1L))
      imputed[[nm]] <- factor(tm$levels[int_vals + 1L],
                              levels = tm$levels, ordered = TRUE)

    } else if (tm$type == "binary") {
      avg_prob <- rowMeans(sapply(decode_results, function(dr) {
        dr$probabilities[[nm]]
      }))
      probs[[nm]] <- avg_prob
      pred_class <- ifelse(avg_prob >= 0.5, tm$levels[2], tm$levels[1])
      imputed[[nm]] <- factor(pred_class, levels = tm$levels)

    } else if (tm$type == "categorical") {
      prob_mats <- lapply(decode_results, function(dr) dr$probabilities[[nm]])
      avg_prob <- Reduce("+", prob_mats) / M
      colnames(avg_prob) <- tm$levels
      probs[[nm]] <- avg_prob
      pred_idx <- apply(avg_prob, 1, which.max)
      imputed[[nm]] <- factor(tm$levels[pred_idx], levels = tm$levels)

    } else if (tm$type == "proportion") {
      vals <- rowMeans(sapply(decode_results, function(dr) dr$imputed[[nm]]))
      imputed[[nm]] <- vals

    } else if (tm$type == "zi_count") {
      # Average the expected values across imputations
      vals <- rowMeans(sapply(decode_results, function(dr) {
        as.numeric(dr$imputed[[nm]])
      }))
      imputed[[nm]] <- as.integer(pmax(round(vals), 0L))
      # Average the P(non-zero) probabilities
      avg_pnz <- rowMeans(sapply(decode_results, function(dr) {
        dr$probabilities[[nm]]
      }))
      probs[[nm]] <- avg_pnz
    }
  }

  list(imputed = imputed, probabilities = probs, latent_mean = latent_mean)
}


# ---- Internal: MC dropout SE computation ----------------------------------------

compute_mc_se <- function(decode_results, trait_map, baseline_se,
                          latent_mean, species_names) {
  # Combines two sources of uncertainty:
  # 1. Baseline (BM) SE -- phylogenetic imputation uncertainty
  # 2. MC dropout between-imputation SD -- GNN correction uncertainty
  # Combined in quadrature: SE_total = sqrt(SE_BM^2 + SE_MC^2)
  M <- length(decode_results)
  n <- nrow(decode_results[[1]]$imputed)
  n_traits <- length(trait_map)

  se_mat <- matrix(NA_real_, nrow = n, ncol = n_traits)
  colnames(se_mat) <- vapply(trait_map, "[[", character(1), "name")
  rownames(se_mat) <- species_names

  for (tm in trait_map) {
    nm <- tm$name
    lc <- tm$latent_cols

    if (tm$type %in% c("continuous", "count")) {
      # MC dropout between-imputation SD in original scale
      vals <- sapply(decode_results, function(dr) as.numeric(dr$imputed[[nm]]))
      se_mc <- apply(vals, 1, stats::sd)

      # Baseline SE in original scale
      se_latent <- baseline_se[, lc[1]]
      se_bm_orig <- se_latent * tm$sd
      if (tm$type == "continuous" && isTRUE(tm$log_transform)) {
        pred_log <- latent_mean[, lc[1]] * tm$sd + tm$mean
        se_bm_orig <- exp(pred_log) * se_bm_orig
      } else if (tm$type == "count") {
        pred_log1p <- latent_mean[, lc[1]] * tm$sd + tm$mean
        se_bm_orig <- exp(pred_log1p) * se_bm_orig
      }

      # Combine in quadrature
      se_mat[, nm] <- sqrt(se_bm_orig^2 + se_mc^2)

    } else if (tm$type == "ordinal") {
      vals <- sapply(decode_results, function(dr) {
        as.integer(dr$imputed[[nm]]) - 1L
      })
      se_mc <- apply(vals, 1, stats::sd)

      se_latent <- baseline_se[, lc[1]]
      se_bm_orig <- se_latent * tm$sd

      se_mat[, nm] <- sqrt(se_bm_orig^2 + se_mc^2)

    } else if (tm$type == "binary") {
      # Uncertainty = min(p, 1-p) of the mean probability across MC runs.
      # Between-imputation SD of probabilities is also available but
      # min(p, 1-p) is more interpretable (probability of being wrong).
      prob_vals <- sapply(decode_results, function(dr) dr$probabilities[[nm]])
      avg_prob  <- rowMeans(prob_vals)
      se_mat[, nm] <- pmin(avg_prob, 1 - avg_prob)

    } else if (tm$type == "categorical") {
      # Uncertainty = 1 - max class probability of the mean prob across runs.
      prob_mats <- lapply(decode_results, function(dr) dr$probabilities[[nm]])
      avg_prob  <- Reduce("+", prob_mats) / M
      se_mat[, nm] <- 1 - apply(avg_prob, 1, max)

    } else if (tm$type == "proportion") {
      # Between-imputation SD on (0,1) scale + baseline SE via delta method
      vals <- sapply(decode_results, function(dr) as.numeric(dr$imputed[[nm]]))
      se_mc <- apply(vals, 1, stats::sd)
      # Baseline SE on logit scale -> (0,1) via delta method: se * p * (1-p)
      se_latent <- baseline_se[, lc[1]]
      se_logit <- se_latent * tm$sd
      pred_logit <- latent_mean[, lc[1]] * tm$sd + tm$mean
      pred_p <- stats::plogis(pred_logit)
      se_bm_orig <- se_logit * pred_p * (1 - pred_p)
      se_mat[, nm] <- sqrt(se_bm_orig^2 + se_mc^2)

    } else if (tm$type == "zi_count") {
      # Between-imputation SD of expected values
      vals <- sapply(decode_results, function(dr) as.numeric(dr$imputed[[nm]]))
      se_mat[, nm] <- apply(vals, 1, stats::sd)
    }
  }
  se_mat
}


# ---- Internal: single-run SE computation ----------------------------------------

compute_single_se <- function(latent_mat, probs, trait_map, baseline_se,
                              species_names) {
  # Propagate baseline (BM) SE to original scale.
  # The baseline SE captures phylogenetic imputation uncertainty.
  # This uncertainty is always present, regardless of GNN correction magnitude.
  n <- nrow(latent_mat)
  n_traits <- length(trait_map)

  se_mat <- matrix(NA_real_, nrow = n, ncol = n_traits)
  colnames(se_mat) <- vapply(trait_map, "[[", character(1), "name")
  rownames(se_mat) <- species_names

  for (tm in trait_map) {
    nm <- tm$name
    lc <- tm$latent_cols

    if (tm$type == "continuous") {
      # Baseline SE in latent (z-score) scale -> original scale
      se_latent <- baseline_se[, lc]
      se_orig <- se_latent * tm$sd
      if (isTRUE(tm$log_transform)) {
        pred_log <- latent_mat[, lc] * tm$sd + tm$mean
        se_orig <- exp(pred_log) * se_orig  # delta method
      }
      se_mat[, nm] <- se_orig

    } else if (tm$type == "count") {
      se_latent <- baseline_se[, lc]
      se_log1p <- se_latent * tm$sd
      pred_log1p <- latent_mat[, lc] * tm$sd + tm$mean
      se_mat[, nm] <- exp(pred_log1p) * se_log1p  # delta method

    } else if (tm$type == "ordinal") {
      se_latent <- baseline_se[, lc]
      se_mat[, nm] <- se_latent * tm$sd

    } else if (tm$type == "binary") {
      # Uncertainty = probability of the non-modal (less likely) class.
      # Ranges 0 (certain) to 0.5 (maximally uncertain).
      # sqrt(p*(1-p)) is the Bernoulli SD of a *single draw*, not the SE of
      # the estimated probability, so we report min(p, 1-p) instead.
      prob <- probs[[nm]]
      se_mat[, nm] <- pmin(prob, 1 - prob)

    } else if (tm$type == "categorical") {
      # Uncertainty = 1 - max class probability (margin from certainty).
      # Ranges 0 (certain) to (K-1)/K (maximally uncertain).
      # Entropy would also work but has different units; 1-max(p) is
      # interpretable as the probability of being wrong.
      prob_mat <- probs[[nm]]
      se_mat[, nm] <- 1 - apply(prob_mat, 1, max)

    } else if (tm$type == "proportion") {
      # Baseline SE on logit scale -> (0,1) via delta method
      se_latent <- baseline_se[, lc]
      se_logit <- se_latent * tm$sd
      pred_logit <- latent_mat[, lc] * tm$sd + tm$mean
      pred_p <- stats::plogis(pred_logit)
      se_mat[, nm] <- se_logit * pred_p * (1 - pred_p)

    } else if (tm$type == "zi_count") {
      # Approximate SE via Bernoulli SD of the gate * count magnitude
      p_nz <- expit(latent_mat[, lc[1]])
      # Gate uncertainty
      se_gate <- sqrt(p_nz * (1 - p_nz))
      # Magnitude SE from baseline
      se_mag_latent <- baseline_se[, lc[2]]
      pred_log1p <- latent_mat[, lc[2]] * tm$sd + tm$mean
      count_hat <- pmax(expm1(pred_log1p), 0)
      se_mag_orig <- exp(pred_log1p) * se_mag_latent * tm$sd
      # Combined SE: delta method for product p_nz * count_hat
      se_mat[, nm] <- sqrt((count_hat * se_gate)^2 + (p_nz * se_mag_orig)^2)
    }
  }
  se_mat
}


# ---- Internal: latent-scale SE for coverage computation --------------------------
#
# Returns an n x p_latent matrix of SE in latent (z-score) scale.
# For continuous/count/ordinal: baseline SE from BM.
# For binary/categorical: NA (coverage is not computed for discrete types).
# When multiple imputations exist, combines baseline SE with between-imputation
# SD in latent scale via quadrature.

compute_latent_se <- function(latent_runs, trait_map, baseline_se,
                              species_names) {
  n <- nrow(latent_runs[[1]])
  p_latent <- ncol(latent_runs[[1]])
  n_imp <- length(latent_runs)

  se_lat <- matrix(NA_real_, nrow = n, ncol = p_latent)
  rownames(se_lat) <- species_names

  for (tm in trait_map) {
    lc <- tm$latent_cols

    if (tm$type %in% c("continuous", "count", "ordinal", "proportion")) {
      for (j in lc) {
        se_bm <- baseline_se[, j]

        if (n_imp > 1L) {
          vals <- sapply(latent_runs, function(lr) lr[, j])
          se_mc <- apply(vals, 1, stats::sd)
          se_lat[, j] <- sqrt(se_bm^2 + se_mc^2)
        } else {
          se_lat[, j] <- se_bm
        }
      }
    } else if (tm$type == "zi_count") {
      # Only magnitude column (col 2) gets latent SE
      j <- lc[2]
      se_bm <- baseline_se[, j]
      if (n_imp > 1L) {
        vals <- sapply(latent_runs, function(lr) lr[, j])
        se_mc <- apply(vals, 1, stats::sd)
        se_lat[, j] <- sqrt(se_bm^2 + se_mc^2)
      } else {
        se_lat[, j] <- se_bm
      }
    }
    # binary/categorical: leave as NA (coverage not applicable)
  }
  se_lat
}


# ---- Internal: reconstruct latent column names from trait_map --------------------

latent_names_from_map <- function(trait_map) {
  nms <- character(0)
  for (tm in trait_map) {
    if (tm$type == "categorical") {
      nms <- c(nms, paste0(tm$name, "=", tm$levels))
    } else if (tm$type == "zi_count") {
      nms <- c(nms, paste0(tm$name, "_gate"), paste0(tm$name, "_mag"))
    } else {
      nms <- c(nms, tm$name)
    }
  }
  nms
}


# ---- Print methods ---------------------------------------------------------------

#' @export
print.pigauto_fit <- function(x, ...) {
  cat("pigauto_fit\n")
  cat("  Species :", length(x$species_names), "\n")
  cat("  Traits  :", length(x$trait_names),
      "--", paste(x$trait_names, collapse = ", "), "\n")

  if (!is.null(x$trait_map)) {
    types <- vapply(x$trait_map, "[[", character(1), "type")
    type_tab <- table(types)
    cat("  Types   :",
        paste(names(type_tab), type_tab, sep = "=", collapse = ", "), "\n")
  }

  cat("  Architecture: hidden_dim =", x$model_config$hidden_dim,
      "| k_eigen =", x$model_config$k_eigen, "\n")

  if (is.finite(x$val_rmse)) {
    cat("  Best val loss :", round(x$val_rmse, 4), "\n")
  }
  if (!is.na(x$test_rmse)) {
    cat("  Test loss     :", round(x$test_rmse, 4), "\n")
  }
  if (!is.null(x$calibrated_gates)) {
    cat("  Gate calibration: yes\n")
  }
  if (!is.null(x$conformal_scores)) {
    cs <- x$conformal_scores[!is.na(x$conformal_scores)]
    cat("  Conformal scores:", length(cs), "traits\n")
  }
  invisible(x)
}


#' @export
print.pigauto_pred <- function(x, ...) {
  cat("pigauto_pred\n")
  cat("  Species :", length(x$species_names), "\n")
  cat("  Traits  :", length(x$trait_names), "\n")
  cat("  Imputations:", x$n_imputations, "\n")

  if (!is.null(x$trait_map)) {
    types <- vapply(x$trait_map, "[[", character(1), "type")
    type_tab <- table(types)
    cat("  Types   :",
        paste(names(type_tab), type_tab, sep = "=", collapse = ", "), "\n")
  }

  if (!is.null(x$probabilities) && length(x$probabilities) > 0) {
    cat("  Probability traits:",
        paste(names(x$probabilities), collapse = ", "), "\n")
  }
  if (!is.null(x$conformal_lower)) {
    cat("  Conformal 95% intervals: yes\n")
  }
  invisible(x)
}
