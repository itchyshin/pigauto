#' Fit a pigauto model for trait imputation
#'
#' Trains a pigauto model: a gated ensemble of a phylogenetic baseline
#' and an attention-based graph neural network correction, implemented
#' as an internal torch module (\code{ResidualPhyloDAE}; "Residual" here
#' refers to the ResNet-style skip connections in the GNN layers, not
#' to a statistical residual). For continuous, count, and ordinal traits
#' the baseline is Brownian motion (phylogenetic correlation matrix); for binary
#' and categorical traits it is phylogenetic label propagation. Supports
#' all five trait types via a unified latent space.
#'
#' @details
#' **Blend formulation:**
#' The prediction is \eqn{\hat{x} = (1-r)\mu + r\delta}, where \eqn{\mu} is
#' the BM baseline, \eqn{\delta} is the model's direct prediction, and
#' \eqn{r = \sigma(\rho) \times \mathrm{cap}} is a per-column learnable gate
#' bounded in \eqn{(0, \mathrm{gate\_cap})}.  When \eqn{r = 0}, the prediction
#' collapses to the baseline.
#' The gate is regularised toward zero via the shrinkage penalty on
#' \eqn{\delta - \mu}, so the model defaults to the baseline unless the
#' GNN's correction demonstrably helps on the validation set.
#'
#' **Training objective** (per epoch):
#' \enumerate{
#'   \item A random subset of observed cells is corrupted with a learnable
#'     mask token.
#'   \item The model predicts \eqn{\delta} from graph context.
#'   \item Loss = type-specific reconstruction on corrupted cells +
#'     \code{lambda_shrink} * MSE(\eqn{\delta - \mu}) +
#'     \code{lambda_gate} * MSE(\eqn{r}).
#' }
#'
#' The gate penalty on \eqn{r} is necessary because when \eqn{\delta =
#' \mu} (the BM-optimal solution for observed cells), the reconstruction
#' and shrinkage losses both equal zero regardless of \eqn{r}, leaving no
#' gradient to close the gate.  The explicit penalty ensures gates default
#' toward zero when the GNN correction provides no benefit.
#'
#' **Type-specific losses:**
#' \describe{
#'   \item{continuous/count/ordinal}{MSE}
#'   \item{binary}{BCE with logits}
#'   \item{categorical}{cross-entropy over K latent columns}
#' }
#'
#' @param data object of class \code{"pigauto_data"}.
#' @param tree object of class \code{"phylo"}.
#' @param splits list (output of \code{\link{make_missing_splits}}) or
#'   \code{NULL}.
#' @param graph list (output of \code{\link{build_phylo_graph}}) or
#'   \code{NULL}.
#' @param baseline list (output of \code{\link{fit_baseline}}) or \code{NULL}.
#' @param hidden_dim integer. Hidden layer width (default \code{64}).
#' @param k_eigen integer. Number of spectral node features (default
#'   \code{8}).
#' @param dropout numeric. Dropout rate (default \code{0.10}).
#' @param lr numeric. AdamW learning rate (default \code{0.003}).
#' @param weight_decay numeric. AdamW weight decay (default \code{1e-4}).
#' @param epochs integer. Maximum training epochs (default \code{3000}).
#' @param n_gnn_layers integer. Number of graph message-passing layers
#'   (default \code{2}).  Each layer has its own learnable alpha gate,
#'   layer normalisation, and ResNet-style skip connection.
#' @param gate_cap numeric. Upper bound for the per-column blend gate
#'   (default \code{0.8}).  Safety comes from regularisation, not the cap.
#' @param use_attention logical. Use attention in the GNN layers (default
#'   \code{TRUE}).
#' @param corruption_rate numeric. Final corruption fraction if
#'   \code{corruption_ramp > 0}; otherwise the fixed corruption rate per
#'   epoch (default \code{0.55}).
#' @param corruption_start numeric. Initial corruption fraction for the
#'   curriculum schedule (default \code{0.20}).  Ignored if
#'   \code{corruption_ramp = 0}.
#' @param corruption_ramp integer. Epochs over which corruption linearly
#'   ramps from \code{corruption_start} to \code{corruption_rate}
#'   (default \code{500}).  Set to \code{0} for fixed corruption.
#' @param refine_steps integer. Iterative refinement steps at inference
#'   (default \code{8}).
#' @param lambda_shrink numeric. Weight on the shrinkage penalty
#'   \code{||delta - baseline||^2} that keeps the GNN correction close
#'   to the phylogenetic baseline (default \code{0.03}).
#' @param lambda_gate numeric. Weight on the gate regularisation penalty
#'   that pushes learnable gates toward zero.  Prevents gates from staying
#'   open when the GNN provides no useful correction (default \code{0.01}).
#' @param warmup_epochs integer. Linear learning-rate warmup over the first
#'   N epochs (default \code{200}).  After warmup, a cosine schedule decays
#'   the LR to \code{1e-5}.
#' @param edge_dropout numeric. Fraction of adjacency edges randomly zeroed
#'   each training epoch for graph regularisation (default \code{0.1}).
#'   Set to \code{0} to disable.
#' @param eval_every integer. Evaluate on val every N epochs (default
#'   \code{100}).
#' @param patience integer. Early-stopping patience in eval cycles (default
#'   \code{10}).
#' @param clip_norm numeric. Gradient clip norm (default \code{1.0}).
#' @param verbose logical. Print training progress (default \code{TRUE}).
#' @param seed integer. Random seed (default \code{1}).
#' @return An object of class \code{"pigauto_fit"}.
#' @importFrom torch torch_tensor torch_float torch_bool torch_long
#' @importFrom torch optim_adam with_no_grad nn_utils_clip_grad_norm_
#' @importFrom torch nnf_mse_loss torch_rand_like torch_where torch_zeros_like
#' @export
fit_pigauto <- function(
    data,
    tree,
    splits            = NULL,
    graph             = NULL,
    baseline          = NULL,
    hidden_dim        = 64L,
    k_eigen           = "auto",
    n_gnn_layers      = 2L,
    gate_cap          = 0.8,
    use_attention     = TRUE,
    use_transformer_blocks = TRUE,
    n_heads           = 4L,
    ffn_mult          = 4L,
    dropout           = 0.10,
    lr                = 0.003,
    weight_decay      = 1e-4,
    epochs            = 3000L,
    corruption_rate   = 0.55,
    corruption_start  = 0.20,
    corruption_ramp   = 500L,
    refine_steps      = 8L,
    lambda_shrink     = 0.03,
    lambda_gate       = 0.01,
    warmup_epochs     = 200L,
    edge_dropout      = 0.1,
    eval_every        = 100L,
    patience          = 10L,
    clip_norm         = 1.0,
    verbose           = TRUE,
    seed              = 1L
) {
  if (!inherits(data, "pigauto_data")) {
    stop("'data' must be a pigauto_data object.")
  }

  set.seed(seed)
  torch::torch_manual_seed(seed)

  device <- get_device()
  if (verbose) message("Using device: ", as.character(device))

  # ---- Graph ----------------------------------------------------------------
  if (is.null(graph)) {
    if (verbose) message("Computing phylogenetic graph...")
    graph <- build_phylo_graph(tree, k_eigen = k_eigen)
  }
  # Use the actual k_eigen from the graph (may differ from input if "auto")
  k_eigen <- ncol(graph$coords)

  # ---- Baseline -------------------------------------------------------------
  if (is.null(baseline)) {
    if (verbose) message("Fitting baseline...")
    # Pass graph through so fit_baseline can reuse graph$D instead of
    # calling ape::cophenetic.phylo() a second time on the same tree.
    baseline <- fit_baseline(data, tree, splits = splits, graph = graph)
  }

  # ---- Trait map ------------------------------------------------------------
  trait_map <- data$trait_map
  has_trait_map <- !is.null(trait_map)

  # ---- Multi-obs handling --------------------------------------------------
  multi_obs  <- isTRUE(data$multi_obs)
  n_obs      <- nrow(data$X_scaled)
  n_species  <- data$n_species %||% n_obs
  obs_to_sp  <- data$obs_to_species  # NULL when single-obs

  # ---- Data preparation ----------------------------------------------------
  n <- n_obs                           # n = n_obs (rows in X)
  p <- ncol(data$X_scaled)

  X_truth <- data$X_scaled             # n_obs x p_latent

  # Baseline mu is at species level (n_species x p).
  # For multi-obs, expand to observation level.
  MU_species <- baseline$mu            # n_species x p_latent
  if (multi_obs) {
    MU <- MU_species[obs_to_sp, , drop = FALSE]
    rownames(MU) <- NULL
  } else {
    MU <- MU_species
  }

  # Fill missing cells with baseline predictions
  X_fill <- X_truth
  X_fill[is.na(X_fill)] <- MU[is.na(X_fill)]

  DELTA_tgt <- X_fill - MU            # observed - baseline (legacy
                                       # tensor t_DELTA below; kept for
                                       # backward compatibility, not
                                       # used as a training target)

  # Observed mask (TRUE = actually observed, not just filled)
  M_obs_mat <- !is.na(X_truth)
  if (!is.null(splits)) M_obs_mat[c(splits$val_idx, splits$test_idx)] <- FALSE

  # Held-out-masked input: like X_fill but with val/test cells replaced
  # with the baseline prediction.  This is what the model SHOULD see at
  # val/test evaluation time and at gate calibration time, so that it
  # does not have direct access to the held-out truth values.  Without
  # this masking, val cells' true values leak into the model input and
  # the calibration signal becomes meaningless.
  X_fill_heldout <- X_fill
  if (!is.null(splits)) {
    hold_idx <- c(splits$val_idx, splits$test_idx)
    X_fill_heldout[hold_idx] <- MU[hold_idx]
  }

  # ---- Trait-level corruption mask expansion --------------------------------
  # For categorical traits: corrupt all K latent columns together
  if (has_trait_map) {
    # Build mapping from latent column -> trait index (for grouped corruption)
    latent_to_trait <- integer(p)
    for (i in seq_along(trait_map)) {
      latent_to_trait[trait_map[[i]]$latent_cols] <- i
    }
    n_orig_traits <- length(trait_map)
  }

  # ---- Tensors ---------------------------------------------------------------
  t_X      <- torch::torch_tensor(X_fill,    dtype = torch::torch_float(),
                                  device = device)
  t_X_eval <- torch::torch_tensor(X_fill_heldout, dtype = torch::torch_float(),
                                  device = device)
  t_MU     <- torch::torch_tensor(MU,        dtype = torch::torch_float(),
                                  device = device)
  t_DELTA  <- torch::torch_tensor(DELTA_tgt, dtype = torch::torch_float(),
                                  device = device)
  t_M_obs  <- torch::torch_tensor(M_obs_mat, dtype = torch::torch_bool(),
                                  device = device)
  t_adj    <- torch::torch_tensor(graph$adj,    dtype = torch::torch_float(),
                                  device = device)
  t_coords <- torch::torch_tensor(graph$coords, dtype = torch::torch_float(),
                                  device = device)

  # Squared cophenetic distances for B2 rate-aware attention.
  # graph$D_sq is set by build_phylo_graph() and is NOT freed in impute.R's
  # cleanup block (only D and R_phy are freed there), so it is available here.
  # Falls back to NULL if the graph was built before B2.1 (old cache).
  t_D_sq <- if (!is.null(graph$D_sq)) {
    torch::torch_tensor(graph$D_sq, dtype = torch::torch_float(),
                        device = device)
  } else {
    NULL
  }

  # Drop the cophenetic distance matrix from fit_pigauto's local view of
  # the graph now that we have the torch tensors we need. This keeps the
  # fit object slim on disk (predict() only reads $adj and $coords), but
  # does NOT free memory at the caller's scope -- R's copy-on-modify
  # means the caller still holds the original list. Callers that want
  # to free memory during training should set graph$D <- NULL before
  # calling fit_pigauto(); impute() already does this, and the
  # bundled benchmark scripts do the same.
  graph$D <- NULL

  # Observation-to-species mapping for multi-obs (1-indexed for torch-R)
  if (multi_obs) {
    t_obs_to_sp <- torch::torch_tensor(
      as.integer(obs_to_sp),
      dtype = torch::torch_long(), device = device
    )
  } else {
    t_obs_to_sp <- NULL
  }

  # Val / test masks
  if (!is.null(splits)) {
    val_mat <- matrix(FALSE, n, p); val_mat[splits$val_idx]  <- TRUE
    t_val   <- torch::torch_tensor(val_mat, dtype = torch::torch_bool(),
                                   device = device)
    t_truth <- torch::torch_tensor(X_truth, dtype = torch::torch_float(),
                                   device = device)
    # Replace NaN in truth with 0 (won't affect masked evaluation)
    t_truth_safe <- t_truth$clone()
    t_truth_safe[t_truth$isnan()] <- 0
    has_val <- TRUE
  } else {
    has_val <- FALSE
  }

  # ---- Covariates (environmental conditioners) ------------------------------
  has_covariates <- !is.null(data$covariates)
  n_cov_cols     <- if (has_covariates) ncol(data$covariates) else 0L
  t_covariates   <- NULL
  if (has_covariates) {
    t_covariates <- torch::torch_tensor(
      data$covariates, dtype = torch::torch_float(), device = device
    )
  }

  # ---- Model ----------------------------------------------------------------
  cov_dim <- p + 1L + n_cov_cols
  # n_user_cov tells the model how many observation-level user covariates
  # exist (the last n_cov_cols columns of the covs tensor).  When multi_obs
  # is TRUE and n_user_cov > 0, the model's obs_refine MLP re-injects
  # these covariates after species-level GNN message passing so that
  # different observations of the same species get distinct predictions.
  n_user_cov <- if (multi_obs) n_cov_cols else 0L
  model <- ResidualPhyloDAE(
    input_dim              = p,
    hidden_dim             = as.integer(hidden_dim),
    coord_dim              = as.integer(k_eigen),
    cov_dim                = as.integer(cov_dim),
    per_column_rs          = TRUE,
    n_gnn_layers           = as.integer(n_gnn_layers),
    gate_cap               = gate_cap,
    use_attention          = use_attention,
    n_user_cov             = as.integer(n_user_cov),
    dropout                = dropout,
    use_transformer_blocks = use_transformer_blocks,
    n_heads                = as.integer(n_heads),
    ffn_mult               = as.integer(ffn_mult)
  )

  # Type-aware gate init: set res_raw so effective gate starts at a
  # non-trivial value for all trait types.  Discrete traits previously
  # started near 0 (disc_target = 0.001), which killed gradient flow
  # through the sigmoid and prevented the GNN from learning discrete
  # corrections.  Starting at 0.10 gives a healthy gradient in both
  # directions.  logit(target / gate_cap) gives the correct res_raw.
  cont_target <- 0.135
  disc_target <- 0.10
  cont_raw <- log(cont_target / gate_cap / (1 - cont_target / gate_cap))
  disc_raw <- log(disc_target / gate_cap / (1 - disc_target / gate_cap))

  if (has_trait_map) {
    init_vals <- rep(cont_raw, p)
    for (tm in trait_map) {
      if (tm$type %in% c("binary", "categorical", "ordinal")) {
        init_vals[tm$latent_cols] <- disc_raw
      } else if (tm$type == "zi_count") {
        # Gate column (col 1): like binary (start closed)
        init_vals[tm$latent_cols[1]] <- disc_raw
        # Magnitude column (col 2): like count (start slightly open)
        init_vals[tm$latent_cols[2]] <- cont_raw
      }
      # proportion: uses cont_raw (default), no special case needed
    }
  } else {
    init_vals <- rep(cont_raw, p)
  }
  torch::with_no_grad({
    model$res_raw$copy_(
      torch::torch_tensor(init_vals, dtype = torch::torch_float())
    )
  })

  model$to(device = device)
  opt <- torch::optim_adamw(model$parameters, lr = lr,
                             weight_decay = weight_decay)

  # ---- Training loop --------------------------------------------------------
  best_val      <- Inf
  best_state    <- NULL
  patience_left <- patience
  history       <- data.frame(
    epoch = integer(0), loss_rec = double(0),
    loss_shrink = double(0), loss_gate = double(0),
    val_loss = double(0), lr = double(0)
  )

  for (epoch in seq_len(epochs)) {
    model$train()
    opt$zero_grad()

    # ---- Learning rate schedule (warmup + cosine decay) --------------------
    scheduled_lr <- cosine_lr(epoch, warmup_epochs, epochs, lr)
    opt$param_groups[[1]]$lr <- scheduled_lr

    # ---- Adaptive corruption rate ------------------------------------------
    if (corruption_ramp > 0L) {
      eff_corruption <- corruption_start +
        (corruption_rate - corruption_start) *
        min(epoch / corruption_ramp, 1.0)
    } else {
      eff_corruption <- corruption_rate
    }

    # ---- Edge dropout (graph regularisation) --------------------------------
    if (edge_dropout > 0 && model$training) {
      edge_mask <- torch::torch_rand_like(t_adj) > edge_dropout
      t_adj_train <- t_adj * edge_mask$to(dtype = torch::torch_float())
      # Renormalise rows to maintain expected message scale
      row_sums <- t_adj_train$sum(dim = 2L, keepdim = TRUE)$clamp(min = 1e-8)
      orig_sums <- t_adj$sum(dim = 2L, keepdim = TRUE)$clamp(min = 1e-8)
      t_adj_train <- t_adj_train * (orig_sums / row_sums)
    } else {
      t_adj_train <- t_adj
    }

    # ---- Corruption masking -------------------------------------------------
    if (has_trait_map) {
      # Corrupt at trait level, then expand to latent columns
      u_trait <- matrix(stats::runif(n * n_orig_traits), n, n_orig_traits)
      corrupt_trait <- u_trait < eff_corruption
      # Only corrupt where observed
      obs_trait <- matrix(FALSE, n, n_orig_traits)
      for (i in seq_along(trait_map)) {
        lc <- trait_map[[i]]$latent_cols
        obs_trait[, i] <- as.logical(M_obs_mat[, lc[1]])
      }
      corrupt_trait <- corrupt_trait & obs_trait
      # Expand to latent columns
      corrupt_latent <- matrix(FALSE, n, p)
      for (i in seq_along(trait_map)) {
        for (lc in trait_map[[i]]$latent_cols) {
          corrupt_latent[, lc] <- corrupt_trait[, i]
        }
      }
      masked_bool <- torch::torch_tensor(corrupt_latent, dtype = torch::torch_bool(),
                                          device = device)
    } else {
      u           <- torch::torch_rand_like(t_X)
      masked_bool <- (u < eff_corruption) & t_M_obs
    }

    if (as.numeric(masked_bool$sum()$cpu()$item()) == 0L) next

    Mask_t  <- masked_bool$to(dtype = torch::torch_float())
    mask_ind <- Mask_t$any(dim = 2L, keepdim = TRUE)$to(
      dtype = torch::torch_float()
    )

    covs_t <- make_covs_tensor(t_MU, mask_ind, t_covariates)

    tok   <- model$mask_token$expand(c(n, p))
    X_in  <- torch::torch_where(masked_bool, tok, t_X)

    out   <- model(X_in, t_coords, covs_t, t_adj_train, t_obs_to_sp,
                   D_sq = t_D_sq)
    delta <- out$delta
    rs    <- out$rs

    # Direct prediction loss: train delta to predict the truth directly.
    # This gives delta a full-strength gradient signal on every corrupted
    # cell, unlike the earlier blend formulation where the effective
    # gradient on delta was scaled by rs (approx 0.1 to 0.2), causing
    # delta to stay close to the baseline and calibration to fall back
    # to gate = 0.
    if (has_trait_map) {
      loss_rec <- compute_mixed_loss(delta, t_X, masked_bool, trait_map)
    } else {
      loss_rec <- torch::nnf_mse_loss(
        delta[masked_bool], t_X[masked_bool]
      )
    }

    # Auxiliary blend loss: (1-rs)*mu + rs*delta against the truth.
    # This gives rs a gradient signal so that the model can learn which
    # columns benefit from the GNN, but with a small weight so it doesn't
    # dominate the direct supervision on delta.  Post-hoc calibration on
    # the validation set still has the final say over the gate.
    pred_blend <- (1 - rs) * t_MU + rs * delta
    if (has_trait_map) {
      loss_blend <- compute_mixed_loss(pred_blend, t_X, masked_bool, trait_map)
    } else {
      loss_blend <- torch::nnf_mse_loss(
        pred_blend[masked_bool], t_X[masked_bool]
      )
    }

    # Shrinkage: keep delta from drifting to extreme values on cells with
    # very weak supervision.  Lighter than before because direct
    # supervision handles most of the regularisation.
    deviation <- (delta - t_MU)[masked_bool]
    loss_shrink <- torch::torch_mean(deviation$pow(2))

    # Gate regularisation pushes rs toward zero by default, so the model
    # only opens a gate when the direct delta signal is clearly helpful.
    loss_gate <- torch::torch_mean(rs$pow(2))

    loss <- loss_rec + 0.1 * loss_blend +
            lambda_shrink * loss_shrink + lambda_gate * loss_gate
    loss$backward()

    if (!all_grads_finite(model$parameters)) {
      opt$zero_grad(); next
    }
    torch::nn_utils_clip_grad_norm_(model$parameters, max_norm = clip_norm)
    opt$step()

    # ---- Evaluation ---------------------------------------------------------
    if (epoch %% eval_every == 0L) {
      val_loss_val <- NA_real_
      if (has_val) {
        model$eval()
        torch::with_no_grad({
          mask_ind0 <- torch::torch_zeros(c(n, 1L), device = device)
          covs0     <- make_covs_tensor(t_MU, mask_ind0, t_covariates)
          # Use t_X_eval (val/test cells replaced with baseline) so
          # that held-out truth does not leak into the model input.
          out0      <- model(t_X_eval, t_coords, covs0, t_adj, t_obs_to_sp,
                             D_sq = t_D_sq)
          pred0     <- (1 - out0$rs) * t_MU + out0$rs * out0$delta
        })
        if (has_trait_map) {
          val_loss_val <- composite_val_loss(pred0, t_truth_safe, t_val,
                                             trait_map)
        } else {
          val_loss_val <- as.numeric(
            rmse_torch(pred0, t_truth_safe, t_val)$cpu()$item()
          )
        }
      }

      history <- rbind(history, data.frame(
        epoch       = epoch,
        loss_rec    = as.numeric(loss_rec$item()),
        loss_shrink = as.numeric(loss_shrink$item()),
        loss_gate   = as.numeric(loss_gate$item()),
        val_loss    = if (is.na(val_loss_val)) NA_real_ else val_loss_val,
        lr          = scheduled_lr
      ))

      if (verbose) {
        msg <- sprintf(
          "epoch %4d | rec %.4f | shrink %.4f | gate %.4f | lr %.1e",
          epoch, as.numeric(loss_rec$item()),
          as.numeric(loss_shrink$item()),
          as.numeric(loss_gate$item()),
          scheduled_lr
        )
        if (!is.na(val_loss_val)) {
          msg <- paste0(msg, sprintf(" | val %.4f | best %.4f | pat %d",
                                     val_loss_val, best_val, patience_left))
        }
        message(msg)
      }

      if (!is.na(val_loss_val)) {
        if (val_loss_val + 1e-6 < best_val) {
          best_val      <- val_loss_val
          best_state    <- model$state_dict()
          patience_left <- patience
        } else {
          patience_left <- patience_left - 1L
        }
        if (patience_left <= 0L) {
          if (verbose) message("Early stopping at epoch ", epoch)
          break
        }
      }
    }
  }

  if (!is.null(best_state)) model$load_state_dict(best_state)

  # ---- Post-training gate calibration (validation set) --------------------
  calibrated_gates <- NULL
  if (has_val && has_trait_map) {
    if (verbose) message("Calibrating gates on validation set...")
    model$eval()

    # Forward pass on val set (use t_X_eval so held-out truths do not leak).
    torch::with_no_grad({
      mask_ind0 <- torch::torch_zeros(c(n, 1L), device = device)
      covs0     <- make_covs_tensor(t_MU, mask_ind0, t_covariates)
      out_cal   <- model(t_X_eval, t_coords, covs0, t_adj, t_obs_to_sp,
                         D_sq = t_D_sq)
      delta_cal <- as.matrix(out_cal$delta$cpu())
      mu_cal    <- as.matrix(t_MU$cpu())
    })

    X_truth_r    <- X_truth
    val_mask_mat <- matrix(FALSE, n, p)
    val_mask_mat[splits$val_idx] <- TRUE

    calibrated_gates <- calibrate_gates(
      trait_map       = trait_map,
      mu_cal          = mu_cal,
      delta_cal       = delta_cal,
      X_truth_r       = X_truth_r,
      val_mask_mat    = val_mask_mat,
      gate_grid       = seq(0, gate_cap, length.out = 9L),
      gate_cap        = gate_cap,
      seed            = seed,
      latent_names    = data$latent_names,
      verbose         = verbose
    )
  }

  # ---- Conformal prediction scores (validation set) -----------------------
  conformal_scores <- NULL
  if (has_val && has_trait_map && !is.null(calibrated_gates)) {
    conformal_scores <- compute_conformal_scores(
      trait_map        = trait_map,
      calibrated_gates = calibrated_gates,
      mu_cal           = mu_cal,
      delta_cal        = delta_cal,
      X_truth_r        = X_truth_r,
      val_mask_mat     = val_mask_mat,
      verbose          = verbose
    )
  }

  # ---- Final test loss -------------------------------------------------------
  test_loss <- NA_real_
  if (!is.null(splits) && length(splits$test_idx) > 0L) {
    test_mat <- matrix(FALSE, n, p); test_mat[splits$test_idx] <- TRUE
    t_test   <- torch::torch_tensor(test_mat, dtype = torch::torch_bool(),
                                    device = device)
    model$eval()
    torch::with_no_grad({
      mask_ind0 <- torch::torch_zeros(c(n, 1L), device = device)
      covs0     <- make_covs_tensor(t_MU, mask_ind0, t_covariates)
      # Use t_X_eval so held-out test cells are not leaked to the model.
      out0      <- model(t_X_eval, t_coords, covs0, t_adj, t_obs_to_sp,
                         D_sq = t_D_sq)
      pred0     <- (1 - out0$rs) * t_MU + out0$rs * out0$delta
      t_truth_f <- torch::torch_tensor(X_truth, dtype = torch::torch_float(),
                                       device = device)
      t_truth_f[t_truth_f$isnan()] <- 0
    })
    if (has_trait_map) {
      test_loss <- composite_val_loss(pred0, t_truth_f, t_test, trait_map)
    } else {
      test_loss <- as.numeric(
        rmse_torch(pred0, t_truth_f, t_test)$cpu()$item()
      )
    }
  }

  # ---- Build result object ---------------------------------------------------
  model_config <- list(
    hidden_dim             = hidden_dim,
    k_eigen                = k_eigen,
    n_gnn_layers           = n_gnn_layers,
    gate_cap               = gate_cap,
    use_attention          = use_attention,
    use_transformer_blocks = use_transformer_blocks,
    n_heads                = as.integer(n_heads),
    ffn_mult               = as.integer(ffn_mult),
    dropout                = dropout,
    refine_steps           = refine_steps,
    cov_dim                = cov_dim,
    input_dim              = p,
    per_column_rs          = TRUE,
    n_user_cov             = n_user_cov
  )

  # Backward-compat: store val_rmse and test_rmse names
  # (graph$D was stripped earlier, right after tensor creation, so the
  # graph stored here is already the slim version without the cophenetic
  # distance matrix.)
  structure(
    list(
      model_state    = model$state_dict(),
      model_config   = model_config,
      graph          = graph,
      baseline       = baseline,
      norm           = list(
        means         = data$means,
        sds           = data$sds,
        log_transform = data$log_transform
      ),
      species_names  = data$species_names,
      obs_species    = data$obs_species,
      obs_to_species = data$obs_to_species,
      n_species      = n_species,
      n_obs          = n_obs,
      multi_obs      = multi_obs,
      trait_names    = data$trait_names,
      latent_names   = data$latent_names,
      trait_map      = trait_map,
      splits         = splits,
      history        = history,
      val_rmse         = best_val,
      test_rmse        = test_loss,
      calibrated_gates = calibrated_gates,
      conformal_scores = conformal_scores,
      covariates       = data$covariates,
      cov_means        = data$cov_means,
      cov_sds          = data$cov_sds,
      cov_names        = data$cov_names
    ),
    class = "pigauto_fit"
  )
}
