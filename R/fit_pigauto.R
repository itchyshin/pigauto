#' Fit a Residual Phylo-DAE model for trait imputation
#'
#' Trains a \code{ResidualPhyloDAE} that blends a phylogenetic baseline with
#' a learned graph-based predictor.  Supports continuous, binary, categorical,
#' ordinal, and count traits via a unified latent space.
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
#'   layer normalisation, and residual connection.
#' @param gate_cap numeric. Upper bound for the per-column residual gate
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
#' @param lambda_shrink numeric. Weight on the residual shrinkage loss
#'   (default \code{0.03}).
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

  DELTA_tgt <- X_fill - MU            # residual target

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

  # ---- Model ----------------------------------------------------------------
  cov_dim <- p + 1L
  model <- ResidualPhyloDAE(
    input_dim      = p,
    hidden_dim     = as.integer(hidden_dim),
    coord_dim      = as.integer(k_eigen),
    cov_dim        = as.integer(cov_dim),
    per_column_rs  = TRUE,
    n_gnn_layers   = as.integer(n_gnn_layers),
    gate_cap       = gate_cap,
    use_attention  = use_attention
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
      }
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

    covs_t   <- torch::torch_cat(list(t_MU, mask_ind), dim = 2L)

    tok   <- model$mask_token$expand(c(n, p))
    X_in  <- torch::torch_where(masked_bool, tok, t_X)

    out   <- model(X_in, t_coords, covs_t, t_adj_train, t_obs_to_sp)
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
          covs0     <- torch::torch_cat(list(t_MU, mask_ind0), dim = 2L)
          # Use t_X_eval (val/test cells replaced with baseline) so
          # that held-out truth does not leak into the model input.
          out0      <- model(t_X_eval, t_coords, covs0, t_adj, t_obs_to_sp)
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

    calibrated_gates <- numeric(p)
    # 9-point grid reduces multiple-testing inflation versus the previous
    # 17-point grid while still covering the gate_cap range.
    gate_grid <- seq(0, gate_cap, length.out = 9L)

    torch::with_no_grad({
      mask_ind0 <- torch::torch_zeros(c(n, 1L), device = device)
      covs0     <- torch::torch_cat(list(t_MU, mask_ind0), dim = 2L)
      # Use t_X_eval so val-cell truths do not leak into the model input
      # during gate calibration (otherwise the GNN trivially reconstructs
      # them and every gate looks "helpful").
      out_cal   <- model(t_X_eval, t_coords, covs0, t_adj, t_obs_to_sp)
      delta_cal <- as.matrix(out_cal$delta$cpu())
      mu_cal    <- as.matrix(t_MU$cpu())
    })

    X_truth_r <- X_truth  # original data with NAs
    val_mask_mat <- matrix(FALSE, n, p)
    val_mask_mat[splits$val_idx] <- TRUE

    # --- Per-trait calibration with split-validation cross-check ----------
    #
    # The validation set is split in half.  Half A is used to pick the
    # best gate by argmin val loss, and half B is used to verify that
    # the chosen gate actually improves over the baseline.  The gate
    # is only accepted if BOTH halves agree that it helps.  This simple
    # sample-splitting guard prevents the calibrator from locking onto
    # gates that happen to look good on the full val set by chance --
    # the leading cause of val to test generalisation failures (e.g.
    # BM degrading by 2% or mixed categorical accuracy dropping 8+ pp).
    cal_min_rel_gain <- 0.02    # 2% relative improvement required

    for (tm in trait_map) {
      lc <- tm$latent_cols
      val_cells <- val_mask_mat[, lc[1]]
      if (length(lc) > 1L) {
        truth_ok <- rowSums(is.na(X_truth_r[, lc, drop = FALSE])) == 0L
      } else {
        truth_ok <- !is.na(X_truth_r[, lc[1]])
      }
      val_cells <- val_cells & truth_ok
      val_row_idx <- which(val_cells)
      n_val <- length(val_row_idx)
      if (n_val == 0) {
        calibrated_gates[lc] <- 0
        next
      }

      # Split val row indices into two halves (deterministic given seed)
      set.seed(seed + 17L)
      perm <- sample(n_val)
      half_a <- val_row_idx[perm[seq_len(floor(n_val / 2))]]
      half_b <- val_row_idx[perm[(floor(n_val / 2) + 1L):n_val]]
      if (length(half_a) == 0L || length(half_b) == 0L) {
        # Tiny val set: fall back to using the full set for both
        half_a <- val_row_idx
        half_b <- val_row_idx
      }

      # Helper: mean loss on a given set of row indices.
      # For binary/categorical we calibrate against 0-1 loss (the metric the
      # user actually sees) rather than cross-entropy.  CE and accuracy can
      # disagree on small validation sets: a gate that improves CE may
      # nudge several borderline probabilities the wrong way and degrade
      # argmax accuracy.  Calibrating against 0-1 loss makes the val signal
      # match the test metric, which is critical when n_val is small.
      cal_mean_loss <- function(g, rows) {
        if (length(rows) == 0L) return(Inf)
        if (tm$type %in% c("continuous", "count", "ordinal")) {
          pred_j <- (1 - g) * mu_cal[rows, lc[1]] +
                    g * delta_cal[rows, lc[1]]
          mean((pred_j - X_truth_r[rows, lc[1]])^2)

        } else if (tm$type == "binary") {
          pred_j <- (1 - g) * mu_cal[rows, lc[1]] +
                    g * delta_cal[rows, lc[1]]
          # 0-1 loss using probability threshold 0.5 (i.e. logit > 0).
          pred_class <- as.numeric(pred_j > 0)
          truth_j    <- X_truth_r[rows, lc[1]]
          mean(pred_class != truth_j)

        } else if (tm$type == "categorical") {
          logits <- (1 - g) * mu_cal[rows, lc, drop = FALSE] +
                    g * delta_cal[rows, lc, drop = FALSE]
          # Argmax-based 0-1 loss against the one-hot truth.
          pred_class  <- max.col(logits, ties.method = "first")
          truth_mat   <- X_truth_r[rows, lc, drop = FALSE]
          truth_class <- max.col(truth_mat, ties.method = "first")
          mean(pred_class != truth_class)
        } else {
          Inf
        }
      }

      # For discrete traits we also require an absolute minimum cell-level
      # improvement, because relative gains on 0-1 loss are deceptively
      # small (a 2% relative gain over baseline 0.35 is only 0.007 ~ 0.2
      # cells out of 30, easily noise).
      is_discrete <- tm$type %in% c("binary", "categorical")
      min_abs_a <- if (is_discrete) 2 / max(length(half_a), 1L) else 0
      min_abs_b <- if (is_discrete) 1 / max(length(half_b), 1L) else 0

      # 1. On half A: find the gate that minimises val loss, subject to
      #    the relative-gain floor and (for discrete) the absolute floor.
      loss_a_0  <- cal_mean_loss(0, half_a)
      best_g    <- 0
      best_la   <- loss_a_0
      for (g in gate_grid) {
        if (g == 0) next
        loss_a_g <- cal_mean_loss(g, half_a)
        if (!is.finite(loss_a_g)) next
        rel <- (loss_a_0 - loss_a_g) / max(loss_a_0, 1e-12)
        abs_gain <- loss_a_0 - loss_a_g
        if (rel >= cal_min_rel_gain && abs_gain >= min_abs_a &&
            loss_a_g < best_la) {
          best_g  <- g
          best_la <- loss_a_g
        }
      }

      # 2. On half B: verify that the chosen gate actually helps.  If it
      #    does not, fall back to gate = 0 (baseline).
      if (best_g > 0) {
        loss_b_0 <- cal_mean_loss(0,     half_b)
        loss_b_g <- cal_mean_loss(best_g, half_b)
        rel_b    <- (loss_b_0 - loss_b_g) / max(loss_b_0, 1e-12)
        abs_b    <- loss_b_0 - loss_b_g
        # The verification threshold is half the calibration threshold:
        # we only require the improvement to persist, not to re-pass the
        # strict bar on the second half.  Discrete traits also need the
        # absolute cell floor.
        if (!is.finite(rel_b) ||
            rel_b < (cal_min_rel_gain / 2) ||
            abs_b < min_abs_b) {
          best_g <- 0
        }
      }

      calibrated_gates[lc] <- best_g
    }

    if (verbose) {
      gate_summary <- round(calibrated_gates, 3)
      names(gate_summary) <- if (!is.null(data$latent_names)) data$latent_names else paste0("col", seq_len(p))
      message("Calibrated gates: ", paste(names(gate_summary), gate_summary, sep = "=", collapse = ", "))
    }
  }

  # ---- Conformal prediction scores (validation set) -----------------------
  conformal_scores <- NULL
  if (has_val && has_trait_map) {
    if (verbose) message("Computing conformal prediction scores...")

    # Use calibrated gates if available, otherwise learned gates
    if (!is.null(calibrated_gates)) {
      gates_to_use <- calibrated_gates
    } else {
      gates_to_use <- as.numeric(out_cal$rs$cpu()$squeeze())
    }

    # Compute calibrated predictions on validation set
    pred_cal <- matrix(0, n, p)
    for (j in seq_len(p)) {
      pred_cal[, j] <- (1 - gates_to_use[j]) * mu_cal[, j] +
                        gates_to_use[j] * delta_cal[, j]
    }

    conformal_scores <- rep(NA_real_, length(trait_map))
    names(conformal_scores) <- vapply(trait_map, "[[", character(1), "name")

    for (tm in trait_map) {
      if (!(tm$type %in% c("continuous", "count", "ordinal"))) next
      lc <- tm$latent_cols
      val_cells <- val_mask_mat[, lc[1]]
      if (sum(val_cells) == 0) next

      residuals <- abs(X_truth_r[val_cells, lc[1]] - pred_cal[val_cells, lc[1]])
      residuals <- residuals[is.finite(residuals)]
      n_val <- length(residuals)
      if (n_val == 0L) next
      # Conformal quantile: ceil((1-alpha)(1+1/n_val))-th value
      alpha <- 0.05
      q_level <- min(ceiling((1 - alpha) * (n_val + 1)) / n_val, 1)
      conformal_scores[tm$name] <- as.numeric(
        stats::quantile(residuals, q_level, na.rm = TRUE)
      )
    }

    if (verbose) {
      cs_print <- round(conformal_scores[!is.na(conformal_scores)], 4)
      message("Conformal scores (latent scale): ",
              paste(names(cs_print), cs_print, sep = "=", collapse = ", "))
    }
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
      covs0     <- torch::torch_cat(list(t_MU, mask_ind0), dim = 2L)
      # Use t_X_eval so held-out test cells are not leaked to the model.
      out0      <- model(t_X_eval, t_coords, covs0, t_adj, t_obs_to_sp)
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
    hidden_dim      = hidden_dim,
    k_eigen         = k_eigen,
    n_gnn_layers    = n_gnn_layers,
    gate_cap        = gate_cap,
    use_attention   = use_attention,
    dropout         = dropout,
    refine_steps    = refine_steps,
    cov_dim         = cov_dim,
    input_dim       = p,
    per_column_rs   = TRUE
  )

  # Backward-compat: store val_rmse and test_rmse names
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
      conformal_scores = conformal_scores
    ),
    class = "pigauto_fit"
  )
}
