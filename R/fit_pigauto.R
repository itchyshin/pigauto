#' Fit a Residual Phylo-DAE model for trait imputation
#'
#' Trains a \code{ResidualPhyloDAE} that learns a bounded residual correction
#' on top of a baseline.  Supports continuous, binary, categorical, ordinal,
#' and count traits via a unified latent space.
#'
#' @details
#' **Training objective** (per epoch):
#' \enumerate{
#'   \item A random subset of observed cells is corrupted with a learnable
#'     mask token.
#'   \item The model predicts the residual \eqn{\delta = X - \mu_{baseline}}.
#'   \item Loss = type-specific reconstruction on corrupted cells +
#'     \code{lambda_shrink} * MSE toward zero.
#' }
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
#' @param corruption_rate numeric. Fraction of observed cells corrupted per
#'   epoch (default \code{0.55}).
#' @param refine_steps integer. Iterative refinement steps at inference
#'   (default \code{8}).
#' @param lambda_shrink numeric. Weight on the residual shrinkage loss
#'   (default \code{0.03}).
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
    splits          = NULL,
    graph           = NULL,
    baseline        = NULL,
    hidden_dim      = 64L,
    k_eigen         = 8L,
    dropout         = 0.10,
    lr              = 0.003,
    weight_decay    = 1e-4,
    epochs          = 3000L,
    corruption_rate = 0.55,
    refine_steps    = 8L,
    lambda_shrink   = 0.03,
    eval_every      = 100L,
    patience        = 10L,
    clip_norm       = 1.0,
    verbose         = TRUE,
    seed            = 1L
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

  # ---- Baseline -------------------------------------------------------------
  if (is.null(baseline)) {
    if (verbose) message("Fitting baseline...")
    baseline <- fit_baseline(data, tree, splits = splits)
  }

  # ---- Trait map ------------------------------------------------------------
  trait_map <- data$trait_map
  has_trait_map <- !is.null(trait_map)

  # ---- Data preparation ----------------------------------------------------
  n <- nrow(data$X_scaled)
  p <- ncol(data$X_scaled)

  X_truth <- data$X_scaled          # n x p_latent
  MU      <- baseline$mu            # n x p_latent

  # Fill missing cells with baseline predictions
  X_fill <- X_truth
  X_fill[is.na(X_fill)] <- MU[is.na(X_fill)]

  DELTA_tgt <- X_fill - MU          # residual target

  # Observed mask (TRUE = actually observed, not just filled)
  M_obs_mat <- !is.na(X_truth)
  if (!is.null(splits)) M_obs_mat[c(splits$val_idx, splits$test_idx)] <- FALSE

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
    input_dim  = p,
    hidden_dim = as.integer(hidden_dim),
    coord_dim  = as.integer(k_eigen),
    cov_dim    = as.integer(cov_dim)
  )
  model$to(device = device)
  opt <- torch::optim_adamw(model$parameters, lr = lr,
                             weight_decay = weight_decay)

  # ---- Training loop --------------------------------------------------------
  best_val      <- Inf
  best_state    <- NULL
  patience_left <- patience
  history       <- data.frame(
    epoch = integer(0), loss_rec = double(0),
    loss_shrink = double(0), val_loss = double(0)
  )

  for (epoch in seq_len(epochs)) {
    model$train()
    opt$zero_grad()

    # ---- Corruption masking -------------------------------------------------
    if (has_trait_map) {
      # Corrupt at trait level, then expand to latent columns
      u_trait <- matrix(stats::runif(n * n_orig_traits), n, n_orig_traits)
      corrupt_trait <- u_trait < corruption_rate
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
      masked_bool <- (u < corruption_rate) & t_M_obs
    }

    if (as.numeric(masked_bool$sum()$cpu()$item()) == 0L) next

    Mask_t  <- masked_bool$to(dtype = torch::torch_float())
    mask_ind <- Mask_t$any(dim = 2L, keepdim = TRUE)$to(
      dtype = torch::torch_float()
    )
    covs_t   <- torch::torch_cat(list(t_MU, mask_ind), dim = 2L)

    tok   <- model$mask_token$expand(c(n, p))
    X_in  <- torch::torch_where(masked_bool, tok, t_X)

    out   <- model(X_in, t_coords, covs_t, t_adj)
    delta <- out$delta
    rs    <- out$rs

    pred_scaled <- t_MU + rs * delta

    # ---- Loss ---------------------------------------------------------------
    if (has_trait_map) {
      loss_rec <- compute_mixed_loss(pred_scaled, t_X, masked_bool, trait_map)
    } else {
      loss_rec <- torch::nnf_mse_loss(
        pred_scaled[masked_bool], t_X[masked_bool]
      )
    }

    loss_shrink <- torch::nnf_mse_loss(
      delta[masked_bool],
      torch::torch_zeros_like(delta[masked_bool])
    )
    loss <- loss_rec + lambda_shrink * loss_shrink
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
          out0      <- model(t_X, t_coords, covs0, t_adj)
          pred0     <- t_MU + out0$rs * out0$delta
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
        val_loss    = if (is.na(val_loss_val)) NA_real_ else val_loss_val
      ))

      if (verbose) {
        msg <- sprintf(
          "epoch %4d | rec %.4f | shrink %.4f",
          epoch, as.numeric(loss_rec$item()),
          as.numeric(loss_shrink$item())
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
      out0      <- model(t_X, t_coords, covs0, t_adj)
      pred0     <- t_MU + out0$rs * out0$delta
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
    dropout         = dropout,
    refine_steps    = refine_steps,
    cov_dim         = cov_dim,
    input_dim       = p
  )

  # Backward-compat: store val_rmse and test_rmse names
  structure(
    list(
      model_state   = model$state_dict(),
      model_config  = model_config,
      graph         = graph,
      baseline      = baseline,
      norm          = list(
        means         = data$means,
        sds           = data$sds,
        log_transform = data$log_transform
      ),
      species_names = data$species_names,
      trait_names   = data$trait_names,
      latent_names  = data$latent_names,
      trait_map     = trait_map,
      splits        = splits,
      history       = history,
      val_rmse      = best_val,
      test_rmse     = test_loss
    ),
    class = "pigauto_fit"
  )
}
