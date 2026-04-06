#' Fit a Residual Phylo-DAE model for trait imputation
#'
#' Trains a \code{ResidualPhyloDAE} that learns a bounded residual correction
#' on top of a Brownian Motion baseline from \pkg{Rphylopars}.  Because the
#' residual scale is initialised near zero and constrained to \eqn{(0, 0.5)},
#' the model cannot perform substantially worse than the BM baseline.
#'
#' @details
#' **Training objective** (per epoch):
#' \enumerate{
#'   \item A random subset of observed cells is corrupted with a learnable
#'     mask token.
#'   \item The model predicts the residual \eqn{\delta = X - \mu_{BM}}.
#'   \item Loss = MSE on corrupted cells + \code{lambda_shrink} * MSE toward
#'     zero (shrinks the residual toward the BM baseline).
#' }
#'
#' **Early stopping** is based on RMSE on the validation cells supplied in
#' \code{splits}.
#'
#' @param data object of class \code{"pigauto_data"}.
#' @param tree object of class \code{"phylo"}.
#' @param splits list (output of \code{\link{make_missing_splits}}) or
#'   \code{NULL}. Needed for validation-based early stopping.
#' @param graph list (output of \code{\link{build_phylo_graph}}) or
#'   \code{NULL} (computed automatically with \code{k_eigen} and
#'   \code{sigma_mult = 0.5}).
#' @param baseline list (output of \code{\link{fit_baseline}}) or \code{NULL}
#'   (computed automatically using BM).
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
#' @examples
#' \dontrun{
#' data(avonet300, tree300, package = "pigauto")
#' traits <- avonet300; rownames(traits) <- traits$Species_Key
#' traits$Species_Key <- NULL
#' pd     <- preprocess_traits(traits, tree300)
#' splits <- make_missing_splits(pd$X_scaled)
#' fit    <- fit_pigauto(pd, tree300, splits, epochs = 500, verbose = FALSE)
#' print(fit)
#' }
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
    if (verbose) message("Fitting Rphylopars BM baseline...")
    baseline <- fit_baseline(data, tree, splits = splits)
  }

  # ---- Data preparation ----------------------------------------------------
  n <- nrow(data$X_scaled)
  p <- ncol(data$X_scaled)

  X_truth <- data$X_scaled          # n x p, z-score scale, NAs possible
  MU      <- baseline$mu            # n x p, BM predictions

  # Fill truly missing cells with BM predictions
  X_fill <- X_truth
  X_fill[is.na(X_fill)] <- MU[is.na(X_fill)]

  DELTA_tgt <- X_fill - MU          # residual target (n x p)

  # Observed mask (TRUE = actually observed, not just filled)
  M_obs_mat <- !is.na(X_truth)
  if (!is.null(splits)) M_obs_mat[c(splits$val_idx, splits$test_idx)] <- FALSE

  # Tensors
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

  # Val / test boolean masks for RMSE in torch
  if (!is.null(splits)) {
    val_mat <- matrix(FALSE, n, p); val_mat[splits$val_idx]  <- TRUE
    t_val   <- torch::torch_tensor(val_mat, dtype = torch::torch_bool(),
                                   device = device)
    t_truth <- torch::torch_tensor(X_truth, dtype = torch::torch_float(),
                                   device = device)
    has_val <- TRUE
  } else {
    has_val <- FALSE
  }

  # ---- Model ----------------------------------------------------------------
  # cov_dim = p (mu_BM per trait) + 1 (per-species mask indicator)
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
    loss_shrink = double(0), val_rmse = double(0)
  )

  for (epoch in seq_len(epochs)) {
    model$train()
    opt$zero_grad()

    # Corrupt only truly observed cells
    u            <- torch::torch_rand_like(t_X)
    masked_bool  <- (u < corruption_rate) & t_M_obs
    if (as.numeric(masked_bool$sum()$cpu()$item()) == 0L) next

    Mask_t  <- masked_bool$to(dtype = torch::torch_float())  # per-cell
    # Per-species mask indicator: 1 if any trait corrupted for that species
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

    loss_rec    <- torch::nnf_mse_loss(
      pred_scaled[masked_bool], t_X[masked_bool]
    )
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
      val_rmse <- NA_real_
      if (has_val) {
        model$eval()
        torch::with_no_grad({
          mask_ind0 <- torch::torch_zeros(c(n, 1L), device = device)
          covs0     <- torch::torch_cat(list(t_MU, mask_ind0), dim = 2L)
          out0      <- model(t_X, t_coords, covs0, t_adj)
          pred0     <- t_MU + out0$rs * out0$delta
          Xcurr     <- t_X * (!t_M_obs)$logical_not()$to(
            dtype = torch::torch_float()
          ) + pred0 * (!t_M_obs)$to(dtype = torch::torch_float())
        })
        val_rmse <- as.numeric(
          rmse_torch(Xcurr, t_truth, t_val)$cpu()$item()
        )
      }

      history <- rbind(history, data.frame(
        epoch       = epoch,
        loss_rec    = loss_rec$item(),
        loss_shrink = loss_shrink$item(),
        val_rmse    = if (is.na(val_rmse)) NA_real_ else val_rmse
      ))

      if (verbose) {
        msg <- sprintf(
          "epoch %4d | rec %.4f | shrink %.4f",
          epoch, loss_rec$item(), loss_shrink$item()
        )
        if (!is.na(val_rmse)) {
          msg <- paste0(msg, sprintf(" | val %.4f | best %.4f | pat %d",
                                     val_rmse, best_val, patience_left))
        }
        message(msg)
      }

      if (!is.na(val_rmse)) {
        if (val_rmse + 1e-6 < best_val) {
          best_val      <- val_rmse
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

  # ---- Final test RMSE -------------------------------------------------------
  test_rmse <- NA_real_
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
      Xcurr     <- t_X * t_M_obs$to(dtype = torch::torch_float()) +
        pred0 * (!t_M_obs)$to(dtype = torch::torch_float())
      t_truth_f <- torch::torch_tensor(X_truth, dtype = torch::torch_float(),
                                       device = device)
    })
    test_rmse <- as.numeric(
      rmse_torch(Xcurr, t_truth_f, t_test)$cpu()$item()
    )
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
      splits        = splits,
      history       = history,
      val_rmse      = best_val,
      test_rmse     = test_rmse
    ),
    class = "pigauto_fit"
  )
}
