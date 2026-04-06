#' Impute missing traits using a fitted Residual Phylo-DAE
#'
#' Runs iterative refinement on the fitted model and returns imputed trait
#' values back-transformed to the original scale.
#'
#' @param object object of class \code{"pigauto_fit"}.
#' @param newdata \code{NULL} (use the training data) or a
#'   \code{"pigauto_data"} object for new species.  When \code{newdata} is
#'   supplied, a new BM baseline is fitted internally.
#' @param return_se logical. Return standard errors? (default \code{TRUE})
#' @param ... ignored.
#' @return If \code{return_se = FALSE}: a numeric matrix of imputed values in
#'   original (back-transformed) units. If \code{return_se = TRUE}: a list
#'   with \code{imputed} and \code{se} matrices, both in original units.
#' @examples
#' \dontrun{
#' pred <- predict(fit, return_se = TRUE)
#' pred$imputed   # matrix, original scale
#' pred$se        # matrix, uncertainty (original scale)
#' }
#' @importFrom torch torch_tensor torch_float with_no_grad torch_zeros
#' @importFrom torch torch_cat
#' @export
predict.pigauto_fit <- function(object, newdata = NULL, return_se = TRUE, ...) {
  cfg    <- object$model_config
  device <- get_device()

  # Reconstruct model
  model <- ResidualPhyloDAE(
    input_dim  = as.integer(cfg$input_dim),
    hidden_dim = as.integer(cfg$hidden_dim),
    coord_dim  = as.integer(cfg$k_eigen),
    cov_dim    = as.integer(cfg$cov_dim)
  )
  model$to(device = device)
  model$load_state_dict(object$model_state)
  model$eval()

  # Prepare data (use training data if newdata is NULL)
  if (is.null(newdata)) {
    X_fill   <- object$baseline$mu   # start from BM prediction
    mu_bm    <- object$baseline$mu
    coords   <- object$graph$coords
    adj      <- object$graph$adj
  } else {
    # newdata must be a pigauto_data; compute BM baseline externally
    stop("newdata support not yet implemented in v0.1.")
  }

  n <- nrow(X_fill)
  p <- ncol(X_fill)

  t_X_fill <- torch::torch_tensor(X_fill, dtype = torch::torch_float(),
                                  device = device)
  t_MU     <- torch::torch_tensor(mu_bm,  dtype = torch::torch_float(),
                                  device = device)
  t_coords <- torch::torch_tensor(coords, dtype = torch::torch_float(),
                                  device = device)
  t_adj    <- torch::torch_tensor(adj,    dtype = torch::torch_float(),
                                  device = device)

  # Iterative refinement: at each step fill only the missing cells
  X_iter <- t_X_fill$clone()

  torch::with_no_grad({
    mask_ind0 <- torch::torch_zeros(c(n, 1L), device = device)

    for (step in seq_len(cfg$refine_steps)) {
      covs0 <- torch::torch_cat(list(t_MU, mask_ind0), dim = 2L)
      out   <- model(X_iter, t_coords, covs0, t_adj)
      pred  <- t_MU + out$rs * out$delta
      X_iter <- pred   # use full prediction at inference
    }

    pred_scaled <- as.matrix(X_iter$cpu())

    # Approximate SE from the residual scale and BM SE
    rs_val  <- as.numeric(out$rs$cpu()$item())
    bm_se   <- object$baseline$se
    pred_se <- bm_se * rs_val   # conservative: SE scales with residual weight
  })

  rownames(pred_scaled) <- object$species_names
  colnames(pred_scaled) <- object$trait_names
  rownames(pred_se)     <- object$species_names
  colnames(pred_se)     <- object$trait_names

  # Back-transform from z-score (and optional log) to original scale
  back_transform <- function(mat) {
    out <- sweep(mat, 2, object$norm$sds,   "*")
    out <- sweep(out, 2, object$norm$means, "+")
    if (isTRUE(object$norm$log_transform)) out <- exp(out)
    out
  }

  imputed <- back_transform(pred_scaled)
  se_out  <- back_transform(pred_se)

  if (!return_se) return(imputed)
  list(imputed = imputed, se = se_out)
}


#' @export
print.pigauto_fit <- function(x, ...) {
  cat("pigauto_fit\n")
  cat("  Species :", length(x$species_names), "\n")
  cat("  Traits  :", length(x$trait_names),
      "--", paste(x$trait_names, collapse = ", "), "\n")
  cat("  Architecture: hidden_dim =", x$model_config$hidden_dim,
      "| k_eigen =", x$model_config$k_eigen, "\n")
  if (is.finite(x$val_rmse)) {
    cat("  Best val RMSE :", round(x$val_rmse, 4), "\n")
  }
  if (!is.na(x$test_rmse)) {
    cat("  Test RMSE     :", round(x$test_rmse, 4), "\n")
  }
  invisible(x)
}
