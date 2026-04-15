#' Evaluate a fitted pigauto model on its test set
#'
#' Computes type-specific performance metrics for each trait.  By default
#' the model is evaluated on the test split stored in the fit object, but
#' an alternative \code{data} / \code{splits} can be supplied.
#'
#' @details
#' **Metrics by trait type:**
#' \describe{
#'   \item{continuous}{RMSE, Pearson r, MAE}
#'   \item{count}{RMSE, MAE, Pearson r}
#'   \item{binary}{Accuracy, Brier score}
#'   \item{categorical}{Accuracy (overall)}
#'   \item{ordinal}{RMSE (on integer scale), Spearman rho}
#' }
#'
#' When conformal scores are present in the fit, conformal coverage at
#' the 95\% nominal level is also reported for continuous, count, and
#' ordinal traits.
#'
#' When the fit includes a baseline, baseline metrics are appended with
#' \code{method = "baseline"} for direct comparison.
#'
#' @param fit pigauto_fit object.
#' @param data pigauto_data object (default: \code{NULL}, uses the
#'   training data stored in the fit via \code{predict()}).
#' @param splits splits object (default: \code{NULL}, uses
#'   \code{fit$splits}).
#' @return A \code{data.frame} with columns: \code{method}, \code{trait},
#'   \code{type}, \code{metric}, \code{value}, \code{n_test}.
#' @examples
#' \dontrun{
#' eval_df <- evaluate(fit)
#' eval_df[eval_df$metric == "rmse", ]
#' }
#' @export
evaluate <- function(fit, data = NULL, splits = NULL) {
  if (!inherits(fit, "pigauto_fit")) {
    stop("'fit' must be a pigauto_fit object.")
  }

  # Resolve splits
  if (is.null(splits)) splits <- fit$splits
  if (is.null(splits) || length(splits$test_idx) == 0L) {
    stop("No test split available. Supply 'splits' or refit with splits.")
  }

  # Get predictions
  pred <- predict(fit, return_se = TRUE)

  # Get truth in latent scale
  trait_map <- fit$trait_map

  # Truth in latent scale.  The fit does not store X_scaled directly,
  # so the original pigauto_data object is required.
  if (!is.null(data)) {
    if (!inherits(data, "pigauto_data")) {
      stop("'data' must be a pigauto_data object.")
    }
    truth_latent <- data$X_scaled
  } else {
    # Attempt to recover truth from the fit object.
    # The fit does not store X_scaled, so we reconstruct what we can.
    # evaluate_imputation() was the previous API for this -- but it
    # required truth.  For the new API, data is strongly recommended.
    stop("'data' is required. Pass the pigauto_data object used for fitting.")
  }

  n <- nrow(truth_latent)
  p <- ncol(truth_latent)

  # Linear indices for test cells
  test_idx <- splits$test_idx

  # Convert linear indices to (row, col)
  row_i <- ((test_idx - 1L) %% n) + 1L
  col_j <- ceiling(test_idx / n)

  # ---- Evaluate pigauto predictions ----------------------------------------
  pigauto_rows <- eval_test_cells(
    pred_latent  = pred$imputed_latent,
    pred_obj     = pred,
    truth_latent = truth_latent,
    row_i        = row_i,
    col_j        = col_j,
    trait_map    = trait_map,
    n            = n,
    method       = "pigauto",
    conformal_scores = fit$conformal_scores
  )

  # ---- Evaluate baseline predictions if available --------------------------
  baseline_rows <- NULL
  if (!is.null(fit$baseline)) {
    mu_species <- fit$baseline$mu
    multi_obs  <- isTRUE(fit$multi_obs)
    if (multi_obs) {
      obs_to_sp <- fit$obs_to_species
      mu <- mu_species[obs_to_sp, , drop = FALSE]
      rownames(mu) <- NULL
    } else {
      mu <- mu_species
    }
    baseline_rows <- eval_test_cells(
      pred_latent  = mu,
      pred_obj     = NULL,
      truth_latent = truth_latent,
      row_i        = row_i,
      col_j        = col_j,
      trait_map    = trait_map,
      n            = n,
      method       = "baseline",
      conformal_scores = NULL
    )
  }

  out <- rbind(pigauto_rows, baseline_rows)
  rownames(out) <- NULL
  out
}


# ---- Internal: evaluate test cells for a single method -----------------------

eval_test_cells <- function(pred_latent, pred_obj, truth_latent,
                            row_i, col_j, trait_map, n,
                            method, conformal_scores = NULL) {

  # Legacy path: no trait_map (all continuous)
  if (is.null(trait_map)) {
    return(eval_test_cells_legacy(
      pred_latent, truth_latent, row_i, col_j, n, method
    ))
  }

  rows <- list()

  for (tm in trait_map) {
    nm <- tm$name
    tp <- tm$type
    lc <- tm$latent_cols

    if (tp %in% c("continuous", "count", "ordinal", "proportion")) {
      # Single latent column per trait
      keep <- which(col_j == lc[1])
      ri   <- row_i[keep]
      if (length(ri) == 0L) next

      p_j <- pred_latent[ri, lc[1]]
      t_j <- truth_latent[ri, lc[1]]
      ok  <- is.finite(p_j) & is.finite(t_j)
      if (sum(ok) == 0L) next

      n_ok <- sum(ok)

      if (tp %in% c("continuous", "proportion")) {
        rmse_val <- rmse_vec(t_j[ok], p_j[ok])
        mae_val  <- mae_vec(t_j[ok], p_j[ok])
        r_val    <- if (n_ok > 2L) stats::cor(t_j[ok], p_j[ok]) else NA_real_

        rows <- c(rows, list(
          data.frame(method = method, trait = nm, type = tp,
                     metric = "rmse", value = rmse_val, n_test = n_ok,
                     stringsAsFactors = FALSE),
          data.frame(method = method, trait = nm, type = tp,
                     metric = "pearson_r", value = r_val, n_test = n_ok,
                     stringsAsFactors = FALSE),
          data.frame(method = method, trait = nm, type = tp,
                     metric = "mae", value = mae_val, n_test = n_ok,
                     stringsAsFactors = FALSE)
        ))

      } else if (tp == "count") {
        rmse_val <- rmse_vec(t_j[ok], p_j[ok])
        mae_val  <- mae_vec(t_j[ok], p_j[ok])
        r_val    <- if (n_ok > 2L) stats::cor(t_j[ok], p_j[ok]) else NA_real_

        rows <- c(rows, list(
          data.frame(method = method, trait = nm, type = tp,
                     metric = "rmse", value = rmse_val, n_test = n_ok,
                     stringsAsFactors = FALSE),
          data.frame(method = method, trait = nm, type = tp,
                     metric = "mae", value = mae_val, n_test = n_ok,
                     stringsAsFactors = FALSE),
          data.frame(method = method, trait = nm, type = tp,
                     metric = "pearson_r", value = r_val, n_test = n_ok,
                     stringsAsFactors = FALSE)
        ))

      } else {
        # ordinal: RMSE on integer scale, Spearman rho
        rmse_val <- rmse_vec(t_j[ok], p_j[ok])
        rho_val  <- if (n_ok > 2L) {
          stats::cor(t_j[ok], p_j[ok], method = "spearman")
        } else {
          NA_real_
        }

        rows <- c(rows, list(
          data.frame(method = method, trait = nm, type = tp,
                     metric = "rmse", value = rmse_val, n_test = n_ok,
                     stringsAsFactors = FALSE),
          data.frame(method = method, trait = nm, type = tp,
                     metric = "spearman_rho", value = rho_val, n_test = n_ok,
                     stringsAsFactors = FALSE)
        ))
      }

      # Conformal coverage for this trait
      if (!is.null(conformal_scores) && !is.na(conformal_scores[nm])) {
        q <- conformal_scores[nm]
        lower <- p_j[ok] - q
        upper <- p_j[ok] + q
        cov <- mean(t_j[ok] >= lower & t_j[ok] <= upper)
        rows <- c(rows, list(
          data.frame(method = method, trait = nm, type = tp,
                     metric = "conformal_coverage_95",
                     value = cov, n_test = n_ok,
                     stringsAsFactors = FALSE)
        ))
      }

    } else if (tp == "binary") {
      keep <- which(col_j == lc[1])
      ri   <- row_i[keep]
      if (length(ri) == 0L) next

      t_j <- truth_latent[ri, lc[1]]
      ok  <- is.finite(t_j)
      if (sum(ok) == 0L) next

      n_ok <- sum(ok)

      # Probabilities: from pred_obj if available, else sigmoid of latent
      if (!is.null(pred_obj) && !is.null(pred_obj$probabilities[[nm]])) {
        prob_j <- pred_obj$probabilities[[nm]][ri]
      } else {
        prob_j <- expit(pred_latent[ri, lc[1]])
      }

      pred_class  <- as.integer(prob_j >= 0.5)
      truth_class <- as.integer(round(t_j))

      acc_val   <- accuracy_vec(truth_class[ok], pred_class[ok])
      brier_val <- brier_vec(t_j[ok], prob_j[ok])

      rows <- c(rows, list(
        data.frame(method = method, trait = nm, type = tp,
                   metric = "accuracy", value = acc_val, n_test = n_ok,
                   stringsAsFactors = FALSE),
        data.frame(method = method, trait = nm, type = tp,
                   metric = "brier", value = brier_val, n_test = n_ok,
                   stringsAsFactors = FALSE)
      ))

    } else if (tp == "categorical") {
      keep <- which(col_j == lc[1])
      ri   <- row_i[keep]
      if (length(ri) == 0L) next

      truth_oh <- truth_latent[ri, lc, drop = FALSE]
      ok       <- stats::complete.cases(truth_oh)
      if (sum(ok) == 0L) next

      n_ok <- sum(ok)
      truth_class <- apply(truth_oh[ok, , drop = FALSE], 1, which.max)
      pred_logits <- pred_latent[ri, lc, drop = FALSE]
      pred_class  <- apply(pred_logits[ok, , drop = FALSE], 1, which.max)

      acc_val <- mean(truth_class == pred_class)

      rows <- c(rows, list(
        data.frame(method = method, trait = nm, type = tp,
                   metric = "accuracy", value = acc_val, n_test = n_ok,
                   stringsAsFactors = FALSE)
      ))

      # Per-class accuracy
      K <- length(tm$levels)
      for (k in seq_len(K)) {
        in_class <- truth_class == k
        if (sum(in_class) == 0L) next
        acc_k <- mean(pred_class[in_class] == k)
        rows <- c(rows, list(
          data.frame(method = method, trait = nm, type = tp,
                     metric = paste0("accuracy_", tm$levels[k]),
                     value = acc_k, n_test = sum(in_class),
                     stringsAsFactors = FALSE)
        ))
      }

    } else if (tp == "zi_count") {
      # ZI count: evaluate on the gate column (col 1)
      keep <- which(col_j == lc[1])
      ri   <- row_i[keep]
      if (length(ri) == 0L) next

      # Gate evaluation: truth column 1 is 0/1
      t_gate <- truth_latent[ri, lc[1]]
      ok_gate <- is.finite(t_gate)
      if (sum(ok_gate) == 0L) next
      n_ok <- sum(ok_gate)

      # Predicted expected value (from pred_obj or latent decode)
      if (!is.null(pred_obj) && !is.null(pred_obj$probabilities[[nm]])) {
        p_nz <- pred_obj$probabilities[[nm]][ri]
      } else {
        p_nz <- expit(pred_latent[ri, lc[1]])
      }

      # Zero-accuracy: predict zero if p_nz < 0.5
      pred_zero <- as.integer(p_nz < 0.5)
      truth_zero <- as.integer(t_gate < 0.5)
      zero_acc <- mean(pred_zero[ok_gate] == truth_zero[ok_gate])

      # Brier score on the gate
      brier_val <- brier_vec(t_gate[ok_gate], p_nz[ok_gate])

      # RMSE on expected value (latent scale is less interpretable for ZI)
      # Use column 1 presence + column 2 magnitude to compute EV
      pred_ev <- rep(0, length(ri))
      truth_ev <- rep(0, length(ri))
      # Truth: reconstruct integer values
      for (idx in seq_along(ri)) {
        if (ok_gate[idx]) {
          if (truth_latent[ri[idx], lc[1]] > 0.5) {
            # Non-zero: reconstruct count
            mag_val <- truth_latent[ri[idx], lc[2]]
            if (is.finite(mag_val)) {
              truth_ev[idx] <- expm1(mag_val * tm$sd + tm$mean)
            }
          }
          # Predicted EV
          mag_pred <- pred_latent[ri[idx], lc[2]]
          if (is.finite(mag_pred)) {
            count_hat <- pmax(expm1(mag_pred * tm$sd + tm$mean), 0)
            pred_ev[idx] <- p_nz[idx] * count_hat
          }
        }
      }
      rmse_val <- rmse_vec(truth_ev[ok_gate], pred_ev[ok_gate])
      mae_val  <- mae_vec(truth_ev[ok_gate], pred_ev[ok_gate])
      r_val    <- if (n_ok > 2L) {
        stats::cor(truth_ev[ok_gate], pred_ev[ok_gate])
      } else {
        NA_real_
      }

      rows <- c(rows, list(
        data.frame(method = method, trait = nm, type = tp,
                   metric = "rmse", value = rmse_val, n_test = n_ok,
                   stringsAsFactors = FALSE),
        data.frame(method = method, trait = nm, type = tp,
                   metric = "mae", value = mae_val, n_test = n_ok,
                   stringsAsFactors = FALSE),
        data.frame(method = method, trait = nm, type = tp,
                   metric = "pearson_r", value = r_val, n_test = n_ok,
                   stringsAsFactors = FALSE),
        data.frame(method = method, trait = nm, type = tp,
                   metric = "zero_accuracy", value = zero_acc, n_test = n_ok,
                   stringsAsFactors = FALSE),
        data.frame(method = method, trait = nm, type = tp,
                   metric = "brier", value = brier_val, n_test = n_ok,
                   stringsAsFactors = FALSE)
      ))

    } else if (tp == "multi_proportion") {
      # Row-level evaluation: a held-out row has ALL K latent cells missing
      # together (group-corruption). We evaluate the full composition.
      keep <- which(col_j == lc[1])
      ri   <- row_i[keep]
      if (length(ri) == 0L) next

      # Truth in CLR space: the first CLR column as a row-selector proxy.
      truth_clr <- truth_latent[ri, lc, drop = FALSE]
      ok <- stats::complete.cases(truth_clr)
      if (sum(ok) == 0L) next
      n_ok <- sum(ok)
      ri_ok <- ri[ok]

      pred_clr  <- pred_latent[ri_ok, lc, drop = FALSE]
      truth_clr <- truth_clr[ok, , drop = FALSE]

      # Aitchison distance = Euclidean distance in CLR space (per cell).
      # We report the average Aitchison distance across held-out rows,
      # calculated on the un-z-scored CLR values.
      pred_clr_un  <- pred_clr
      truth_clr_un <- truth_clr
      for (k in seq_len(tm$n_latent)) {
        pred_clr_un[, k]  <- pred_clr[, k]  * tm$sd[k] + tm$mean[k]
        truth_clr_un[, k] <- truth_clr[, k] * tm$sd[k] + tm$mean[k]
      }
      aitch <- sqrt(rowSums((pred_clr_un - truth_clr_un)^2))
      aitch_mean <- mean(aitch)

      # RMSE on z-scored CLR (directly comparable to continuous traits)
      rmse_val <- rmse_vec(as.vector(truth_clr), as.vector(pred_clr))

      # MAE on the simplex after softmax decode
      #   pred probabilities:
      pred_clr_un <- pred_clr_un - rowMeans(pred_clr_un)
      truth_clr_un <- truth_clr_un - rowMeans(truth_clr_un)
      ex_p <- exp(pred_clr_un  - apply(pred_clr_un,  1, max))
      ex_t <- exp(truth_clr_un - apply(truth_clr_un, 1, max))
      prop_pred  <- ex_p / rowSums(ex_p)
      prop_truth <- ex_t / rowSums(ex_t)
      simplex_mae <- mean(abs(prop_pred - prop_truth))

      rows <- c(rows, list(
        data.frame(method = method, trait = nm, type = tp,
                   metric = "aitchison", value = aitch_mean, n_test = n_ok,
                   stringsAsFactors = FALSE),
        data.frame(method = method, trait = nm, type = tp,
                   metric = "rmse_clr", value = rmse_val, n_test = n_ok,
                   stringsAsFactors = FALSE),
        data.frame(method = method, trait = nm, type = tp,
                   metric = "simplex_mae", value = simplex_mae, n_test = n_ok,
                   stringsAsFactors = FALSE)
      ))
    }
  }

  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}


# ---- Internal: legacy all-continuous evaluation ------------------------------

eval_test_cells_legacy <- function(pred_latent, truth_latent,
                                   row_i, col_j, n, method) {
  trait_names <- colnames(truth_latent)
  if (is.null(trait_names)) {
    trait_names <- paste0("trait", seq_len(ncol(truth_latent)))
  }

  rows <- list()
  for (j in seq_len(ncol(truth_latent))) {
    keep <- which(col_j == j)
    ri   <- row_i[keep]
    if (length(ri) == 0L) next

    p_j <- pred_latent[ri, j]
    t_j <- truth_latent[ri, j]
    ok  <- is.finite(p_j) & is.finite(t_j)
    if (sum(ok) == 0L) next

    n_ok     <- sum(ok)
    rmse_val <- rmse_vec(t_j[ok], p_j[ok])
    mae_val  <- mae_vec(t_j[ok], p_j[ok])
    r_val    <- if (n_ok > 2L) stats::cor(t_j[ok], p_j[ok]) else NA_real_

    rows <- c(rows, list(
      data.frame(method = method, trait = trait_names[j], type = "continuous",
                 metric = "rmse", value = rmse_val, n_test = n_ok,
                 stringsAsFactors = FALSE),
      data.frame(method = method, trait = trait_names[j], type = "continuous",
                 metric = "pearson_r", value = r_val, n_test = n_ok,
                 stringsAsFactors = FALSE),
      data.frame(method = method, trait = trait_names[j], type = "continuous",
                 metric = "mae", value = mae_val, n_test = n_ok,
                 stringsAsFactors = FALSE)
    ))
  }

  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}


# ===========================================================================
# compare_methods
# ===========================================================================

#' Compare BM baseline and pigauto methods across replicates
#'
#' Runs a full comparison of the BM baseline and pigauto (with attention +
#' calibration) over multiple random seeds.  Returns a tidy data.frame
#' suitable for plotting or downstream analysis.
#'
#' @details
#' For each seed the function:
#' \enumerate{
#'   \item Creates train/val/test splits.
#'   \item Fits the phylogenetic BM baseline.
#'   \item Fits pigauto (with attention and calibration enabled by default).
#'   \item Evaluates both methods on the test split.
#'   \item Collects results into a single data.frame.
#' }
#'
#' @param data pigauto_data object.
#' @param tree phylo object.
#' @param splits pre-computed splits (applied to all reps) or \code{NULL}
#'   to create fresh splits per seed.
#' @param seeds integer vector of random seeds for replication.
#' @param epochs number of training epochs.
#' @param verbose logical.
#' @param ... additional arguments passed to \code{\link{fit_pigauto}}.
#' @return A \code{data.frame} with columns: \code{method}, \code{trait},
#'   \code{type}, \code{metric}, \code{value}, \code{rep}.
#' @examples
#' \dontrun{
#' cmp <- compare_methods(pd, tree300, seeds = 1:3, epochs = 500)
#' # Summarise across reps
#' aggregate(value ~ method + trait + metric, data = cmp, FUN = mean)
#' }
#' @export
compare_methods <- function(data, tree, splits = NULL, seeds = 1:3,
                            epochs = 500L, verbose = TRUE, ...) {
  if (!inherits(data, "pigauto_data")) {
    stop("'data' must be a pigauto_data object.")
  }
  if (!inherits(tree, "phylo")) {
    stop("'tree' must be a phylo object.")
  }

  all_results <- vector("list", length(seeds))

  for (i in seq_along(seeds)) {
    s <- seeds[i]
    if (verbose) message("=== Replicate ", i, "/", length(seeds),
                         " (seed ", s, ") ===")

    # ---- Splits -----------------------------------------------------------
    if (is.null(splits)) {
      sp <- make_missing_splits(data$X_scaled, seed = s,
                                trait_map = data$trait_map)
    } else {
      sp <- splits
    }

    # ---- Baseline ---------------------------------------------------------
    if (verbose) message("Fitting baseline...")
    bl <- fit_baseline(data, tree, splits = sp)

    # ---- Pigauto (attention + calibration) --------------------------------
    if (verbose) message("Fitting pigauto...")
    fit <- fit_pigauto(
      data    = data,
      tree    = tree,
      splits  = sp,
      baseline = bl,
      epochs  = as.integer(epochs),
      verbose = verbose,
      seed    = s,
      ...
    )

    # ---- Evaluate ---------------------------------------------------------
    if (verbose) message("Evaluating...")
    eval_df <- evaluate(fit, data = data, splits = sp)
    eval_df$rep <- i

    all_results[[i]] <- eval_df
  }

  out <- do.call(rbind, all_results)
  rownames(out) <- NULL
  out
}


# ===========================================================================
# summary.pigauto_fit
# ===========================================================================

#' Summary method for pigauto_fit objects
#'
#' Prints a formatted evaluation table including per-trait metrics,
#' gate calibration, and conformal coverage.  Requires the original
#' \code{pigauto_data} object to compute test-set performance.
#'
#' @param object pigauto_fit object.
#' @param data pigauto_data object used for fitting (optional; when
#'   \code{NULL} the summary skips per-trait test metrics).
#' @param ... ignored.
#' @return Invisibly returns the evaluation data.frame (or \code{NULL} if
#'   \code{data} is not supplied).
#' @examples
#' \dontrun{
#' summary(fit, data = pd)
#' }
#' @export
summary.pigauto_fit <- function(object, ..., data = NULL) {
  cfg       <- object$model_config
  trait_map <- object$trait_map
  n_species <- object$n_species %||% length(object$species_names)
  n_traits  <- if (!is.null(trait_map)) length(trait_map) else length(object$trait_names)

  # ---- Header -------------------------------------------------------------
  rule <- strrep("-", 56)
  cat("pigauto_fit summary\n")
  cat(rule, "\n")
  n_epochs <- if (nrow(object$history) > 0L) max(object$history$epoch) else 0L
  cat(sprintf("Species: %d | Traits: %d | Epochs: %d\n",
              n_species, n_traits, n_epochs))
  cat(sprintf("Architecture: hidden_dim=%d, k_eigen=%d, n_gnn=%d, attention=%s\n",
              cfg$hidden_dim, cfg$k_eigen,
              cfg$n_gnn_layers %||% 1L,
              if (isTRUE(cfg$use_attention)) "TRUE" else "FALSE"))

  if (!is.null(trait_map)) {
    types <- vapply(trait_map, "[[", character(1), "type")
    type_tab <- table(types)
    cat("Trait types: ",
        paste(names(type_tab), type_tab, sep = "=", collapse = ", "), "\n")
  }
  cat("\n")

  # ---- Validation / test loss ---------------------------------------------
  if (is.finite(object$val_rmse)) {
    cat(sprintf("Best validation loss: %.4f\n", object$val_rmse))
  }
  if (!is.na(object$test_rmse)) {
    cat(sprintf("Test loss:            %.4f\n", object$test_rmse))
  }
  cat("\n")

  # ---- Per-trait performance (if data available) ---------------------------
  eval_df <- NULL
  if (!is.null(data) && !is.null(object$splits) &&
      length(object$splits$test_idx) > 0L) {

    eval_df <- tryCatch(
      evaluate(object, data = data, splits = object$splits),
      error = function(e) {
        message("Could not compute test metrics: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(eval_df)) {
      cat("Trait Performance (test set):\n")
      cat(strrep("-", 56), "\n")

      # Pivot to wide format for display: one row per (method, trait)
      pigauto_df <- eval_df[eval_df$method == "pigauto", ]
      trait_names_unique <- unique(pigauto_df$trait)

      # Determine which metric columns to show
      has_rmse     <- any(pigauto_df$metric == "rmse")
      has_pearson  <- any(pigauto_df$metric == "pearson_r")
      has_accuracy <- any(pigauto_df$metric == "accuracy")
      has_brier    <- any(pigauto_df$metric == "brier")
      has_spearman <- any(pigauto_df$metric == "spearman_rho")
      has_mae      <- any(pigauto_df$metric == "mae")

      # Header
      hdr <- sprintf("  %-25s %-13s %6s", "Trait", "Type", "n")
      if (has_rmse)     hdr <- paste0(hdr, sprintf(" %8s", "RMSE"))
      if (has_pearson)  hdr <- paste0(hdr, sprintf(" %8s", "r"))
      if (has_spearman) hdr <- paste0(hdr, sprintf(" %8s", "rho"))
      if (has_accuracy) hdr <- paste0(hdr, sprintf(" %8s", "Acc"))
      if (has_brier)    hdr <- paste0(hdr, sprintf(" %8s", "Brier"))
      if (has_mae)      hdr <- paste0(hdr, sprintf(" %8s", "MAE"))
      cat(hdr, "\n")
      cat("  ", strrep("-", nchar(hdr) - 2), "\n", sep = "")

      for (tnm in trait_names_unique) {
        sub <- pigauto_df[pigauto_df$trait == tnm &
                            !grepl("^accuracy_", pigauto_df$metric) &
                            pigauto_df$metric != "conformal_coverage_95", ]
        if (nrow(sub) == 0L) next

        tp   <- sub$type[1]
        n_t  <- sub$n_test[1]

        get_val <- function(m) {
          v <- sub$value[sub$metric == m]
          if (length(v) == 0L) NA_real_ else v[1]
        }

        line <- sprintf("  %-25s %-13s %6d", tnm, tp, n_t)
        if (has_rmse)     line <- paste0(line, format_metric(get_val("rmse")))
        if (has_pearson)  line <- paste0(line, format_metric(get_val("pearson_r")))
        if (has_spearman) line <- paste0(line, format_metric(get_val("spearman_rho")))
        if (has_accuracy) line <- paste0(line, format_metric(get_val("accuracy")))
        if (has_brier)    line <- paste0(line, format_metric(get_val("brier")))
        if (has_mae)      line <- paste0(line, format_metric(get_val("mae")))
        cat(line, "\n")
      }

      # Baseline comparison row (if available)
      bl_df <- eval_df[eval_df$method == "baseline", ]
      if (nrow(bl_df) > 0L) {
        cat("\n  Baseline (BM) for comparison:\n")
        cat("  ", strrep("-", nchar(hdr) - 2), "\n", sep = "")

        for (tnm in trait_names_unique) {
          sub <- bl_df[bl_df$trait == tnm &
                         !grepl("^accuracy_", bl_df$metric) &
                         bl_df$metric != "conformal_coverage_95", ]
          if (nrow(sub) == 0L) next

          tp   <- sub$type[1]
          n_t  <- sub$n_test[1]

          get_val <- function(m) {
            v <- sub$value[sub$metric == m]
            if (length(v) == 0L) NA_real_ else v[1]
          }

          line <- sprintf("  %-25s %-13s %6d", tnm, tp, n_t)
          if (has_rmse)     line <- paste0(line, format_metric(get_val("rmse")))
          if (has_pearson)  line <- paste0(line, format_metric(get_val("pearson_r")))
          if (has_spearman) line <- paste0(line, format_metric(get_val("spearman_rho")))
          if (has_accuracy) line <- paste0(line, format_metric(get_val("accuracy")))
          if (has_brier)    line <- paste0(line, format_metric(get_val("brier")))
          if (has_mae)      line <- paste0(line, format_metric(get_val("mae")))
          cat(line, "\n")
        }
      }
      cat("\n")
    }
  }

  # ---- Gate calibration ---------------------------------------------------
  if (!is.null(object$calibrated_gates)) {
    cat("Gate Calibration:\n")

    gates <- object$calibrated_gates
    if (!is.null(trait_map)) {
      # Show one gate value per original trait (first latent col)
      gate_vals <- vapply(trait_map, function(tm) {
        gates[tm$latent_cols[1]]
      }, numeric(1))
      names(gate_vals) <- vapply(trait_map, "[[", character(1), "name")
    } else {
      gate_vals <- gates
      names(gate_vals) <- object$trait_names
    }

    # Print in rows of up to 5
    nms <- names(gate_vals)
    n_gates <- length(nms)
    chunk_size <- 5L
    for (start in seq(1L, n_gates, by = chunk_size)) {
      end <- min(start + chunk_size - 1L, n_gates)
      items <- vapply(start:end, function(k) {
        sprintf("%s: %.3f", nms[k], gate_vals[k])
      }, character(1))
      cat("  ", paste(items, collapse = "  "), "\n")
    }
    cat("\n")
  }

  # ---- Conformal coverage -------------------------------------------------
  if (!is.null(object$conformal_scores)) {
    cs <- object$conformal_scores[!is.na(object$conformal_scores)]
    if (length(cs) > 0L) {
      cat("Conformal Scores (latent scale, 95% nominal):\n")
      items <- vapply(names(cs), function(nm) {
        sprintf("%s: %.4f", nm, cs[nm])
      }, character(1))
      cat("  ", paste(items, collapse = "  "), "\n")

      # If we computed eval_df, show actual coverage
      if (!is.null(eval_df)) {
        cov_rows <- eval_df[eval_df$method == "pigauto" &
                              eval_df$metric == "conformal_coverage_95", ]
        if (nrow(cov_rows) > 0L) {
          cat("Conformal Coverage (test set):\n")
          items <- vapply(seq_len(nrow(cov_rows)), function(k) {
            sprintf("%s: %.1f%%", cov_rows$trait[k], cov_rows$value[k] * 100)
          }, character(1))
          cat("  ", paste(items, collapse = "  "), "\n")
        }
      }
      cat("\n")
    }
  }

  invisible(eval_df)
}


# ---- Internal: format a metric value for display ----------------------------

format_metric <- function(x) {
  if (is.na(x)) {
    sprintf(" %8s", "-")
  } else {
    sprintf(" %8.3f", x)
  }
}
