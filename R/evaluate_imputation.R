#' Evaluate imputation performance against known values
#'
#' Computes type-specific metrics for each trait on the validation and test
#' splits.  When a \code{trait_map} is supplied, metrics are dispatched per
#' trait type; otherwise the function falls back to continuous-only metrics
#' (RMSE, Pearson r, 95% coverage).
#'
#' @details
#' **Metrics by trait type:**
#' \describe{
#'   \item{continuous}{RMSE, Pearson r, 95\% interval coverage (if SE
#'     supplied)}
#'   \item{proportion}{RMSE, Pearson r, 95\% interval coverage (if SE
#'     supplied)}
#'   \item{count}{RMSE, MAE, Pearson r}
#'   \item{ordinal}{RMSE, Spearman rho}
#'   \item{binary}{Accuracy, Brier score}
#'   \item{categorical}{Accuracy}
#'   \item{zi_count}{RMSE, MAE, Pearson r, zero-accuracy, Brier score
#'     (on gate)}
#' }
#'
#' For binary and categorical traits the function accepts either a
#' \code{pigauto_pred} object (preferred, gives access to probabilities) or
#' raw matrices (latent scale).
#'
#' @param pred predicted values: either a numeric matrix in latent scale
#'   (same dimensions as \code{truth}), or a \code{"pigauto_pred"} object
#'   from \code{\link{predict.pigauto_fit}}.
#' @param truth numeric matrix of true values in latent scale (from
#'   \code{pigauto_data$X_scaled}).
#' @param splits list (output of \code{\link{make_missing_splits}}).
#' @param pred_se numeric matrix of prediction SEs (same scale as
#'   \code{pred}).  Used for 95\% coverage of continuous traits.  Ignored
#'   when \code{pred} is a \code{pigauto_pred} (uses \code{pred$se}).
#' @param trait_map list of trait descriptors (from \code{pigauto_data}).
#'   If \code{NULL} and \code{pred} is not a \code{pigauto_pred}, the
#'   v0.1 all-continuous evaluation is used.
#' @return A \code{data.frame} with columns \code{split}, \code{trait},
#'   \code{type}, \code{n}, and type-specific metric columns.
#' @examples
#' \dontrun{
#' # From pigauto_pred object
#' eval_df <- evaluate_imputation(pred_obj, pd$X_scaled, splits)
#'
#' # From raw latent matrix
#' eval_df <- evaluate_imputation(bl$mu, pd$X_scaled, splits,
#'                                 trait_map = pd$trait_map)
#' }
#' @importFrom stats complete.cases
#' @export
evaluate_imputation <- function(pred, truth, splits, pred_se = NULL,
                                trait_map = NULL) {
  # Extract from pigauto_pred if supplied
  if (inherits(pred, "pigauto_pred")) {
    trait_map <- pred$trait_map
    probs     <- pred$probabilities
    pred_mat  <- pred$imputed_latent
    # Use latent-scale SE for coverage computation (matches truth scale)
    if (is.null(pred_se) && !is.null(pred$se_latent)) {
      pred_se <- pred$se_latent
    }
  } else {
    pred_mat <- pred
    probs    <- NULL
  }

  if (!is.matrix(pred_mat) || !is.matrix(truth)) {
    stop("'pred' and 'truth' must be matrices (or pred a pigauto_pred).")
  }
  if (!identical(dim(pred_mat), dim(truth))) {
    stop("Prediction and truth matrices must have identical dimensions.")
  }

  n <- nrow(truth)
  p <- ncol(truth)

  # ---- Legacy: no trait_map, all-continuous --------------------------------
  if (is.null(trait_map)) {
    return(evaluate_continuous_legacy(pred_mat, truth, splits, pred_se))
  }

  # ---- Mixed-type evaluation -----------------------------------------------
  eval_split <- function(idx_latent, label) {
    # Convert latent linear indices to (row, col)
    row_i <- ((idx_latent - 1L) %% n) + 1L
    col_j <- ceiling(idx_latent / n)

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

        p_j <- pred_mat[ri, lc[1]]
        t_j <- truth[ri, lc[1]]
        ok  <- is.finite(p_j) & is.finite(t_j)
        if (sum(ok) == 0L) next

        rmse_val <- rmse_vec(t_j[ok], p_j[ok])

        if (tp %in% c("continuous", "proportion")) {
          r_val <- if (sum(ok) > 1L) stats::cor(t_j[ok], p_j[ok]) else NA_real_

          cov95 <- NA_real_
          if (!is.null(pred_se)) {
            se_j  <- pred_se[ri, lc[1]]
            lower <- p_j - 1.96 * se_j
            upper <- p_j + 1.96 * se_j
            cov95 <- mean((t_j[ok] >= lower[ok]) & (t_j[ok] <= upper[ok]))
          }

          rows[[length(rows) + 1L]] <- data.frame(
            split = label, trait = nm, type = tp, n = sum(ok),
            rmse = rmse_val, pearson_r = r_val, coverage_95 = cov95,
            mae = NA_real_, spearman_rho = NA_real_,
            accuracy = NA_real_, brier = NA_real_,
            zero_accuracy = NA_real_,
            aitchison = NA_real_, rmse_clr = NA_real_, simplex_mae = NA_real_,
          stringsAsFactors = FALSE
          )

        } else if (tp == "count") {
          r_val   <- if (sum(ok) > 1L) stats::cor(t_j[ok], p_j[ok]) else NA_real_
          mae_val <- mae_vec(t_j[ok], p_j[ok])

          rows[[length(rows) + 1L]] <- data.frame(
            split = label, trait = nm, type = tp, n = sum(ok),
            rmse = rmse_val, pearson_r = r_val, coverage_95 = NA_real_,
            mae = mae_val, spearman_rho = NA_real_,
            accuracy = NA_real_, brier = NA_real_,
            zero_accuracy = NA_real_,
            aitchison = NA_real_, rmse_clr = NA_real_, simplex_mae = NA_real_,
          stringsAsFactors = FALSE
          )

        } else {
          # ordinal
          rho_val <- if (sum(ok) > 1L) {
            stats::cor(t_j[ok], p_j[ok], method = "spearman")
          } else {
            NA_real_
          }

          rows[[length(rows) + 1L]] <- data.frame(
            split = label, trait = nm, type = tp, n = sum(ok),
            rmse = rmse_val, pearson_r = NA_real_, coverage_95 = NA_real_,
            mae = NA_real_, spearman_rho = rho_val,
            accuracy = NA_real_, brier = NA_real_,
            zero_accuracy = NA_real_,
            aitchison = NA_real_, rmse_clr = NA_real_, simplex_mae = NA_real_,
          stringsAsFactors = FALSE
          )
        }

      } else if (tp == "binary") {
        keep <- which(col_j == lc[1])
        ri   <- row_i[keep]
        if (length(ri) == 0L) next

        t_j <- truth[ri, lc[1]]
        ok  <- is.finite(t_j)
        if (sum(ok) == 0L) next

        # Binary predictions: threshold at 0.5 after sigmoid
        prob_j <- expit(pred_mat[ri, lc[1]])
        pred_class <- as.integer(prob_j >= 0.5)
        truth_class <- as.integer(round(t_j))

        acc_val   <- accuracy_vec(truth_class[ok], pred_class[ok])
        brier_val <- brier_vec(t_j[ok], prob_j[ok])

        rows[[length(rows) + 1L]] <- data.frame(
          split = label, trait = nm, type = tp, n = sum(ok),
          rmse = NA_real_, pearson_r = NA_real_, coverage_95 = NA_real_,
          mae = NA_real_, spearman_rho = NA_real_,
          accuracy = acc_val, brier = brier_val,
          zero_accuracy = NA_real_,
          aitchison = NA_real_, rmse_clr = NA_real_, simplex_mae = NA_real_,
          stringsAsFactors = FALSE
        )

      } else if (tp == "categorical") {
        # Species are masked if ALL K latent cols are masked
        keep <- which(col_j == lc[1])
        ri   <- row_i[keep]
        if (length(ri) == 0L) next

        # Truth: one-hot in latent columns -> argmax
        truth_oh <- truth[ri, lc, drop = FALSE]
        ok       <- complete.cases(truth_oh)
        if (sum(ok) == 0L) next

        truth_class <- apply(truth_oh[ok, , drop = FALSE], 1, which.max)

        # Pred: softmax of latent logits -> argmax
        pred_logits <- pred_mat[ri, lc, drop = FALSE]
        pred_class  <- apply(pred_logits[ok, , drop = FALSE], 1, which.max)

        acc_val <- mean(truth_class == pred_class)

        rows[[length(rows) + 1L]] <- data.frame(
          split = label, trait = nm, type = tp, n = sum(ok),
          rmse = NA_real_, pearson_r = NA_real_, coverage_95 = NA_real_,
          mae = NA_real_, spearman_rho = NA_real_,
          accuracy = acc_val, brier = NA_real_,
          zero_accuracy = NA_real_,
          aitchison = NA_real_, rmse_clr = NA_real_, simplex_mae = NA_real_,
          stringsAsFactors = FALSE
        )

      } else if (tp == "zi_count") {
        # ZI count: 2 latent columns (gate + magnitude)
        keep <- which(col_j == lc[1])
        ri   <- row_i[keep]
        if (length(ri) == 0L) next

        # Gate column truth (0/1)
        t_gate <- truth[ri, lc[1]]
        ok_gate <- is.finite(t_gate)
        if (sum(ok_gate) == 0L) next
        n_ok <- sum(ok_gate)

        # Gate probability: sigmoid of latent
        if (!is.null(probs) && !is.null(probs[[nm]])) {
          p_nz <- probs[[nm]][ri]
        } else {
          p_nz <- expit(pred_mat[ri, lc[1]])
        }

        # Zero-accuracy: predict zero if p_nz < 0.5
        pred_zero  <- as.integer(p_nz < 0.5)
        truth_zero <- as.integer(t_gate < 0.5)
        zero_acc   <- mean(pred_zero[ok_gate] == truth_zero[ok_gate])

        # Brier score on the gate
        brier_val <- brier_vec(t_gate[ok_gate], p_nz[ok_gate])

        # RMSE / MAE / Pearson r on expected value (original count scale)
        pred_ev  <- rep(0, length(ri))
        truth_ev <- rep(0, length(ri))
        for (idx in seq_along(ri)) {
          if (ok_gate[idx]) {
            # Truth EV
            if (truth[ri[idx], lc[1]] > 0.5) {
              mag_val <- truth[ri[idx], lc[2]]
              if (is.finite(mag_val)) {
                truth_ev[idx] <- expm1(mag_val * tm$sd + tm$mean)
              }
            }
            # Predicted EV = p(non-zero) * count_hat
            mag_pred <- pred_mat[ri[idx], lc[2]]
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

        rows[[length(rows) + 1L]] <- data.frame(
          split = label, trait = nm, type = tp, n = n_ok,
          rmse = rmse_val, pearson_r = r_val, coverage_95 = NA_real_,
          mae = mae_val, spearman_rho = NA_real_,
          accuracy = NA_real_, brier = brier_val,
          zero_accuracy = zero_acc,
          aitchison = NA_real_, rmse_clr = NA_real_, simplex_mae = NA_real_,
          stringsAsFactors = FALSE
        )

      } else if (tp == "multi_proportion") {
        # Row-level: all K CLR cells of a held-out row are held out together.
        # Use column lc[1] as the proxy row selector.
        keep <- which(col_j == lc[1])
        ri   <- row_i[keep]
        if (length(ri) == 0L) next

        truth_clr <- truth[ri, lc, drop = FALSE]
        ok <- stats::complete.cases(truth_clr)
        if (sum(ok) == 0L) next
        n_ok <- sum(ok)
        ri_ok <- ri[ok]
        pred_clr  <- pred_mat[ri_ok, lc, drop = FALSE]
        truth_clr <- truth_clr[ok, , drop = FALSE]

        # Aitchison distance on UN-z-scored CLR values
        pred_un  <- pred_clr
        truth_un <- truth_clr
        for (k in seq_len(tm$n_latent)) {
          pred_un[, k]  <- pred_clr[, k]  * tm$sd[k] + tm$mean[k]
          truth_un[, k] <- truth_clr[, k] * tm$sd[k] + tm$mean[k]
        }
        aitch_mean <- mean(sqrt(rowSums((pred_un - truth_un)^2)))

        # RMSE on z-scored CLR
        rmse_val <- rmse_vec(as.vector(truth_clr), as.vector(pred_clr))

        # Simplex MAE after softmax decode
        pred_un  <- pred_un  - rowMeans(pred_un)
        truth_un <- truth_un - rowMeans(truth_un)
        ex_p <- exp(pred_un  - apply(pred_un,  1, max))
        ex_t <- exp(truth_un - apply(truth_un, 1, max))
        prop_p <- ex_p / rowSums(ex_p)
        prop_t <- ex_t / rowSums(ex_t)
        simplex_mae <- mean(abs(prop_p - prop_t))

        rows[[length(rows) + 1L]] <- data.frame(
          split = label, trait = nm, type = tp, n = n_ok,
          rmse = NA_real_, pearson_r = NA_real_, coverage_95 = NA_real_,
          mae = NA_real_, spearman_rho = NA_real_,
          accuracy = NA_real_, brier = NA_real_, zero_accuracy = NA_real_,
          aitchison = aitch_mean, rmse_clr = rmse_val,
          simplex_mae = simplex_mae,
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(rows) == 0L) return(NULL)
    do.call(rbind, rows)
  }

  result <- list()
  if (length(splits$val_idx) > 0L) {
    result$val <- eval_split(splits$val_idx, "val")
  }
  if (length(splits$test_idx) > 0L) {
    result$test <- eval_split(splits$test_idx, "test")
  }

  out <- do.call(rbind, result[!sapply(result, is.null)])
  if (is.null(out)) {
    out <- data.frame(
      split = character(0), trait = character(0), type = character(0),
      n = integer(0), rmse = double(0), pearson_r = double(0),
      coverage_95 = double(0), mae = double(0), spearman_rho = double(0),
      accuracy = double(0), brier = double(0), zero_accuracy = double(0),
      stringsAsFactors = FALSE
    )
  }
  rownames(out) <- NULL
  out
}


# ---- Legacy continuous-only evaluation ----------------------------------------

evaluate_continuous_legacy <- function(pred, truth, splits, pred_se) {
  trait_names <- colnames(truth)
  if (is.null(trait_names)) trait_names <- paste0("trait", seq_len(ncol(truth)))
  n <- nrow(truth)

  eval_split <- function(idx, label) {
    row_i <- ((idx - 1L) %% n) + 1L
    col_j <- ceiling(idx / n)

    lapply(seq_len(ncol(truth)), function(j) {
      keep <- which(col_j == j)
      ri   <- row_i[keep]
      if (length(ri) == 0L) return(NULL)

      p_j <- pred[ri, j]
      t_j <- truth[ri, j]
      ok  <- is.finite(p_j) & is.finite(t_j)
      if (sum(ok) == 0L) return(NULL)

      rmse_j <- rmse_vec(t_j[ok], p_j[ok])
      r_j    <- if (sum(ok) > 1L) stats::cor(t_j[ok], p_j[ok]) else NA_real_

      cov95 <- NA_real_
      if (!is.null(pred_se) && sum(ok) > 0L) {
        se_j  <- pred_se[ri, j]
        lower <- p_j - 1.96 * se_j
        upper <- p_j + 1.96 * se_j
        cov95 <- mean((t_j[ok] >= lower[ok]) & (t_j[ok] <= upper[ok]))
      }

      data.frame(
        split = label, trait = trait_names[j], type = "continuous",
        n = sum(ok), rmse = rmse_j, pearson_r = r_j,
        coverage_95 = cov95, mae = NA_real_, spearman_rho = NA_real_,
        accuracy = NA_real_, brier = NA_real_, zero_accuracy = NA_real_,
        aitchison = NA_real_, rmse_clr = NA_real_, simplex_mae = NA_real_,
          stringsAsFactors = FALSE
      )
    }) |>
      (\(lst) do.call(rbind, lst[!sapply(lst, is.null)]))()
  }

  rows <- list()
  if (length(splits$val_idx) > 0L) rows$val <- eval_split(splits$val_idx, "val")
  if (length(splits$test_idx) > 0L) rows$test <- eval_split(splits$test_idx, "test")

  do.call(rbind, rows[!sapply(rows, is.null)])
}
