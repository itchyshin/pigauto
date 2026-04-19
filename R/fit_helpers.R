# R/fit_helpers.R
#
# Internal helpers extracted from fit_pigauto() to keep the training
# function below ~600 lines and each logical block independently readable.
# None of these functions are exported.

# ---------------------------------------------------------------------------
# make_covs_tensor()
# ---------------------------------------------------------------------------
#
# Concatenate the baseline-mean matrix (t_MU), the corruption/inference
# indicator (mask_ind), and any user-supplied covariates into the single
# covs tensor expected by ResidualPhyloDAE$forward().
#
# Signature:
#   covs_t = [t_MU | mask_ind | t_covariates]   (n_obs × (p + 1 + k))
#
# Called during the training loop (mask_ind = actual corruption mask) and
# during evaluation/calibration/test passes (mask_ind = zeros).
#
# @param t_MU         torch tensor (n_obs × p) — baseline mean predictions
# @param mask_ind     torch tensor (n_obs × 1) — 1 if cell was corrupted
# @param t_covariates torch tensor or NULL — user covariates (n_obs × k)
# @return torch tensor (n_obs × (p + 1 + k))
make_covs_tensor <- function(t_MU, mask_ind, t_covariates) {
  cov_parts <- list(t_MU, mask_ind)
  if (!is.null(t_covariates))
    cov_parts[[length(cov_parts) + 1L]] <- t_covariates
  torch::torch_cat(cov_parts, dim = 2L)
}

# ---------------------------------------------------------------------------
# calibrate_gates()
# ---------------------------------------------------------------------------
#
# Post-training per-trait gate calibration on a held-out validation set.
#
# Uses a 9-point grid search with a half-A / half-B split-validation
# cross-check to prevent the GNN from degrading traits the BM baseline
# already handles well.  For discrete traits the calibration metric is
# 0-1 loss (not cross-entropy) so the val signal matches the test metric.
#
# The gate is accepted only if the improvement passes BOTH a relative-gain
# floor (cal_min_rel_gain) AND, for discrete traits, an absolute cell-level
# floor (2 correctly-flipped cells on half A, 1 on half B).
#
# @param trait_map        list of trait descriptors from preprocess_traits()
# @param mu_cal           numeric matrix (n_obs × p) — baseline predictions
# @param delta_cal        numeric matrix (n_obs × p) — GNN predictions
# @param X_truth_r        numeric matrix (n_obs × p) — original data with NAs
# @param val_mask_mat     logical matrix (n_obs × p) — TRUE = validation cell
# @param gate_grid        numeric vector — candidate gate values (9-pt grid)
# @param gate_cap         numeric scalar — maximum allowed gate value
# @param cal_min_rel_gain numeric — minimum relative improvement (default 0.02)
# @param seed             integer — random seed for half-split (seed + 17L)
# @param latent_names     character vector or NULL — for verbose output
# @param verbose          logical
# @return named numeric vector of length p; values in [0, gate_cap]
calibrate_gates <- function(trait_map, mu_cal, delta_cal,
                            X_truth_r, val_mask_mat,
                            gate_grid, gate_cap,
                            cal_min_rel_gain = 0.02,
                            gate_method = c("single_split", "median_splits"),
                            gate_splits_B = 31L,
                            min_val_cells = 10L,
                            seed = 1L,
                            latent_names = NULL,
                            verbose = FALSE) {
  gate_method <- match.arg(gate_method)
  p                <- ncol(mu_cal)
  calibrated_gates <- numeric(p)
  low_val_traits   <- character(0)

  for (tm in trait_map) {
    lc <- tm$latent_cols

    # Restrict to validation rows that have finite truth values
    val_cells <- val_mask_mat[, lc[1]]
    if (length(lc) > 1L) {
      truth_ok <- rowSums(is.na(X_truth_r[, lc, drop = FALSE])) == 0L
    } else {
      truth_ok <- !is.na(X_truth_r[, lc[1]])
    }
    val_cells   <- val_cells & truth_ok
    val_row_idx <- which(val_cells)
    n_val       <- length(val_row_idx)

    if (n_val == 0) {
      calibrated_gates[lc] <- 0
      next
    }

    if (n_val < min_val_cells) {
      low_val_traits <- c(low_val_traits,
                          sprintf("%s (n=%d)", tm$name, n_val))
    }

    # ------------------------------------------------------------------
    # Select the half-A / half-B split(s) used for grid search.
    #
    # "single_split" (default, backward-compat): one deterministic split
    # keyed off `seed + 17L`. Fast; unstable when n_val is small (e.g.
    # 7 val cells -> 3-4 in half_a -> gate picks are essentially random
    # between 0 and the max grid value).
    #
    # "median_splits": repeat the whole grid + half-B-verify procedure
    # for each of `gate_splits_B` random splits, take the MEDIAN gate
    # across the B runs. Shrinks the bimodal-gate flip-flop observed at
    # small n_val without any retraining. Default B = 31 gives a cheap
    # variance reduction and is odd so the median is well-defined.
    # ------------------------------------------------------------------
    split_seeds <- if (gate_method == "single_split") {
      (seed + 17L)
    } else {
      as.integer(seed + 17L + seq_len(gate_splits_B) - 1L)
    }

    resolve_best_g_one_split <- function(ss) {
      set.seed(ss)
      perm_i <- sample(n_val)
      half_a_i <- val_row_idx[perm_i[seq_len(floor(n_val / 2))]]
      half_b_i <- val_row_idx[perm_i[(floor(n_val / 2) + 1L):n_val]]
      if (length(half_a_i) == 0L || length(half_b_i) == 0L) {
        half_a_i <- val_row_idx
        half_b_i <- val_row_idx
      }
      list(half_a = half_a_i, half_b = half_b_i)
    }

    split_pairs <- lapply(split_seeds, resolve_best_g_one_split)
    # For the single-split path we keep half_a / half_b as before; for
    # the median-splits path we rebuild them inside the loop below.
    half_a <- split_pairs[[1]]$half_a
    half_b <- split_pairs[[1]]$half_b

    # Helper: mean loss for gate g on a set of row indices.
    # Uses 0-1 loss for binary/categorical so that the val signal matches
    # the test metric (argmax accuracy), not cross-entropy.
    cal_mean_loss <- function(g, rows) {
      if (length(rows) == 0L) return(Inf)
      if (tm$type %in% c("continuous", "count", "ordinal", "proportion")) {
        pred_j <- (1 - g) * mu_cal[rows, lc[1]] +
                  g       * delta_cal[rows, lc[1]]
        mean((pred_j - X_truth_r[rows, lc[1]])^2)

      } else if (tm$type == "multi_proportion") {
        # MSE in CLR space, averaged across K components.
        pred_mat  <- (1 - g) * mu_cal[rows, lc, drop = FALSE] +
                     g       * delta_cal[rows, lc, drop = FALSE]
        truth_mat <- X_truth_r[rows, lc, drop = FALSE]
        mean((pred_mat - truth_mat)^2)

      } else if (tm$type == "binary") {
        pred_j     <- (1 - g) * mu_cal[rows, lc[1]] +
                      g       * delta_cal[rows, lc[1]]
        pred_class <- as.numeric(pred_j > 0)
        truth_j    <- X_truth_r[rows, lc[1]]
        mean(pred_class != truth_j)

      } else if (tm$type == "categorical") {
        logits      <- (1 - g) * mu_cal[rows, lc, drop = FALSE] +
                       g       * delta_cal[rows, lc, drop = FALSE]
        pred_class  <- max.col(logits,    ties.method = "first")
        truth_mat   <- X_truth_r[rows, lc, drop = FALSE]
        truth_class <- max.col(truth_mat, ties.method = "first")
        mean(pred_class != truth_class)

      } else if (tm$type == "zi_count") {
        # Coupled gate: single g applied to both gate and magnitude columns.
        gate_pred  <- (1 - g) * mu_cal[rows, lc[1]] +
                      g       * delta_cal[rows, lc[1]]
        mag_pred   <- (1 - g) * mu_cal[rows, lc[2]] +
                      g       * delta_cal[rows, lc[2]]
        p_nz       <- expit(gate_pred)
        count_hat  <- pmax(expm1(mag_pred * tm$sd + tm$mean), 0)
        pred_ev    <- p_nz * count_hat
        truth_gate <- X_truth_r[rows, lc[1]]
        truth_mag  <- X_truth_r[rows, lc[2]]
        truth_ev   <- rep(0, length(rows))
        nz         <- which(truth_gate > 0.5 & is.finite(truth_mag))
        truth_ev[nz] <- expm1(truth_mag[nz] * tm$sd + tm$mean)
        mean((pred_ev - truth_ev)^2)

      } else {
        Inf
      }
    }

    # Absolute minimum cell-level improvement floor for discrete traits
    # (relative gains on 0-1 loss are small: a 2% relative gain over 0.35
    # baseline is only ~0.007, easily noise on small val sets).
    is_discrete <- tm$type %in% c("binary", "categorical")

    # Inner: grid search on half_a + half_b verification, returns best_g.
    resolve_one_split <- function(ha, hb) {
      min_abs_a <- if (is_discrete) 2 / max(length(ha), 1L) else 0
      min_abs_b <- if (is_discrete) 1 / max(length(hb), 1L) else 0

      loss_a_0 <- cal_mean_loss(0, ha)
      best_g   <- 0
      best_la  <- loss_a_0
      for (g in gate_grid) {
        if (g == 0) next
        loss_a_g <- cal_mean_loss(g, ha)
        if (!is.finite(loss_a_g)) next
        rel      <- (loss_a_0 - loss_a_g) / max(loss_a_0, 1e-12)
        abs_gain <- loss_a_0 - loss_a_g
        if (rel >= cal_min_rel_gain && abs_gain >= min_abs_a &&
            loss_a_g < best_la) {
          best_g  <- g
          best_la <- loss_a_g
        }
      }

      if (best_g > 0) {
        loss_b_0 <- cal_mean_loss(0,      hb)
        loss_b_g <- cal_mean_loss(best_g, hb)
        rel_b    <- (loss_b_0 - loss_b_g) / max(loss_b_0, 1e-12)
        abs_b    <- loss_b_0 - loss_b_g
        if (!is.finite(rel_b) ||
            rel_b < (cal_min_rel_gain / 2) ||
            abs_b < min_abs_b) {
          best_g <- 0
        }
      }
      best_g
    }

    if (gate_method == "single_split") {
      best_g <- resolve_one_split(half_a, half_b)
    } else {
      # median_splits: run grid + verify for each of B splits, take median
      best_g_vec <- vapply(split_pairs,
                           function(sp) resolve_one_split(sp$half_a, sp$half_b),
                           numeric(1))
      best_g <- as.numeric(stats::median(best_g_vec))
    }

    calibrated_gates[lc] <- best_g
  }

  if (verbose) {
    gate_summary <- round(calibrated_gates, 3)
    names(gate_summary) <- if (!is.null(latent_names))
                             latent_names
                           else
                             paste0("col", seq_len(p))
    message("Calibrated gates: ",
            paste(names(gate_summary), gate_summary, sep = "=", collapse = ", "))
  }

  # Small-val warning: if any trait had fewer than `min_val_cells`
  # validation cells, the gate calibration AND the conformal score will
  # be noisy for that trait. Users see under-/over-coverage in that
  # regime; surface this up-front so they know to (a) increase the
  # fraction of held-out data, (b) use gate_method = "median_splits" +
  # conformal_method = "bootstrap" to smooth the per-split variance,
  # or (c) accept that intervals are approximate.
  if (length(low_val_traits) > 0L) {
    warning("Small validation set for ", length(low_val_traits),
            " trait(s): ", paste(low_val_traits, collapse = ", "),
            ". Calibrated gate and conformal scores will be noisy; ",
            "coverage may deviate from the 95%% target. See ",
            "`?fit_pigauto` under 'Calibration at small n' for ",
            "smoothing options.",
            call. = FALSE)
  }

  calibrated_gates
}

# ---------------------------------------------------------------------------
# compute_conformal_scores()
# ---------------------------------------------------------------------------
#
# Compute per-trait conformal prediction scores from validation residuals.
#
# Produces the (1 - alpha) conformal quantile of |truth - pred| on held-out
# validation cells for continuous, count, ordinal, and proportion traits.
# Used at prediction time to construct marginal conformal intervals with
# guaranteed (1 - alpha) coverage.
#
# Two estimators are supported:
#
# - `"split"` (default, backward compat): the classical Vovk conformal
#   quantile at level `ceil((1 - alpha)(n + 1))/n` on the val residuals.
#   Fast, exact finite-sample coverage guarantee under exchangeability,
#   but extremely noisy when `n_val` is small (~10 cells): at `n_val = 7`
#   the estimator is essentially `max(residuals)`, and a single seed
#   reshuffle can move it by an order of magnitude.
#
# - `"bootstrap"` (new): compute `B` bootstrap resamples of the val
#   residuals, take the split-conformal quantile on each resample, then
#   average (mean). This reduces the score's Monte-Carlo variance ~√B×
#   at the cost of a single `B`-sample loop per trait.  It loses the
#   exact exchangeability guarantee from split conformal but stays in the
#   same asymptotic regime.  Recommended when `n_val` per trait is < 30.
#
# @param trait_map         list of trait descriptors
# @param calibrated_gates  named numeric vector (length p)
# @param mu_cal            numeric matrix (n × p) — baseline predictions
# @param delta_cal         numeric matrix (n × p) — GNN predictions
# @param X_truth_r         numeric matrix (n × p) — truth with NAs
# @param val_mask_mat      logical matrix (n × p) — TRUE = validation cell
# @param alpha             numeric — miscoverage level; default 0.05 (→ 95%)
# @param method            `"split"` or `"bootstrap"`
# @param bootstrap_B       integer — bootstrap resamples when `method = "bootstrap"`
# @param verbose           logical
# @return named numeric vector; NA for discrete traits or empty val sets
compute_conformal_scores <- function(trait_map, calibrated_gates,
                                     mu_cal, delta_cal,
                                     X_truth_r, val_mask_mat,
                                     alpha = 0.05,
                                     method = c("split", "bootstrap"),
                                     bootstrap_B = 500L,
                                     verbose = FALSE) {
  method <- match.arg(method)
  p <- ncol(mu_cal)
  n <- nrow(mu_cal)

  # Blended predictions using calibrated gates
  pred_cal <- matrix(0, n, p)
  for (j in seq_len(p)) {
    pred_cal[, j] <- (1 - calibrated_gates[j]) * mu_cal[, j] +
                      calibrated_gates[j]       * delta_cal[, j]
  }

  conformal_scores <- rep(NA_real_, length(trait_map))
  names(conformal_scores) <- vapply(trait_map, "[[", character(1), "name")

  for (tm in trait_map) {
    if (!(tm$type %in% c("continuous", "count", "ordinal", "proportion"))) next
    lc        <- tm$latent_cols
    val_cells <- val_mask_mat[, lc[1]]
    if (sum(val_cells) == 0) next

    residuals <- abs(X_truth_r[val_cells, lc[1]] - pred_cal[val_cells, lc[1]])
    residuals <- residuals[is.finite(residuals)]
    n_val     <- length(residuals)
    if (n_val == 0L) next

    q_level <- min(ceiling((1 - alpha) * (n_val + 1)) / n_val, 1)

    if (method == "split") {
      conformal_scores[tm$name] <- as.numeric(
        stats::quantile(residuals, q_level, na.rm = TRUE)
      )
    } else {
      # bootstrap: average the split quantile over B resamples
      B <- as.integer(bootstrap_B)
      qs <- vapply(seq_len(B), function(b) {
        idx <- sample.int(n_val, n_val, replace = TRUE)
        as.numeric(stats::quantile(residuals[idx], q_level, na.rm = TRUE))
      }, numeric(1))
      conformal_scores[tm$name] <- mean(qs, na.rm = TRUE)
    }
  }

  if (verbose) {
    cs_print <- round(conformal_scores[!is.na(conformal_scores)], 4)
    message(sprintf("Conformal scores (latent scale, method=%s): %s",
                    method,
                    paste(names(cs_print), cs_print, sep = "=", collapse = ", ")))
  }

  conformal_scores
}
