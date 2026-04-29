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
                            safety_floor = FALSE,
                            mean_baseline_per_col = NULL,
                            simplex_step = 0.05,
                            min_val_cells = 10L,
                            seed = 1L,
                            latent_names = NULL,
                            verbose = FALSE) {
  gate_method <- match.arg(gate_method)
  if (safety_floor) {
    if (is.null(mean_baseline_per_col)) {
      stop("safety_floor = TRUE requires mean_baseline_per_col (numeric, length = ncol(mu_cal))")
    }
    stopifnot(length(mean_baseline_per_col) == ncol(mu_cal))
    simplex <- simplex_grid(step = simplex_step)
  } else {
    simplex <- NULL
  }
  p                <- ncol(mu_cal)
  calibrated_gates <- numeric(p)
  r_cal_bm_vec    <- numeric(p)
  r_cal_gnn_vec   <- numeric(p)
  r_cal_mean_vec  <- numeric(p)
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
      calibrated_gates[lc]    <- 0
      r_cal_bm_vec[lc]        <- 1   # pure BM fallback
      r_cal_gnn_vec[lc]       <- 0
      r_cal_mean_vec[lc]      <- 0
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

    # Helper: mean loss for gate g on a set of row indices.
    # Uses 0-1 loss for binary/categorical so that the val signal matches
    # the test metric (argmax accuracy), not cross-entropy.
    cal_mean_loss <- function(g, rows) {
      if (length(rows) == 0L) return(Inf)
      mean_j <- if (safety_floor) mean_baseline_per_col[lc[1]] else NA_real_
      # When safety_floor = TRUE, g is a length-3 vector (r_bm, r_gnn, r_mean).
      # When FALSE, g is a scalar r_gnn and the legacy (1-g)*mu + g*delta blend
      # applies.
      blend1 <- function(x_bm, x_delta, mean_scalar = mean_j) {
        if (safety_floor) {
          g[1] * x_bm + g[2] * x_delta + g[3] * mean_scalar
        } else {
          # g may be a scalar (legacy caller) OR a length-3 promoted vector
          # from cand_grid where g[1]=1-r_gnn, g[2]=r_gnn, g[3]=0.
          # Extract the GNN weight unambiguously.
          r_gnn <- if (length(g) == 1L) g else g[2L]
          (1 - r_gnn) * x_bm + r_gnn * x_delta
        }
      }
      if (tm$type %in% c("continuous", "count", "ordinal", "proportion")) {
        pred_j <- blend1(mu_cal[rows, lc[1]], delta_cal[rows, lc[1]])
        mean((pred_j - X_truth_r[rows, lc[1]])^2)

      } else if (tm$type == "multi_proportion") {
        # multi_proportion NOT supported in safety_floor (requires per-component
        # mean vector; out of scope for this spec).
        if (safety_floor) return(Inf)
        # Legacy path: g may be scalar OR length-3 (c(1-r, r, 0)) from cand_grid.
        # Extract r_gnn unambiguously to avoid broadcasting mismatch against n×K mat.
        r_gnn     <- if (length(g) == 1L) g else g[2L]
        pred_mat  <- (1 - r_gnn) * mu_cal[rows, lc, drop = FALSE] +
                     r_gnn       * delta_cal[rows, lc, drop = FALSE]
        truth_mat <- X_truth_r[rows, lc, drop = FALSE]
        mean((pred_mat - truth_mat)^2)

      } else if (tm$type == "binary") {
        pred_j     <- blend1(mu_cal[rows, lc[1]], delta_cal[rows, lc[1]])
        pred_class <- as.numeric(pred_j > 0)
        truth_j    <- X_truth_r[rows, lc[1]]
        mean(pred_class != truth_j)

      } else if (tm$type == "categorical") {
        logits <- matrix(0, nrow = length(rows), ncol = length(lc))
        for (kk in seq_along(lc)) {
          mean_k <- if (safety_floor) mean_baseline_per_col[lc[kk]] else NA_real_
          logits[, kk] <- blend1(mu_cal[rows, lc[kk]], delta_cal[rows, lc[kk]],
                                 mean_scalar = mean_k)
        }
        pred_class  <- max.col(logits,    ties.method = "first")
        truth_mat   <- X_truth_r[rows, lc, drop = FALSE]
        truth_class <- max.col(truth_mat, ties.method = "first")
        mean(pred_class != truth_class)

      } else if (tm$type == "zi_count") {
        gate_pred  <- blend1(mu_cal[rows, lc[1]], delta_cal[rows, lc[1]])
        mag_pred   <- blend1(mu_cal[rows, lc[2]], delta_cal[rows, lc[2]])
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

    # Inner: grid search on half_a + half_b verification, returns best weight vector.
    resolve_one_split <- function(ha, hb) {
      min_abs_a <- if (is_discrete) 2 / max(length(ha), 1L) else 0
      min_abs_b <- if (is_discrete) 1 / max(length(hb), 1L) else 0

      # Candidate grid is 3-column: for safety_floor use the simplex,
      # for legacy 1-D use the scalar gate_grid promoted to (1-g, g, 0).
      cand_grid <- if (safety_floor) simplex else {
        matrix(c(1 - gate_grid, gate_grid, rep(0, length(gate_grid))),
               ncol = 3L,
               dimnames = list(NULL, c("r_bm", "r_gnn", "r_mean")))
      }

      # "Pure BM" reference: (1,0,0) in simplex form when safety_floor = TRUE;
      # scalar 0 (i.e. g=0, no GNN) in legacy form when safety_floor = FALSE.
      # Both cal_mean_loss calls use the same dispatch path so no recycling occurs.
      ref_w          <- if (safety_floor) c(1, 0, 0) else 0
      loss_a_pure_bm <- cal_mean_loss(ref_w, ha)
      best_w         <- ref_w
      # Defensive: if loss_a_pure_bm is non-finite (e.g. multi_proportion
      # CLR latents with too few half-A cells for some component), seed
      # best_la with Inf so any finite candidate la beats it.  Same NA-
      # safety pattern as the half-B fix in commit 573decd.  Without
      # this, `la < best_la` errors with "missing value where TRUE/FALSE
      # needed" when best_la is NA.
      best_la        <- if (is.finite(loss_a_pure_bm)) loss_a_pure_bm else Inf
      pure_bm_is_finite <- is.finite(loss_a_pure_bm)
      for (ci in seq_len(nrow(cand_grid))) {
        w_try <- cand_grid[ci, ]
        la    <- cal_mean_loss(w_try, ha)
        if (!is.finite(la)) next
        if (!safety_floor) {
          # Legacy half-A filter (bit-identical to Task 3): skip candidates
          # that do not beat pure-BM by cal_min_rel_gain AND min_abs_a.
          r_gnn_try <- w_try[2L]
          if (r_gnn_try == 0) next
          # If pure-BM loss is non-finite, the legacy filter cannot
          # evaluate; admit the candidate (the half-B verification step
          # is the safety net).
          if (pure_bm_is_finite) {
            rel      <- (loss_a_pure_bm - la) / max(loss_a_pure_bm, 1e-12)
            abs_gain <- loss_a_pure_bm - la
            if (rel < cal_min_rel_gain || abs_gain < min_abs_a) next
          }
        }
        if (la < best_la) {
          best_la <- la
          best_w  <- w_try
        }
      }

      # Half-B verification: require the winning w's half-B loss to also
      # beat the pure-BM half-B loss by cal_min_rel_gain (legacy: /2) AND
      # by the discrete absolute cell floor.
      loss_b_pure_bm      <- cal_mean_loss(ref_w,    hb)
      loss_b_best         <- cal_mean_loss(best_w,   hb)
      rel_gain_b          <- (loss_b_pure_bm - loss_b_best) / max(loss_b_pure_bm, 1e-12)
      abs_gain_b          <- loss_b_pure_bm - loss_b_best
      rel_gain_threshold  <- if (safety_floor) cal_min_rel_gain else cal_min_rel_gain / 2
      # Defensive: if either half-B loss is NA (e.g. a multi_proportion CLR
      # column with too few cells in half-B for some component), the gain
      # cannot be evaluated.  Treat as "verification failed" -> fall back
      # to the safe reference, preserving the safety guarantee.
      verification_failed <-
        !is.finite(rel_gain_b) ||
        rel_gain_b < rel_gain_threshold ||
        (is_discrete && (!is.finite(abs_gain_b) || abs_gain_b < min_abs_b))
      if (verification_failed) {
        # Revert to best safe fallback.  When safety_floor = TRUE the pure-mean
        # point (0,0,1) is always a valid option — pick whichever of BM and mean
        # is better on half-B so the safety guarantee is preserved.
        if (safety_floor) {
          mean_ref_w     <- c(0, 0, 1)
          loss_b_mean    <- cal_mean_loss(mean_ref_w, hb)
          # If either loss is NA, prefer the BM ref (safer than mean
          # for continuous-like latents); otherwise pick whichever is
          # smaller on half-B.
          best_w <- if (is.finite(loss_b_mean) && is.finite(loss_b_pure_bm)) {
            if (loss_b_mean <= loss_b_pure_bm) mean_ref_w else ref_w
          } else {
            ref_w
          }
        } else {
          best_w <- ref_w
        }
      }

      # Always return a length-3 vector so vapply(., numeric(3L)) works uniformly.
      if (length(best_w) == 1L) c(1 - best_w, best_w, 0) else best_w
    }

    best_w_across_splits <- vapply(split_pairs, function(pair)
      resolve_one_split(pair$half_a, pair$half_b),
      numeric(3L))    # 3 rows x B cols
    w_final <- apply(best_w_across_splits, 1L, stats::median)
    # Guard against degenerate median = c(0, 0, 0), which can occur when
    # adversarial splits alternate between corners (e.g. >50% of splits pick
    # (1,0,0), >50% pick (0,0,1), and each coordinate's median is 0). Pure-mean
    # (0,0,1) is the safe fallback matching the invariant guarantee.
    if (!is.finite(sum(w_final)) || sum(w_final) < 1e-10) {
      w_final <- if (safety_floor) c(0, 0, 1) else c(1, 0, 0)
    } else {
      w_final <- w_final / sum(w_final)   # renormalise (medians don't preserve sum)
    }

    # safety_floor guarantee (continuous family, byte-for-byte
    # equivalent to the pre-2026-04-29 implementation): the final blend
    # on the full val set must not be worse than the pure-mean baseline.
    # If median-of-splits picked a blend that loses to mean on the full
    # val set, fall back to pure mean. Pure-BM is trusted via the
    # half-A/half-B + median cross-check above; checking blend vs
    # pure-BM on the full val set would be too noisy at typical val
    # sizes (~5 % sampling SD on 30 %-MAR / n=1500) and spuriously
    # overrides the blend on high-phylo-signal continuous traits.
    if (safety_floor && length(val_row_idx) > 0L &&
        tm$type %in% c("continuous", "count", "ordinal", "proportion")) {
      j1 <- lc[1]
      mean_j1 <- mean_baseline_per_col[j1]
      blend_full <- w_final[1] * mu_cal[val_row_idx, j1] +
                    w_final[2] * delta_cal[val_row_idx, j1] +
                    w_final[3] * mean_j1
      mean_full  <- rep(mean_j1, length(val_row_idx))
      truth_full <- X_truth_r[val_row_idx, j1]
      ok         <- is.finite(truth_full)
      if (any(ok)) {
        mse_blend <- mean((blend_full[ok] - truth_full[ok])^2)
        mse_mean  <- mean((mean_full[ok]  - truth_full[ok])^2)
        if (mse_blend > mse_mean + 1e-12) {
          w_final <- c(0, 0, 1)
        }
      }
    }

    # Strict val-floor (DISCRETE only — Tier-1 fix 2026-04-29):
    # binary, categorical, zi_count types previously bypassed any final
    # invariant check, which let through non-corner blends that lost to
    # pure-BM at the cell level (binary -5..-9 pp, categorical -5..-12 pp,
    # zi_count +5..+23 % RMSE in the 2026-04-28 benches). For these
    # types we apply a strict full-val check using `cal_mean_loss`
    # (0-1 loss for binary/categorical, gate+magnitude composite for
    # zi_count): if the blend's val-loss exceeds either pure-BM or
    # pure-MEAN by 1e-12, override to whichever pure corner is better.
    # Discrete losses are step functions over a small set of integers,
    # so the strict check is well-conditioned and the val→test
    # extrapolation noise that hurts the same check on continuous
    # types (see e.g. AVONET Mass) is much smaller here.
    if (length(val_row_idx) > 0L &&
        tm$type %in% c("binary", "categorical", "zi_count")) {
      lb <- cal_mean_loss(w_final, val_row_idx)
      bm_corner <- if (safety_floor) c(1, 0, 0) else 0
      lm_bm     <- cal_mean_loss(bm_corner, val_row_idx)
      lm_mean   <- if (safety_floor) cal_mean_loss(c(0, 0, 1), val_row_idx) else Inf
      coerce_loss <- function(x) if (is.finite(x)) x else Inf
      lb     <- coerce_loss(lb)
      lm_bm  <- coerce_loss(lm_bm)
      lm_mean <- coerce_loss(lm_mean)
      tol <- 1e-12
      if (lb > lm_bm + tol || lb > lm_mean + tol) {
        w_final <- if (lm_bm <= lm_mean) c(1, 0, 0) else c(0, 0, 1)
      }
    }

    calibrated_gates[lc] <- w_final[2]   # legacy scalar = r_gnn
    r_cal_bm_vec[lc]     <- w_final[1]
    r_cal_gnn_vec[lc]    <- w_final[2]
    r_cal_mean_vec[lc]   <- w_final[3]
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

  if (!is.null(latent_names)) {
    if (length(latent_names) == length(calibrated_gates)) {
      names(calibrated_gates) <- latent_names
    }
    if (length(latent_names) == length(r_cal_bm_vec)) {
      names(r_cal_bm_vec)   <- latent_names
      names(r_cal_gnn_vec)  <- latent_names
      names(r_cal_mean_vec) <- latent_names
    }
  }
  list(r_cal_bm = r_cal_bm_vec,
        r_cal_gnn = r_cal_gnn_vec,
        r_cal_mean = r_cal_mean_vec)
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
