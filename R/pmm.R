# Phase G' (2026-05-01): Predictive Mean Matching (PMM) for pigauto.
#
# Standard imputation technique (Little 1988; Buuren & Groothuis-Oudshoorn
# 2011 mice).  After the BM + GNN + blend pipeline produces a continuous
# prediction for each cell, PMM replaces the imputed value with an
# observed value from a donor selected by predicted-value proximity.
# Imputations are guaranteed to lie within the observed data range
# (no extrapolation possible by construction).
#
# Targets the AVONET Mass tail-extrapolation mode documented in
# useful/MEMO_2026-05-01_avonet_mass_diag.md and continues from Phase G's
# clamp_outliers band-aid.  Whereas clamp_outliers caps at obs_max *
# clamp_factor, PMM picks an actual observed value -- principled, no
# free heuristic parameter.

# pmm_impute_one_trait
# --------------------
# Per-trait, per-draw PMM.  For each missing cell (NA in `truth`), find
# the K observed cells with predicted values closest to the missing
# cell's predicted value, sample one, and return that donor's observed
# truth value.  Observed cells are passed through unchanged.
#
# Phase G'' (2026-05-01): added `when` argument to scope PMM to
# extrapolating predictions only.  When predictions are within
# [obs_min, obs_max], retaining the prediction is more accurate than
# snapping to a donor (the AVONET seed-2032 Mass regression of +805 %
# in PR #61 documents this).  Default `when = "outside_observed"`
# applies PMM only as a safety net for extrapolation; opt in to
# `when = "always"` for the mice-style "PMM for all imputations"
# behaviour (used in MI workflows where between-draw donor variance
# matters).
#
# @param predictions numeric vector length n.  Per-cell predictions
#   (post-back-transform, original units).
# @param truth numeric vector length n.  Original-units observed values
#   with NA at cells to impute.
# @param K integer >= 1.  Donor pool size.  Reduced to length(observed)
#   when fewer than K observed cells exist.
# @param seed integer.  Per-call random seed for the within-K donor
#   choice.  Use distinct seeds across multi-imputation draws to get
#   different donors (preserving between-draw variance).
# @param when character.  When to apply PMM.  \code{"outside_observed"}
#   (default) only swaps cells whose prediction is < min(observed)
#   or > max(observed); cells whose prediction is in-range pass through
#   unchanged.  \code{"always"} swaps every missing cell (the original
#   mice convention).
# @return numeric vector length n with NA cells filled by donor values
#   (or by passthrough predictions when \code{when = "outside_observed"}
#   and the prediction is in-range).  Whenever PMM triggers, the donor
#   is an observed value from \code{truth}.
# @keywords internal
pmm_impute_one_trait <- function(predictions, truth, K = 5L, seed = 1L,
                                  when = c("outside_observed", "always")) {
  stopifnot(length(predictions) == length(truth))
  when <- match.arg(when)
  K <- as.integer(K)
  if (K < 1L) stop("'K' must be >= 1.", call. = FALSE)

  obs_idx  <- which(!is.na(truth))
  miss_idx <- which(is.na(truth))

  # No donors or no missing cells -> no-op
  if (length(obs_idx) == 0L || length(miss_idx) == 0L) {
    return(truth)
  }

  K_eff <- min(K, length(obs_idx))
  obs_pred  <- predictions[obs_idx]
  obs_truth <- truth[obs_idx]

  # Phase G'' triggers PMM only when prediction is outside the
  # OBSERVED range.  obs_truth is the ground-truth donor pool;
  # obs_min / obs_max define the in-range window.
  obs_min <- min(obs_truth, na.rm = TRUE)
  obs_max <- max(obs_truth, na.rm = TRUE)

  result <- truth
  rng_state <- .Random.seed_or_null()
  on.exit(restore_rng_state(rng_state), add = TRUE)
  set.seed(seed)

  for (i in miss_idx) {
    pred_i <- predictions[i]
    if (!is.finite(pred_i)) {
      # Pred is NA/NaN/Inf: PMM cannot match.  Leave as NA -- caller can
      # decide whether to fall back.
      next
    }
    if (identical(when, "outside_observed") &&
        pred_i >= obs_min && pred_i <= obs_max) {
      # Phase G'': prediction is in observed range; trust it.
      result[i] <- pred_i
      next
    }
    # K-nearest by predicted-value distance.  order() gives stable
    # tie-breaking (first-occurrence index wins).
    dists   <- abs(obs_pred - pred_i)
    nearest <- order(dists)[seq_len(K_eff)]
    # Pick one uniformly at random from the K nearest
    if (K_eff == 1L) {
      donor <- nearest
    } else {
      donor <- nearest[sample.int(K_eff, 1L)]
    }
    result[i] <- obs_truth[donor]
  }
  result
}

# pmm_eligible_types
# ------------------
# Returns the set of trait types where PMM should be applied.
# Currently: log-transformed continuous, count, and zi_count magnitude
# (where the back-transform amplifies tail overshoots), plus proportion
# (mild benefit -- preserves observed proportion shape).  Discrete-class
# types (binary, categorical, ordinal) are no-op since their decode is
# already from observed levels by construction.  multi_proportion: no-op
# (softmax is bounded; donor matching across the K-vector is awkward).
#
# @keywords internal
pmm_eligible_types <- function() {
  c("continuous", "count", "proportion", "zi_count")
}

pmm_is_eligible <- function(tm) {
  if (identical(tm$type, "continuous")) {
    return(isTRUE(tm$log_transform))
  }
  tm$type %in% pmm_eligible_types()
}

# .Random.seed_or_null / restore_rng_state
# ----------------------------------------
# Save and restore the RNG state so PMM's seed-setting doesn't leak
# upstream.  R's set.seed has global effect; saving + restoring around
# our per-trait calls keeps the user's RNG state untouched.
# @keywords internal
.Random.seed_or_null <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv())
  } else {
    NULL
  }
}

restore_rng_state <- function(state) {
  if (is.null(state)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(list = ".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
  invisible(NULL)
}

# apply_pmm_to_decoded
# --------------------
# Apply PMM per-trait per-draw to a list of decoded predictions.
# Modifies decoded$imputed[[trait]] in place for each eligible trait.
# Discrete/non-eligible traits are passed through.
#
# @param decode_results list of M decoded results, each with $imputed
#   data.frame.
# @param trait_map list of trait_map entries.
# @param X_orig n x n_traits original-units data.frame.  NA at
#   originally-missing cells (PMM imputes those; observed cells pass
#   through).  This is recovered from object$X_scaled at the call site.
# @param K integer donor pool size.
# @param base_seed integer.  Per-draw seed = base_seed + draw_idx.
# @return list of M decoded results, with PMM applied.
# @keywords internal
apply_pmm_to_decoded <- function(decode_results, trait_map, X_orig,
                                  K = 5L, base_seed = 1L,
                                  when = c("outside_observed", "always")) {
  when <- match.arg(when)
  K <- as.integer(K)
  M <- length(decode_results)
  for (m in seq_len(M)) {
    for (tm in trait_map) {
      if (!pmm_is_eligible(tm)) next
      nm <- tm$name
      preds <- as.numeric(decode_results[[m]]$imputed[[nm]])
      truth <- as.numeric(X_orig[[nm]])
      if (length(preds) != length(truth)) {
        next  # shape mismatch -- skip rather than error (multi-obs edge)
      }
      decode_results[[m]]$imputed[[nm]] <-
        pmm_impute_one_trait(preds, truth,
                              K = K, seed = base_seed + m, when = when)
      # Coerce back to integer for count and zi_count
      if (tm$type == "count") {
        decode_results[[m]]$imputed[[nm]] <-
          as.integer(decode_results[[m]]$imputed[[nm]])
      } else if (tm$type == "zi_count") {
        decode_results[[m]]$imputed[[nm]] <-
          as.integer(decode_results[[m]]$imputed[[nm]])
      }
    }
  }
  decode_results
}

# recover_X_orig
# --------------
# Reconstruct original-units observed data from object$X_scaled (z-scored,
# possibly log-transformed) using the trait_map normalisation params.
# Returns a data.frame n_obs x n_traits.  Observed cells get original
# values; missing cells stay NA.  Used by predict.pigauto_fit() to give
# PMM access to the donor pool.
#
# @keywords internal
recover_X_orig <- function(X_scaled, trait_map) {
  if (is.null(X_scaled)) return(NULL)
  n <- nrow(X_scaled)
  out <- vector("list", length(trait_map))
  out_names <- character(length(trait_map))
  for (i in seq_along(trait_map)) {
    tm <- trait_map[[i]]
    out_names[i] <- tm$name
    lc <- tm$latent_cols
    if (tm$type == "continuous") {
      vals <- X_scaled[, lc[1]] * tm$sd + tm$mean
      if (isTRUE(tm$log_transform)) vals <- exp(vals)
      out[[i]] <- vals
    } else if (tm$type == "count") {
      vals <- X_scaled[, lc[1]] * tm$sd + tm$mean
      vals <- expm1(vals)
      vals[!is.na(vals)] <- pmax(round(vals[!is.na(vals)]), 0)
      out[[i]] <- as.integer(vals)
    } else if (tm$type == "proportion") {
      vals <- X_scaled[, lc[1]] * tm$sd + tm$mean
      out[[i]] <- stats::plogis(vals)
    } else if (tm$type == "zi_count") {
      gate <- X_scaled[, lc[1]]
      mag  <- X_scaled[, lc[2]] * tm$sd + tm$mean
      mag  <- expm1(mag)
      vals <- ifelse(is.na(gate) | gate == 0, 0,
                     pmax(round(mag), 0))
      # NA only when gate is NA AND magnitude is NA (preserve ambiguity)
      vals[is.na(gate)] <- NA_integer_
      out[[i]] <- as.integer(vals)
    } else {
      # Other types (binary, categorical, ordinal, multi_proportion) are
      # not PMM-eligible; we don't need to recover originals for them.
      out[[i]] <- rep(NA, n)
    }
  }
  df <- as.data.frame(out, col.names = out_names, stringsAsFactors = FALSE)
  rownames(df) <- rownames(X_scaled)
  df
}
