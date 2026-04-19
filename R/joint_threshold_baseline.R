#' Aggregate a multi-obs X_scaled matrix to species-level with per-type semantics
#'
#' For multi-obs pigauto_data, the Level-C helpers need species-level input.
#' Collapses obs rows per species using the appropriate aggregator per trait
#' type: mean for continuous-family, threshold-at-0.5 for binary/ZI-gate,
#' argmax-one-hot for categorical. Single-obs data passes through unchanged.
#'
#' Splits are masked in obs space first, then re-derived in species space:
#' a species cell is held-out iff any of its obs cells was held-out.
#'
#' @param data pigauto_data (single- or multi-obs).
#' @param splits optional output of make_missing_splits (obs-level).
#' @return list(X_species, splits_species, species_names).
#' @keywords internal
#' @noRd
aggregate_to_species <- function(data, splits = NULL, soft_aggregate = FALSE) {
  if (!isTRUE(data$multi_obs)) {
    # Single-obs: still need to apply split masking to match prior behaviour
    # of build_liability_matrix(), which masked val/test cells up front.
    X_single <- data$X_scaled
    if (!is.null(splits)) {
      X_single[splits$val_idx]  <- NA
      X_single[splits$test_idx] <- NA
    }
    return(list(X_species        = X_single,
                splits_species   = splits,
                species_names    = data$species_names %||% rownames(data$X_scaled),
                is_proportion_col = rep(FALSE, ncol(X_single))))
  }

  trait_map <- data$trait_map
  X         <- data$X_scaled
  n_obs     <- data$n_obs
  n_species <- data$n_species
  spp       <- data$species_names
  obs_to_sp <- data$obs_to_species

  # Apply split mask first (in obs space) so aggregation propagates NA correctly
  obs_original <- !is.na(X)
  if (!is.null(splits)) {
    X[splits$val_idx]  <- NA
    X[splits$test_idx] <- NA
  }

  # Per-column trait-type map (zi_count expands: gate=binary, mag=continuous)
  type_by_col <- rep(NA_character_, ncol(X))
  for (tm in trait_map) {
    if (tm$type == "zi_count") {
      type_by_col[tm$latent_cols[1]] <- "binary"
      type_by_col[tm$latent_cols[2]] <- "continuous"
    } else {
      type_by_col[tm$latent_cols] <- tm$type
    }
  }

  X_species <- matrix(NA_real_, nrow = n_species, ncol = ncol(X),
                      dimnames = list(spp, colnames(X)))

  is_proportion_col <- rep(FALSE, ncol(X))

  for (j in seq_len(ncol(X))) {
    tp <- type_by_col[j]
    sp_vals <- tapply(X[, j], obs_to_sp, function(v) {
      v <- v[!is.na(v)]
      if (length(v) == 0L) NA_real_ else mean(v)
    })
    # Reorder to spp order via integer index
    sp_vals <- sp_vals[as.character(seq_len(n_species))]
    if (tp == "binary") {
      if (soft_aggregate) {
        X_species[, j] <- unname(sp_vals)  # keep proportion
        is_proportion_col[j] <- TRUE
      } else {
        X_species[, j] <- ifelse(is.na(sp_vals), NA_real_,
                                  as.numeric(sp_vals >= 0.5))
      }
    } else {
      X_species[, j] <- unname(sp_vals)
    }
  }

  # Re-enforce one-hot for categorical via argmax (or preserve proportions in soft mode)
  for (tm in trait_map) {
    if (tm$type != "categorical") next
    k_cols <- tm$latent_cols
    if (soft_aggregate) {
      for (kc in k_cols) is_proportion_col[kc] <- TRUE
      next  # preserve proportions
    }
    for (s in seq_len(n_species)) {
      props <- X_species[s, k_cols]
      if (all(is.na(props))) next
      winner <- which.max(props)
      X_species[s, k_cols] <- 0
      X_species[s, k_cols[winner]] <- 1
    }
  }

  # Re-derive splits in species space
  splits_species <- NULL
  if (!is.null(splits)) {
    sp_originally_observed <- matrix(FALSE, nrow = n_species, ncol = ncol(X))
    for (s in seq_len(n_species)) {
      rows <- which(obs_to_sp == s)
      sp_originally_observed[s, ] <- apply(obs_original[rows, , drop = FALSE], 2, any)
    }
    sp_held  <- sp_originally_observed & is.na(X_species)
    held_lin <- which(sp_held)
    n_val  <- length(splits$val_idx)
    n_test <- length(splits$test_idx)
    if (length(held_lin) > 0L && (n_val + n_test) > 0L) {
      n_val_sp <- round(length(held_lin) * n_val / (n_val + n_test))
      val_sp  <- held_lin[seq_len(n_val_sp)]
      test_sp <- held_lin[setdiff(seq_along(held_lin), seq_len(n_val_sp))]
    } else {
      val_sp  <- integer(0)
      test_sp <- integer(0)
    }
    mask_species <- !is.na(X_species)
    splits_species <- list(val_idx  = val_sp,
                           test_idx = test_sp,
                           mask     = mask_species)
  }

  list(X_species        = X_species,
       splits_species   = splits_species,
       species_names    = spp,
       is_proportion_col = is_proportion_col)
}

#' Assemble liability-scale input matrix for the threshold-model joint baseline
#'
#' For each trait type:
#'   - continuous/count/ordinal/proportion: pass through the z-scored latent value
#'     directly (it IS the liability on the preprocessed scale).
#'   - binary: apply truncated-Gaussian E-step with a vague N(0,1) prior so
#'     observed y=1 cells get a positive posterior mean and y=0 cells a
#'     negative one.
#'   - zi_count gate: handled as binary (col 1 of the zi_count latent pair).
#'
#' Cells that are missing (NA in X_scaled) or held out by `splits` are set to
#' NA so Rphylopars marginalises over them.
#'
#' Categorical traits handled separately by fit_ovr_categorical_fits().
#' Skipped here — not included in the joint liability matrix.
#'
#' multi_proportion columns are also skipped — they have the same
#' rank-(K-1) phylopars instability as categorical K-col groups.
#'
#' Multi-obs datasets are aggregated to species-level via
#' `aggregate_to_species()` at the entry point; downstream logic is
#' single-obs from there on.
#'
#' @param data pigauto_data from preprocess_traits.
#' @param splits output of make_missing_splits() or NULL.
#' @return list(X_liab, liab_cols, liab_types) where X_liab is
#'   n_species x length(liab_cols), liab_cols is the integer index into the
#'   full p_latent space, and liab_types is a character vector aligned with
#'   liab_cols.
#' @keywords internal
#' @noRd
build_liability_matrix <- function(data, splits = NULL, soft_aggregate = FALSE,
                                   sd_prior_vec = NULL,
                                   mu_prior_mat = NULL,
                                   sd_prior_mat = NULL) {
  trait_map <- data$trait_map
  # Prior lookup precedence (highest first):
  #   1. mu_prior_mat / sd_prior_mat (Phase 7, per-cell, n x K_liab)
  #   2. sd_prior_vec (Phase 6, per-column, length K_liab)
  #   3. NULL default — plug-in N(0, 1) byte-identical to v0.9.1.
  use_mat <- !is.null(mu_prior_mat) || !is.null(sd_prior_mat)
  mu_i_j_lookup <- function(i, j) {
    if (!is.null(mu_prior_mat)) mu_prior_mat[i, j] else 0
  }
  sd_i_j_lookup <- function(i, j) {
    if (!is.null(sd_prior_mat)) return(sd_prior_mat[i, j])
    if (!is.null(sd_prior_vec)) return(sd_prior_vec[j])
    1
  }
  # Collapse multi-obs to species level (no-op for single-obs). Splits are
  # masked internally in obs space before aggregation, so no separate
  # masking step is needed below.
  agg    <- aggregate_to_species(data, splits = splits, soft_aggregate = soft_aggregate)
  X      <- agg$X_species
  splits <- agg$splits_species
  n      <- nrow(X)

  # is_proportion_col: TRUE for cols where the aggregated value is a proportion
  # (soft binary or soft categorical), not a hard 0/1 or one-hot.
  is_prop <- agg$is_proportion_col
  # Overlay explicit proportion_cols if present on data (used by OVR path)
  if (!is.null(data$proportion_cols)) {
    is_prop[data$proportion_cols] <- TRUE
  }

  liab_cols  <- integer(0)
  liab_types <- character(0)
  liab_tms   <- list()
  for (tm in trait_map) {
    if (tm$type %in% c("continuous", "count", "ordinal", "proportion")) {
      liab_cols  <- c(liab_cols, tm$latent_cols)
      liab_types <- c(liab_types, tm$type)
      liab_tms   <- c(liab_tms, list(tm))
    } else if (tm$type == "binary") {
      liab_cols  <- c(liab_cols, tm$latent_cols)
      liab_types <- c(liab_types, "binary")
      liab_tms   <- c(liab_tms, list(tm))
    } else if (tm$type == "zi_count") {
      # Phase 5: ZI gate col joins the joint as a binary-like liability.
      # The magnitude col (col 2) is already in bm_cols at the caller side.
      liab_cols  <- c(liab_cols, tm$latent_cols[1])
      liab_types <- c(liab_types, "binary")         # treat as binary for E-step
      liab_tms   <- c(liab_tms, list(list(type = "binary")))
    }
    # Categorical traits handled separately by fit_ovr_categorical_fits().
    # Skip here — do not include categorical cols in the joint liability
    # matrix.
    # multi_proportion: also skipped (same phylopars rank-(K-1) instability
    # as categorical K-col groups).
  }

  X_liab <- matrix(NA_real_, nrow = n, ncol = length(liab_cols))
  colnames(X_liab) <- colnames(X)[liab_cols]
  rownames(X_liab) <- rownames(X)

  for (j in seq_along(liab_cols)) {
    col_idx <- liab_cols[j]
    tp      <- liab_types[j]
    tm_j    <- liab_tms[[j]]
    src_col <- X[, col_idx]

    if (tp == "binary") {
      # Observed value is either hard 0/1 (hard aggregation or single-obs) or
      # a species-level proportion in [0,1] (soft aggregation for multi-obs).
      # Dispatch to soft E-step when is_prop indicates a proportion value.
      for (i in seq_len(n)) {
        v <- src_col[i]
        if (is.na(v)) next
        mu_i <- mu_i_j_lookup(i, j); sd_i <- sd_i_j_lookup(i, j)
        if (is_prop[col_idx]) {
          post <- estep_liability_binary_soft(p = as.numeric(v),
                                               mu_prior = mu_i, sd_prior = sd_i)
        } else {
          post <- estep_liability(tm_j, observed = v,
                                   mu_prior = mu_i, sd_prior = sd_i)
        }
        X_liab[i, j] <- post$mean
      }

    } else if (tp == "ordinal") {
      # Ordinal: interval-truncated Gaussian E-step via the dispatcher.
      # estep_liability() handles z-score roundtrip internally.
      for (i in seq_len(n)) {
        v <- src_col[i]
        if (is.na(v)) next
        mu_i <- mu_i_j_lookup(i, j); sd_i <- sd_i_j_lookup(i, j)
        post <- estep_liability(tm_j, observed = v,
                                  mu_prior = mu_i, sd_prior = sd_i)
        X_liab[i, j] <- post$mean
      }

    } else {
      # Continuous-family (continuous/count/proportion): pass through
      X_liab[, j] <- src_col
    }
  }

  list(X_liab = X_liab, liab_cols = liab_cols, liab_types = liab_types)
}

#' Fit a joint MVN BM on liability-scale columns (Phase 3)
#'
#' Runs `build_liability_matrix()` then delegates to `Rphylopars::phylopars()`.
#' Returns raw liability-scale posterior. Decoding to per-type output scale
#' happens upstream in the caller / in `fit_baseline()` glue.
#'
#' @inheritParams fit_joint_mvn_baseline
#' @return list(mu_liab, se_liab, liab_cols, liab_types).
#' @keywords internal
#' @noRd
fit_joint_threshold_baseline <- function(data, tree, splits, graph = NULL,
                                        soft_aggregate = FALSE,
                                        sd_prior_vec = NULL) {
  stopifnot(joint_mvn_available())

  built <- build_liability_matrix(data, splits = splits,
                                  soft_aggregate = soft_aggregate,
                                  sd_prior_vec = sd_prior_vec)
  X_liab     <- built$X_liab
  liab_cols  <- built$liab_cols
  liab_types <- built$liab_types

  spp <- if (!is.null(data$species_names)) data$species_names else rownames(data$X_scaled)

  # Phylopars needs >= 2 non-NA observations per column; treat anything fewer
  # the same as all-NA: excluded from the fit, left as NA in the output.
  has_obs <- apply(X_liab, 2, function(col) sum(!is.na(col)) >= 2L)
  fit_cols <- which(has_obs)

  n_species  <- nrow(X_liab)
  n_all_cols <- ncol(X_liab)
  mu_liab    <- matrix(NA_real_, nrow = n_species, ncol = n_all_cols,
                        dimnames = list(spp, colnames(X_liab)))
  se_liab    <- matrix(NA_real_, nrow = n_species, ncol = n_all_cols,
                        dimnames = list(spp, colnames(X_liab)))

  phylopars_fit <- NULL
  if (length(fit_cols) >= 1L) {
    X_fit <- X_liab[, fit_cols, drop = FALSE]

    # Rphylopars rejects column names with '=' (e.g. "z=A"). Sanitise to
    # safe names for the Rphylopars call and restore originals afterwards.
    orig_colnames   <- colnames(X_fit)
    safe_colnames   <- make.names(orig_colnames, unique = TRUE)
    colnames(X_fit) <- safe_colnames

    df_in          <- as.data.frame(X_fit)
    df_in$species  <- spp
    df_in          <- df_in[, c("species", safe_colnames), drop = FALSE]

    fit <- Rphylopars::phylopars(
      trait_data = df_in,
      tree       = tree,
      model      = "BM"
    )

    tip_rows <- match(spp, rownames(fit$anc_recon))
    mu_fit   <- fit$anc_recon[tip_rows, , drop = FALSE]
    se_fit   <- sqrt(fit$anc_var[tip_rows, , drop = FALSE])

    # Validate shape
    if (ncol(mu_fit) != length(fit_cols)) {
      stop("fit_joint_threshold_baseline: Rphylopars returned ",
           ncol(mu_fit), " columns but ", length(fit_cols),
           " were passed. Column alignment would be ambiguous.",
           call. = FALSE)
    }

    # Restore original column names before writing back to output matrices
    colnames(mu_fit) <- orig_colnames
    colnames(se_fit) <- orig_colnames

    mu_liab[, fit_cols] <- mu_fit
    se_liab[, fit_cols] <- se_fit
    phylopars_fit <- fit
  }

  # phylopars_fit and fit_cols_idx are additive fields used by Phase 6 EM to
  # extract the per-column BM rate for the next iteration's prior. Existing
  # callers that don't look for them are unaffected.
  list(mu_liab      = mu_liab,
       se_liab      = se_liab,
       liab_cols    = liab_cols,
       liab_types   = liab_types,
       phylopars_fit = phylopars_fit,
       fit_cols_idx = fit_cols)
}

#' Decode liability-scale posterior to logit(P(y=1)) for binary traits
#'
#' Given posterior N(mu_liab, se_liab^2) for a latent liability L and the
#' threshold model y = 1 iff L > 0, the marginal probability is
#'   P(y=1) = pnorm(mu_liab / sqrt(1 + se_liab^2)).
#' We return the logit of this clipped to [0.01, 0.99] so downstream code
#' (GNN blending, BCE loss) stays numerically stable. Matches the clip used
#' by the label-propagation path in `fit_baseline()`.
#'
#' @param mu_liab numeric (possibly vectorised) posterior mean on liability scale.
#' @param se_liab numeric (same length) posterior SD on liability scale.
#' @return list(p, mu_logit) each the same length as the inputs.
#' @keywords internal
#' @noRd
decode_binary_liability <- function(mu_liab, se_liab) {
  sigma_marg <- sqrt(1 + se_liab^2)
  p <- stats::pnorm(mu_liab / sigma_marg)
  p <- pmin(pmax(p, 0.01), 0.99)
  list(p = p, mu_logit = stats::qlogis(p))
}

#' Decode ordinal liability posterior to z-scored integer class
#'
#' Given liability posterior N(mu_liab, se_liab^2) and ordinal thresholds,
#' find the most likely class and return its z-scored integer representation.
#' Matches the existing z-scored integer convention in X_scaled so the GNN's
#' MSE loss target is preserved.
#'
#' @param mu_liab numeric vector, posterior mean on liability scale.
#' @param se_liab numeric vector (same length), posterior SD.
#' @param tm trait_map entry for this ordinal trait.
#' @return list(mu_z, se_z) z-scored integer values.
#' @keywords internal
#' @noRd
decode_ordinal_liability <- function(mu_liab, se_liab, tm) {
  info <- liability_info(tm)
  thresholds <- info$thresholds
  K <- length(thresholds) + 1L
  # Most likely class: which interval does the posterior mean fall in?
  class_idx <- findInterval(mu_liab, thresholds) + 1L
  class_idx <- pmax(1L, pmin(K, class_idx))
  # Convert to z-scored integer: (class_0indexed - mean) / sd
  mu_z <- (class_idx - 1 - tm$mean) / tm$sd
  list(mu_z = mu_z, se_z = rep(0, length(mu_liab)))
}

# Extract per-column BM rate variance from an Rphylopars fit.
# Rphylopars stores the BM rate matrix in `$pars$phylocov`; some older
# builds used `$pars$sigma2`. Both are K x K matrices (or 1x1 for a single
# column). For a K-col joint fit we return the K diagonal elements
# (off-diagonals are the phylogenetic covariances, which Phase 6 does NOT
# use). For scalar sigma2 back-compat, handle that too.
# Phase 6 EM uses this to build the next iteration's sd_prior_vec via
# `sqrt(pmax(extract_liability_variances(...), tiny))`.
#' @keywords internal
#' @noRd
extract_liability_variances <- function(phylopars_fit) {
  Sigma <- phylopars_fit$pars$phylocov
  if (is.null(Sigma)) Sigma <- phylopars_fit$pars$sigma2  # older Rphylopars
  if (is.matrix(Sigma)) {
    as.numeric(diag(Sigma))
  } else {
    as.numeric(Sigma)
  }
}

# Build per-cell conditional-MVN prior for Phase 7 EM.
#
# Given K x K Σ and n x K posterior-mean liability matrix L_hat from the
# previous iteration, return (n x K) matrices mu_prior and sd_prior where
# for each cell (i, j):
#   mu[i, j] = Σ[j, -j] %*% solve(Σ[-j, -j]) %*% L_hat[i, -j observed]
#   sd[i, j] = sqrt(Σ[j, j] - Σ[j, -j] %*% solve(Σ[-j, -j]) %*% Σ[-j, j])
# Missing-pattern handling: if L_hat[i, -j] has NAs, restrict the condition
# to the observed subset with a smaller solve(); if ALL are NA, fall back
# to the unconditional prior (mu = 0, sd = sqrt(Σ[j, j])) for that cell.
#
# K is typically small (2–10) so per-cell matrix inversion with restricted
# subsets is cheap relative to the phylopars fit that produced Σ.
#' @keywords internal
#' @noRd
build_conditional_prior <- function(Sigma, L_hat, eps = 1e-8) {
  if (!is.matrix(Sigma) || nrow(Sigma) != ncol(Sigma)) {
    stop("build_conditional_prior: Sigma must be a square matrix.",
         call. = FALSE)
  }
  K <- ncol(Sigma)
  if (K != ncol(L_hat)) {
    stop("build_conditional_prior: ncol(L_hat) = ", ncol(L_hat),
         " != ncol(Sigma) = ", K, ".", call. = FALSE)
  }
  n <- nrow(L_hat)

  mu_mat <- matrix(0,  nrow = n, ncol = K)
  sd_mat <- matrix(1,  nrow = n, ncol = K)

  sd_plain <- sqrt(pmax(diag(Sigma), eps))  # fallback prior when all-NA

  for (j in seq_len(K)) {
    sigma_jj  <- Sigma[j, j]
    sigma_jnj <- Sigma[j, -j, drop = FALSE]         # 1 x (K-1)
    sigma_njn <- Sigma[-j, -j, drop = FALSE]        # (K-1) x (K-1)

    # All-observed cells share one precomputed inverse. Cells with NAs use
    # a per-cell inverse on the observed subset.
    inv_full  <- NULL

    for (i in seq_len(n)) {
      li <- L_hat[i, -j, drop = TRUE]
      na <- is.na(li)

      if (all(na)) {
        # No conditioning info — fall back to unconditional prior.
        mu_mat[i, j] <- 0
        sd_mat[i, j] <- sd_plain[j]
        next
      }

      if (!any(na)) {
        if (is.null(inv_full)) {
          inv_full <- solve(sigma_njn + diag(eps, K - 1))
        }
        mu_mat[i, j] <- as.numeric(sigma_jnj %*% inv_full %*% li)
        var_ij <- sigma_jj - as.numeric(sigma_jnj %*% inv_full %*% t(sigma_jnj))
        sd_mat[i, j] <- sqrt(max(var_ij, eps))
      } else {
        # Restricted conditional on observed subset of -j
        obs <- which(!na)
        s_jnj_obs <- sigma_jnj[, obs, drop = FALSE]
        s_njn_obs <- sigma_njn[obs, obs, drop = FALSE]
        inv_obs <- solve(s_njn_obs + diag(eps, length(obs)))
        mu_mat[i, j] <- as.numeric(s_jnj_obs %*% inv_obs %*% li[obs])
        var_ij <- sigma_jj -
          as.numeric(s_jnj_obs %*% inv_obs %*% t(s_jnj_obs))
        sd_mat[i, j] <- sqrt(max(var_ij, eps))
      }
    }
  }

  list(mu_prior = mu_mat, sd_prior = sd_mat)
}

# Relative Frobenius distance between two variance vectors (or matrices).
# Used as the Phase 6 EM early-stop criterion:
#   delta = ||v_new - v_old||_F / ||v_old||_F
# Returns +Inf if `v_old` has zero norm (first real iter).
#' @keywords internal
#' @noRd
rel_frobenius <- function(v_new, v_old) {
  den <- sqrt(sum(v_old^2))
  if (!is.finite(den) || den == 0) return(Inf)
  num <- sqrt(sum((v_new - v_old)^2))
  num / den
}

# Phase 6 EM wrapper around fit_joint_threshold_baseline().
# Iterates:
#   1. Run fit_joint_threshold_baseline() with sd_prior_vec = sqrt(diag(Sigma_prev)).
#   2. Extract new Sigma from the fit's phylopars_fit slot.
#   3. Check relative Frobenius convergence; break if delta < em_tol.
# Returns the final fit (shape identical to fit_joint_threshold_baseline's
# return) augmented with $em_state. $phylopars_fit and $fit_cols_idx are
# returned unchanged — Phase 6 doesn't hide them.
#' @keywords internal
#' @noRd
fit_joint_threshold_baseline_em <- function(data, tree, splits,
                                             graph = NULL,
                                             soft_aggregate = FALSE,
                                             em_iterations = 5L,
                                             em_tol = 1e-3) {
  stopifnot(joint_mvn_available(), em_iterations >= 1L)

  sd_prior_vec <- NULL  # iter 1 = plug-in N(0, 1)
  prev_vars    <- NULL
  delta        <- Inf
  base         <- NULL
  iter_run     <- 0L
  n_liab_cols  <- NULL   # set on iter 1 after first build

  for (iter in seq_len(em_iterations)) {
    new_base <- tryCatch(
      fit_joint_threshold_baseline(data, tree, splits = splits, graph = graph,
                                    soft_aggregate = soft_aggregate,
                                    sd_prior_vec = sd_prior_vec),
      error = function(e) NULL
    )

    if (is.null(new_base)) {
      if (is.null(base)) {
        stop("fit_joint_threshold_baseline_em: phylopars failed on iter 1; ",
             "no fallback available. Try em_iterations = 0L.", call. = FALSE)
      }
      warning("fit_joint_threshold_baseline_em: phylopars failed at iter ",
              iter, ". Returning previous iteration's baseline.",
              call. = FALSE)
      break
    }

    iter_run <- iter
    base     <- new_base
    if (is.null(n_liab_cols)) n_liab_cols <- length(new_base$liab_cols)

    if (is.null(new_base$phylopars_fit)) {
      # All cols < 2 obs; no fit possible. Degenerate case — bail.
      break
    }
    new_vars <- extract_liability_variances(new_base$phylopars_fit)

    if (!is.null(prev_vars)) {
      delta <- rel_frobenius(new_vars, prev_vars)
      if (delta < em_tol) break
    }

    prev_vars    <- new_vars
    # Align sd_prior_vec to the FULL liab_cols. Unfitted cols stay at 1
    # (plug-in fallback).
    sd_prior_vec <- rep(1, n_liab_cols)
    sd_prior_vec[new_base$fit_cols_idx] <- sqrt(pmax(new_vars, 1e-8))
  }

  base$em_state <- list(
    iterations_run = iter_run,
    converged      = is.finite(delta) && delta < em_tol,
    final_delta    = delta
  )
  base
}

