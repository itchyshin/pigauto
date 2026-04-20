#' Fit K independent threshold-joint baselines for a categorical trait (OVR)
#'
#' For each class k in 1..K, constructs a synthetic binary trait
#' "is_class_k vs rest" and runs `fit_joint_threshold_baseline()` on a
#' modified pigauto_data where the categorical is replaced by this
#' synthetic binary. Returns the K per-species P(class_k) probabilities.
#'
#' This is the OVR (one-vs-rest) strategy BACE uses. Each individual fit
#' has only 1 categorical-related column so Rphylopars stays
#' well-conditioned.
#'
#' @param data pigauto_data from preprocess_traits.
#' @param tree phylo.
#' @param trait_name character, name of the categorical trait.
#' @param splits output of make_missing_splits() or NULL.
#' @param graph optional graph (unused).
#' @return numeric matrix (n_species x K) of P(class_k) per species; NA
#'   for classes whose fit failed.
#' @keywords internal
#' @noRd
fit_ovr_categorical_fits <- function(data, tree, trait_name,
                                      splits = NULL, graph = NULL,
                                      soft_aggregate = FALSE,
                                      sd_prior_k = NULL) {
  stopifnot(joint_mvn_available())
  # sd_prior_k: optional length-K vector of per-class prior SDs for Phase 6
  # EM. NULL default (plug-in 1) is the v0.9.1 byte-identical path.

  # If multi-obs: aggregate X to species level + build a single-obs pigauto_data
  # shell. The K-loop's reduced pd_k objects are then naturally single-obs.
  if (isTRUE(data$multi_obs)) {
    agg <- aggregate_to_species(data, splits = splits,
                                soft_aggregate = soft_aggregate)
    data_single <- data
    data_single$X_scaled       <- agg$X_species
    data_single$multi_obs      <- FALSE
    data_single$n_obs          <- data$n_species
    data_single$obs_to_species <- NULL
    data_single$obs_species    <- NULL
    data   <- data_single
    splits <- agg$splits_species
  }

  # Locate the categorical trait_map entry by matching column name prefix
  tm_cat <- NULL
  for (tm in data$trait_map) {
    if (!identical(tm$type, "categorical")) next
    col_names_k <- colnames(data$X_scaled)[tm$latent_cols]
    if (any(startsWith(col_names_k, paste0(trait_name, "=")))) {
      tm_cat <- tm
      break
    }
  }
  if (is.null(tm_cat)) {
    stop("fit_ovr_categorical_fits: no categorical trait '",
         trait_name, "' found in trait_map.", call. = FALSE)
  }

  K      <- tm_cat$n_latent
  k_cols <- tm_cat$latent_cols
  spp    <- if (!is.null(data$species_names)) data$species_names else rownames(data$X_scaled)
  n      <- length(spp)

  probs <- matrix(NA_real_, nrow = n, ncol = K,
                   dimnames = list(spp, colnames(data$X_scaled)[k_cols]))
  # per_class_fits holds each OVR k's phylopars_fit for Phase 6 EM variance
  # extraction. Attached as an attribute on `probs` so the function's main
  # return shape (n x K matrix) is unchanged for existing callers.
  per_class_fits <- vector("list", K)

  for (k in seq_len(K)) {
    pd_k <- build_ovr_pd(data, tm_cat, k, soft_proportion = soft_aggregate)

    # Reindex splits for the reduced X_scaled
    splits_k <- if (is.null(splits)) NULL else reindex_splits(
      splits, n_rows = n,
      p_old = ncol(data$X_scaled),
      kept_cols = attr(pd_k, "kept_cols")
    )

    # Build per-k sd_prior_vec: only the synthetic binary column uses
    # sd_prior_k[k]; all other liability columns in pd_k (continuous, etc.)
    # stay at 1. Sized on the fly once we see the liab_cols count on iter 1.
    # Simplest: pass sd_prior_k[k] scalar; fit_joint_threshold_baseline will
    # broadcast via a K-element vector below.
    sd_k <- if (is.null(sd_prior_k)) 1 else sd_prior_k[k]

    # To keep the prior targeted at the synthetic binary only, we need to
    # know the position of that column in pd_k's liab_cols. Peek via a dry
    # build. Cheap: build_liability_matrix doesn't call phylopars.
    built_k <- build_liability_matrix(pd_k, splits = splits_k,
                                       soft_aggregate = soft_aggregate)
    sdv <- rep(1, length(built_k$liab_cols))
    synth_col  <- attr(pd_k, "synthetic_bin_col")
    which_synth <- which(built_k$liab_cols == synth_col)
    if (length(which_synth) == 1L) sdv[which_synth] <- sd_k

    jt <- tryCatch(
      fit_joint_threshold_baseline(pd_k, tree, splits = splits_k,
                                    graph = graph,
                                    soft_aggregate = soft_aggregate,
                                    sd_prior_vec = sdv),
      error = function(e) NULL
    )
    if (is.null(jt)) next

    # Find the synthetic binary's col index in the reduced fit output
    our_pos <- attr(pd_k, "synthetic_bin_col")
    bin_idx <- which(jt$liab_types == "binary" &
                       jt$liab_cols == our_pos)
    if (length(bin_idx) == 0L) next

    dec <- decode_binary_liability(mu_liab = jt$mu_liab[, bin_idx],
                                    se_liab = jt$se_liab[, bin_idx])
    probs[, k] <- dec$p
    per_class_fits[[k]] <- list(
      jt           = jt,
      liab_col_idx = bin_idx,             # position in jt's liab_cols
      synth_col    = synth_col             # original pd_k col index
    )
  }

  attr(probs, "per_class_fits") <- per_class_fits
  probs
}

# Helper: construct a modified pigauto_data for class k's OVR fit
# @keywords internal
# @noRd
build_ovr_pd <- function(data, tm_cat, k, soft_proportion = FALSE) {
  K      <- tm_cat$n_latent
  k_cols <- tm_cat$latent_cols
  drop_cols <- setdiff(k_cols, k_cols[k])
  kept_cols <- setdiff(seq_len(ncol(data$X_scaled)), drop_cols)
  col_map   <- setNames(seq_along(kept_cols), as.character(kept_cols))

  pd_k <- data
  pd_k$X_scaled <- data$X_scaled[, kept_cols, drop = FALSE]

  # Rebuild trait_map: replace categorical with synthetic binary, reindex
  # other traits' latent_cols.
  new_tm <- list()
  for (tm in data$trait_map) {
    if (identical(tm$type, "categorical") &&
        identical(tm$latent_cols, k_cols)) {
      new_tm[[length(new_tm) + 1L]] <- list(
        name        = paste0(tm_cat$name %||% "cat", "_is_",
                              tm_cat$levels[k] %||% as.character(k)),
        type        = "binary",
        n_latent    = 1L,
        latent_cols = unname(col_map[as.character(k_cols[k])]),
        levels      = c("no", "yes"),
        mean        = 0,
        sd          = 1
      )
    } else {
      tm2 <- tm
      tm2$latent_cols <- unname(col_map[as.character(tm$latent_cols)])
      new_tm[[length(new_tm) + 1L]] <- tm2
    }
  }
  pd_k$trait_map <- new_tm

  attr(pd_k, "kept_cols") <- kept_cols
  synth_col <- unname(col_map[as.character(k_cols[k])])
  attr(pd_k, "synthetic_bin_col") <- synth_col
  if (soft_proportion) {
    pd_k$proportion_cols <- synth_col
  }
  pd_k
}

# Helper: reindex splits linear indices into a reduced column set
# @keywords internal
# @noRd
reindex_splits <- function(splits, n_rows, p_old, kept_cols) {
  col_map <- setNames(seq_along(kept_cols), as.character(kept_cols))
  reindex <- function(lin_idx) {
    if (length(lin_idx) == 0L) return(integer(0))
    row_i <- ((lin_idx - 1L) %% n_rows) + 1L
    col_j <- ((lin_idx - 1L) %/% n_rows) + 1L
    keep  <- col_j %in% kept_cols
    new_col <- col_map[as.character(col_j[keep])]
    (unname(new_col) - 1L) * n_rows + row_i[keep]
  }
  sp_k <- splits
  sp_k$val_idx  <- reindex(splits$val_idx)
  sp_k$test_idx <- reindex(splits$test_idx)
  if (!is.null(splits$mask)) {
    sp_k$mask <- splits$mask[, kept_cols, drop = FALSE]
  }
  sp_k
}

#' Phase 6 EM wrapper around fit_ovr_categorical_fits().
#'
#' Iterates a length-K per-class variance vector `sigma_cat`. Each iter runs
#' K OVR fits with `sd_prior_k = sqrt(sigma_cat)`; after all K complete,
#' extract the synthetic-binary liability variance from each fit's
#' phylopars_fit, update `sigma_cat`, and check relative Frobenius
#' convergence.
#'
#' Returns the same shape as fit_ovr_categorical_fits() (n x K probs matrix)
#' with `em_state` attached as an attribute. Phase 6 diagonal-only: the
#' off-diagonal Σ_cat[j, k] is not used, consistent with the design spec.
#'
#' @inheritParams fit_ovr_categorical_fits
#' @param em_iterations integer, max number of EM iters (>= 1).
#' @param em_tol numeric, relative Frobenius early-stop threshold.
#' @keywords internal
#' @noRd
fit_ovr_categorical_fits_em <- function(data, tree, trait_name,
                                         splits = NULL, graph = NULL,
                                         soft_aggregate = FALSE,
                                         em_iterations = 5L,
                                         em_tol = 1e-3) {
  stopifnot(em_iterations >= 1L)

  # Discover K by running a single plug-in iter first.
  iter1 <- fit_ovr_categorical_fits(data, tree, trait_name = trait_name,
                                     splits = splits, graph = graph,
                                     soft_aggregate = soft_aggregate,
                                     sd_prior_k = NULL)
  K <- ncol(iter1)
  sigma_cat <- rep(1, K)
  probs     <- iter1
  iter_run  <- 1L
  delta     <- Inf

  prev_sigma <- NULL
  # Extract iter-1 variances from the per_class_fits attribute.
  pcf <- attr(iter1, "per_class_fits")
  new_sigma <- vapply(seq_len(K), function(k) {
    cf <- pcf[[k]]
    if (is.null(cf) || is.null(cf$jt$phylopars_fit)) return(NA_real_)
    v  <- extract_liability_variances(cf$jt$phylopars_fit)
    if (length(v) == 0L) return(NA_real_)
    # Which fit col is the synthetic binary? Match synth_col against
    # fit_cols_idx of the jt (positions in the full liab_cols vector).
    # cf$liab_col_idx is the bin_idx in jt$liab_cols; we need its position
    # in jt$fit_cols_idx.
    fit_idx_pos <- match(cf$liab_col_idx, cf$jt$fit_cols_idx)
    if (is.na(fit_idx_pos)) return(NA_real_)
    v[fit_idx_pos]
  }, numeric(1))
  new_sigma[is.na(new_sigma)] <- sigma_cat[is.na(new_sigma)]
  prev_sigma <- new_sigma
  sigma_cat  <- new_sigma

  if (em_iterations >= 2L) {
    for (iter in 2:em_iterations) {
      sd_prior_k <- sqrt(pmax(sigma_cat, 1e-8))
      iter_probs <- tryCatch(
        fit_ovr_categorical_fits(data, tree, trait_name = trait_name,
                                  splits = splits, graph = graph,
                                  soft_aggregate = soft_aggregate,
                                  sd_prior_k = sd_prior_k),
        error = function(e) NULL
      )
      if (is.null(iter_probs)) {
        warning("fit_ovr_categorical_fits_em: OVR K-loop failed at iter ",
                iter, "; returning previous iter's fits.", call. = FALSE)
        break
      }
      iter_run <- iter
      probs    <- iter_probs

      pcf <- attr(iter_probs, "per_class_fits")
      new_sigma <- vapply(seq_len(K), function(k) {
        cf <- pcf[[k]]
        if (is.null(cf) || is.null(cf$jt$phylopars_fit)) return(NA_real_)
        v  <- extract_liability_variances(cf$jt$phylopars_fit)
        fit_idx_pos <- match(cf$liab_col_idx, cf$jt$fit_cols_idx)
        if (is.na(fit_idx_pos) || length(v) == 0L) return(NA_real_)
        v[fit_idx_pos]
      }, numeric(1))
      new_sigma[is.na(new_sigma)] <- sigma_cat[is.na(new_sigma)]

      delta <- rel_frobenius(new_sigma, prev_sigma)
      if (delta < em_tol) {
        sigma_cat <- new_sigma
        break
      }
      prev_sigma <- new_sigma
      sigma_cat  <- new_sigma
    }
  }

  attr(probs, "em_state") <- list(
    iterations_run = iter_run,
    converged      = is.finite(delta) && delta < em_tol,
    final_delta    = delta,
    sigma_cat      = sigma_cat
  )
  probs
}

#' Normalise K independent OVR probabilities into log-probabilities
#'
#' @param probs numeric matrix (n x K), unnormalised P(class_k).
#' @return numeric matrix (n x K), log-probabilities summing to 1 per row.
#' @keywords internal
#' @noRd
decode_ovr_categorical <- function(probs) {
  K <- ncol(probs)
  out <- probs
  for (i in seq_len(nrow(probs))) {
    row <- probs[i, ]
    na_mask <- is.na(row)
    if (any(!na_mask)) row[!na_mask] <- pmax(row[!na_mask], 0.01)
    if (all(na_mask)) {
      # All NA: uniform
      row <- rep(1/K, K)
    } else if (any(na_mask)) {
      s <- sum(row[!na_mask])
      cap <- 1 - 0.01 * sum(na_mask)
      if (s > cap) row[!na_mask] <- row[!na_mask] * cap / s
      row[na_mask] <- (1 - sum(row[!na_mask])) / sum(na_mask)
    } else {
      row <- row / sum(row)
    }
    row <- pmax(row, 0.01)
    row <- row / sum(row)
    out[i, ] <- log(row)
  }
  out
}
