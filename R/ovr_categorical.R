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
                                      soft_aggregate = FALSE) {
  stopifnot(joint_mvn_available())

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

  for (k in seq_len(K)) {
    pd_k <- build_ovr_pd(data, tm_cat, k, soft_proportion = soft_aggregate)

    # Reindex splits for the reduced X_scaled
    splits_k <- if (is.null(splits)) NULL else reindex_splits(
      splits, n_rows = n,
      p_old = ncol(data$X_scaled),
      kept_cols = attr(pd_k, "kept_cols")
    )

    jt <- tryCatch(
      fit_joint_threshold_baseline(pd_k, tree, splits = splits_k,
                                    graph = graph,
                                    soft_aggregate = soft_aggregate),
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
  }

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
