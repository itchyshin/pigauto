#' Assemble liability-scale input matrix for the threshold-model joint baseline
#'
#' For each trait type:
#'   - continuous/count/ordinal/proportion: pass through the z-scored latent value
#'     directly (it IS the liability on the preprocessed scale).
#'   - binary: apply truncated-Gaussian E-step with a vague N(0,1) prior so
#'     observed y=1 cells get a positive posterior mean and y=0 cells a
#'     negative one.
#'
#' Cells that are missing (NA in X_scaled) or held out by `splits` are set to
#' NA so Rphylopars marginalises over them.
#'
#' Categorical, multi_proportion, and zi_count columns are NOT included in
#' Phase 3 — they are handled by the existing label-propagation / BM paths
#' and skipped here.
#'
#' Single-observation mode only. Multi-obs datasets should take the
#' per-column baseline path; the dispatcher in `fit_baseline()` already
#' enforces this.
#'
#' @param data pigauto_data from preprocess_traits.
#' @param splits output of make_missing_splits() or NULL.
#' @return list(X_liab, liab_cols, liab_types) where X_liab is
#'   n_species x length(liab_cols), liab_cols is the integer index into the
#'   full p_latent space, and liab_types is a character vector aligned with
#'   liab_cols.
#' @keywords internal
#' @noRd
build_liability_matrix <- function(data, splits = NULL) {
  stopifnot(!isTRUE(data$multi_obs))
  trait_map <- data$trait_map
  X         <- data$X_scaled
  n         <- nrow(X)

  # Apply split masking up front (same leakage protection as Phase 2)
  if (!is.null(splits)) {
    X[splits$val_idx]  <- NA
    X[splits$test_idx] <- NA
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
    }
    # categorical / multi_proportion / zi_count: skipped in Phase 3
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
      # Observed 0/1 -> truncated-Gaussian posterior mean with N(0,1) prior.
      # Route through the estep_liability() dispatcher so the contract stays
      # consistent across trait types (CLAUDE.md mandates dispatcher use).
      for (i in seq_len(n)) {
        v <- src_col[i]
        if (is.na(v)) next
        post <- estep_liability(tm_j, observed = v, mu_prior = 0, sd_prior = 1)
        X_liab[i, j] <- post$mean
      }
    } else {
      # Continuous-family: pass through
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
fit_joint_threshold_baseline <- function(data, tree, splits, graph = NULL) {
  stopifnot(joint_mvn_available())

  built <- build_liability_matrix(data, splits = splits)
  X_liab     <- built$X_liab
  liab_cols  <- built$liab_cols
  liab_types <- built$liab_types

  spp <- if (!is.null(data$species_names)) data$species_names else rownames(data$X_scaled)

  df_in          <- as.data.frame(X_liab)
  df_in$species  <- spp
  df_in          <- df_in[, c("species", colnames(X_liab)), drop = FALSE]

  fit <- Rphylopars::phylopars(
    trait_data = df_in,
    tree       = tree,
    model      = "BM"
  )

  tip_rows <- match(spp, rownames(fit$anc_recon))
  mu_liab  <- fit$anc_recon[tip_rows, , drop = FALSE]
  se_liab  <- sqrt(fit$anc_var[tip_rows, , drop = FALSE])
  rownames(mu_liab) <- spp
  rownames(se_liab) <- spp
  colnames(mu_liab) <- colnames(X_liab)
  colnames(se_liab) <- colnames(X_liab)

  list(mu_liab    = mu_liab,
       se_liab    = se_liab,
       liab_cols  = liab_cols,
       liab_types = liab_types)
}
