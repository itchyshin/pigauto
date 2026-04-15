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
#' @param data pigauto_data from preprocess_traits.
#' @param splits output of make_missing_splits() or NULL.
#' @return list(X_liab, liab_cols, liab_types) where X_liab is
#'   n_species x length(liab_cols), liab_cols is the integer index into the
#'   full p_latent space, and liab_types is a character vector aligned with
#'   liab_cols.
#' @keywords internal
#' @noRd
build_liability_matrix <- function(data, splits = NULL) {
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
  for (tm in trait_map) {
    if (tm$type %in% c("continuous", "count", "ordinal", "proportion")) {
      liab_cols  <- c(liab_cols, tm$latent_cols)
      liab_types <- c(liab_types, tm$type)
    } else if (tm$type == "binary") {
      liab_cols  <- c(liab_cols, tm$latent_cols)
      liab_types <- c(liab_types, "binary")
    }
    # categorical / multi_proportion / zi_count: skipped in Phase 3
  }

  X_liab <- matrix(NA_real_, nrow = n, ncol = length(liab_cols))
  colnames(X_liab) <- colnames(X)[liab_cols]
  rownames(X_liab) <- rownames(X)

  for (j in seq_along(liab_cols)) {
    col_idx <- liab_cols[j]
    tp      <- liab_types[j]
    src_col <- X[, col_idx]

    if (tp == "binary") {
      # Observed 0/1 -> truncated-Gaussian posterior mean with N(0,1) prior
      for (i in seq_len(n)) {
        v <- src_col[i]
        if (is.na(v)) next
        post <- estep_liability_binary(y = as.integer(v),
                                        mu_prior = 0, sd_prior = 1)
        X_liab[i, j] <- post$mean
      }
    } else {
      # Continuous-family: pass through
      X_liab[, j] <- src_col
    }
  }

  list(X_liab = X_liab, liab_cols = liab_cols, liab_types = liab_types)
}
