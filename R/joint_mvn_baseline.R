#' Is the Rphylopars delegate available?
#'
#' Used by `fit_baseline()` to decide whether to dispatch to the
#' joint multivariate-BM baseline (Level C Phase 2). When FALSE,
#' `fit_baseline()` falls back to the per-column BM path.
#'
#' @keywords internal
#' @noRd
joint_mvn_available <- function() {
  requireNamespace("Rphylopars", quietly = TRUE)
}

#' Joint multivariate BM baseline for continuous-family latent columns
#'
#' Delegates to `Rphylopars::phylopars()`. Returns per-cell posterior mean
#' and SE in the same shape as the per-column path in [fit_baseline()].
#' Non-BM columns (binary, categorical) are untouched; they stay at zero
#' in the returned matrices (caller handles them separately).
#'
#' @param data pigauto_data (from `preprocess_traits`).
#' @param tree phylo.
#' @param splits output of `make_missing_splits()`; val/test cells are
#'   masked to NA before the joint fit (no leakage).
#' @param graph output of `build_phylo_graph()` (unused here but kept for
#'   interface parity with `fit_baseline()`).
#' @return list(mu, se), each `n_species x p_latent`.
#' @keywords internal
#' @noRd
fit_joint_mvn_baseline <- function(data, tree, splits, graph = NULL) {
  stopifnot(joint_mvn_available())

  trait_map <- data$trait_map
  X         <- data$X_scaled
  n         <- nrow(X)
  p         <- ncol(X)

  # Species names: use species_names if present (single-obs path), else rownames
  spp <- if (!is.null(data$species_names)) data$species_names else rownames(X)

  # Identify BM-eligible columns — same rule as fit_baseline.R lines 107-128
  bm_cols <- integer(0)
  for (tm in trait_map) {
    if (tm$type %in% c("continuous", "count", "ordinal", "proportion",
                       "multi_proportion")) {
      bm_cols <- c(bm_cols, tm$latent_cols)
    } else if (tm$type == "zi_count") {
      # Magnitude column only (col 2); gate column (col 1) is handled by LP
      bm_cols <- c(bm_cols, tm$latent_cols[2])
    }
  }

  # Subset to BM-eligible columns
  X_bm <- X[, bm_cols, drop = FALSE]

  # Mask val/test cells before fitting (no leakage) — same linear-index approach
  # as fit_baseline.R lines 93-95, but restricted to the BM-column subset.
  if (!is.null(splits)) {
    hold_idx <- c(splits$val_idx, splits$test_idx)
    # Decompose linear indices from full p-column space to BM-column space
    row_i <- ((hold_idx - 1L) %% n) + 1L
    col_j <- ((hold_idx - 1L) %/% n) + 1L
    bm_match <- match(col_j, bm_cols)          # which held-out cols are BM cols
    keep <- !is.na(bm_match)
    if (any(keep)) {
      local_col <- bm_match[keep]
      local_row <- row_i[keep]
      for (k in seq_along(local_row)) {
        X_bm[local_row[k], local_col[k]] <- NA_real_
      }
    }
  }

  # Build trait_data for Rphylopars: species column first, then trait columns
  df_in          <- as.data.frame(X_bm)
  df_in$species  <- spp
  # Rphylopars requires species as first column
  df_in          <- df_in[, c("species", colnames(X_bm)), drop = FALSE]

  fit <- Rphylopars::phylopars(
    trait_data = df_in,
    tree       = tree,
    model      = "BM"
  )

  # anc_recon and anc_var are (n_tips + n_internal) x q matrices.
  # Tip rows come first and their rownames match tip labels.
  # We extract the tip block by matching species names.
  tip_rows <- match(spp, rownames(fit$anc_recon))
  mu_bm    <- fit$anc_recon[tip_rows, , drop = FALSE]
  se_bm    <- sqrt(fit$anc_var[tip_rows, , drop = FALSE])

  # Assemble full p_latent output matrices (non-BM cols stay at 0)
  mu <- matrix(0, n, p, dimnames = list(spp, colnames(X)))
  se <- matrix(0, n, p, dimnames = list(spp, colnames(X)))
  mu[, bm_cols] <- mu_bm
  se[, bm_cols] <- se_bm

  list(mu = mu, se = se)
}
