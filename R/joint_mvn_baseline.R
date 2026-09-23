#' Is the joint MVN baseline available?
#'
#' Historically gated on `requireNamespace("Rphylopars")`. As of the
#' in-house Sigma solver (R/joint_mvn_solver.R) the joint MVN baseline
#' has no external dependency and is always available. Kept as a
#' function (rather than a constant) so call sites in `fit_baseline()`
#' don't need to change.
#'
#' @keywords internal
#' @noRd
joint_mvn_available <- function() {
  TRUE
}

#' Joint multivariate BM baseline for continuous-family latent columns
#'
#' Fits a matrix-normal BM model jointly across BM-eligible latent
#' columns via the internal `fit_mvn_bm_inhouse()` solver
#' (R/joint_mvn_solver.R). Returns
#' per-cell posterior mean and SE in the same shape as the per-column
#' path in [fit_baseline()]. Non-BM columns (binary, categorical) are
#' untouched; they stay at zero in the returned matrices (caller
#' handles them separately).
#'
#' @param data pigauto_data (from `preprocess_traits`).
#' @param tree phylo.
#' @param splits output of `make_missing_splits()`; val/test cells are
#'   masked to NA before the joint fit (no leakage).
#' @param graph output of `build_phylo_graph()` (unused here but kept for
#'   interface parity with `fit_baseline()`).
#' @param joint_solver character, `"inhouse"` (default) or `"rphylopars"`.
#'   See `fit_joint_solver()` in R/joint_mvn_solver.R.
#' @param joint_refine_iter integer, default `0L`. See `fit_joint_solver()`
#'   in R/joint_mvn_solver.R.
#' @param lambda_mode character, `"fixed_1"` (default) or `"estimate"`.
#'   S4 (dispatcher, feat/joint-lambda-default): every column this function
#'   fits IS already continuous-family (this path only ever runs when no
#'   binary/ordinal columns are present -- see `fit_baseline.R`'s
#'   `use_continuous_joint` guard), so under `"estimate"` every column is
#'   eligible for its own Pagel's lambda (`lambda_cols = NULL` downstream).
#'   `"cv"` / `"bayes"` have no joint analogue; callers must translate those
#'   to `"fixed_1"` before calling this function (see `fit_baseline.R`).
#' @param lambda_fixed optional named numeric vector (names = latent column
#'   names of `data$X_scaled`) giving a fixed lambda per column, overriding
#'   `lambda_mode` entirely (spec 4.5 predict-time rebuild: typically a
#'   previous fit's own `$lambda_per_trait`, replayed rather than
#'   re-estimated).
#' @return list(mu, se, lambda_per_trait, lambda_block), mu/se each
#'   `n_species x p_latent`; `lambda_per_trait` is named by the BM-eligible
#'   columns actually fit (a subset of `colnames(data$X_scaled)`).
#' @keywords internal
#' @noRd
fit_joint_mvn_baseline <- function(data, tree, splits, graph = NULL,
                                   soft_aggregate = FALSE,
                                   joint_solver = "inhouse",
                       predict_method = "per_column",
                                   joint_refine_iter = 0L,
                                   lambda_mode = "fixed_1",
                                   lambda_fixed = NULL) {
  stopifnot(joint_mvn_available())

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

  # Dispatch to the in-house solver (default) or Rphylopars (joint_solver
  # = "rphylopars", with fallback to in-house on failure). Returns a list
  # with $anc_recon and $anc_var on n_tips x q (no internal nodes).
  L_in <- X_bm
  rownames(L_in) <- spp

  # S4: every column of L_in is already continuous-family (this function
  # only ever runs when fit_baseline()'s use_continuous_joint fires, which
  # requires zero binary/ordinal columns), so there is no liability-column
  # exclusion to compute here -- lambda_cols = NULL (all columns eligible)
  # covers it. lambda_fixed, when supplied, takes priority over lambda_mode
  # entirely (predict-time rebuild).
  lambda_arg <- if (!is.null(lambda_fixed)) {
    unname(lambda_fixed[colnames(L_in)])
  } else if (identical(lambda_mode, "estimate")) {
    "estimate"
  } else {
    "fixed_1"
  }
  fit <- fit_joint_solver(L = L_in, tree = tree, joint_solver = joint_solver,
                          predict_method = predict_method,
                          joint_refine_iter = joint_refine_iter,
                          lambda = lambda_arg)

  tip_rows <- match(spp, rownames(fit$anc_recon))
  mu_bm    <- fit$anc_recon[tip_rows, , drop = FALSE]
  se_bm    <- sqrt(fit$anc_var[tip_rows, , drop = FALSE])

  # Assemble full p_latent output matrices (non-BM cols stay at 0)
  mu <- matrix(0, n, p, dimnames = list(spp, colnames(X)))
  se <- matrix(0, n, p, dimnames = list(spp, colnames(X)))
  mu[, bm_cols] <- mu_bm
  se[, bm_cols] <- se_bm

  # fit$lambda_per_trait is a scalar NA_real_ on the rphylopars path
  # (fit_joint_solver() does not read phylopars' own lambda back in this
  # slice) -- broadcast to a named length-K vector so the caller always
  # gets one entry per column, honestly reporting "not tracked" as NA
  # rather than defaulting those columns to a possibly-wrong 1.
  lambda_per_trait <- fit$lambda_per_trait
  if (!is.null(lambda_per_trait)) {
    if (length(lambda_per_trait) != ncol(L_in)) {
      lambda_per_trait <- rep(lambda_per_trait[1], ncol(L_in))
    }
    names(lambda_per_trait) <- colnames(L_in)
  }

  list(mu = mu, se = se, lambda_per_trait = lambda_per_trait,
       lambda_block = fit$lambda_block)
}
