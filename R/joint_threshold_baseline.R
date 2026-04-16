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
  }

  list(mu_liab    = mu_liab,
       se_liab    = se_liab,
       liab_cols  = liab_cols,
       liab_types = liab_types)
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

