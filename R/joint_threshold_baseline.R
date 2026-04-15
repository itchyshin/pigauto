#' Assemble liability-scale input matrix for the threshold-model joint baseline
#'
#' For each trait type:
#'   - continuous/count/ordinal/proportion: pass through the z-scored latent value
#'     directly (it IS the liability on the preprocessed scale).
#'   - binary: apply truncated-Gaussian E-step with a vague N(0,1) prior so
#'     observed y=1 cells get a positive posterior mean and y=0 cells a
#'     negative one.
#'   - categorical: K latent cols per trait; observed rows go through the
#'     K-dim plug-in E-step (`estep_liability_categorical`) producing a
#'     sum-zero K-vector posterior mean. Missing rows stay NA in all K
#'     cols (row-level missingness: the whole one-hot is observed or
#'     none of it is).
#'
#' Cells that are missing (NA in X_scaled) or held out by `splits` are set to
#' NA so Rphylopars marginalises over them.
#'
#' multi_proportion and zi_count columns are NOT included in Phase 4 — they
#' are handled by the existing label-propagation / BM paths and skipped here.
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
build_liability_matrix <- function(data, splits = NULL,
                                    cat_encoding = c("joint_K", "ovr")) {
  stopifnot(!isTRUE(data$multi_obs))
  cat_encoding <- match.arg(cat_encoding)
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
    } else if (tm$type == "categorical") {
      # K latent cols per categorical trait; ALL K cols join the joint
      # liability. Use the same `tm` entry for every one so downstream
      # decoding knows which cols group together.
      K <- tm$n_latent
      for (kk in seq_len(K)) {
        liab_cols  <- c(liab_cols, tm$latent_cols[kk])
        liab_types <- c(liab_types, "categorical")
        liab_tms   <- c(liab_tms, list(tm))
      }
    }
    # multi_proportion / zi_count: skipped in Phase 4
  }

  X_liab <- matrix(NA_real_, nrow = n, ncol = length(liab_cols))
  colnames(X_liab) <- colnames(X)[liab_cols]
  rownames(X_liab) <- rownames(X)

  # For joint_K categorical: populate K cols at once per trait — not per
  # column. Track which categorical trait cols we've already handled so we
  # skip them on subsequent iterations of the per-column loop.
  cat_cols_done <- integer(0)

  for (j in seq_along(liab_cols)) {
    col_idx <- liab_cols[j]
    tp      <- liab_types[j]
    tm_j    <- liab_tms[[j]]
    src_col <- X[, col_idx]

    if (tp == "categorical" && cat_encoding == "joint_K") {
      if (col_idx %in% cat_cols_done) next
      K       <- tm_j$n_latent
      k_cols  <- tm_j$latent_cols
      k_local <- match(k_cols, liab_cols)
      for (i in seq_len(n)) {
        obs_K <- X[i, k_cols]
        if (any(is.na(obs_K))) next
        post <- estep_liability(tm_j,
                                 observed = obs_K,
                                 mu_prior = rep(0, K),
                                 sd_prior = rep(1, K))
        X_liab[i, k_local] <- post$mean
      }
      cat_cols_done <- c(cat_cols_done, k_cols)

    } else if (tp == "categorical" && cat_encoding == "ovr") {
      for (i in seq_len(n)) {
        v <- src_col[i]
        if (is.na(v)) next
        post <- estep_liability_binary(y = as.integer(round(v)),
                                        mu_prior = 0, sd_prior = 1)
        X_liab[i, j] <- post$mean
      }

    } else if (tp == "binary") {
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
fit_joint_threshold_baseline <- function(data, tree, splits, graph = NULL,
                                          cat_encoding = c("joint_K", "ovr")) {
  stopifnot(joint_mvn_available())
  cat_encoding <- match.arg(cat_encoding)

  built <- build_liability_matrix(data, splits = splits, cat_encoding = cat_encoding)
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

    df_in          <- as.data.frame(X_fit)
    df_in$species  <- spp
    df_in          <- df_in[, c("species", colnames(X_fit)), drop = FALSE]

    fit <- Rphylopars::phylopars(
      trait_data = df_in,
      tree       = tree,
      model      = "BM"
    )

    tip_rows <- match(spp, rownames(fit$anc_recon))
    mu_fit   <- fit$anc_recon[tip_rows, , drop = FALSE]
    se_fit   <- sqrt(fit$anc_var[tip_rows, , drop = FALSE])

    # Validate shape: if Rphylopars still drops or duplicates columns, fail
    # loudly rather than mis-align. has_obs already filters columns with < 2
    # observations (phylopars errors on those internally).
    if (ncol(mu_fit) != length(fit_cols)) {
      stop("fit_joint_threshold_baseline: Rphylopars returned ",
           ncol(mu_fit), " columns but ", length(fit_cols),
           " were passed. Column alignment would be ambiguous.",
           call. = FALSE)
    }

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

#' Decode K-dim liability posterior to log-probabilities for categorical traits
#'
#' Two encodings:
#'   joint_K: K-dim sum-zero liability -> log-softmax over posterior means.
#'   ovr:     K independent probit probabilities -> renormalised across K.
#'
#' Both paths clip each class probability at 0.01 and renormalise so probs
#' sum to 1. Output matches the LP baseline's log-probability convention.
#'
#' @param mu_K numeric length K posterior mean.
#' @param se_K numeric length K posterior SD (optional).
#' @param cat_encoding "joint_K" or "ovr".
#' @return list(log_probs = K-length numeric vector).
#' @keywords internal
#' @noRd
decode_categorical_liability <- function(mu_K, se_K = NULL,
                                          cat_encoding = c("joint_K", "ovr")) {
  cat_encoding <- match.arg(cat_encoding)
  if (cat_encoding == "joint_K") {
    mu_max    <- max(mu_K)
    log_denom <- mu_max + log(sum(exp(mu_K - mu_max)))
    probs     <- exp(mu_K - log_denom)
  } else {
    if (is.null(se_K)) se_K <- rep(0, length(mu_K))
    probs <- stats::pnorm(mu_K / sqrt(1 + se_K^2))
    probs <- probs / sum(probs)
  }
  probs <- pmax(probs, 0.01)
  probs <- probs / sum(probs)
  list(log_probs = log(probs))
}
