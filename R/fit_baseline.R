#' Fit a phylogenetic BM baseline
#'
#' Fits an internal univariate Brownian Motion baseline using the
#' phylogenetic correlation matrix \eqn{R = \mathrm{cov2cor}(\mathrm{vcv}(\mathrm{tree}))}
#' and returns imputed means and standard errors for every species.
#'
#' @details
#' When \code{splits} is supplied the val and test cells are masked to
#' \code{NA} before fitting, so the baseline is evaluated under the same
#' conditions as \code{\link{fit_pigauto}}.
#'
#' For continuous, count, ordinal, and proportion traits (which are
#' continuous in latent space), each column is imputed independently
#' via conditional multivariate normal on the phylogenetic correlation
#' matrix: GLS phylogenetic mean, REML variance, and conditional
#' \eqn{E[y_m | y_o]}.  Binary and categorical traits use phylogenetic
#' label propagation: each species receives a personalised baseline
#' computed as the phylo-similarity-weighted average of observed values,
#' using a Gaussian kernel on cophenetic distances.
#'
#' @param data object of class \code{"pigauto_data"}.
#' @param tree object of class \code{"phylo"}.
#' @param splits list (output of \code{\link{make_missing_splits}}) or
#'   \code{NULL}.
#' @param model character. Evolutionary model: \code{"BM"} (default) or
#'   \code{"OU"}.
#' @param graph optional list returned by \code{\link{build_phylo_graph}}.
#'   When supplied, \code{graph$D} (cophenetic distances) is reused for
#'   label propagation and \code{graph$R_phy} (phylogenetic correlation
#'   matrix) is reused for BM imputation, avoiding duplicate \eqn{O(n^2)}
#'   allocations. When \code{NULL} (default), both matrices are computed
#'   here.
#' @return A list with:
#'   \describe{
#'     \item{mu}{Numeric matrix (n_species x p_latent), baseline means in
#'       latent scale.}
#'     \item{se}{Numeric matrix (n_species x p_latent), standard errors.}
#'   }
#' @examples
#' \dontrun{
#' data(avonet300, tree300, package = "pigauto")
#' traits <- avonet300; rownames(traits) <- traits$Species_Key
#' traits$Species_Key <- NULL
#' pd     <- preprocess_traits(traits, tree300)
#' splits <- make_missing_splits(pd$X_scaled, trait_map = pd$trait_map)
#' bl     <- fit_baseline(pd, tree300, splits)
#' }
#' @importFrom stats complete.cases rnorm rbinom
#' @export
fit_baseline <- function(data, tree, splits = NULL, model = "BM",
                         graph = NULL,
                         cat_encoding = c("joint_K", "ovr")) {
  if (!inherits(data, "pigauto_data")) {
    stop("'data' must be a pigauto_data object (output of preprocess_traits).")
  }
  if (!inherits(tree, "phylo")) stop("'tree' must be a phylo object.")
  cat_encoding <- match.arg(cat_encoding)

  X   <- data$X_scaled
  p   <- ncol(X)

  multi_obs <- isTRUE(data$multi_obs)
  if (multi_obs) {
    # Multi-obs: X is n_obs x p, species_names is n_species
    n_obs     <- data$n_obs
    n_species <- data$n_species
    spp       <- data$species_names          # unique species (n_species)
    obs_spp   <- data$obs_species            # species per obs (n_obs)
    obs_to_sp <- data$obs_to_species         # integer mapping (n_obs)
  } else {
    n_obs     <- nrow(X)
    n_species <- n_obs
    spp       <- data$species_names
    obs_spp   <- spp
    obs_to_sp <- NULL
  }

  # ---- Phylogenetic similarity for discrete-trait label propagation ------
  # Reuse build_phylo_graph()'s cached D if supplied; otherwise compute.
  # At n = 10,000 each cophenetic() call is ~15 seconds and ~800 MB of
  # allocation, so caching through `graph` is a meaningful speedup even
  # though this stage is not the dominant scaling bottleneck.
  if (!is.null(graph) && !is.null(graph$D)) {
    D_phylo <- graph$D
  } else {
    D_phylo <- ape::cophenetic.phylo(tree)
  }
  # Reorder to match species order
  D_phylo  <- D_phylo[spp, spp]
  sigma_lp <- stats::median(D_phylo) * 0.5
  sim_phylo <- exp(-(D_phylo^2) / (2 * sigma_lp^2))
  diag(sim_phylo) <- 0  # exclude self for label propagation

  # Mask val + test cells before fitting
  if (!is.null(splits)) {
    X[splits$val_idx]  <- NA
    X[splits$test_idx] <- NA
  }

  trait_map <- data$trait_map

  # Output at species level (n_species x p)
  mu <- matrix(0, nrow = n_species, ncol = p)
  se <- matrix(0, nrow = n_species, ncol = p)
  dimnames(mu) <- list(spp, data$latent_names)
  dimnames(se) <- list(spp, data$latent_names)

  # ---- Identify BM-eligible columns (continuous in latent space) -----------
  bm_cols <- integer(0)
  has_multi_proportion <- FALSE  # track for joint-dispatch guard
  zi_mag_fallback <- integer(0)  # ZI magnitude cols with too few non-zero obs
  for (tm in trait_map) {
    if (tm$type %in% c("continuous", "count", "ordinal", "proportion",
                       "multi_proportion")) {
      # multi_proportion: K independent BM fits, one per CLR column
      if (tm$type == "multi_proportion") has_multi_proportion <- TRUE
      bm_cols <- c(bm_cols, tm$latent_cols)
    } else if (tm$type == "zi_count") {
      # Magnitude column (col 2) is BM-eligible if enough non-zero obs
      mag_col <- tm$latent_cols[2]
      n_finite <- sum(is.finite(X[, mag_col]))
      if (n_finite >= 5L) {
        bm_cols <- c(bm_cols, mag_col)
      } else {
        # Fallback: constant imputation (global mean of non-zero values)
        zi_mag_fallback <- c(zi_mag_fallback, mag_col)
        finite_vals <- X[is.finite(X[, mag_col]), mag_col]
        mu[, mag_col] <- if (length(finite_vals) > 0) mean(finite_vals) else 0
        se[, mag_col] <- if (length(finite_vals) > 1) stats::sd(finite_vals) else 0
      }
    }
  }

  # ---- Level-C Phase 2, 3, 4 & 5: joint baseline dispatch -------------------
  # Binary and ZI-count gate cols both take the threshold-joint path (both
  # are binary-like observations with a truncated-Gaussian E-step).
  binary_cols  <- integer(0)
  zi_gate_cols <- integer(0)
  cat_cols     <- integer(0)
  for (tm in trait_map) {
    if (tm$type == "binary") {
      binary_cols <- c(binary_cols, tm$latent_cols)
    } else if (tm$type == "zi_count") {
      zi_gate_cols <- c(zi_gate_cols, tm$latent_cols[1])
    } else if (tm$type == "categorical") {
      cat_cols <- c(cat_cols, tm$latent_cols)
    }
  }
  # ZI gates join binary for dispatch purposes.
  binary_cols <- c(binary_cols, zi_gate_cols)

  # Phase 4 scope: threshold-joint fires when (binary) + BM are present.
  # Categorical-only (no binary) datasets fall back to Phase 2 MVN + LP:
  # Rphylopars has numerical instability with multi-categorical liability
  # matrices (the rank-(K-1) drop + multiple cat groups combine badly).
  # Phase 6 EM will refine this once Sigma is estimated stably.
  use_threshold_joint <- length(binary_cols) >= 1L &&
    length(bm_cols) >= 1L &&
    !has_multi_proportion && !multi_obs &&
    joint_mvn_available()

  use_continuous_joint <- !use_threshold_joint &&
    length(bm_cols) >= 2L &&
    !has_multi_proportion && !multi_obs &&
    joint_mvn_available()

  if (use_threshold_joint) {
    jt <- fit_joint_threshold_baseline(data, tree, splits = splits,
                                        graph = graph,
                                        cat_encoding = cat_encoding)

    populated_cols <- integer(0)

    # Continuous-family (mu_liab on z-score scale)
    cont_idx <- which(jt$liab_types != "binary" & jt$liab_types != "categorical")
    for (idx in cont_idx) {
      col <- jt$liab_cols[idx]
      if (any(!is.na(jt$mu_liab[, idx]))) {
        mu[, col] <- jt$mu_liab[, idx]
        se[, col] <- jt$se_liab[, idx]
        populated_cols <- c(populated_cols, col)
      }
    }

    # Binary -> logit(P)
    bin_idx <- which(jt$liab_types == "binary")
    for (idx in bin_idx) {
      col <- jt$liab_cols[idx]
      if (all(is.na(jt$mu_liab[, idx]))) next
      dec <- decode_binary_liability(mu_liab = jt$mu_liab[, idx],
                                      se_liab = jt$se_liab[, idx])
      mu[, col] <- dec$mu_logit
      se[, col] <- 0
      populated_cols <- c(populated_cols, col)
    }

    # Categorical -> K log-probs (both encodings, per trait)
    for (tm in trait_map) {
      if (tm$type != "categorical") next
      k_cols    <- tm$latent_cols
      local_idx <- match(k_cols, jt$liab_cols)
      if (any(is.na(local_idx))) next
      # If ALL K cols are all-NA, trait was filtered -> fall through to LP
      if (all(apply(jt$mu_liab[, local_idx, drop = FALSE], 2,
                    function(col) all(is.na(col))))) next
      for (i in seq_len(nrow(jt$mu_liab))) {
        mu_K <- jt$mu_liab[i, local_idx]
        se_K <- jt$se_liab[i, local_idx]
        dec  <- decode_categorical_liability(mu_K = mu_K, se_K = se_K,
                                              cat_encoding = cat_encoding)
        mu[i, k_cols] <- dec$log_probs
      }
      se[, k_cols]   <- 0
      populated_cols <- c(populated_cols, k_cols)
    }

    bm_cols     <- setdiff(bm_cols,     populated_cols)
    binary_cols <- setdiff(binary_cols, populated_cols)
    cat_cols    <- setdiff(cat_cols,    populated_cols)

  } else if (use_continuous_joint) {
    joint <- fit_joint_mvn_baseline(data, tree, splits = splits, graph = graph)
    mu[, bm_cols] <- joint$mu[, bm_cols]
    se[, bm_cols] <- joint$se[, bm_cols]
    bm_cols <- integer(0)
  }

  # ---- Internal BM imputation on BM-eligible columns -----------------------
  if (model == "OU") {
    message("OU not yet supported by the internal BM baseline; using BM. ",
            "Install Rphylopars for OU support.")
  }

  if (length(bm_cols) > 0) {
    # Retrieve or compute the phylogenetic correlation matrix
    if (!is.null(graph) && !is.null(graph$R_phy)) {
      R_phy <- graph$R_phy
    } else {
      R_phy <- phylo_cor_matrix(tree)
    }
    R_phy <- R_phy[spp, spp]

    # Aggregate multi-obs to species-level means for BM imputation
    if (multi_obs) {
      X_sp <- matrix(NA_real_, n_species, length(bm_cols))
      colnames(X_sp) <- colnames(X)[bm_cols]
      for (j in seq_along(bm_cols)) {
        col_vals <- X[, bm_cols[j]]
        sp_means <- tapply(col_vals, obs_spp, function(v) {
          v <- v[!is.na(v)]
          if (length(v) == 0L) NA_real_ else mean(v)
        })
        X_sp[match(names(sp_means), spp), j] <- as.numeric(sp_means)
      }
    } else {
      X_sp <- X[, bm_cols, drop = FALSE]
    }

    # Impute each BM-eligible column independently
    for (j in seq_along(bm_cols)) {
      res_j <- bm_impute_col(X_sp[, j], R_phy)
      mu[, bm_cols[j]] <- res_j$mu
      se[, bm_cols[j]] <- res_j$se
    }
  }

  # ---- Binary baseline: phylogenetic label propagation -------------------
  for (tm in trait_map) {
    if (tm$type != "binary") next
    lc   <- tm$latent_cols
    # If the Phase 3 threshold-joint populated this col, skip LP.
    # binary_cols is the set of UNpopulated binary latent cols after the
    # threshold dispatch; if our `lc` is not in binary_cols, the joint
    # fit already handled it.
    if (!all(lc %in% binary_cols)) next

    # Get species-level observations
    if (multi_obs) {
      sp_vals <- tapply(X[, lc], obs_spp, function(v) {
        v <- v[!is.na(v)]
        if (length(v) == 0) NA_real_ else mean(v)
      })
      vals_species <- rep(NA_real_, n_species)
      names(vals_species) <- spp
      vals_species[names(sp_vals)] <- as.numeric(sp_vals)
    } else {
      vals_species <- X[, lc]
      names(vals_species) <- spp
    }

    observed <- !is.na(vals_species)
    if (sum(observed) == 0) {
      mu[, lc] <- logit(0.5)
      se[, lc] <- 0
      next
    }

    # Phylo-weighted probability for each species
    sim_obs <- sim_phylo[, observed, drop = FALSE]
    row_weights <- rowSums(sim_obs)
    row_weights[row_weights < 1e-10] <- 1e-10
    probs <- as.numeric(sim_obs %*% vals_species[observed]) / row_weights
    probs <- pmin(pmax(probs, 0.01), 0.99)  # clip for stability

    mu[, lc] <- logit(probs)
    se[, lc] <- 0
  }

  # ---- Categorical baseline: phylogenetic label propagation ---------------
  for (tm in trait_map) {
    if (tm$type != "categorical") next
    K    <- tm$n_latent
    lc   <- tm$latent_cols
    if (!all(lc %in% cat_cols)) next  # handled by threshold-joint path
    oh   <- X[, lc, drop = FALSE]  # n_obs x K one-hot (with NAs)

    # Get species-level one-hot observations
    if (multi_obs) {
      # Average one-hot within species (handles multiple obs)
      oh_species <- matrix(NA_real_, n_species, K)
      for (s in seq_len(n_species)) {
        rows <- which(obs_to_sp == s)
        obs_rows <- rows[complete.cases(oh[rows, , drop = FALSE])]
        if (length(obs_rows) > 0) {
          oh_species[s, ] <- colMeans(oh[obs_rows, , drop = FALSE])
        }
      }
    } else {
      oh_species <- oh
    }

    # Which species have observed values
    observed <- complete.cases(oh_species)
    if (sum(observed) == 0) {
      freqs <- rep(1 / K, K)
      log_freqs <- log(freqs)
      for (k in seq_len(K)) mu[, lc[k]] <- log_freqs[k]
      se[, lc] <- 0
      next
    }

    # Phylo-weighted category probabilities per species
    sim_obs <- sim_phylo[, observed, drop = FALSE]
    row_weights <- rowSums(sim_obs)
    row_weights[row_weights < 1e-10] <- 1e-10

    # Weighted category probs: (n_species x K)
    weighted_probs <- (sim_obs %*% oh_species[observed, , drop = FALSE]) /
      row_weights
    # Add small floor and renormalise
    weighted_probs <- pmax(weighted_probs, 1e-6)
    weighted_probs <- weighted_probs / rowSums(weighted_probs)

    for (k in seq_len(K)) {
      mu[, lc[k]] <- log(weighted_probs[, k])
    }
    se[, lc] <- 0
  }

  # ---- ZI count gate baseline: phylogenetic label propagation ---------------
  for (tm in trait_map) {
    if (tm$type != "zi_count") next
    lc_gate <- tm$latent_cols[1]
    # Phase 5: if threshold-joint handled this gate, skip LP
    if (!(lc_gate %in% binary_cols)) next

    # Get species-level gate values (0 = zero, 1 = non-zero)
    if (multi_obs) {
      sp_vals <- tapply(X[, lc_gate], obs_spp, function(v) {
        v <- v[!is.na(v)]
        if (length(v) == 0) NA_real_ else mean(v)
      })
      vals_species <- rep(NA_real_, n_species)
      names(vals_species) <- spp
      vals_species[names(sp_vals)] <- as.numeric(sp_vals)
    } else {
      vals_species <- X[, lc_gate]
      names(vals_species) <- spp
    }

    observed <- !is.na(vals_species)
    if (sum(observed) == 0) {
      mu[, lc_gate] <- logit(0.5)
      se[, lc_gate] <- 0
      next
    }

    # Phylo-weighted non-zero probability
    sim_obs <- sim_phylo[, observed, drop = FALSE]
    row_weights <- rowSums(sim_obs)
    row_weights[row_weights < 1e-10] <- 1e-10
    probs <- as.numeric(sim_obs %*% vals_species[observed]) / row_weights
    probs <- pmin(pmax(probs, 0.01), 0.99)

    mu[, lc_gate] <- logit(probs)
    se[, lc_gate] <- 0
  }

  list(mu = mu, se = se)
}
