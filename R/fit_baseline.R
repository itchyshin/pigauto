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
#' @param multi_obs_aggregation character. How to aggregate multiple
#'   observations per species before the Level-C (Rphylopars) baseline:
#'   \code{"hard"} (default) thresholds binary proportions at 0.5 and uses
#'   argmax for categorical, matching Phase 10 behaviour.  \code{"soft"}
#'   preserves species-level proportions and dispatches the truncated-Gaussian
#'   soft E-step (\code{estep_liability_binary_soft}) so that intermediate
#'   class frequencies contribute fractional liability evidence.  Only
#'   relevant for multi-obs data with binary or categorical traits when the
#'   Level-C joint baseline is active.
#' @param em_iterations integer. Number of Phase 6 EM iterations for the
#'   threshold-joint baseline (binary + ordinal + OVR categorical). Default
#'   \code{0L} disables the EM loop and preserves v0.9.1 output byte-for-byte.
#'   When \code{>= 1}, the BM rate \eqn{\Sigma} learned by
#'   \code{Rphylopars::phylopars()} at iteration \eqn{k} is fed back as the
#'   per-trait prior SD at iteration \eqn{k+1}, up to \code{em_iterations}
#'   times or until \code{em_tol} convergence.  \code{em_iterations = 1L} is
#'   a degenerate single-pass run and produces the same baseline output as
#'   \code{0L}; \code{>= 2L} is needed for actual iteration. Only affects
#'   the threshold-joint path (continuous-only traits pass through the
#'   existing joint MVN path unchanged).
#' @param em_tol numeric. Relative-Frobenius convergence tolerance for the
#'   Phase 6 / 7 EM loop. Early-stops when
#'   \eqn{||\Sigma_k - \Sigma_{k-1}||_F / ||\Sigma_{k-1}||_F < }
#'   \code{em_tol}.  Default \code{1e-3}.
#' @param em_offdiag logical. Phase 7 opt-in: when \code{TRUE} AND
#'   \code{em_iterations >= 2L}, each liability cell's prior at iteration
#'   \eqn{k+1} is the conditional-MVN \eqn{(\mu, sd)} given the posterior
#'   liability of other traits at iteration \eqn{k}, using the full off-
#'   diagonal entries of \eqn{\Sigma}. Binary + ordinal only (OVR categorical
#'   stays on Phase 6 diagonal). Default \code{FALSE} preserves Phase 6
#'   behaviour.
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
                         multi_obs_aggregation = c("hard", "soft"),
                         em_iterations = 0L,
                         em_tol = 1e-3,
                         em_offdiag = FALSE) {
  multi_obs_aggregation <- match.arg(multi_obs_aggregation)
  soft_aggregate <- identical(multi_obs_aggregation, "soft")
  em_iterations <- as.integer(em_iterations)
  em_offdiag    <- isTRUE(em_offdiag)
  if (!is.finite(em_iterations) || em_iterations < 0L) {
    stop("'em_iterations' must be a non-negative integer.", call. = FALSE)
  }
  if (em_offdiag && em_iterations < 2L) {
    # Silent: em_offdiag has no effect at em=0 (no EM at all) or em=1
    # (plug-in path; no previous Σ to condition on).
    em_offdiag <- FALSE
  }
  if (!inherits(data, "pigauto_data")) {
    stop("'data' must be a pigauto_data object (output of preprocess_traits).")
  }
  if (!inherits(tree, "phylo")) stop("'tree' must be a phylo object.")

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
  ordinal_cols <- integer(0)
  for (tm in trait_map) {
    if (tm$type == "binary") {
      binary_cols <- c(binary_cols, tm$latent_cols)
    } else if (tm$type == "zi_count") {
      zi_gate_cols <- c(zi_gate_cols, tm$latent_cols[1])
    } else if (tm$type == "categorical") {
      cat_cols <- c(cat_cols, tm$latent_cols)
    } else if (tm$type == "ordinal") {
      ordinal_cols <- c(ordinal_cols, tm$latent_cols)
    }
  }
  # ZI gates join binary for dispatch purposes.
  binary_cols <- c(binary_cols, zi_gate_cols)

  # Phase 4 scope: threshold-joint fires when (binary) + BM are present.
  # Categorical-only (no binary) datasets fall back to Phase 2 MVN + LP:
  # Rphylopars has numerical instability with multi-categorical liability
  # matrices (the rank-(K-1) drop + multiple cat groups combine badly).
  # Phase 6 EM will refine this once Sigma is estimated stably.
  use_threshold_joint <- (length(binary_cols) + length(ordinal_cols)) >= 1L &&
    length(bm_cols) >= 1L &&
    !has_multi_proportion &&
    joint_mvn_available()

  use_continuous_joint <- !use_threshold_joint &&
    length(bm_cols) >= 2L &&
    !has_multi_proportion &&
    joint_mvn_available()


  if (use_threshold_joint) {
    jt <- if (em_iterations >= 1L) {
      fit_joint_threshold_baseline_em(data, tree, splits = splits,
                                       graph = graph,
                                       soft_aggregate = soft_aggregate,
                                       em_iterations = em_iterations,
                                       em_tol = em_tol,
                                       em_offdiag = em_offdiag)
    } else {
      fit_joint_threshold_baseline(data, tree, splits = splits,
                                    graph = graph,
                                    soft_aggregate = soft_aggregate)
    }

    populated_cols <- integer(0)

    # Continuous-family passthrough (mu_liab on z-score scale).
    # Excludes binary (needs logit decode) and ordinal (needs threshold decode).
    cont_idx <- which(!(jt$liab_types %in% c("binary", "categorical", "ordinal")))
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

    # Ordinal -> z-scored integer class via threshold decode
    ord_idx <- which(jt$liab_types == "ordinal")
    ordinal_threshold_populated <- integer(0)   # subset for path-selection below
    for (idx in ord_idx) {
      col <- jt$liab_cols[idx]
      if (all(is.na(jt$mu_liab[, idx]))) next
      # Find the trait_map entry for this ordinal col
      tm_ord <- NULL
      for (tm in trait_map) {
        if (tm$type == "ordinal" && col %in% tm$latent_cols) {
          tm_ord <- tm; break
        }
      }
      if (is.null(tm_ord)) next
      dec <- decode_ordinal_liability(mu_liab = jt$mu_liab[, idx],
                                        se_liab = jt$se_liab[, idx],
                                        tm = tm_ord)
      mu[, col] <- dec$mu_z
      se[, col] <- 0
      populated_cols <- c(populated_cols, col)
      ordinal_threshold_populated <- c(ordinal_threshold_populated, col)
    }

    bm_cols      <- setdiff(bm_cols,      populated_cols)
    binary_cols  <- setdiff(binary_cols,  populated_cols)
    ordinal_cols <- setdiff(ordinal_cols, populated_cols)

    # ---- Per-trait ordinal path selection (Opus #6, 2026-04-30) -----------
    # The threshold-joint path is theoretically more flexible than per-
    # column BM-via-MVN on z-scored integer class for ordinal traits, but
    # at small K (especially K=3, e.g. AVONET Migration) the K-1 thresholds
    # are pinned to a narrow band by phylopars EM and produce systematically
    # worse predictions than a Gaussian conditional MVN on z-scored
    # integers.  See `useful/MEMO_2026-04-29_phase6_migration_bisect.md`
    # for the bisect localising the regression to commit a541dbd.
    #
    # Rather than ship a `K <= 3 -> LP` heuristic, we compute BOTH paths
    # for each populated ordinal trait and pick the lower-val-MSE one
    # against the held-out val cells in `data$X_scaled`.  Single-obs only
    # for now (multi-obs would require species aggregation of the
    # alternative path; out of scope for this fix).
    ordinal_path_chosen <- character(0)
    if (length(ordinal_threshold_populated) > 0L &&
        !is.null(splits) && !multi_obs) {
      R_phy_local <- if (!is.null(graph) && !is.null(graph$R_phy)) {
        graph$R_phy[spp, spp]
      } else {
        phylo_cor_matrix(tree)[spp, spp]
      }
      # Linear-index decode helpers for splits$val_idx (integer indices
      # into the original n_obs x p_latent matrix).
      n_rows_sp <- nrow(data$X_scaled)
      val_idx   <- splits$val_idx
      val_col   <- ((val_idx - 1L) %/% n_rows_sp) + 1L
      val_row   <- ((val_idx - 1L) %% n_rows_sp) + 1L
      truth_full <- data$X_scaled

      for (col in ordinal_threshold_populated) {
        val_rows_j <- val_row[val_col == col]
        if (length(val_rows_j) == 0L) {
          ordinal_path_chosen[as.character(col)] <- "threshold_joint"
          next
        }
        # Threshold-joint prediction (currently in mu, possibly
        # species-level; for single-obs n_species == n_obs).
        tj_pred <- mu[, col]
        # BM-via-MVN alternative on the masked z-scored ordinal column.
        bm_res <- bm_impute_col(X[, col], R_phy_local)
        # Val MSE for both paths.
        truth_j  <- truth_full[val_rows_j, col]
        finite_t <- is.finite(truth_j)
        if (!any(finite_t)) {
          ordinal_path_chosen[as.character(col)] <- "threshold_joint"
          next
        }
        tj_diff  <- tj_pred[val_rows_j[finite_t]] - truth_j[finite_t]
        bm_diff  <- bm_res$mu[val_rows_j[finite_t]] - truth_j[finite_t]
        tj_mse   <- if (any(is.finite(tj_diff))) {
                      mean(tj_diff[is.finite(tj_diff)]^2)
                    } else NA_real_
        bm_mse   <- if (any(is.finite(bm_diff))) {
                      mean(bm_diff[is.finite(bm_diff)]^2)
                    } else NA_real_

        if (is.finite(bm_mse) && is.finite(tj_mse) && bm_mse < tj_mse) {
          mu[, col] <- bm_res$mu
          se[, col] <- bm_res$se
          ordinal_path_chosen[as.character(col)] <- "bm_mvn"
        } else {
          ordinal_path_chosen[as.character(col)] <- "threshold_joint"
        }
      }
    }

  } else if (use_continuous_joint) {
    joint <- fit_joint_mvn_baseline(data, tree, splits = splits, graph = graph,
                                     soft_aggregate = soft_aggregate)
    mu[, bm_cols] <- joint$mu[, bm_cols]
    se[, bm_cols] <- joint$se[, bm_cols]
    bm_cols <- integer(0)
  }

  # ---- Categorical -> K independent OVR fits (Phase 6) -------------------
  # Each categorical trait gets K separate threshold-joint fits, one per
  # class. This sidesteps the rank-(K-1) phylopars instability of the
  # single-fit approach and fires regardless of whether the threshold-joint
  # / continuous-joint dispatchers above ran. If phylopars is unavailable
  # OR a fit fails for any reason, the per-trait result falls through to LP
  # below.
  if (length(cat_cols) > 0L && joint_mvn_available()) {
    for (tm in trait_map) {
      if (tm$type != "categorical") next
      k_cols <- tm$latent_cols
      if (!all(k_cols %in% cat_cols)) next  # already handled
      # Extract the trait name from the "<name>=<level>" column names
      col_name_1 <- colnames(data$X_scaled)[k_cols[1]]
      trait_name <- sub("=.*$", "", col_name_1)
      probs <- tryCatch(
        if (em_iterations >= 1L) {
          fit_ovr_categorical_fits_em(data, tree, trait_name = trait_name,
                                       splits = splits, graph = graph,
                                       soft_aggregate = soft_aggregate,
                                       em_iterations = em_iterations,
                                       em_tol = em_tol)
        } else {
          fit_ovr_categorical_fits(data, tree, trait_name = trait_name,
                                    splits = splits, graph = graph,
                                    soft_aggregate = soft_aggregate)
        },
        error = function(e) NULL
      )
      if (is.null(probs)) next
      # If OVR came back all-NA (every class's fit failed), leave for LP.
      if (all(is.na(probs))) next
      log_probs <- decode_ovr_categorical(probs)
      mu[, k_cols] <- log_probs
      se[, k_cols] <- 0
      cat_cols <- setdiff(cat_cols, k_cols)
    }
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

    # ---- Covariate-aware design matrix (Fix G, 2026-04-25) ------------------
    # When `data$covariates` is non-NULL, the BM baseline becomes a GLS
    # regression on covariates: y = X*beta + u, u ~ MVN(0, sigma^2 * R).
    # This puts linear covariate effects into the BASELINE so the GNN's
    # delta only has to learn the nonlinear / interactive residuals.
    #
    # Without this, the GNN had to re-derive linear cov effects from
    # scratch through several non-linear layers + a regularised gate.
    # Empirically that converged to a much worse solution than direct
    # GLS regression (see useful/GNN_ARCHITECTURE_EXPLAINED.md).
    cov_design <- NULL
    if (!is.null(data$covariates)) {
      cov_mat <- as.matrix(data$covariates)
      if (multi_obs) {
        # Aggregate covariates to species level (mean across obs per species)
        cov_sp <- matrix(NA_real_, n_species, ncol(cov_mat))
        colnames(cov_sp) <- colnames(cov_mat)
        for (j in seq_len(ncol(cov_mat))) {
          sp_means <- tapply(cov_mat[, j], obs_spp, mean, na.rm = TRUE)
          cov_sp[match(names(sp_means), spp), j] <- as.numeric(sp_means)
        }
        cov_design <- cbind(intercept = 1, cov_sp)
      } else {
        cov_design <- cbind(intercept = 1, cov_mat)
      }
      # Replace any residual NAs (defensive): mean impute by column
      for (j in seq_len(ncol(cov_design))) {
        bad <- !is.finite(cov_design[, j])
        if (any(bad)) cov_design[bad, j] <- mean(cov_design[!bad, j])
      }
    }

    # Impute each BM-eligible column (covariate-aware when cov_design supplied)
    for (j in seq_along(bm_cols)) {
      if (is.null(cov_design)) {
        res_j <- bm_impute_col(X_sp[, j], R_phy)
      } else {
        res_j <- bm_impute_col_with_cov(X_sp[, j], cov_design, R_phy)
      }
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

  out <- list(mu = mu, se = se)
  if (exists("ordinal_path_chosen", inherits = FALSE) &&
      length(ordinal_path_chosen) > 0L) {
    out$ordinal_path_chosen <- ordinal_path_chosen
  }
  out
}
