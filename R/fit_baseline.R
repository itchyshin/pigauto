#' Fit the phylogenetic baseline
#'
#' Dispatches to pigauto's phylogenetic baseline machinery and returns imputed
#' latent-scale means and standard errors for every species.
#'
#' @details
#' When \code{splits} is supplied the val and test cells are masked to
#' \code{NA} before fitting, so the baseline is evaluated under the same
#' conditions as \code{\link{fit_pigauto}}.
#'
#' Continuous-family columns use Brownian-motion conditional MVN baselines on
#' the phylogenetic correlation matrix, either independently or through the
#' joint MVN path when the data and optional dependencies support it. Binary,
#' ordinal, categorical, and zero-inflated gate columns use the appropriate
#' label-propagation or threshold/liability baseline candidates, with
#' per-column fallbacks when a joint path is not available.
#'
#' \strong{Covariates and the joint baseline (P1-8)}: \code{data$covariates}
#' is only used by the per-column BM path (\code{bm_impute_col_with_cov()}).
#' The joint MVN and threshold-joint (Rphylopars) baselines do not accept a
#' covariate design matrix, so when a joint path is selected (BM-eligible
#' columns >= 2, or binary/ordinal cols present, with Rphylopars available)
#' any supplied covariates are ignored for the BASELINE and a warning is
#' emitted; covariates still reach the GNN correction via
#' \code{\link{fit_pigauto}} regardless of which baseline path fires.
#'
#' \strong{Per-type lambda dispatch (arc/lambda-per-type)}: \code{lambda_mode}
#' only ever governs the baseline for CONTINUOUS-FAMILY columns (continuous,
#' count, ordinal, proportion, zi_count magnitude). Binary, ordinal,
#' categorical, and zero-inflated gate columns keep the threshold-joint /
#' OVR-categorical joint baseline (always fit at lambda = 1) regardless of
#' \code{lambda_mode} -- there is no discrete-trait analogue of Pagel's
#' lambda, and previously forcing these columns onto label propagation any
#' time \code{lambda_mode != "fixed_1"} cost 19pp of Trophic.Level accuracy
#' on AVONET (0.789 -> 0.600; see
#' \code{docs/dev-log/2026-08-16-external-comparison-results.md}). When
#' \code{lambda_mode != "fixed_1"} and the threshold-joint baseline fires
#' for a dataset with binary/ordinal AND continuous-family columns, the
#' joint liability fit still uses the continuous-family columns internally
#' to inform the joint Sigma (and hence the binary/ordinal posteriors); only
#' its continuous-column baseline OUTPUT is discarded in favour of the
#' lambda-aware per-column BM fit.
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
#'   observations per species before the Level-C joint baseline:
#'   \code{"hard"} (default) thresholds binary proportions at 0.5 and uses
#'   argmax for categorical, matching Phase 10 behaviour.  \code{"soft"}
#'   preserves species-level proportions and dispatches the truncated-Gaussian
#'   soft E-step (\code{estep_liability_binary_soft}) so that intermediate
#'   class frequencies contribute fractional liability evidence.  Only
#'   relevant for multi-obs data with binary or categorical traits when the
#'   Level-C joint baseline is active.
#' @param lambda_mode character. Pagel-lambda mode for the CONTINUOUS-FAMILY
#'   baseline (continuous, count, ordinal, proportion, zi_count magnitude
#'   columns only -- see \dQuote{Per-type lambda dispatch} in Details).
#'   \code{"fixed_1"} preserves the default Brownian correlation matrix;
#'   \code{"estimate"}, \code{"cv"}, and \code{"bayes"} delegate lambda
#'   handling to the per-column BM path. Binary/categorical/zi_gate columns
#'   are unaffected by \code{lambda_mode} and keep the threshold-joint /
#'   OVR-categorical baseline. \strong{Covariate caveat}: when
#'   \code{data$covariates} is supplied, the per-column path switches to
#'   \code{bm_impute_col_with_cov()}, which has no lambda argument and
#'   always fits at lambda = 1; \code{lambda_mode != "fixed_1"} is then
#'   silently ignored for BM-eligible columns and a warning is emitted.
#' @param em_iterations integer. Number of Phase 6 EM iterations for the
#'   threshold-joint baseline (binary + ordinal + OVR categorical). Default
#'   \code{0L} disables the EM loop and preserves v0.9.1 output byte-for-byte.
#'   When \code{>= 1}, the BM rate \eqn{\Sigma} learned by the in-house
#'   joint solver (\code{R/joint_mvn_solver.R}) at iteration \eqn{k} is fed back as the
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
#' @param joint_solver character. Which solver estimates the joint
#'   Sigma / posterior for the joint MVN, threshold-joint, and OVR
#'   categorical baselines. \code{"inhouse"} (default) uses the
#'   single-pass in-house solver (\code{R/joint_mvn_solver.R}) and is
#'   byte-identical to prior releases. \code{"rphylopars"} delegates to
#'   \code{Rphylopars::phylopars()}'s converged REML fit, which measured
#'   0.14-1.27 lower z-RMSE than the in-house solver on AVONET300 (see
#'   \code{docs/dev-log/2026-08-16-continuous-gap-diagnosis.md}); on
#'   failure or non-finite output it falls back to \code{"inhouse"} with
#'   a warning. Only affects the joint MVN / threshold-joint / OVR
#'   categorical paths above; ignored when those paths don't fire. Note
#'   that \code{lambda_mode != "fixed_1"} disables the continuous-only
#'   joint MVN path (it has no lambda argument) but no longer disables
#'   the threshold-joint / OVR-categorical paths -- see
#'   \dQuote{Per-type lambda dispatch} in Details.
#' @param joint_refine_iter integer, default \code{0L}. Enables
#'   cross-trait refinement of the joint baseline's cell imputations
#'   using the estimated Sigma (the in-house solver's \code{max_iter}
#'   EM cell-refinement; \code{R/joint_mvn_solver.R}). \code{0L} preserves
#'   current behaviour byte-for-byte. The refinement is guarded: the
#'   Sigma step must shrink each iteration, or the loop rolls back to the
#'   last good iterate and sets \code{$diverged}.
#'
#'   \strong{Where it helps, measured}
#'   (\code{docs/dev-log/2026-08-17-refinement-results.md}): at Pagel's
#'   lambda around 0.2, \code{joint_refine_iter = 3L} lowered simulated
#'   RMSE by ~13\% and moved SE coverage from 0.943 to 0.964. At lambda
#'   near 1 it changes essentially nothing: on AVONET300, whose four
#'   continuous traits sit at lambda 0.993-0.998, every effect was within
#'   one Monte Carlo standard error. Note the reason is that refinement
#'   has nothing to add at high signal -- \emph{not} that the guard
#'   intervenes; on that data the guard does not fire. Use it when
#'   phylogenetic signal is moderate or low; expect no benefit on
#'   strongly conserved traits.
#' @return A list with:
#'   \describe{
#'     \item{mu}{Numeric matrix (n_species x p_latent), baseline means in
#'       latent scale.}
#'     \item{se}{Numeric matrix (n_species x p_latent), standard errors.}
#'   }
#' @examples
#' \donttest{
#' data(avonet300, tree300, package = "pigauto")
#' tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
#' traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
#'                      c("Mass", "Wing.Length"), drop = FALSE]
#' rownames(traits) <- tree$tip.label
#' pd     <- preprocess_traits(traits, tree)
#' splits <- make_missing_splits(pd$X_scaled, trait_map = pd$trait_map)
#' bl     <- fit_baseline(pd, tree, splits)
#' }
#' @importFrom stats complete.cases rnorm rbinom
#' @export
fit_baseline <- function(data, tree, splits = NULL, model = "BM",
                         graph = NULL,
                         multi_obs_aggregation = c("hard", "soft"),
                         lambda_mode = c("fixed_1", "estimate", "cv", "bayes"),
                         em_iterations = 0L,
                         em_tol = 1e-3,
                         em_offdiag = FALSE,
                         joint_solver = c("inhouse", "rphylopars"),
                         joint_refine_iter = 0L) {
  multi_obs_aggregation <- match.arg(multi_obs_aggregation)
  soft_aggregate <- identical(multi_obs_aggregation, "soft")
  lambda_mode <- match.arg(lambda_mode)
  joint_solver <- match.arg(joint_solver)
  if (!is.numeric(joint_refine_iter) || length(joint_refine_iter) != 1L ||
      !is.finite(joint_refine_iter) ||
      joint_refine_iter != as.integer(joint_refine_iter) ||
      joint_refine_iter < 0L) {
    stop("'joint_refine_iter' must be a non-negative integer scalar.",
         call. = FALSE)
  }
  joint_refine_iter <- as.integer(joint_refine_iter)
  # Translate dispatcher mode to the kernel-layer lambda argument.
  # "fixed_1"  -> lambda = 1.0   (back-compat; bit-identical to v0.9.x)
  # "estimate" -> lambda = "estimate" (per-column ML)
  # "cv"       -> per-column CV-selected lambda, computed below
  bm_lambda <- switch(lambda_mode,
    "fixed_1"  = 1.0,
    "estimate" = "estimate",
    "cv"       = "cv",
    "bayes"    = "bayes",
    1.0
  )
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
  #
  # arc/lambda-per-type (2026-08): `lambda_mode` now governs ONLY where
  # CONTINUOUS-FAMILY columns get their baseline mu/se, not whether the
  # discrete-trait joint machinery runs at all. Previously
  # `lambda_mode != "fixed_1"` set `force_per_column <- TRUE`, which also
  # disabled `use_threshold_joint` and the OVR-categorical loop below --
  # i.e. binary/ordinal/categorical traits were pushed onto plain label
  # propagation any time a user asked for lambda estimation on their
  # continuous traits. Measured cost: on AVONET data this dropped
  # Trophic.Level (categorical) accuracy from 0.789 to 0.600 (19pp) while
  # lambda_mode="bayes" only ever improves continuous traits -- discrete
  # traits have no lambda concept and were always fit at lambda = 1
  # regardless (see docs/dev-log/2026-08-16-external-comparison-results.md
  # and NEWS.md). `force_per_column` below therefore now gates ONLY
  # `use_continuous_joint` (the continuous-only joint MVN path, which has
  # no lambda argument); `use_threshold_joint` and the OVR-categorical
  # dispatch further down are unconditionally eligible whenever their
  # other preconditions hold, independent of `lambda_mode`.
  #
  # This creates a genuine hybrid inside `use_threshold_joint`: the joint
  # liability fit still INCLUDES continuous-family columns (they inform
  # the joint Sigma that the binary/ordinal posteriors condition on), but
  # when `lambda_mode != "fixed_1"` we discard the joint fit's continuous-
  # column OUTPUT and let those columns fall through to the per-column
  # lambda-aware `bm_impute_col(..., lambda = bm_lambda)` path below (see
  # the `cont_idx` block a few lines down). Binary/ordinal/categorical/
  # zi_gate columns keep the joint/OVR baseline, fit at lambda = 1 as they
  # always have -- lambda_mode never touches them.
  force_per_column <- lambda_mode %in% c("estimate", "cv", "bayes")
  use_threshold_joint <- (length(binary_cols) + length(ordinal_cols)) >= 1L &&
    length(bm_cols) >= 1L &&
    !has_multi_proportion &&
    joint_mvn_available()

  use_continuous_joint <- !use_threshold_joint &&
    length(bm_cols) >= 2L &&
    !has_multi_proportion &&
    !force_per_column &&
    joint_mvn_available()

  # P1-8: the joint MVN / threshold-joint (Rphylopars) paths have no covariate
  # design matrix, so user covariates never reach the BASELINE when a joint
  # path fires. They still reach the GNN correction in fit_pigauto(). Warn
  # once rather than silently ignoring them.
  if (!is.null(data$covariates) && (use_threshold_joint || use_continuous_joint)) {
    warning(
      "'covariates' are ignored by the joint ",
      if (use_threshold_joint) "threshold-joint" else "MVN",
      " baseline (Rphylopars has no covariate design). The baseline is ",
      "covariate-free; covariates still reach the GNN correction. To get a ",
      "covariate-aware baseline, force the per-column path (e.g. supply a ",
      "single BM-eligible trait, or uninstall/mask Rphylopars).",
      call. = FALSE
    )
  }

  if (use_threshold_joint) {
    jt <- if (em_iterations >= 1L) {
      fit_joint_threshold_baseline_em(data, tree, splits = splits,
                                       graph = graph,
                                       soft_aggregate = soft_aggregate,
                                       em_iterations = em_iterations,
                                       em_tol = em_tol,
                                       em_offdiag = em_offdiag,
                                       joint_solver = joint_solver,
                                       joint_refine_iter = joint_refine_iter)
    } else {
      fit_joint_threshold_baseline(data, tree, splits = splits,
                                    graph = graph,
                                    soft_aggregate = soft_aggregate,
                                    joint_solver = joint_solver,
                                    joint_refine_iter = joint_refine_iter)
    }

    populated_cols <- integer(0)

    # Continuous-family passthrough (mu_liab on z-score scale).
    # Excludes binary (needs logit decode) and ordinal (needs threshold decode).
    #
    # arc/lambda-per-type hybrid: the joint liability fit above (`jt`) still
    # USED these continuous-family columns internally to estimate the joint
    # Sigma that the binary/ordinal posteriors condition on -- that's
    # unavoidable and desirable (it's the whole point of the joint model).
    # But its continuous-column mu/se are always fit at lambda = 1
    # (Rphylopars / the in-house solver have no lambda argument). When the
    # caller asked for lambda estimation (`lambda_mode != "fixed_1"`), we
    # discard that lambda=1 continuous OUTPUT here and leave these columns
    # OFF `populated_cols`, so they fall through to `bm_cols` and get fit by
    # the lambda-aware per-column path a few hundred lines down instead
    # (`bm_impute_col(..., lambda = bm_lambda)`). Binary/ordinal columns are
    # unaffected -- they are populated from `jt` below regardless of
    # `lambda_mode`.
    cont_idx <- which(!(jt$liab_types %in% c("binary", "categorical", "ordinal")))
    if (identical(lambda_mode, "fixed_1")) {
      for (idx in cont_idx) {
        col <- jt$liab_cols[idx]
        if (any(!is.na(jt$mu_liab[, idx]))) {
          mu[, col] <- jt$mu_liab[, idx]
          se[, col] <- jt$se_liab[, idx]
          populated_cols <- c(populated_cols, col)
        }
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
        bm_res <- bm_impute_col(X[, col], R_phy_local, lambda = bm_lambda)
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

        # ---- Phase F (2026-05-01): LP via K-class OVR --------------------
        # Run K binary LPs (one per ordinal class), normalise across K,
        # then compute E[class] = sum_k p_k * k on the integer-class scale
        # and z-score back to match threshold_joint / bm_mvn output scale.
        # Targets the K=3 ordinal regime (e.g. AVONET Migration) where
        # threshold_joint is systematically misspecified; see
        # `useful/MEMO_2026-04-29_phase6_migration_bisect.md` and
        # `specs/2026-05-01-phase-f-lp-corner-ordinal-design.md`.
        lp_pred_z <- NULL
        lp_se_z   <- NULL
        lp_mse    <- NA_real_
        # Locate the trait_map entry for this latent col.
        tm_ord <- NULL
        for (tm_x in trait_map) {
          if (identical(tm_x$type, "ordinal") && col %in% tm_x$latent_cols) {
            tm_ord <- tm_x; break
          }
        }
        if (!is.null(tm_ord) && !is.null(sim_phylo) &&
            length(tm_ord$levels) >= 2L &&
            is.finite(tm_ord$sd) && tm_ord$sd > 0) {
          K_ord  <- length(tm_ord$levels)
          z_vals <- X[, col]
          cls_int <- round(z_vals * tm_ord$sd + tm_ord$mean)
          obs_lp  <- !is.na(cls_int) &
                     cls_int >= 1L & cls_int <= K_ord
          if (sum(obs_lp) >= K_ord) {
            P_lp <- matrix(0.0, nrow = length(spp), ncol = K_ord)
            sim_obs <- sim_phylo[, obs_lp, drop = FALSE]
            rw      <- rowSums(sim_obs)
            rw[rw < 1e-10] <- 1e-10
            for (kk in seq_len(K_ord)) {
              y_k <- as.numeric(cls_int[obs_lp] == kk)
              P_lp[, kk] <- as.numeric(sim_obs %*% y_k) / rw
            }
            P_lp <- pmax(P_lp, 1e-6)
            P_lp <- P_lp / rowSums(P_lp)
            e_cls <- as.numeric(P_lp %*% seq_len(K_ord))
            var_cls <- vapply(seq_len(nrow(P_lp)), function(i) {
              sum(P_lp[i, ] * (seq_len(K_ord) - e_cls[i])^2)
            }, numeric(1L))
            lp_pred_z <- (e_cls - tm_ord$mean) / tm_ord$sd
            lp_se_z   <- sqrt(pmax(var_cls, 0)) / tm_ord$sd
            lp_diff <- lp_pred_z[val_rows_j[finite_t]] - truth_j[finite_t]
            lp_mse <- if (any(is.finite(lp_diff))) {
                        mean(lp_diff[is.finite(lp_diff)]^2)
                      } else NA_real_
          }
        }

        # Pick the lowest finite val MSE among the three options.
        mses <- c(threshold_joint = tj_mse,
                  bm_mvn          = bm_mse,
                  lp              = lp_mse)
        mses <- mses[is.finite(mses)]
        if (length(mses) == 0L) {
          ordinal_path_chosen[as.character(col)] <- "threshold_joint"
        } else {
          chosen <- names(mses)[which.min(mses)]
          if (chosen == "bm_mvn") {
            mu[, col] <- bm_res$mu
            se[, col] <- bm_res$se
          } else if (chosen == "lp" && !is.null(lp_pred_z)) {
            mu[, col] <- lp_pred_z
            se[, col] <- lp_se_z
          } # else: threshold_joint already in mu/se
          ordinal_path_chosen[as.character(col)] <- chosen
        }
      }
    }

  } else if (use_continuous_joint) {
    joint <- fit_joint_mvn_baseline(data, tree, splits = splits, graph = graph,
                                     soft_aggregate = soft_aggregate,
                                     joint_solver = joint_solver,
                                     joint_refine_iter = joint_refine_iter)
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
  #
  # arc/lambda-per-type: NOT gated on `force_per_column` -- categorical
  # traits have no lambda concept (OVR fits are threshold-joint calls,
  # always at lambda = 1), so `lambda_mode` never forces them onto label
  # propagation. See the dispatch comment above `use_threshold_joint`.
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
                                       em_tol = em_tol,
                                       joint_solver = joint_solver,
                                       joint_refine_iter = joint_refine_iter)
        } else {
          fit_ovr_categorical_fits(data, tree, trait_name = trait_name,
                                    splits = splits, graph = graph,
                                    soft_aggregate = soft_aggregate,
                                    joint_solver = joint_solver,
                                    joint_refine_iter = joint_refine_iter)
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
      # `bm_impute_col_with_cov()` has no lambda argument -- the
      # covariate-aware GLS path always runs at lambda = 1, so
      # lambda_mode = "estimate" / "cv" / "bayes" is silently ignored
      # for BM-eligible columns when covariates are supplied. Warn once
      # (not per column).
      if (lambda_mode != "fixed_1") {
        warning(
          "lambda_mode = '", lambda_mode, "' is not supported by the ",
          "covariate-aware BM baseline; bm_impute_col_with_cov() always ",
          "runs at lambda = 1 when covariates are supplied. Pagel's ",
          "lambda is ignored for BM-eligible columns in this fit.",
          call. = FALSE
        )
      }
    }

    # Impute each BM-eligible column (covariate-aware when cov_design supplied)
    for (j in seq_along(bm_cols)) {
      if (is.null(cov_design)) {
        res_j <- bm_impute_col(X_sp[, j], R_phy, lambda = bm_lambda)
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
