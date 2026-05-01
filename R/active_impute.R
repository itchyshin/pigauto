# R/active_impute.R
#
# Active imputation: suggest which currently-missing cell, if observed
# next, would most reduce TOTAL predictive variance across all
# remaining missing cells.  Built on the BM conditional-MVN
# variance-update formula derived via a Sherman-Morrison rank-1
# inverse update.
#
# Inferential target: phylogenetically-informed sampling-design
# guidance for partially-observed datasets.  When you can measure k
# more species, which ones?  Higher delta_var_total means observing
# that cell most reduces uncertainty across the rest of your tree.
#
# Scope (v1, 2026-04-30):
#   * Continuous-family traits only (continuous, count, ordinal,
#     proportion).  These share the BM-via-conditional-MVN baseline
#     where the variance update has a closed form.
#   * Single-obs only.  Multi-obs would require species-level
#     aggregation of the variance update; deferred.
#   * Discrete traits (binary, categorical, zi_count) skipped silently
#     in the public API.  Discrete uncertainty reduction would use
#     expected-entropy-reduction (defer to v2).
#
# References:
#   Hansen, T.F. (1997). Stabilizing selection and the comparative
#     analysis of adaptation. Evolution 51(5):1341-1351.  Sec 2.3
#     gives the conditional-MVN BM formulas pigauto inherits.
#   Cohn, D.A. et al. (1996). Active learning with statistical models.
#     JAIR 4:129-145.  General framework for variance-reduction-based
#     active sampling.

# ---- bm_variance_reduction --------------------------------------------------
#
# Closed-form variance reduction (single-obs continuous BM):
#
# Given y (length n_species, NA at miss_idx) and R (n x n correlation
# matrix), for each currently-missing cell s_new in miss_idx, return:
#     delta_V(s_new) := sum_{i in miss_idx} ( var_i_current - var_i_after )
# where var_i_current = sigma2 * (1 - h_i) is the current BM posterior
# variance at i and var_i_after is the new variance after additionally
# observing s_new (at any value -- the variance update is value-free).
#
# Math (Sherman-Morrison rank-1 inverse update on adding row+col to
# R_oo):
#   R_oo' = [R_oo, r; r', 1]      where r = R[obs, s_new]
#   R_oo'^{-1}  = [R_oo^{-1} + uu'/alpha, -u/alpha; -u'/alpha, 1/alpha]
#   where u = R_oo^{-1} r, alpha = 1 - r' u = 1 - h_{s_new}.
# For another miss cell i, with v = R[i, obs]:
#   h_i' = (v, R[i, s_new]) R_oo'^{-1} (v, R[i, s_new])'
#        = h_i + (R[i, s_new] - v' u)^2 / alpha
# So Delta_h_i = (R[i, s_new] - h_i_to_s_new)^2 / alpha
#    Delta_var_i = sigma2 * Delta_h_i.
# Plus at s_new itself, current variance sigma2 * alpha drops to 0.
#
# Total over ALL currently-missing cells (including s_new):
#   delta_V(s_new) = sigma2 * (alpha + sum_{i != s_new} Delta_h_i)
#                   = sigma2 * sum_{i in miss} D[i, k]^2 / alpha
# where D[i, k] = R[miss[i], miss[k]] - h_to_full[i, k] is the residual
# matrix and D[k, k] = alpha by construction.
#
# Vectorised: U = R_oo^{-1} %*% R[obs, miss], h_to_full = t(U) %*% R[obs, miss],
# D = R[miss, miss] - h_to_full, output = sigma2 * colSums(D^2) / diag(D).
#
# Cost: O(n_o^2 n_m) for U and h_to_full, O(n_m^2) per query thereafter.
# Total complexity O(n_m^2 + n_o^2 n_m) per trait, dominated by the
# matrix solve at small n_m.
#
# @param y      numeric vector (length n_species).  NA marks missing.
# @param R      phylogenetic correlation matrix (n x n).
# @param nugget numeric scalar.  Added to diag(R[obs, obs]) before
#               solve(); same default as bm_impute_col().
# @return Named numeric vector of length(miss_idx) with total variance
#   reduction at each candidate cell.  Empty if n_obs < 5 or
#   length(miss_idx) == 0.
# @keywords internal
bm_variance_reduction <- function(y, R, nugget = 1e-6) {
  obs_idx  <- which(!is.na(y))
  miss_idx <- which(is.na(y))
  n_o <- length(obs_idx)
  n_m <- length(miss_idx)
  if (n_o < 5L || n_m == 0L) return(numeric(0L))

  R_oo  <- R[obs_idx, obs_idx, drop = FALSE]
  R_om  <- R[obs_idx, miss_idx, drop = FALSE]
  R_mm  <- R[miss_idx, miss_idx, drop = FALSE]

  # Cholesky with nugget back-off (matches bm_impute_col)
  L <- NULL
  nug <- nugget
  for (attempt in seq_len(6L)) {
    R_oo_reg <- R_oo + diag(nug, n_o)
    L <- tryCatch(chol(R_oo_reg), error = function(e) NULL)
    if (!is.null(L)) break
    nug <- nug * 2
  }
  if (is.null(L)) {
    warning("bm_variance_reduction: Cholesky failed after nugget back-off; ",
            "returning zeros", call. = FALSE)
    out <- rep(0, n_m); names(out) <- names(y)[miss_idx]
    return(out)
  }
  chol_solve <- function(b) backsolve(L, forwardsolve(t(L), b))

  # GLS phylogenetic mean + REML variance (same as bm_impute_col)
  ones    <- rep(1, n_o)
  y_o     <- y[obs_idx]
  a       <- chol_solve(ones)
  b       <- chol_solve(y_o)
  mu_hat  <- sum(b) / sum(a)
  e       <- y_o - mu_hat
  e_solve <- chol_solve(e)
  sigma2  <- as.numeric(crossprod(e, e_solve)) / max(n_o - 1L, 1L)

  # Closed-form variance update
  U        <- chol_solve(R_om)              # n_o x n_m, U = R_oo^{-1} R_om
  h_full   <- crossprod(R_om, U)            # n_m x n_m, t(U) R_om = R_mo R_oo^{-1} R_om
  D        <- R_mm - h_full                 # residual matrix
  alpha    <- diag(D)                       # = 1 - h_self for each miss cell
  alpha    <- pmax(alpha, .Machine$double.eps)  # guard against tiny / zero

  total_reduction <- as.numeric(sigma2 * colSums(D^2) / alpha)
  total_reduction <- pmax(total_reduction, 0)   # numerical floor

  if (!is.null(names(y))) {
    names(total_reduction) <- names(y)[miss_idx]
  }
  total_reduction
}

# ---- lp_entropy_reduction_binary --------------------------------------------
#
# Closed-form expected entropy reduction (single-obs binary trait, LP).
#
# Given y in {0, 1, NA} and sim (n x n similarity matrix, diag = 0 by
# convention -- self-exclusion for label propagation), for each
# currently-missing cell s_new in miss_idx, return the expected total
# entropy reduction across all currently-missing cells if s_new were
# observed next.
#
# The current LP probability at species i is:
#     p_i = (sim[i, obs] %*% y_obs) / sum(sim[i, obs])
# Entropy at i is H(p_i) = -p log p - (1 - p) log(1 - p) (binary).
#
# If we imagine observing s_new with unknown value y_new, the new LP
# probability is:
#     p_i_new = (sim[i, obs] %*% y_obs + sim[i, s_new] * y_new) /
#                (row_w[i] + sim[i, s_new])
# We don't know y_new ahead of time, so we use the current LP estimate
# at s_new (call it q = p_s_new) as P(y_new = 1).  The expected new
# entropy at i is:
#     E[H_i | obs s_new] = (1 - q) * H(p_i_y0) + q * H(p_i_y1)
# where p_i_y0, p_i_y1 are the LP updates with y_new = 0, 1 respectively.
#
# Total expected entropy reduction = current_total_H - expected_new_total_H,
# accounting for the fact that observing s_new drops its own entropy to 0.
#
# Cost: O(n) per candidate, O(n * n_m) total.
#
# @param y    numeric vector (length n).  NA = missing; non-NA values must
#             be 0 or 1.
# @param sim  n x n similarity matrix.  diag = 0 by LP convention.
# @return Named numeric vector of length(miss_idx) with total entropy
#   reduction at each candidate cell.  Empty if no observed or no
#   missing.  Non-negative (numerical floor at 0).
# @keywords internal
lp_entropy_reduction_binary <- function(y, sim) {
  obs_idx  <- which(!is.na(y))
  miss_idx <- which(is.na(y))
  n_o <- length(obs_idx); n_m <- length(miss_idx)
  if (n_o == 0L || n_m == 0L) return(numeric(0L))

  H <- function(p) {
    # Vectorised binary entropy with H(0) = H(1) = 0
    out <- rep(0, length(p))
    keep <- p > 0 & p < 1
    out[keep] <- -p[keep] * log(p[keep]) - (1 - p[keep]) * log(1 - p[keep])
    out
  }

  # Current LP probabilities at every species
  sim_obs <- sim[, obs_idx, drop = FALSE]                    # n x n_o
  numerator <- as.numeric(sim_obs %*% y[obs_idx])           # length n
  row_w   <- rowSums(sim_obs)                                # length n
  row_w_safe <- pmax(row_w, 1e-12)
  p_all <- numerator / row_w_safe
  p_all <- pmin(pmax(p_all, 0.01), 0.99)   # match LP path's clamp
  H_cur <- H(p_all[miss_idx])

  total_H_cur <- sum(H_cur)

  out <- numeric(n_m)
  for (k in seq_len(n_m)) {
    s_new <- miss_idx[k]
    q <- p_all[s_new]                          # current LP estimate
    sim_to_s_new <- sim[, s_new]               # length n
    denom_new <- row_w + sim_to_s_new
    denom_new <- pmax(denom_new, 1e-12)
    # LP after observing s_new with y_new = 0 / 1
    p_y0 <- numerator / denom_new
    p_y1 <- (numerator + sim_to_s_new) / denom_new
    p_y0 <- pmin(pmax(p_y0, 0.01), 0.99)
    p_y1 <- pmin(pmax(p_y1, 0.01), 0.99)
    # Expected entropy at miss cells (other than s_new -- s_new drops to 0)
    miss_other <- miss_idx[miss_idx != s_new]
    if (length(miss_other) == 0L) {
      sum_E <- 0
    } else {
      E_H <- (1 - q) * H(p_y0[miss_other]) + q * H(p_y1[miss_other])
      sum_E <- sum(E_H)
    }
    out[k] <- total_H_cur - sum_E
  }
  out <- pmax(out, 0)   # numerical floor
  if (!is.null(names(y))) names(out) <- names(y)[miss_idx]
  out
}

# ---- lp_entropy_reduction_categorical ---------------------------------------
#
# K-ary version of the binary entropy-reduction formula.  Given a
# K-class one-hot matrix `oh` (n x K, NA in entire row = missing) and
# similarity matrix `sim`, return expected total entropy reduction
# across all currently-missing rows for each candidate.
#
# Current K-class LP at species i:
#     p_i_k = (sim[i, obs] %*% oh_obs[, k]) / row_w[i]
#     where oh_obs is the one-hot matrix at observed rows.
# Entropy: H(p_i) = -sum_k p_i_k log p_i_k (with floors).
#
# After observing s_new with class k_new:
#     p_i_k_new = (numerator[i, k] + (k == k_new) * sim[i, s_new]) /
#                  (row_w[i] + sim[i, s_new])
# Expected entropy: sum_k_new q_k_new * H(p_i | y_new = k_new), where
# q_k_new = current LP estimate at s_new.
#
# Cost: O(n * K) per candidate, O(n * n_m * K) total.
#
# @param oh   n x K one-hot matrix (numeric 0/1, NA in entire row = missing)
# @param sim  n x n similarity matrix
# @return Named numeric vector of length(miss_rows) with total entropy
#   reduction.
# @keywords internal
lp_entropy_reduction_categorical <- function(oh, sim) {
  K <- ncol(oh)
  n <- nrow(oh)
  obs_idx  <- which(stats::complete.cases(oh))
  miss_idx <- setdiff(seq_len(n), obs_idx)
  n_o <- length(obs_idx); n_m <- length(miss_idx)
  if (n_o == 0L || n_m == 0L) return(numeric(0L))

  H_cat <- function(P) {
    # Row-wise K-ary entropy: P is m x K with rows summing to 1
    P <- pmax(P, 1e-12)
    P <- P / rowSums(P)
    -rowSums(P * log(P))
  }

  sim_obs <- sim[, obs_idx, drop = FALSE]            # n x n_o
  numerator <- sim_obs %*% oh[obs_idx, , drop = FALSE]  # n x K
  row_w <- rowSums(sim_obs)                          # length n
  row_w_safe <- pmax(row_w, 1e-12)
  p_all <- numerator / row_w_safe
  p_all <- pmax(p_all, 1e-12)
  p_all <- p_all / rowSums(p_all)                    # row-stochastic

  H_cur <- H_cat(p_all[miss_idx, , drop = FALSE])
  total_H_cur <- sum(H_cur)

  out <- numeric(n_m)
  for (k in seq_len(n_m)) {
    s_new <- miss_idx[k]
    q <- p_all[s_new, ]                              # K-vector
    sim_to_s_new <- sim[, s_new]                     # length n
    denom_new <- row_w + sim_to_s_new
    denom_new <- pmax(denom_new, 1e-12)
    miss_other <- miss_idx[miss_idx != s_new]
    if (length(miss_other) == 0L) {
      out[k] <- total_H_cur
      next
    }
    sum_E <- 0
    for (kc in seq_len(K)) {
      # Observing s_new in class kc
      num_kc <- numerator
      num_kc[, kc] <- num_kc[, kc] + sim_to_s_new
      P_kc <- num_kc / denom_new
      P_kc <- pmax(P_kc, 1e-12)
      P_kc <- P_kc / rowSums(P_kc)
      H_kc <- H_cat(P_kc[miss_other, , drop = FALSE])
      sum_E <- sum_E + q[kc] * sum(H_kc)
    }
    out[k] <- total_H_cur - sum_E
  }
  out <- pmax(out, 0)
  if (!is.null(rownames(oh))) names(out) <- rownames(oh)[miss_idx]
  out
}

# ---- suggest_next_observation -----------------------------------------------

#' Suggest which cell to observe next to maximise imputation precision
#'
#' For a fitted \code{pigauto_result} (returned by \code{\link{impute}}),
#' compute the closed-form expected reduction in total predictive
#' variance across all currently-missing cells if each candidate cell
#' were observed next.  Useful for sampling-design guidance: when you
#' have time/budget to measure \emph{k} more species, this function
#' tells you which ones contribute most to imputation precision.
#'
#' @details
#' Two metrics are supported, dispatched by trait type:
#'
#' \strong{Continuous-family traits} (continuous, count, ordinal,
#' proportion) use the BM variance-reduction formula.  The variance-
#' reduction formula is derived from a Sherman-Morrison rank-1 inverse
#' update on the BM conditional MVN: adding species \eqn{s} to the
#' observed set updates the inverse correlation matrix by a known
#' closed form.  For each candidate cell \eqn{(s, t)},
#' \deqn{
#'   \Delta V(s, t) = \sigma_t^2
#'     \sum_{i \in \mathrm{miss}_t} \frac{D_{ik}^2}{\alpha_k}
#' }
#' where \eqn{D = R_{mm} - R_{mo} R_{oo}^{-1} R_{om}} is the residual
#' matrix at currently-missing cells, \eqn{\alpha_k = D_{kk}} is the
#' current relative leverage of cell \eqn{k}, and \eqn{\sigma_t^2} is
#' the REML BM variance for trait \eqn{t}.
#'
#' \strong{Discrete traits} (binary, categorical) use a label-
#' propagation expected-entropy-reduction formula.  The current LP
#' probability at species \eqn{i} is
#' \eqn{p_i = \mathrm{sim}[i, \mathrm{obs}] y_{\mathrm{obs}} /
#' \sum \mathrm{sim}[i, \mathrm{obs}]}, with entropy
#' \eqn{H(p_i) = -\sum_k p_{i,k} \log p_{i,k}}.  After observing
#' \eqn{s_{\mathrm{new}}} with unknown class \eqn{y_{\mathrm{new}}},
#' the new LP probability has a closed form, and the expected entropy
#' is averaged over \eqn{P(y_{\mathrm{new}})} = current LP estimate at
#' \eqn{s_{\mathrm{new}}}.  Total expected entropy reduction sums
#' across all currently-missing cells (the entropy at
#' \eqn{s_{\mathrm{new}}} itself drops to 0).
#'
#' \strong{Variance and entropy are NOT directly comparable.}  The
#' output sorts within each metric and the cross-metric ordering by
#' \code{delta} is approximate.  When you want a strict ranking,
#' filter by \code{metric} first.
#'
#' Reductions are summed across the included traits for each species
#' when \code{by = "species"}, supporting the typical use case where
#' measuring a species observes all of its currently-missing traits
#' at once.  At \code{by = "species"}, the per-trait variance and
#' entropy reductions are summed separately into
#' \code{delta_var_total} and \code{delta_entropy_total} columns; the
#' \code{delta} column is whichever is non-NA (or
#' \code{delta_var_total} when both are populated).  Cross-type
#' species-level ranking is approximate -- see the variance-vs-
#' entropy caveat above.
#'
#' \code{zi_count} (v2): observing a missing zi_count cell reveals
#' the gate value (entropy reduction at the gate column, computed via
#' the LP binary formula) AND, with probability \eqn{p_{\mathrm{gate}}}
#' (current LP estimate at \eqn{s_{\mathrm{new}}}), reveals a
#' magnitude (variance reduction at the magnitude column, computed
#' via the BM Sherman-Morrison formula on the gate=1 subset).  Output
#' rows for zi_count populate BOTH \code{delta_var_total} (= expected
#' magnitude variance reduction = \eqn{p_{\mathrm{gate}} \times \Delta
#' V_{\mathrm{mag}}}) AND \code{delta_entropy_total} (= gate entropy
#' reduction).  \code{metric} is set to \code{"variance"} so the row
#' sorts on the magnitude scale; \code{delta_entropy_total} is
#' available for users who care about gate-uncertainty separately.
#'
#' \code{multi_proportion} (v2): observing a row reveals all K
#' simplex components simultaneously.  Per-component variance
#' reductions are computed via BM Sherman-Morrison on each CLR-z
#' latent column, summed across components.  \code{metric} is
#' \code{"variance"}; \code{delta_var_total} is the K-component sum.
#'
#' @param result A \code{pigauto_result} object returned by
#'   \code{\link{impute}}.  Must have been produced from single-obs
#'   data; multi-obs inputs error with a clear message.
#' @param top_n integer, default \code{10L}.  Number of suggestions to
#'   return (descending by \code{delta}).
#' @param by character, one of \code{"cell"} (default) or
#'   \code{"species"}.  \code{"cell"} returns individual
#'   \code{(species, trait)} pairs.  \code{"species"} aggregates by
#'   species (summing reductions across the species' currently-missing
#'   traits).
#' @param types character vector of pigauto trait types to include.
#'   Default includes all eight supported types: \code{continuous},
#'   \code{count}, \code{ordinal}, \code{proportion}, \code{binary},
#'   \code{categorical}, \code{zi_count} (added v2, 2026-05-01),
#'   \code{multi_proportion} (added v2).
#' @return A data.frame of class \code{"pigauto_active"}.  Columns
#'   when \code{by = "cell"}: \code{species}, \code{trait},
#'   \code{type}, \code{metric} (\code{"variance"} or \code{"entropy"}),
#'   \code{delta}, \code{delta_var_total} (NA for discrete rows), and
#'   \code{delta_entropy_total} (NA for continuous rows), sorted by
#'   \code{delta} descending.  When \code{by = "species"}:
#'   \code{species}, \code{delta_var_total}, \code{delta_entropy_total},
#'   \code{n_traits_missing}, sorted by the SUM of available metrics.
#' @examples
#' \dontrun{
#' data(avonet300, tree300, package = "pigauto")
#' res <- impute(avonet300, tree300)
#' suggest_next_observation(res, top_n = 5)              # top-5 cells
#' suggest_next_observation(res, top_n = 10, by = "species")  # top-10 species
#'
#' # Continuous only:
#' suggest_next_observation(res, top_n = 10,
#'   types = c("continuous", "count", "ordinal", "proportion"))
#' }
#' @seealso \code{\link{impute}}
#' @export
suggest_next_observation <- function(result, top_n = 10L,
                                       by = c("cell", "species"),
                                       types = c("continuous", "count",
                                                 "ordinal", "proportion",
                                                 "binary", "categorical",
                                                 "zi_count",
                                                 "multi_proportion")) {
  by <- match.arg(by)
  if (!inherits(result, "pigauto_result")) {
    stop("'result' must be a pigauto_result object (from impute()).",
         call. = FALSE)
  }
  if (!is.numeric(top_n) || length(top_n) != 1L || top_n < 1L) {
    stop("'top_n' must be a single positive integer.", call. = FALSE)
  }
  top_n <- as.integer(top_n)

  data <- result$data
  fit  <- result$fit

  if (isTRUE(data$multi_obs)) {
    stop("suggest_next_observation: multi_obs input not yet supported. ",
         "v1 ships single-obs only; multi-obs deferred to v2.",
         call. = FALSE)
  }

  X <- data$X_scaled
  trait_map <- data$trait_map
  spp <- data$species_names

  cont_types <- c("continuous", "count", "ordinal", "proportion")
  disc_types <- c("binary", "categorical")
  # v2 (2026-05-01): zi_count needs BOTH R_phy (magnitude) and sim_phylo
  # (gate); multi_proportion needs R_phy only (per-component BM).
  has_cont <- any(types %in% cont_types) &&
              any(vapply(trait_map,
                          function(tm) tm$type %in% intersect(types, cont_types),
                          logical(1)))
  has_disc <- any(types %in% disc_types) &&
              any(vapply(trait_map,
                          function(tm) tm$type %in% intersect(types, disc_types),
                          logical(1)))
  has_zi   <- "zi_count" %in% types &&
              any(vapply(trait_map, function(tm) identical(tm$type, "zi_count"),
                          logical(1)))
  has_mp   <- "multi_proportion" %in% types &&
              any(vapply(trait_map,
                          function(tm) identical(tm$type, "multi_proportion"),
                          logical(1)))
  needs_R_phy    <- has_cont || has_zi || has_mp
  needs_sim_phylo <- has_disc || has_zi

  # Recover R_phy: needed for continuous variance reduction (continuous
  # family + zi_count magnitude + multi_proportion per-component).
  R_phy <- NULL
  if (needs_R_phy) {
    R_phy <- fit$graph$R_phy
    if (is.null(R_phy) && !is.null(result$tree)) {
      R_phy <- phylo_cor_matrix(result$tree)
    }
    if (is.null(R_phy)) {
      stop("suggest_next_observation: phylogenetic correlation matrix not ",
           "available in fit$graph$R_phy and tree not in result. Refit via ",
           "impute() to populate.", call. = FALSE)
    }
    R_phy <- R_phy[spp, spp, drop = FALSE]
  }

  # Recover sim_phylo (Gaussian kernel on cophenetic distances): needed
  # for discrete entropy reduction (binary + categorical + zi_count gate).
  # Reconstruct from result$tree (the similarity matrix is not cached
  # in fit$graph after impute()'s cleanup block).
  sim_phylo <- NULL
  if (needs_sim_phylo) {
    if (is.null(result$tree)) {
      stop("suggest_next_observation: result$tree is missing; cannot ",
           "reconstruct similarity matrix for discrete traits. Refit ",
           "via impute() (>=v0.9.1.9008 stores tree in result).",
           call. = FALSE)
    }
    D_phylo <- ape::cophenetic.phylo(result$tree)
    D_phylo <- D_phylo[spp, spp]
    sigma_lp <- stats::median(D_phylo) * 0.5
    sim_phylo <- exp(-(D_phylo^2) / (2 * sigma_lp^2))
    diag(sim_phylo) <- 0   # LP convention: exclude self
  }

  # Per-trait loop: dispatch by type
  cell_dfs <- list()
  for (tm in trait_map) {
    if (!(tm$type %in% types)) next

    if (tm$type %in% cont_types) {
      j <- tm$latent_cols[1L]
      y <- X[, j]
      names(y) <- spp
      if (sum(!is.na(y)) < 5L || sum(is.na(y)) == 0L) next

      delta_v <- bm_variance_reduction(y, R_phy)
      if (length(delta_v) == 0L) next

      miss_idx <- which(is.na(y))
      cell_dfs[[length(cell_dfs) + 1L]] <- data.frame(
        species              = spp[miss_idx],
        trait                = tm$name,
        type                 = tm$type,
        metric               = "variance",
        delta                = as.numeric(delta_v),
        delta_var_total      = as.numeric(delta_v),
        delta_entropy_total  = NA_real_,
        stringsAsFactors = FALSE)

    } else if (identical(tm$type, "binary")) {
      j <- tm$latent_cols[1L]
      y <- X[, j]
      names(y) <- spp
      if (sum(!is.na(y)) < 1L || sum(is.na(y)) == 0L) next

      delta_e <- lp_entropy_reduction_binary(y, sim_phylo)
      if (length(delta_e) == 0L) next

      miss_idx <- which(is.na(y))
      cell_dfs[[length(cell_dfs) + 1L]] <- data.frame(
        species              = spp[miss_idx],
        trait                = tm$name,
        type                 = tm$type,
        metric               = "entropy",
        delta                = as.numeric(delta_e),
        delta_var_total      = NA_real_,
        delta_entropy_total  = as.numeric(delta_e),
        stringsAsFactors = FALSE)

    } else if (identical(tm$type, "categorical")) {
      lc <- tm$latent_cols
      oh <- X[, lc, drop = FALSE]
      rownames(oh) <- spp
      # rows where any latent col is NA -> entire row missing
      row_obs <- stats::complete.cases(oh)
      if (sum(row_obs) < 1L || sum(!row_obs) == 0L) next

      delta_e <- lp_entropy_reduction_categorical(oh, sim_phylo)
      if (length(delta_e) == 0L) next

      miss_idx <- which(!row_obs)
      cell_dfs[[length(cell_dfs) + 1L]] <- data.frame(
        species              = spp[miss_idx],
        trait                = tm$name,
        type                 = tm$type,
        metric               = "entropy",
        delta                = as.numeric(delta_e),
        delta_var_total      = NA_real_,
        delta_entropy_total  = as.numeric(delta_e),
        stringsAsFactors = FALSE)

    } else if (identical(tm$type, "zi_count")) {
      # v2 (2026-05-01): hybrid variance + entropy.
      # Layout: lc[1] = gate (0/1, NA at user-missing), lc[2] = magnitude
      # (z-scored log1p, NA when gate=0 OR user-missing).
      lc_gate <- tm$latent_cols[1L]
      lc_mag  <- tm$latent_cols[2L]
      y_gate  <- X[, lc_gate]
      y_mag   <- X[, lc_mag]
      names(y_gate) <- spp
      names(y_mag)  <- spp

      miss_gate <- which(is.na(y_gate))    # user-missing cells
      if (length(miss_gate) == 0L)        next
      if (sum(!is.na(y_gate)) < 1L)       next

      # Gate entropy reduction (binary LP)
      delta_e <- lp_entropy_reduction_binary(y_gate, sim_phylo)
      if (length(delta_e) != length(miss_gate)) next  # defensive

      # Probability gate is 1 at user-missing cells (current LP estimate)
      obs_gate_idx <- which(!is.na(y_gate))
      sim_obs_g    <- sim_phylo[, obs_gate_idx, drop = FALSE]
      row_w_g      <- rowSums(sim_obs_g)
      row_w_g_safe <- pmax(row_w_g, 1e-12)
      p_gate_all   <- as.numeric(sim_obs_g %*% y_gate[obs_gate_idx]) / row_w_g_safe
      p_gate_all   <- pmin(pmax(p_gate_all, 0.01), 0.99)
      p_gate_user  <- p_gate_all[miss_gate]

      # Magnitude variance reduction (BM Sherman-Morrison) on the gate=1
      # subset.  bm_variance_reduction returns one value per
      # magnitude-missing cell; we pick the values for user-missing cells
      # (subset of magnitude-missing).
      delta_v_full <- bm_variance_reduction(y_mag, R_phy)
      miss_mag <- which(is.na(y_mag))
      if (length(delta_v_full) != length(miss_mag)) {
        # bm_variance_reduction returned empty (e.g. n_obs < 5); skip
        # variance contribution but still report entropy
        var_red_user <- rep(0, length(miss_gate))
      } else {
        var_red_user <- delta_v_full[match(miss_gate, miss_mag)]
        var_red_user[is.na(var_red_user)] <- 0   # not in mag-miss (impossible)
      }
      delta_v_expected <- p_gate_user * var_red_user

      cell_dfs[[length(cell_dfs) + 1L]] <- data.frame(
        species              = spp[miss_gate],
        trait                = tm$name,
        type                 = tm$type,
        metric               = "variance",
        delta                = as.numeric(delta_v_expected),
        delta_var_total      = as.numeric(delta_v_expected),
        delta_entropy_total  = as.numeric(delta_e),
        stringsAsFactors = FALSE)

    } else if (identical(tm$type, "multi_proportion")) {
      # v2 (2026-05-01): per-component BM variance reduction summed
      # across K components.  Per CLAUDE.md, multi_proportion has K
      # latent CLR-z-scored cols and rows are EITHER fully observed OR
      # fully missing (no partial observations of K-component
      # composition).  So miss_idx is well-defined as rows with any
      # latent col missing (= all latent cols missing).
      lc <- tm$latent_cols
      n_components <- length(lc)
      X_mp <- X[, lc, drop = FALSE]
      rownames(X_mp) <- spp
      row_obs <- stats::complete.cases(X_mp)
      if (sum(row_obs) < 5L || sum(!row_obs) == 0L) next

      delta_v_per_component <- numeric(0)
      for (k in seq_along(lc)) {
        y_k <- X_mp[, k]
        names(y_k) <- spp
        delta_k <- bm_variance_reduction(y_k, R_phy)
        if (length(delta_k) == 0L) {
          delta_v_per_component <- numeric(0); break
        }
        if (length(delta_v_per_component) == 0L) {
          delta_v_per_component <- delta_k
        } else {
          # Component-wise sum (each component has the same miss_idx by
          # the row-level masking convention)
          if (length(delta_k) != length(delta_v_per_component)) {
            delta_v_per_component <- numeric(0); break
          }
          delta_v_per_component <- delta_v_per_component + delta_k
        }
      }
      if (length(delta_v_per_component) == 0L) next

      miss_idx <- which(!row_obs)
      cell_dfs[[length(cell_dfs) + 1L]] <- data.frame(
        species              = spp[miss_idx],
        trait                = tm$name,
        type                 = tm$type,
        metric               = "variance",
        delta                = as.numeric(delta_v_per_component),
        delta_var_total      = as.numeric(delta_v_per_component),
        delta_entropy_total  = NA_real_,
        stringsAsFactors = FALSE)
    }
  }

  if (length(cell_dfs) == 0L) {
    out <- data.frame(species = character(0L),
                      trait = character(0L),
                      type = character(0L),
                      metric = character(0L),
                      delta = numeric(0L),
                      delta_var_total = numeric(0L),
                      delta_entropy_total = numeric(0L),
                      stringsAsFactors = FALSE)
    return(structure(out, class = c("pigauto_active", "data.frame"),
                     by = by))
  }

  all_cells <- do.call(rbind, cell_dfs)

  out <- if (identical(by, "cell")) {
    all_cells <- all_cells[order(-all_cells$delta), , drop = FALSE]
    rownames(all_cells) <- NULL
    utils::head(all_cells, top_n)
  } else {
    # Aggregate to species: sum of variance + entropy per species
    agg_var <- stats::aggregate(delta_var_total ~ species,
                                  data = all_cells,
                                  FUN = function(x) sum(x, na.rm = TRUE))
    agg_ent <- stats::aggregate(delta_entropy_total ~ species,
                                  data = all_cells,
                                  FUN = function(x) sum(x, na.rm = TRUE))
    agg_n   <- stats::aggregate(trait ~ species,
                                  data = all_cells, FUN = length)
    agg <- merge(merge(agg_var, agg_ent, by = "species", all = TRUE),
                  agg_n, by = "species", all = TRUE)
    names(agg)[names(agg) == "trait"] <- "n_traits_missing"
    # Convert NaN/NA to 0 for sums
    agg$delta_var_total[is.na(agg$delta_var_total) |
                        is.nan(agg$delta_var_total)]     <- 0
    agg$delta_entropy_total[is.na(agg$delta_entropy_total) |
                            is.nan(agg$delta_entropy_total)] <- 0
    # For ranking: use delta_var_total when present, else delta_entropy_total
    rank_score <- ifelse(agg$delta_var_total > 0,
                         agg$delta_var_total,
                         agg$delta_entropy_total)
    agg <- agg[order(-rank_score), , drop = FALSE]
    rownames(agg) <- NULL
    utils::head(agg, top_n)
  }

  structure(out, class = c("pigauto_active", "data.frame"), by = by)
}

#' @export
print.pigauto_active <- function(x, ...) {
  by <- attr(x, "by")
  if (is.null(by)) by <- "cell"
  cat("Active-imputation suggestions (",
      sprintf("by = %s, n = %d", by, nrow(x)), ")\n", sep = "")
  if (nrow(x) == 0L) {
    cat("  No candidates -- dataset is fully observed on continuous-family traits.\n")
    return(invisible(x))
  }
  attrs <- attributes(x)
  attrs$class <- "data.frame"
  attrs$by <- NULL
  attributes(x) <- attrs
  print(x, ...)
  cat("\nObserve the cell(s) with the largest delta_var_total to maximally ",
      "reduce total imputation variance.\n", sep = "")
  invisible(x)
}
