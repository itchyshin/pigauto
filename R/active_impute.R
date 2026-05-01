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
#' The variance-reduction formula is derived from a Sherman-Morrison
#' rank-1 inverse update on the BM conditional MVN: adding species
#' \eqn{s} to the observed set updates the inverse correlation matrix
#' by a known closed form.  For each candidate cell \eqn{(s, t)},
#' \deqn{
#'   \Delta V(s, t) = \sigma_t^2
#'     \sum_{i \in \mathrm{miss}_t} \frac{D_{ik}^2}{\alpha_k}
#' }
#' where \eqn{D = R_{mm} - R_{mo} R_{oo}^{-1} R_{om}} is the residual
#' matrix at currently-missing cells, \eqn{\alpha_k = D_{kk}} is the
#' current relative leverage of cell \eqn{k}, and \eqn{\sigma_t^2} is
#' the REML BM variance for trait \eqn{t}.  See \code{R/active_impute.R}
#' source for the full derivation.
#'
#' Reductions are summed across the BM-eligible traits for each species
#' when \code{by = "species"}, supporting the typical use case where
#' measuring a species observes all of its currently-missing
#' continuous-family traits at once.
#'
#' Discrete traits (binary, categorical, zi_count) are silently
#' skipped: their uncertainty is captured by class-probability entropy,
#' not by the BM variance, and the variance-reduction formula above
#' does not apply.  An expected-entropy-reduction extension is queued
#' for v2.
#'
#' @param result A \code{pigauto_result} object returned by
#'   \code{\link{impute}}.  Must have been produced from single-obs
#'   data; multi-obs inputs error with a clear message.
#' @param top_n integer, default \code{10L}.  Number of suggestions to
#'   return (descending by \code{delta_var_total}).
#' @param by character, one of \code{"cell"} (default) or
#'   \code{"species"}.  \code{"cell"} returns individual
#'   \code{(species, trait)} pairs.  \code{"species"} aggregates by
#'   species (summing variance reductions across the species'
#'   currently-missing continuous-family traits).
#' @param types character vector of pigauto trait types to include.
#'   Default includes all continuous-family types
#'   (\code{continuous}, \code{count}, \code{ordinal}, \code{proportion}).
#' @return A data.frame of class \code{"pigauto_active"}.  Columns
#'   when \code{by = "cell"}: \code{species}, \code{trait}, \code{type},
#'   \code{delta_var_total} (descending).  When \code{by = "species"}:
#'   \code{species}, \code{delta_var_total}, \code{n_traits_missing}.
#' @examples
#' \dontrun{
#' data(avonet300, tree300, package = "pigauto")
#' res <- impute(avonet300, tree300)
#' suggest_next_observation(res, top_n = 5)              # top-5 cells
#' suggest_next_observation(res, top_n = 10, by = "species")  # top-10 species
#' }
#' @seealso \code{\link{impute}}
#' @export
suggest_next_observation <- function(result, top_n = 10L,
                                       by = c("cell", "species"),
                                       types = c("continuous", "count",
                                                 "ordinal", "proportion")) {
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

  # Recover R_phy: prefer fit$graph$R_phy (cached), fall back to
  # phylo_cor_matrix(tree) if absent (older fits, defensive)
  R_phy <- fit$graph$R_phy
  if (is.null(R_phy) && !is.null(result$tree)) {
    R_phy <- phylo_cor_matrix(result$tree)
  }
  if (is.null(R_phy)) {
    stop("suggest_next_observation: phylogenetic correlation matrix not ",
         "available in fit$graph$R_phy. Refit with build_phylo_graph() ",
         "to populate it.", call. = FALSE)
  }
  R_phy <- R_phy[spp, spp, drop = FALSE]

  # Per-trait loop: variance reduction on continuous-family columns
  cell_dfs <- list()
  for (tm in trait_map) {
    if (!(tm$type %in% types)) next
    j <- tm$latent_cols[1L]
    y <- X[, j]
    names(y) <- spp
    if (sum(!is.na(y)) < 5L || sum(is.na(y)) == 0L) next

    delta_v <- bm_variance_reduction(y, R_phy)
    if (length(delta_v) == 0L) next

    miss_idx <- which(is.na(y))
    cell_dfs[[length(cell_dfs) + 1L]] <- data.frame(
      species         = spp[miss_idx],
      trait           = tm$name,
      type            = tm$type,
      delta_var_total = as.numeric(delta_v),
      stringsAsFactors = FALSE
    )
  }

  if (length(cell_dfs) == 0L) {
    out <- data.frame(species = character(0L),
                      trait = character(0L),
                      type = character(0L),
                      delta_var_total = numeric(0L),
                      stringsAsFactors = FALSE)
    return(structure(out, class = c("pigauto_active", "data.frame"),
                     by = by))
  }

  all_cells <- do.call(rbind, cell_dfs)

  out <- if (identical(by, "cell")) {
    all_cells <- all_cells[order(-all_cells$delta_var_total), , drop = FALSE]
    rownames(all_cells) <- NULL
    utils::head(all_cells, top_n)
  } else {
    agg_total <- stats::aggregate(delta_var_total ~ species,
                                    data = all_cells, FUN = sum)
    agg_n     <- stats::aggregate(trait ~ species,
                                    data = all_cells, FUN = length)
    agg <- merge(agg_total, agg_n, by = "species")
    names(agg)[names(agg) == "trait"] <- "n_traits_missing"
    agg <- agg[order(-agg$delta_var_total), , drop = FALSE]
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
