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
