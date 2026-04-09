#' Generate M complete datasets for multiple imputation
#'
#' Run pigauto's full imputation pipeline and return `M` stochastic
#' completions of the trait matrix instead of a single point estimate.
#' The `M` datasets are the input needed for the classical multiple
#' imputation workflow: fit a downstream model on each dataset, then
#' pool the results with Rubin's rules via [pool_mi()]. This is the
#' standard way to propagate imputation uncertainty into phylogenetic
#' comparative analyses (PGLS, PGLMM, etc.) rather than treating
#' imputed cells as if they were observed.
#'
#' @param traits data.frame with species as rownames and trait columns.
#'   Same input format as [impute()]. Supported column types are
#'   numeric, integer, factor, ordered factor, and logical.
#' @param tree object of class `phylo` aligned with `traits`.
#' @param m integer. Number of imputation datasets to generate
#'   (default `100`). Internally draws are produced by Monte Carlo
#'   dropout of the GNN, so each dataset differs only at cells that
#'   were missing in the input; observed cells are identical across
#'   all `M` datasets.
#' @param species_col character or `NULL`. If set, marks the column
#'   in `traits` containing species identifiers and enables multiple
#'   observations per species. See [impute()] for details.
#' @param log_transform logical. Auto-log positive continuous columns
#'   (default `TRUE`).
#' @param missing_frac numeric. Fraction of observed cells held out for
#'   validation/test during training (default `0.25`). Passed through
#'   to [impute()].
#' @param epochs integer. Maximum GNN training epochs (default `2000`).
#' @param verbose logical. Print progress (default `TRUE`).
#' @param seed integer. Random seed (default `1`).
#' @param ... additional arguments forwarded to [fit_pigauto()] via
#'   [impute()].
#'
#' @return An object of class `"pigauto_mi"` with components:
#'   \describe{
#'     \item{`datasets`}{A list of length `m`. Each element is a
#'       data.frame with the same shape and column types as the input
#'       `traits`; observed cells are preserved and missing cells are
#'       filled with the corresponding imputation draw. Pass this list
#'       to [with_imputations()] to fit downstream models.}
#'     \item{`m`}{Number of imputations.}
#'     \item{`pooled_point`}{A single data.frame whose missing cells
#'       are replaced by the MC-averaged point estimate. Convenient for
#'       reporting but does *not* propagate imputation uncertainty —
#'       use `datasets` + [pool_mi()] for inference.}
#'     \item{`se`}{Matrix of per-cell standard errors combining the
#'       baseline SE and the between-imputation standard deviation.}
#'     \item{`imputed_mask`}{Logical matrix; `TRUE` where a cell was
#'       originally missing.}
#'     \item{`fit`}{The underlying [`pigauto_fit`][fit_pigauto()]
#'       object, retained for diagnostics and for calls to [predict()]
#'       on new data.}
#'     \item{`data`}{The [`pigauto_data`][preprocess_traits()] object.}
#'     \item{`tree`}{The input phylogeny.}
#'     \item{`species_col`}{Passed-through species-column name or
#'       `NULL`.}
#'   }
#'
#' @details
#' Multiple imputation is a method for doing *downstream analysis*
#' under missing data, not an end in itself. Plugging a single
#' point-estimate imputation into a regression underestimates standard
#' errors because it treats imputed cells as if they were observed.
#' The standard remedy, due to Rubin (1987), is to generate `M`
#' stochastic completions, fit the downstream model on each, and pool
#' the results. `multi_impute()` + [with_imputations()] + [pool_mi()]
#' implement this workflow end to end.
#'
#' Nakagawa & Freckleton (2008, 2011) review the consequences of
#' ignoring missing data in ecological and comparative analyses and
#' argue for multiple imputation as the default.
#'
#' @references
#' Rubin DB (1987). *Multiple Imputation for Nonresponse in Surveys.*
#' Wiley.
#'
#' Nakagawa S, Freckleton RP (2008). "Missing inaction: the dangers of
#' ignoring missing data." *Trends in Ecology & Evolution* 23(11):
#' 592–596.
#'
#' Nakagawa S, Freckleton RP (2011). "Model averaging, missing data and
#' multiple imputation: a case study for behavioural ecology."
#' *Behavioral Ecology and Sociobiology* 65(1): 103–116.
#'
#' @seealso [impute()] for single-point imputation, [with_imputations()]
#'   for applying a model-fitting function across the `M` datasets,
#'   [pool_mi()] for Rubin's rules pooling of the resulting fits.
#'
#' @examples
#' \dontrun{
#' library(pigauto)
#' data(avonet300, tree300)
#' df <- avonet300; rownames(df) <- df$Species_Key; df$Species_Key <- NULL
#'
#' # Generate 100 complete datasets
#' mi <- multi_impute(df, tree300, m = 100)
#' print(mi)
#'
#' # Downstream analysis: phylogenetic GLS via nlme, pooled with Rubin's rules
#' fits <- with_imputations(mi, function(d) {
#'   d$species <- rownames(d)
#'   nlme::gls(
#'     log(Mass) ~ log(Wing.Length),
#'     correlation = ape::corBrownian(phy = tree300, form = ~species),
#'     data = d, method = "ML"
#'   )
#' })
#' pool_mi(fits)
#' }
#'
#' @export
multi_impute <- function(traits, tree, m = 100L,
                         species_col = NULL,
                         log_transform = TRUE,
                         missing_frac = 0.25,
                         epochs = 2000L, verbose = TRUE, seed = 1L, ...) {

  m <- as.integer(m)
  if (!is.finite(m) || m < 2L) {
    stop("`m` must be an integer >= 2 (multiple imputation needs at least ",
         "two draws). Got m = ", m, ".", call. = FALSE)
  }

  # 1. Run the full pipeline with MC dropout sampling.
  res <- impute(
    traits        = traits,
    tree          = tree,
    species_col   = species_col,
    log_transform = log_transform,
    missing_frac  = missing_frac,
    n_imputations = m,
    epochs        = as.integer(epochs),
    verbose       = verbose,
    seed          = as.integer(seed),
    ...
  )

  pred <- res$prediction
  if (is.null(pred$imputed_datasets) || length(pred$imputed_datasets) != m) {
    stop("predict.pigauto_fit() did not return ", m,
         " imputed datasets. This is an internal error — please report.",
         call. = FALSE)
  }

  # 2. Merge each raw imputation with observed cells via build_completed().
  #    build_completed() lives in R/impute.R and is unexported; we reach
  #    it directly because we are in the same package.
  datasets <- lapply(pred$imputed_datasets, function(imp_df) {
    build_completed(traits, imp_df, species_col)$completed
  })

  structure(
    list(
      datasets     = datasets,
      m            = m,
      pooled_point = res$completed,
      se           = pred$se,
      imputed_mask = res$imputed_mask,
      fit          = res$fit,
      data         = res$data,
      tree         = tree,
      species_col  = species_col,
      evaluation   = res$evaluation
    ),
    class = c("pigauto_mi", "list")
  )
}


#' @export
print.pigauto_mi <- function(x, ...) {
  n_sp  <- length(x$data$species_names)
  traits <- vapply(x$data$trait_map, "[[", character(1), "name")
  p <- length(traits)

  total_cells <- length(x$imputed_mask)
  n_imp_cells <- sum(x$imputed_mask)
  pct <- if (total_cells > 0) 100 * n_imp_cells / total_cells else 0

  cat("pigauto multiple imputation\n")
  cat(sprintf("  M        : %d imputations (MC dropout)\n", x$m))
  cat(sprintf("  Species  : %d\n", n_sp))
  cat(sprintf("  Traits   : %d -- %s\n", p,
              paste(traits, collapse = ", ")))
  cat(sprintf("  Cells    : %d imputed / %d total (%.1f%%)\n",
              n_imp_cells, total_cells, pct))

  cat("\n  Access imputation draws:  mi$datasets[[i]]\n")
  cat("  Fit downstream models:    with_imputations(mi, fit_fun)\n")
  cat("  Pool with Rubin's rules:  pool_mi(fits)\n")
  invisible(x)
}
