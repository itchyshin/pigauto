#' Fit a tabular foundation model baseline (no phylogeny)
#'
#' Uses a pretrained tabular foundation model (TabImpute; Feitelberg
#' et al. 2025) to impute missing trait values **without** phylogenetic
#' information.
#'
#' This provides a non-phylogenetic baseline for comparison with
#' \code{\link{fit_baseline}} (BM/OU) and \code{\link{fit_pigauto}}
#' (graph autoencoder).  The output has the same structure as
#' \code{\link{fit_baseline}}, so it can be passed directly to
#' \code{\link{fit_pigauto}} via its \code{baseline} argument or
#' evaluated with \code{\link{evaluate_imputation}}.
#'
#' @details
#' The function calls the Python \pkg{tabimpute} package via
#' \pkg{reticulate}.
#'
#' Run \code{\link{setup_tabpfn}} once to install the required Python
#' packages into a virtual environment.
#'
#' TabImpute has quadratic attention cost and is practical for datasets
#' up to a few hundred species.  For larger datasets consider chunking
#' or using the column-wise TabPFN approach (not yet implemented).
#'
#' Because TabImpute is deterministic, \code{se} is estimated by
#' column-wise residual standard deviation on training cells (a rough
#' surrogate).  Set \code{se_method = "none"} to return \code{NA}
#' instead.
#'
#' @param data Object of class \code{"pigauto_data"} (output of
#'   \code{\link{preprocess_traits}}).
#' @param splits List (output of \code{\link{make_missing_splits}}) or
#'   \code{NULL}.  When provided, val and test cells are masked before
#'   imputation.
#' @param device Character: \code{"cpu"} (default) or \code{"cuda"}.
#' @param se_method Character.  How to approximate standard errors:
#'   \code{"residual"} (default) computes per-trait residual SD on
#'   training cells; \code{"none"} returns \code{NA}.
#' @param envname Python virtualenv name (default \code{"r-tabpfn"}).
#'   Set to \code{NULL} to use the current Python environment.
#' @return A list with:
#'   \describe{
#'     \item{mu}{Numeric matrix (n_species x n_traits), imputed values
#'       in z-score scale.}
#'     \item{se}{Numeric matrix (n_species x n_traits), approximate
#'       standard errors or \code{NA}.}
#'   }
#' @seealso \code{\link{setup_tabpfn}}, \code{\link{fit_baseline}},
#'   \code{\link{evaluate_imputation}}
#' @examples
#' \dontrun{
#' setup_tabpfn()                        # once
#' data(avonet300, tree300, package = "pigauto")
#' traits <- avonet300
#' rownames(traits) <- traits$Species_Key
#' traits$Species_Key <- NULL
#' pd     <- preprocess_traits(traits, tree300)
#' splits <- make_missing_splits(pd$X_scaled)
#'
#' # --- TabPFN alone (no tree) ---
#' bl_tab <- fit_baseline_tabpfn(pd, splits)
#' evaluate_imputation(bl_tab$mu, pd$X_scaled, splits)
#'
#' # --- pigauto with TabPFN as baseline (tree + FM) ---
#' graph <- build_phylo_graph(tree300)
#' fit   <- fit_pigauto(pd, tree300, splits,
#'                      graph = graph, baseline = bl_tab)
#' }
#' @export
fit_baseline_tabpfn <- function(data,
                                splits    = NULL,
                                device    = "cpu",
                                se_method = c("residual", "none"),
                                envname   = "r-tabpfn") {

  if (!inherits(data, "pigauto_data")) {
    stop("'data' must be a pigauto_data object (output of preprocess_traits).")
  }
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for TabPFN imputation.\n",
         "Install with: install.packages('reticulate')")
  }
  se_method <- match.arg(se_method)

  # Activate virtualenv (silently skip if NULL or unavailable)
  if (!is.null(envname)) {
    reticulate::use_virtualenv(envname, required = FALSE)
  }

  X <- data$X_scaled

  # Track which cells are observed training data (before masking)
  obs_train <- !is.na(X)

  # Mask val + test cells
  if (!is.null(splits)) {
    X[splits$val_idx]  <- NA
    X[splits$test_idx] <- NA
    obs_train[splits$val_idx]  <- FALSE
    obs_train[splits$test_idx] <- FALSE
  }

  # Convert NA -> NaN for numpy compatibility
  X[is.na(X)] <- NaN

  # ------ Python call -------------------------------------------------------
  tabimpute <- tryCatch(
    reticulate::import("tabimpute.interface"),
    error = function(e) {
      stop("Cannot import 'tabimpute'. ",
           "Run pigauto::setup_tabpfn() to install Python dependencies.\n",
           "Original error: ", conditionMessage(e))
    }
  )

  imputer    <- tabimpute$ImputePFN(device = device)
  X_py       <- reticulate::r_to_py(X)
  result_py  <- imputer$impute(X_py)
  mu         <- as.matrix(reticulate::py_to_r(result_py))

  # ------ Dimension & name checks -------------------------------------------
  if (!identical(dim(mu), dim(data$X_scaled))) {
    stop("TabImpute returned matrix of unexpected dimensions: ",
         paste(dim(mu), collapse = " x "), " (expected ",
         paste(dim(data$X_scaled), collapse = " x "), ").")
  }
  dimnames(mu) <- list(data$species_names, data$trait_names)

  # ------ Approximate SE ----------------------------------------------------
  se <- matrix(NA_real_, nrow(mu), ncol(mu), dimnames = dimnames(mu))

  if (se_method == "residual") {
    # Per-trait residual SD on observed training cells
    for (j in seq_len(ncol(mu))) {
      keep <- obs_train[, j] & is.finite(data$X_scaled[, j]) & is.finite(mu[, j])
      if (sum(keep) > 1L) {
        se[, j] <- stats::sd(data$X_scaled[keep, j] - mu[keep, j])
      }
    }
  }

  list(mu = mu, se = se)
}
