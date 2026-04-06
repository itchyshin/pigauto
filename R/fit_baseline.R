#' Fit a phylogenetic BM baseline using Rphylopars
#'
#' Fits a Brownian Motion (or OU) model with \pkg{Rphylopars} and returns
#' imputed means and standard errors for every species, including those with
#' missing trait values.
#'
#' @details
#' When \code{splits} is supplied the val and test cells are masked to
#' \code{NA} before fitting, so the baseline is evaluated under the same
#' conditions as \code{\link{fit_pigauto}}.
#'
#' If \pkg{Rphylopars} throws an error (typically a singular covariance
#' matrix from near-zero branch lengths), a small jitter of \code{1e-6} is
#' added to all branch lengths and the fit is retried once.
#'
#' @param data object of class \code{"pigauto_data"} (output of
#'   \code{\link{preprocess_traits}}).
#' @param tree object of class \code{"phylo"}.
#' @param splits list (output of \code{\link{make_missing_splits}}) or
#'   \code{NULL}. When provided, val and test cells are masked before
#'   fitting.
#' @param model character. Evolutionary model: \code{"BM"} (default) or
#'   \code{"OU"}.
#' @return A list with:
#'   \describe{
#'     \item{mu}{Numeric matrix (n_species x n_traits), imputed means in
#'       z-score scale.}
#'     \item{se}{Numeric matrix (n_species x n_traits), standard errors
#'       (square root of ancestral variance).}
#'   }
#' @examples
#' \dontrun{
#' data(avonet300, tree300, package = "pigauto")
#' traits <- avonet300; rownames(traits) <- traits$Species_Key
#' traits$Species_Key <- NULL
#' pd     <- preprocess_traits(traits, tree300)
#' splits <- make_missing_splits(pd$X_scaled)
#' bl     <- fit_baseline(pd, tree300, splits)
#' dim(bl$mu)   # 300 x 4
#' }
#' @export
fit_baseline <- function(data, tree, splits = NULL, model = "BM") {
  if (!inherits(data, "pigauto_data")) {
    stop("'data' must be a pigauto_data object (output of preprocess_traits).")
  }
  if (!inherits(tree, "phylo")) stop("'tree' must be a phylo object.")

  X <- data$X_scaled

  # Mask val + test cells before fitting
  if (!is.null(splits)) {
    X[splits$val_idx]  <- NA
    X[splits$test_idx] <- NA
  }

  n   <- nrow(X)
  p   <- ncol(X)
  spp <- data$species_names

  # Build Rphylopars input data.frame
  df_ph <- as.data.frame(X)
  df_ph <- cbind(species = spp, df_ph)

  fit_rphylo <- function(tr) {
    suppressWarnings(
      Rphylopars::phylopars(
        trait_data  = df_ph,
        tree        = tr,
        model       = model,
        pheno_error = FALSE
      )
    )
  }

  fit <- tryCatch(
    fit_rphylo(tree),
    error = function(e) {
      message("Rphylopars failed (", conditionMessage(e),
              "). Retrying with branch-length jitter.")
      tree_j <- tree
      tree_j$edge.length <- tree_j$edge.length + 1e-6
      tryCatch(
        fit_rphylo(tree_j),
        error = function(e2) {
          stop("Rphylopars failed after jitter retry: ", conditionMessage(e2))
        }
      )
    }
  )

  # Extract tip reconstructions (rows 1:n correspond to tips)
  mu_raw  <- as.matrix(fit$anc_recon[spp, data$trait_names, drop = FALSE])
  var_raw <- as.matrix(fit$anc_var[spp,  data$trait_names, drop = FALSE])

  se_raw  <- sqrt(pmax(var_raw, 0))  # guard tiny negative variances
  rownames(mu_raw) <- spp
  rownames(se_raw) <- spp

  list(mu = mu_raw, se = se_raw)
}
