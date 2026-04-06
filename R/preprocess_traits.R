#' Preprocess trait data: align to tree, optionally log-transform, z-score
#'
#' Aligns species in the trait data frame to the tree, drops species missing
#' from the tree (with a warning), optionally log-transforms and z-scores each
#' trait column. The output rows are reordered to match \code{tree$tip.label}
#' order, which is required for downstream graph operations.
#'
#' @param traits \code{data.frame} with species as row names and numeric traits
#'   as columns (output of \code{\link{read_traits}}).
#' @param tree object of class \code{"phylo"}.
#' @param log_transform logical. Apply \code{log()} before z-scoring? Default
#'   \code{TRUE}; appropriate for most morphological measurements.
#' @param center logical. Subtract column means? Default \code{TRUE}.
#' @param scale logical. Divide by column standard deviations? Default
#'   \code{TRUE}.
#' @return A list of class \code{"pigauto_data"} with components:
#'   \describe{
#'     \item{X_scaled}{Numeric matrix (n_species x n_traits), z-scored.}
#'     \item{X_raw}{Numeric matrix (n_species x n_traits), after optional log
#'       but before z-scoring.}
#'     \item{means}{Named numeric vector of column means used for
#'       standardisation (on the log scale if \code{log_transform = TRUE}).}
#'     \item{sds}{Named numeric vector of column SDs.}
#'     \item{species_names}{Character vector matching \code{tree$tip.label}
#'       order.}
#'     \item{trait_names}{Character vector of trait column names.}
#'     \item{log_transform}{Logical, whether log was applied.}
#'   }
#' @examples
#' data(avonet300, tree300, package = "pigauto")
#' traits <- avonet300
#' rownames(traits) <- traits$Species_Key
#' traits$Species_Key <- NULL
#' pd <- preprocess_traits(traits, tree300)
#' dim(pd$X_scaled)   # 300 x 4
#' @importFrom ape keep.tip
#' @export
preprocess_traits <- function(traits, tree, log_transform = TRUE,
                              center = TRUE, scale = TRUE) {
  if (!is.data.frame(traits)) stop("'traits' must be a data.frame.")
  if (!inherits(tree, "phylo")) stop("'tree' must be a phylo object.")
  if (is.null(rownames(traits))) stop("'traits' must have species as row names.")

  # Align species
  in_tree  <- rownames(traits) %in% tree$tip.label
  in_data  <- tree$tip.label %in% rownames(traits)

  n_dropped <- sum(!in_tree)
  if (n_dropped > 0) {
    warning(n_dropped, " species in traits not found in tree -- dropped.")
    traits <- traits[in_tree, , drop = FALSE]
  }
  n_missing_from_data <- sum(!in_data)
  if (n_missing_from_data > 0) {
    message(n_missing_from_data,
            " tree tip(s) have no trait data and will have all-NA rows.")
  }

  # Reorder to tree$tip.label order (critical for graph operations)
  idx <- match(tree$tip.label, rownames(traits))
  traits <- traits[idx, , drop = FALSE]
  rownames(traits) <- tree$tip.label

  X <- as.matrix(traits)
  storage.mode(X) <- "double"

  # Log transform
  if (log_transform) {
    if (any(X <= 0, na.rm = TRUE)) {
      stop("log_transform = TRUE requires all trait values to be positive. ",
           "Found non-positive value(s). Set log_transform = FALSE or ",
           "remove/transform the offending values first.")
    }
    X <- log(X)
  }

  X_raw <- X

  # Z-score (compute on observed values only)
  col_means <- colMeans(X, na.rm = TRUE)
  col_sds   <- apply(X, 2, stats::sd, na.rm = TRUE)
  col_sds[col_sds == 0] <- 1  # avoid division by zero for constant columns

  if (center) X <- sweep(X, 2, col_means, "-")
  if (scale)  X <- sweep(X, 2, col_sds,   "/")

  structure(
    list(
      X_scaled      = X,
      X_raw         = X_raw,
      means         = col_means,
      sds           = col_sds,
      species_names = tree$tip.label,
      trait_names   = colnames(X),
      log_transform = log_transform
    ),
    class = "pigauto_data"
  )
}

#' @export
print.pigauto_data <- function(x, ...) {
  cat("pigauto_data\n")
  cat("  Species:", length(x$species_names), "\n")
  cat("  Traits:", length(x$trait_names), ":",
      paste(x$trait_names, collapse = ", "), "\n")
  cat("  Log-transformed:", x$log_transform, "\n")
  missing_frac <- mean(is.na(x$X_scaled))
  if (missing_frac > 0) {
    cat("  Missing values:", round(100 * missing_frac, 1), "%\n")
  }
  invisible(x)
}
