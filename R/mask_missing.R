#' Create an observed/missing mask matrix
#'
#' Returns a logical matrix of the same dimensions as \code{X}, with
#' \code{TRUE} where values are observed (not \code{NA}).
#'
#' @param X numeric matrix (species x traits).
#' @return Logical matrix, \code{TRUE} = observed.
#' @examples
#' X <- matrix(c(1, NA, 3, 4), nrow = 2)
#' mask_missing(X)
#' @export
mask_missing <- function(X) {
  !is.na(X)
}

#' Split all matrix cells into train/val/test for imputation evaluation
#'
#' Randomly designates a fraction of cells as "missing" and splits them into
#' validation and test sets. The remaining cells are treated as observed.
#' Masking is applied to \emph{all} cells uniformly (not conditioned on whether
#' values are already \code{NA} in \code{X}).
#'
#' The returned index vectors use linear (column-major) indexing into the
#' flattened matrix, consistent with \code{which()}, \code{[]} single-index
#' subsetting, and \code{torch_tensor} flattening order.
#'
#' @param X numeric matrix (species x traits). Used only for dimensions.
#' @param missing_frac numeric. Fraction of all cells to designate as missing
#'   (default \code{0.25}).
#' @param val_frac numeric. Fraction of missing cells to use for validation
#'   (default \code{0.25}); the rest become the test set.
#' @param seed integer. Random seed for reproducibility (default \code{555}).
#' @return A list with:
#'   \describe{
#'     \item{val_idx}{Integer vector of linear cell indices for validation.}
#'     \item{test_idx}{Integer vector of linear cell indices for test.}
#'     \item{n}{Integer. Number of species (rows).}
#'     \item{p}{Integer. Number of traits (columns).}
#'     \item{mask}{Logical matrix (n x p). \code{TRUE} = observed (not in
#'       val or test).}
#'   }
#' @examples
#' X <- matrix(rnorm(100), nrow = 20)
#' splits <- make_missing_splits(X, missing_frac = 0.25, seed = 1)
#' length(splits$val_idx)   # ~25% of 20 cells designated as val
#' @export
make_missing_splits <- function(X, missing_frac = 0.25, val_frac = 0.25,
                                seed = 555) {
  if (!is.matrix(X)) stop("'X' must be a matrix.")
  if (missing_frac <= 0 || missing_frac >= 1)
    stop("'missing_frac' must be in (0, 1).")
  if (val_frac <= 0 || val_frac >= 1)
    stop("'val_frac' must be in (0, 1).")

  n <- nrow(X)
  p <- ncol(X)
  all_idx <- seq_len(n * p)

  set.seed(seed)
  n_miss  <- floor(missing_frac * length(all_idx))
  miss    <- sample(all_idx, n_miss)

  n_val   <- floor(val_frac * length(miss))
  val_idx  <- miss[seq_len(n_val)]
  test_idx <- miss[(n_val + 1L):length(miss)]

  mask <- matrix(TRUE, nrow = n, ncol = p)
  mask[c(val_idx, test_idx)] <- FALSE

  list(
    val_idx  = val_idx,
    test_idx = test_idx,
    n        = n,
    p        = p,
    mask     = mask
  )
}
