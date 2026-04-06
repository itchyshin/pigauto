#' Create an observed/missing mask matrix
#'
#' Returns a logical matrix of the same dimensions as \code{X}, with
#' \code{TRUE} where values are observed (not \code{NA}).
#'
#' @param X numeric matrix (species x traits or species x latent columns).
#' @return Logical matrix, \code{TRUE} = observed.
#' @examples
#' X <- matrix(c(1, NA, 3, 4), nrow = 2)
#' mask_missing(X)
#' @export
mask_missing <- function(X) {
  !is.na(X)
}

#' Split cells into train/val/test for imputation evaluation
#'
#' Randomly designates a fraction of cells as "missing" and splits them into
#' validation and test sets.  When a \code{trait_map} is supplied, masking
#' operates at the **original trait level** -- all latent columns belonging
#' to one trait are held out together (important for categorical traits).
#'
#' The returned index vectors use linear (column-major) indexing.  Both
#' original-trait-space and latent-space indices are returned when a
#' \code{trait_map} is present.
#'
#' @param X numeric matrix (species x latent columns from
#'   \code{\link{preprocess_traits}}).  Used only for dimensions.
#' @param missing_frac numeric. Fraction of all (species, trait) cells to
#'   designate as missing (default \code{0.25}).
#' @param val_frac numeric. Fraction of missing cells for validation
#'   (default \code{0.25}); the rest become the test set.
#' @param seed integer. Random seed for reproducibility (default \code{555}).
#' @param trait_map list of trait descriptors (from \code{pigauto_data}).
#'   If \code{NULL}, masking is applied per latent column (v0.1 behaviour).
#' @return A list with:
#'   \describe{
#'     \item{val_idx}{Integer vector of linear indices (latent space).}
#'     \item{test_idx}{Integer vector of linear indices (latent space).}
#'     \item{val_idx_trait}{Integer vector in original-trait space (if
#'       \code{trait_map} supplied).}
#'     \item{test_idx_trait}{Integer vector in original-trait space (if
#'       \code{trait_map} supplied).}
#'     \item{n}{Number of species (rows).}
#'     \item{p}{Number of latent columns.}
#'     \item{n_traits}{Number of original traits.}
#'     \item{mask}{Logical matrix (n x p_latent). \code{TRUE} = observed.}
#'   }
#' @examples
#' X <- matrix(rnorm(100), nrow = 20)
#' splits <- make_missing_splits(X, missing_frac = 0.25, seed = 1)
#' length(splits$val_idx)
#' @export
make_missing_splits <- function(X, missing_frac = 0.25, val_frac = 0.25,
                                seed = 555, trait_map = NULL) {
  if (!is.matrix(X)) stop("'X' must be a matrix.")
  if (missing_frac <= 0 || missing_frac >= 1)
    stop("'missing_frac' must be in (0, 1).")
  if (val_frac <= 0 || val_frac >= 1)
    stop("'val_frac' must be in (0, 1).")

  n <- nrow(X)
  p <- ncol(X)

  if (is.null(trait_map)) {
    # ---- v0.1 behaviour: mask per latent column ----------------------------
    all_idx <- seq_len(n * p)
    set.seed(seed)
    n_miss   <- floor(missing_frac * length(all_idx))
    miss     <- sample(all_idx, n_miss)
    n_val    <- floor(val_frac * length(miss))
    val_idx  <- miss[seq_len(n_val)]
    test_idx <- miss[(n_val + 1L):length(miss)]

    mask <- matrix(TRUE, nrow = n, ncol = p)
    mask[c(val_idx, test_idx)] <- FALSE

    return(list(
      val_idx  = val_idx,
      test_idx = test_idx,
      n        = n,
      p        = p,
      n_traits = p,
      mask     = mask
    ))
  }

  # ---- Trait-level masking (mixed types) -----------------------------------
  n_traits <- length(trait_map)
  n_cells  <- n * n_traits  # one cell per (species, original trait)

  set.seed(seed)
  n_miss   <- floor(missing_frac * n_cells)
  miss     <- sample(seq_len(n_cells), n_miss)
  n_val    <- floor(val_frac * length(miss))
  val_trait   <- miss[seq_len(n_val)]
  test_trait  <- miss[(n_val + 1L):length(miss)]

  # Expand original-trait indices to latent-column indices
  val_latent  <- expand_trait_idx_to_latent(val_trait, n, trait_map)
  test_latent <- expand_trait_idx_to_latent(test_trait, n, trait_map)

  mask <- matrix(TRUE, nrow = n, ncol = p)
  mask[c(val_latent, test_latent)] <- FALSE

  list(
    val_idx       = val_latent,
    test_idx      = test_latent,
    val_idx_trait  = val_trait,
    test_idx_trait = test_trait,
    n             = n,
    p             = p,
    n_traits      = n_traits,
    mask          = mask
  )
}


#' Expand original-trait-space linear indices to latent-space indices
#'
#' @param idx integer vector of linear indices in (n x n_traits) space.
#' @param n integer, number of species.
#' @param trait_map list of trait descriptors.
#' @return integer vector of linear indices in (n x p_latent) space.
#' @keywords internal
expand_trait_idx_to_latent <- function(idx, n, trait_map) {
  # Convert linear index to (row, trait_number)
  row_i   <- ((idx - 1L) %% n) + 1L
  trait_j <- ((idx - 1L) %/% n) + 1L

  latent_idx <- integer(0)
  for (i in seq_along(idx)) {
    tm <- trait_map[[trait_j[i]]]
    for (lc in tm$latent_cols) {
      latent_idx <- c(latent_idx, row_i[i] + (lc - 1L) * n)
    }
  }
  latent_idx
}
