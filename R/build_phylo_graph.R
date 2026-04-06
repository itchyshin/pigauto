# Internal helpers -------------------------------------------------------

# Correct spectral features: smallest non-zero Laplacian eigenvectors.
# Returns plain R matrix (n x k); tensors are created in fit_pigauto().
# Source: script/test_new6.R lines 34-56 (fixed ordering).
get_spectral_features_internal <- function(tree, k, sigma_mult) {
  D     <- ape::cophenetic.phylo(tree)
  sigma <- stats::median(D) * sigma_mult
  A     <- exp(-(D^2) / (2 * sigma^2))
  diag(A) <- 0

  L   <- diag(rowSums(A)) - A
  eig <- eigen(L, symmetric = TRUE)

  # eigen() returns eigenvalues in DECREASING order for symmetric matrices;
  # re-sort ascending to get smallest first.
  ord  <- order(eig$values, decreasing = FALSE)
  vecs <- eig$vectors[, ord, drop = FALSE]

  n <- nrow(A)
  if (k + 1L > n) {
    stop("k_eigen (", k, ") is too large for a tree with ", n, " tips. ",
         "Use k_eigen < ", n - 1L, ".")
  }
  # Skip vector 1 (constant, eigenvalue ~0); take next k.
  vecs[, 2L:(k + 1L), drop = FALSE]
}

# Symmetric-normalised Gaussian kernel adjacency + self-loops.
# Source: script/test_new6.R lines 58-72.
get_adj_symnorm_internal <- function(tree, sigma_mult) {
  D     <- ape::cophenetic.phylo(tree)
  sigma <- stats::median(D) * sigma_mult
  A     <- exp(-(D^2) / (2 * sigma^2))
  diag(A) <- 1  # self-loops for stability

  rs        <- rowSums(A) + 1e-8
  Dinv_sqrt <- diag(1 / sqrt(rs))
  Dinv_sqrt %*% A %*% Dinv_sqrt
}


# Exported function -------------------------------------------------------

#' Build a phylogenetic graph representation from a tree
#'
#' Computes a symmetric-normalised Gaussian kernel adjacency matrix and
#' Laplacian spectral node features from the cophenetic distance matrix.
#' Results are optionally cached to an \code{.rds} file.
#'
#' @details
#' The graph is built in three steps: (1) cophenetic distances between all
#' pairs of tips; (2) Gaussian kernel \eqn{A_{ij} = \exp(-d_{ij}^2 /
#' (2\sigma^2))} with \eqn{\sigma = \text{median}(D) \times
#' \code{sigma_mult}}; (3) symmetric normalisation
#' \eqn{\tilde{A} = D^{-1/2} A D^{-1/2}} with self-loops added before
#' normalisation.
#'
#' Spectral node features are the \code{k_eigen} smallest non-zero
#' eigenvectors of the unnormalised Laplacian, which encode the broad
#' cluster structure of the phylogeny.
#'
#' @param tree object of class \code{"phylo"}.
#' @param k_eigen integer. Number of Laplacian eigenvectors to use as node
#'   features (default \code{8}).
#' @param sigma_mult numeric. Bandwidth multiplier:
#'   \eqn{\sigma = \mathrm{median}(D) \times \code{sigma_mult}}
#'   (default \code{0.5}).
#' @param cache_path character or \code{NULL}. Path to an \code{.rds} cache
#'   file. If the file exists and dimensions match, it is loaded instead of
#'   recomputing. \code{NULL} disables caching (default).
#' @return A list with:
#'   \describe{
#'     \item{adj}{Numeric matrix (n x n). Symmetric-normalised adjacency.}
#'     \item{coords}{Numeric matrix (n x k_eigen). Spectral node features.}
#'     \item{n}{Integer. Number of tips.}
#'     \item{sigma}{Numeric. Bandwidth used.}
#'   }
#' @examples
#' set.seed(1)
#' tree <- ape::rtree(30)
#' g <- build_phylo_graph(tree, k_eigen = 4)
#' dim(g$adj)    # 30 x 30
#' dim(g$coords) # 30 x 4
#' @importFrom ape cophenetic.phylo
#' @export
build_phylo_graph <- function(tree, k_eigen = 8L, sigma_mult = 0.5,
                               cache_path = NULL) {
  if (!inherits(tree, "phylo")) stop("'tree' must be a phylo object.")
  n <- length(tree$tip.label)

  if (n > 2000L) {
    warning(
      "build_phylo_graph: tree has ", n, " tips. ",
      "cophenetic() requires O(n^2) memory (~",
      round(n^2 * 8 / 1e6), " MB). ",
      "Consider subsetting your tree for large datasets."
    )
  }

  # Load from cache if available and dimensions match
  if (!is.null(cache_path) && file.exists(cache_path)) {
    cache <- readRDS(cache_path)
    if (isTRUE(nrow(cache$adj) == n) && isTRUE(ncol(cache$coords) == k_eigen)) {
      return(cache)
    }
    message("Cache dimensions mismatch -- recomputing.")
  }

  adj    <- get_adj_symnorm_internal(tree, sigma_mult)
  coords <- get_spectral_features_internal(tree, k_eigen, sigma_mult)
  sigma  <- stats::median(ape::cophenetic.phylo(tree)) * sigma_mult

  result <- list(adj = adj, coords = coords, n = n, sigma = sigma)

  if (!is.null(cache_path)) {
    saveRDS(result, cache_path)
  }

  result
}
