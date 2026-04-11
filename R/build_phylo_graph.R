# Internal helpers -------------------------------------------------------

# Smallest non-zero Laplacian eigenvectors from a precomputed distance matrix.
#
# Given a cophenetic distance matrix `D`, this builds the Gaussian-kernel
# graph Laplacian and returns the (k+1) smallest eigenvectors. The first
# eigenvector (~constant, eigenvalue ~0) is dropped and the remaining `k`
# are returned as an (n x k) matrix.
#
# Uses RSpectra::eigs_sym (sparse Lanczos) when the package is available and
# n > `dense_threshold`, which is O(n * k * iters) instead of O(n^3).
# Falls back to base::eigen() for small trees or when RSpectra is absent,
# so the package works without RSpectra installed.
#
# Why the threshold is 7500 and not smaller: on highly symmetric ultrametric
# trees (e.g. those produced by ape::rcoal) the Laplacian spectrum has tight
# degenerate clusters, and Lanczos needs many iterations to separate
# eigenvalues inside those clusters. Empirically, dense eigen() is faster
# than sparse eigs_sym() on rcoal trees up to ~7500 tips, and only clearly
# wins above that. On realistic asymmetric phylogenies (variable branch
# lengths, cherries of different ages) the crossover is much earlier -- the
# AVONET BirdTree Stage2 Hackett MCC tree (9,993 tips) solves in ~40s
# sparse vs >400s dense, a 10x speedup. The conservative threshold ensures
# the sparse path is only used where it reliably wins even on worst-case
# ultrametric trees, while still unlocking the scaling path for real 10k+
# phylogenies.
spectral_features_from_D <- function(D, k, sigma_mult,
                                     dense_threshold = 7500L) {
  sigma <- stats::median(D) * sigma_mult
  A     <- exp(-(D^2) / (2 * sigma^2))
  diag(A) <- 0

  n <- nrow(A)
  if (k + 1L > n) {
    stop("k_eigen (", k, ") is too large for a tree with ", n, " tips. ",
         "Use k_eigen < ", n - 1L, ".")
  }

  L <- diag(rowSums(A)) - A

  use_rspectra <- n > dense_threshold &&
    requireNamespace("RSpectra", quietly = TRUE)

  if (use_rspectra) {
    # Direct Lanczos via `which = "SA"` (smallest algebraic). This is
    # robust on positive-semidefinite matrices with a zero eigenvalue,
    # unlike shift-and-invert with sigma near zero (which makes the
    # shifted operator nearly singular) or `which = "SM"` (which uses
    # shift-and-invert internally and suffers from the same problem).
    #
    # On highly symmetric ultrametric trees (e.g. ape::rcoal) the
    # Laplacian spectrum contains tight degenerate clusters; Lanczos
    # needs a larger Krylov subspace to resolve enough eigenvalues from
    # inside a cluster. We use `ncv = max(4*(k+1) + 1, 20)` which is
    # conservative enough to converge on rcoal trees from n = 600 to
    # n = 10,000, empirically in a few seconds.
    ncv <- as.integer(min(n, max(4L * (k + 1L) + 1L, 20L)))
    eig <- tryCatch(
      RSpectra::eigs_sym(L, k = k + 1L, which = "SA",
                         opts = list(maxitr = 10000L, tol = 1e-9,
                                     ncv = ncv)),
      warning = function(w) NULL,
      error   = function(e) NULL
    )
    if (!is.null(eig) && length(eig$values) == k + 1L) {
      ord  <- order(eig$values, decreasing = FALSE)
      vecs <- eig$vectors[, ord, drop = FALSE]
      return(vecs[, 2L:(k + 1L), drop = FALSE])
    }
    # If RSpectra failed or returned fewer eigenvalues than requested,
    # fall through to the dense path.
    warning("RSpectra::eigs_sym did not converge with k = ", k + 1L,
            " eigenvalues; falling back to dense eigen().",
            call. = FALSE)
  }

  # Dense fallback: base::eigen() is O(n^3) but rock-solid for small n
  # and is the reference implementation we test against.
  eig  <- eigen(L, symmetric = TRUE)
  ord  <- order(eig$values, decreasing = FALSE)
  vecs <- eig$vectors[, ord, drop = FALSE]
  vecs[, 2L:(k + 1L), drop = FALSE]
}

# Symmetric-normalised Gaussian kernel adjacency + self-loops,
# computed from a precomputed cophenetic distance matrix.
adj_symnorm_from_D <- function(D, sigma_mult) {
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
#' @section Scaling:
#' For trees with more than 7500 tips, \code{build_phylo_graph()} uses
#' \pkg{RSpectra}'s sparse Lanczos eigensolver (\code{eigs_sym}) instead
#' of the dense \code{base::eigen()}. This reduces the spectral step
#' from \eqn{O(n^3)} (hours at \eqn{n = 10{,}000}) to \eqn{O(n \cdot k
#' \cdot \mathrm{iters})} (seconds to a minute). The threshold is
#' conservative because the dense eigensolver is competitive on small
#' and mid-size trees, and on pathological ultrametric simulations
#' (\code{ape::rcoal}) the Laplacian spectrum has tight degenerate
#' clusters that slow down Lanczos until the tree is fairly large. On
#' realistic asymmetric phylogenies with variable branch lengths
#' (e.g. the AVONET BirdTree Stage2 Hackett tree) the sparse path wins
#' much earlier, and at \eqn{n = 9{,}993} it delivers a ~10x speedup
#' over dense. The dense path is kept as the default for smaller trees
#' and as a fallback when \pkg{RSpectra} is not installed, or if
#' Lanczos fails to converge. The cophenetic distance matrix is
#' computed once per call and reused for both the adjacency and the
#' spectral features; it is also returned in the result so that
#' downstream consumers (notably \code{\link{fit_baseline}}) can skip
#' recomputation.
#'
#' @param tree object of class \code{"phylo"}.
#' @param k_eigen integer or \code{"auto"}. Number of Laplacian eigenvectors
#'   to use as node features.  When \code{"auto"} (default), scales with tree
#'   size: \code{min(max(ceiling(n/20), 4), 32)}, giving 4 for very small
#'   trees, 8 for 100-160 tips, 15 for 300 tips, and 32 for 640+ tips.
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
#'     \item{D}{Numeric matrix (n x n). Cophenetic (patristic) distances
#'       between tips, row/column-ordered by \code{tree$tip.label}.
#'       Returned so downstream functions can reuse it instead of calling
#'       \code{ape::cophenetic.phylo()} again.}
#'     \item{R_phy}{Numeric matrix (n x n). Phylogenetic correlation matrix
#'       \code{cov2cor(ape::vcv(tree))} (diagonal = 1). Used by
#'       \code{\link{fit_baseline}} for BM conditional imputation.}
#'   }
#' @examples
#' set.seed(1)
#' tree <- ape::rtree(30)
#' g <- build_phylo_graph(tree, k_eigen = 4)
#' dim(g$adj)    # 30 x 30
#' dim(g$coords) # 30 x 4
#' dim(g$D)      # 30 x 30
#' @importFrom ape cophenetic.phylo
#' @export
build_phylo_graph <- function(tree, k_eigen = "auto", sigma_mult = 0.5,
                               cache_path = NULL) {
  if (!inherits(tree, "phylo")) stop("'tree' must be a phylo object.")
  n <- length(tree$tip.label)

  # Adaptive k_eigen: scales with tree size

  if (identical(k_eigen, "auto") || identical(k_eigen, "Auto")) {
    k_eigen <- as.integer(min(max(ceiling(n / 20), 4L), 32L))
  }
  k_eigen <- as.integer(k_eigen)

  if (n > 10000L) {
    warning(
      "build_phylo_graph: tree has ", n, " tips. ",
      "cophenetic() requires O(n^2) memory (~",
      round(n^2 * 8 / 1e6), " MB). ",
      "Consider subsetting your tree for very large datasets."
    )
  }

  # Load from cache if available and dimensions match
  if (!is.null(cache_path) && file.exists(cache_path)) {
    cache <- readRDS(cache_path)
    if (isTRUE(nrow(cache$adj) == n) && isTRUE(ncol(cache$coords) == k_eigen)) {
      # Older caches may not have $D; backfill so downstream code can rely
      # on its presence.
      if (is.null(cache$D)) {
        cache$D <- ape::cophenetic.phylo(tree)
      }
      if (is.null(cache$R_phy)) {
        cache$R_phy <- phylo_cor_matrix(tree)
      }
      return(cache)
    }
    message("Cache dimensions mismatch -- recomputing.")
  }

  # Compute the cophenetic distance matrix exactly once and reuse it for
  # the adjacency, the spectral features, and sigma. The previous
  # implementation computed cophenetic() three times inside this function
  # and a fourth time inside fit_baseline(); caching D here collapses all
  # four calls into one.
  D <- ape::cophenetic.phylo(tree)

  # Phylogenetic correlation matrix: R = cov2cor(vcv(tree)).

  # Used by fit_baseline() for the internal BM imputation.
  # Cached here alongside D to avoid a second O(n^2) computation.
  R_phy <- phylo_cor_matrix(tree)

  adj    <- adj_symnorm_from_D(D, sigma_mult)
  coords <- spectral_features_from_D(D, k_eigen, sigma_mult)
  sigma  <- stats::median(D) * sigma_mult

  result <- list(adj = adj, coords = coords, n = n, sigma = sigma,
                 D = D, R_phy = R_phy)

  if (!is.null(cache_path)) {
    saveRDS(result, cache_path)
  }

  result
}
