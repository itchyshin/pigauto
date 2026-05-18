# Pagel's lambda BM baseline (v0.10).
#
# Spec: specs/2026-05-18-pagel-lambda-baseline-design.md.
#
# Pagel's lambda is a single-parameter generalization of Brownian motion
# that shrinks the off-diagonal of the phylogenetic correlation matrix
# toward zero:
#
#   R(lambda) = lambda * R + (1 - lambda) * I
#
# This is equivalent to scaling all INTERNAL branch lengths by lambda
# while leaving terminal (tip-incident) edges unchanged: under that
# transform, off-diagonal entries of vcv(tree) scale by lambda while the
# diagonal (root-to-tip total) is preserved. The equivalence is exact for
# ultrametric trees and a close approximation for non-ultrametric trees;
# we rely on the tree-transform form throughout because it keeps the
# downstream sparse Hadfield-Nakagawa machinery working unchanged.
#
# Reference: Pagel (1999) Nature 401: 877-884.

# Scale internal edges of a phylo tree by lambda, with compensation
# on terminal edges that preserves root-to-tip distances for
# ultrametric trees.
#
# Concretely: each terminal edge t (parent height h_p, child = tip)
# is replaced by
#   t' = t + (1 - lambda) * h_p
# and each internal edge i is scaled by lambda. The result is that
# cov2cor(vcv(out)) equals lambda * cov2cor(vcv(tree)) + (1-lambda)*I
# exactly (when the tree is ultrametric).
#
# For non-ultrametric trees the formula is a close approximation but
# not exact; the diagonal of the new correlation matrix still equals
# one, while off-diagonals are scaled approximately by lambda. Users
# with strongly non-ultrametric trees should be aware.
#
# @param tree object of class "phylo".
# @param lambda numeric scalar in [0, 1].
# @return a phylo with edge lengths transformed.
# @noRd
transform_tree_pagel <- function(tree, lambda) {
  if (!inherits(tree, "phylo")) {
    stop("'tree' must be a phylo object.", call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1L ||
      !is.finite(lambda) || lambda < 0 || lambda > 1) {
    stop("'lambda' must be a numeric scalar in [0, 1]; got: ",
         paste(lambda, collapse = ", "), call. = FALSE)
  }
  n_tip <- ape::Ntip(tree)
  # node.depth.edgelength: distance from root for every node, ordered
  # so node ids 1..n_tip are tips and n_tip+1..n_tip+n_int are internal.
  node_depths <- ape::node.depth.edgelength(tree)
  parents <- tree$edge[, 1L]
  children <- tree$edge[, 2L]
  is_terminal <- children <= n_tip

  out <- tree
  # Internal edges: scale by lambda.
  out$edge.length[!is_terminal] <-
    lambda * tree$edge.length[!is_terminal]
  # Terminal edges: t' = t + (1-lambda) * parent_depth, where
  # parent_depth = root-to-parent distance under the original tree.
  parent_depth_term <- node_depths[parents[is_terminal]]
  out$edge.length[is_terminal] <-
    tree$edge.length[is_terminal] + (1 - lambda) * parent_depth_term
  out
}
