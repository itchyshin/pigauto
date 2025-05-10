
#' Build adjacency from phylogenetic tree
#' @noRd
phylo_to_adj <- function(tree) {
  D <- ape::cophenetic.phylo(tree)
  A <- 1/(1 + D)
  diag(A) <- 0
  A
}

