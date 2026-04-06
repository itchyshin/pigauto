#' Read a phylogenetic tree from a file
#'
#' Reads a Newick or NEXUS tree file using \pkg{ape}. Tries
#' \code{ape::read.tree} first; falls back to \code{ape::read.nexus} if that
#' fails.
#'
#' @param path character. Path to the tree file.
#' @return An object of class \code{"phylo"}.
#' @examples
#' \dontrun{
#' tree <- read_tree("path/to/tree.tre")
#' }
#' @importFrom ape read.tree read.nexus
#' @export
read_tree <- function(path) {
  if (!file.exists(path)) {
    stop("Tree file not found: ", path)
  }
  tree <- tryCatch(
    ape::read.tree(path),
    error = function(e) NULL
  )
  if (is.null(tree)) {
    tree <- tryCatch(
      ape::read.nexus(path),
      error = function(e) {
        stop("Could not read tree file as Newick or NEXUS: ", path,
             "\nOriginal error: ", conditionMessage(e))
      }
    )
  }
  tree
}
