#' Read a phylogenetic tree from a file
#'
#' Reads a Newick or NEXUS tree file using \pkg{ape}. Tries
#' \code{ape::read.tree} first; falls back to \code{ape::read.nexus} if that
#' fails.
#'
#' For bird phylogenies, BirdTree provides a species-subset download tool for
#' up to 2,500 species and full-tree distributions for larger sets. Download
#' and unzip the trees, select a Newick tree file, and read it here. BirdTree
#' documents its distributions as Newick trees; this function also accepts
#' NEXUS files.
#'
#' pigauto does not download BirdTree data automatically. Pass the resulting
#' object to \code{impute(traits, tree)}.
#'
#' BirdTree provides tree distributions. Its FAQ recommends using more than
#' 100 draws for full-tree analyses rather than relying on a consensus tree.
#' This function reads one tree file at a time. pigauto's
#' \code{multi_impute_trees()} is experimental and supports descriptive checks
#' of prediction sensitivity only; its stochastic datasets are not validated for
#' downstream inference. See \code{vignette("tree-uncertainty")}.
#'
#' BirdTree asks researchers using full or partial tree data to cite Jetz et
#' al. (2012). If you use the BirdTree web tool, also credit BirdTree.org.
#'
#' @param path character. Path to the tree file.
#' @return An object of class \code{"phylo"}.
#' @seealso \code{\link{multi_impute_trees}}
#' @references Jetz W, Thomas GH, Joy JB, Hartmann K, Mooers AO (2012). The
#'   global diversity of birds in space and time. \emph{Nature}, 491, 444-448.
#'   \doi{10.1038/nature11631}.
#'   BirdTree.org. Downloads, phylogeny subsets, and FAQ.
#'   \url{https://birdtree.org/downloads/}
#'   \url{https://birdtree.org/subsets/}
#'   \url{https://birdtree.org/faq/}.
#' @examples
#' \donttest{
#' path <- tempfile(fileext = ".tre")
#' ape::write.tree(ape::rtree(10L), path)
#' tree <- read_tree(path)
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
