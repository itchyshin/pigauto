#' Read a phylogenetic tree from a file
#'
#' Reads a Newick or NEXUS tree file using \pkg{ape}. Tries
#' \code{ape::read.tree} first; falls back to \code{ape::read.nexus} if that
#' fails.
#'
#' For bird phylogenies, BirdTree provides a species-subset download tool for
#' up to 2,500 species and full-tree distributions for larger sets. Download
#' and unzip the trees, then read a Newick tree file here. BirdTree documents
#' its distributions as Newick trees; this function also accepts NEXUS files.
#' A file with one tree returns a \code{"phylo"} object; a file with several
#' trees returns a \code{"multiPhylo"} object.
#'
#' pigauto does not download BirdTree data automatically. Pass a single-tree
#' result to \code{impute(traits, tree)}. A multi-tree result can be passed to
#' \code{multi_impute_trees(traits, trees)} for experimental prediction-
#' sensitivity checks only; this path has not been validated for downstream
#' inference.
#'
#' BirdTree's subset tool defaults to at least 100 draws and recommends a
#' reasonable sample (more than 100) for analyses using its full-tree
#' distributions. pigauto's \code{multi_impute_trees()} remains experimental:
#' it supports descriptive checks of prediction sensitivity only, and its
#' stochastic datasets are not validated for downstream inference. See
#' \code{vignette("tree-uncertainty")}.
#'
#' BirdTree asks researchers using full or partial tree data to cite Jetz et
#' al. (2012). If you use the BirdTree web tool, also credit BirdTree.org.
#'
#' @param path character. Path to the tree file.
#' @return An object of class \code{"phylo"} when the file contains one tree,
#'   or \code{"multiPhylo"} when it contains multiple trees.
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
