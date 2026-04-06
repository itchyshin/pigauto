#' AVONET morphological trait data for 300 bird species
#'
#' A data frame with 300 rows and 5 columns: Species_Key (character), Mass,
#' Beak.Length_Culmen, Tarsus.Length, Wing.Length (all numeric, in original
#' units). Log-transformed and z-scored values are created by
#' \code{\link{preprocess_traits}}.
#'
#' @format A data frame with 300 rows and 5 variables:
#' \describe{
#'   \item{Species_Key}{Character. Species name in BirdTree format (spaces
#'     replaced by underscores).}
#'   \item{Mass}{Numeric. Body mass (grams).}
#'   \item{Beak.Length_Culmen}{Numeric. Beak length from culmen (mm).}
#'   \item{Tarsus.Length}{Numeric. Tarsus length (mm).}
#'   \item{Wing.Length}{Numeric. Wing length (mm).}
#' }
#' @source Tobias et al. (2022) AVONET: morphological, ecological and
#'   geographical data for all birds. \emph{Ecology Letters}, 25, 581-597.
#'   BirdTree backbone: Hackett et al. MCC tree via BirdTree.org.
"avonet300"


#' Pruned BirdTree phylogeny for the 300 species in \code{avonet300}
#'
#' An object of class \code{'phylo'} from the \pkg{ape} package. The tree is a
#' Maximum Clade Credibility (MCC) tree from the Hackett et al. backbone
#' (Stage2_Hackett_MCC_no_neg.tre), pruned to the 300 species present in
#' \code{\link{avonet300}}.
#'
#' @format An object of class \code{phylo} with 300 tips.
#' @source BirdTree.org (Jetz et al. 2012, Hackett et al. backbone).
"tree300"
