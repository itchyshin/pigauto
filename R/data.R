#' AVONET morphological and ecological trait data for 300 bird species
#'
#' A data frame with 300 rows and 8 columns: 4 continuous morphometric traits
#' and 3 ecological traits (2 categorical, 1 ordinal) for demonstrating
#' mixed-type imputation.  Species names are in the \code{Species_Key}
#' column.  Set \code{rownames(df) <- df$Species_Key; df$Species_Key <- NULL}
#' before passing to \code{\link{preprocess_traits}}.
#'
#' @format A data frame with 300 rows and 8 variables:
#' \describe{
#'   \item{Species_Key}{Character. Species name in BirdTree format (spaces
#'     replaced by underscores).}
#'   \item{Mass}{Numeric. Body mass (grams).}
#'   \item{Beak.Length_Culmen}{Numeric. Beak length from culmen (mm).}
#'   \item{Tarsus.Length}{Numeric. Tarsus length (mm).}
#'   \item{Wing.Length}{Numeric. Wing length (mm).}
#'   \item{Trophic.Level}{Factor. Dietary category: Carnivore, Herbivore,
#'     Omnivore, or Scavenger.}
#'   \item{Primary.Lifestyle}{Factor. Lifestyle category: Aerial, Aquatic,
#'     Generalist, Insessorial, or Terrestrial.}
#'   \item{Migration}{Ordered factor. Migration strategy: Resident < Partial
#'     < Full.}
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


#' Full AVONET morphological and ecological trait data for 9,993 bird species
#'
#' The full-scale counterpart to \code{\link{avonet300}}: every bird species
#' for which AVONET3 and the BirdTree Stage2 Hackett MCC phylogeny agree on
#' both a species label and a complete set of continuous morphometric
#' measurements. The schema is identical to \code{avonet300} (same trait
#' columns, same factor encodings, same \code{Species_Key} column) so any
#' code that runs on \code{avonet300} runs on \code{avonet_full} with no
#' modification.
#'
#' Unlike the 300-species subset, native AVONET missingness is PRESERVED in
#' the two ecological columns: 4 NAs in \code{Trophic.Level} and 20 NAs in
#' \code{Migration}. The four continuous columns are complete by construction.
#'
#' Use \code{avonet_full} + \code{\link{tree_full}} to exercise pigauto at
#' real-world scale (the AVONET missingness sweep benchmark uses this
#' dataset). For quick examples or unit tests, prefer the 300-species
#' \code{\link{avonet300}} subset -- it is ~30x smaller and loads instantly.
#'
#' See \code{\link{avonet300}} for the full column schema (same traits and
#' encodings).
#'
#' @format A data frame with 9,993 rows and 8 variables. Columns match
#'   \code{\link{avonet300}}: \code{Species_Key}, \code{Mass},
#'   \code{Beak.Length_Culmen}, \code{Tarsus.Length}, \code{Wing.Length},
#'   \code{Trophic.Level}, \code{Primary.Lifestyle}, \code{Migration}.
#' @source Tobias et al. (2022) AVONET: morphological, ecological and
#'   geographical data for all birds. \emph{Ecology Letters}, 25, 581-597.
#'   BirdTree backbone: Hackett et al. MCC tree via BirdTree.org.
#' @seealso \code{\link{avonet300}}, \code{\link{tree_full}}
"avonet_full"


#' Pruned BirdTree phylogeny for the 9,993 species in \code{avonet_full}
#'
#' An object of class \code{'phylo'} from the \pkg{ape} package. The tree is
#' the same Stage2_Hackett_MCC Maximum Clade Credibility tree used for
#' \code{\link{tree300}}, but pruned to the 9,993 species present in
#' \code{\link{avonet_full}} rather than a 300-species random subset.
#'
#' Row order in \code{avonet_full} matches tip order in \code{tree_full}:
#' \code{all(avonet_full$Species_Key == tree_full$tip.label)} returns
#' \code{TRUE}.
#'
#' @format An object of class \code{phylo} with 9,993 tips.
#' @source BirdTree.org (Jetz et al. 2012, Hackett et al. backbone).
#' @seealso \code{\link{tree300}}, \code{\link{avonet_full}}
"tree_full"
