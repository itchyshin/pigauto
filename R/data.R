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
#'   Species labels are aligned to the bundled example phylogenies.
"avonet300"


#' Example bird phylogeny for the 300 species in \code{avonet300}
#'
#' An object of class \code{'phylo'} from the \pkg{ape} package. The tree is a
#' posterior Hackett-backbone sample distributed by the MIT-licensed
#' \pkg{megatrees} package, pruned to the 300 species present in
#' \code{\link{avonet300}}. It is intended for examples, not as a consensus
#' or recommended tree for substantive analysis.
#'
#' @format An object of class \code{phylo} with 300 tips.
#' @source Li (2026), \pkg{megatrees} 1.0.0, MIT licence; underlying tree
#'   sample from Jetz et al. (2012), Hackett et al. backbone.
"tree300"


#' Full AVONET morphological and ecological trait data for 9,993 bird species
#'
#' The full-scale counterpart to \code{\link{avonet300}}: every bird species
#' for which AVONET3 and the bundled example phylogeny agree on
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
#'   Species labels are aligned to the bundled example phylogeny.
#' @seealso \code{\link{avonet300}}, \code{\link{tree_full}}
"avonet_full"


#' Example bird phylogeny for the species in \code{avonet_full}
#'
#' An object of class \code{'phylo'} from the \pkg{ape} package. The tree is
#' the same posterior Hackett-backbone sample used for \code{\link{tree300}},
#' but pruned to the species present in \code{\link{avonet_full}} rather than
#' a 300-species random subset.
#'
#' Row order in \code{avonet_full} matches tip order in \code{tree_full}:
#' \code{all(avonet_full$Species_Key == tree_full$tip.label)} returns
#' \code{TRUE}.
#'
#' @format An object of class \code{phylo} with 9,993 tips.
#' @source Li (2026), \pkg{megatrees} 1.0.0, MIT licence; underlying tree
#'   sample from Jetz et al. (2012), Hackett et al. backbone.
#' @seealso \code{\link{tree300}}, \code{\link{avonet_full}}
"tree_full"


#' 50 posterior phylogenies for the 300 species in \code{avonet300}
#'
#' A \code{multiPhylo} list of 50 phylogenetic trees randomly sampled from the
#' BirdTree Hackett backbone posterior (Jetz et al. 2012), each pruned to the
#' 300 species in \code{\link{avonet300}}.  These trees capture phylogenetic
#' uncertainty: topologies and branch lengths vary across the posterior sample.
#'
#' Use with \code{\link{multi_impute_trees}} for experimental sensitivity of
#' point imputations to posterior-tree choice. Tree uncertainty is not
#' supported by the analysis-aware MI backend, and these stochastic datasets
#' are not validated for downstream inference.
#'
#' @format An object of class \code{multiPhylo} containing 50 \code{phylo}
#'   objects, each with 300 tips.
#' @source Li (2026), \pkg{megatrees} 1.0.0, MIT licence; trees pruned from
#'   \code{megatrees::get_tree_bird_n100()}. Underlying posterior from Jetz
#'   et al. (2012), Hackett et al. backbone.
#' @seealso \code{\link{tree300}}, \code{\link{avonet300}},
#'   \code{\link{multi_impute_trees}}
"trees300"


#' Simulated multi-observation-per-species CTmax data
#'
#' A simulated dataset mimicking a thermal tolerance study where each species
#' has multiple measurements of critical thermal maximum (CTmax) taken at
#' different acclimation temperatures.  The data-generating process is
#'
#' \deqn{CTmax_{ij} = 38 + phylo_i + 0.50 \times acclim\_temp_{ij} + \epsilon_{ij}}
#'
#' where \eqn{phylo_i} follows Brownian motion on \code{\link{tree300}},
#' the within-species acclimation response ratio is 0.50, and
#' \eqn{\epsilon_{ij} \sim N(0, 1.5)}.  Thirty percent of species are entirely
#' unobserved (all CTmax values are NA), and an additional 15\% of remaining
#' observations are missing at random.
#'
#' This dataset demonstrates pigauto's multi-observation and observation-level
#' covariate support.  Use with \code{\link{tree300}} and set
#' \code{species_col = "species"} in \code{\link{impute}}.
#'
#' @format A data frame with 1,464 rows and 3 variables:
#' \describe{
#'   \item{species}{Character. Species name matching \code{\link{tree300}} tips.}
#'   \item{acclim_temp}{Numeric. Acclimation temperature (degrees C).
#'     This is an observation-level covariate that varies within species.}
#'   \item{CTmax}{Numeric. Critical thermal maximum (degrees C). Contains
#'     NA values for unobserved species and MCAR missingness.}
#' }
#' @source Simulated.  See \code{data-raw/make_ctmax_sim.R} for the
#'   generation script.
#' @seealso \code{\link{tree300}}, \code{\link{impute}}
"ctmax_sim"
