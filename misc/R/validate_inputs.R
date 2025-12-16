
#' @noRd
validate_inputs <- function(traits, phylo, env_data, species_id) {
  if (! (is.data.frame(traits) || is.matrix(traits)))
    stop('traits must be data.frame or matrix.')
  if (!inherits(phylo, 'phylo'))
    stop('phylo must be ape::phylo.')
  if (missing(env_data) || missing(species_id))
    stop('env_data and species_id must be provided.')
  if (length(species_id)!=nrow(traits))
    stop('species_id length must equal nrow(traits).')
  if (!all(species_id %in% phylo$tip.label))
    stop('species_id values must match phylogeny tip labels.')
}

