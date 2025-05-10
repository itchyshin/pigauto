
#' Impute missing traits
#'
#' High-level wrapper with placeholder backend. Replace TODO sections with
#' full torch training for production.
#' @param traits data.frame of traits (rows = observations/species)
#' @param phylo ape::phylo
#' @param env_data data.frame of environmental predictors
#' @param species_id vector mapping each row to phylo tip label
#' @param method 'vgae','mc_dropout','bootstrap'
#' @param latent_dim integer
#' @param epochs integer training epochs
#' @return phyloimpute_result object
#' @export
impute_phylo <- function(traits, phylo, env_data, species_id,
                         method=c('vgae','mc_dropout','bootstrap'),
                         latent_dim=8, epochs=100) {
  method <- match.arg(method)
  # validate
  validate_inputs(traits, phylo, env_data, species_id)
  types <- detect_trait_types(traits)
  # For placeholder we only numeric
  cont_cols <- types$continuous
  X <- as.matrix(traits[, cont_cols, drop=FALSE])
  mask <- mask_matrix(X)
  X[!mask] <- 0
  if (method=='vgae') {
    model <- build_vgae(ncol(X), latent_dim)
  } else if (method=='mc_dropout') {
    model <- build_dropout_ae(ncol(X), latent_dim)
  } else {
    model <- build_ensemble(ncol(X), latent_dim, n=5)
  }
  # Placeholder impute with column means
  col_means <- colMeans(X)
  imputed <- X
  imputed[!mask] <- col_means[col(imputed)[!mask]]
  res <- list(completed_data = as.data.frame(imputed),
              uncertainty = NULL,
              model=model,
              note='Placeholder functions; add torch training.')
  class(res) <- 'phyloimpute_result'
  res
}

