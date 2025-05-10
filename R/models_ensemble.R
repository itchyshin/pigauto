
#' Build ensemble of dropout models
#' @noRd
build_ensemble <- function(input_dim, latent_dim, n=5) {
  lapply(seq_len(n), function(i) build_dropout_ae(input_dim, latent_dim))
}

