
#' Build placeholder dropout autoencoder
#' @noRd
build_dropout_ae <- function(input_dim, latent_dim) {
  torch::nn_module(
    initialize = function() {
      self$enc <- torch::nn_linear(input_dim, latent_dim)
      self$drop<- torch::nn_dropout(p=0.2)
      self$dec <- torch::nn_linear(latent_dim, input_dim)
    },
    forward = function(x) {
      z <- torch::torch_relu(self$enc(x))
      z <- self$drop(z)
      self$dec(z)
    }
  )()
}

