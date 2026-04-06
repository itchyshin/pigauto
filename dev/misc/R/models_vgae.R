
#' Build placeholder VGAE models
#' @noRd
build_vgae <- function(input_dim, latent_dim) {
  encoder <- torch::nn_module(
    initialize = function() {
      self$fc1 <- torch::nn_linear(input_dim, latent_dim)
    },
    forward = function(x, A) {
      h <- torch::torch_relu(self$fc1(x))
      list(mu=h, logvar=torch::torch_zeros_like(h))
    }
  )
  decoder <- torch::nn_module(
    initialize = function() {
      self$fc <- torch::nn_linear(latent_dim, input_dim)
    },
    forward = function(z) self$fc(z)
  )
  list(encoder=encoder(), decoder=decoder())
}

