# ResidualPhyloDAE — internal nn_module, not exported.
#
# Architecture:
#   Input  : x (n x p), coords (n x k), covs (n x cov_dim)
#   Encoder: two linear layers with ReLU + dropout
#   Message: one-step graph message passing with learnable strength alpha
#   Decoder: two linear layers -> delta (n x p)
#   Output : BM_baseline + sigmoid(res_raw) * cap * delta
#
# Source: script/test_new7.R lines 83-117 (generalised to p > 1 traits).

ResidualPhyloDAE <- torch::nn_module(
  "ResidualPhyloDAE",
  initialize = function(input_dim, hidden_dim, coord_dim, cov_dim) {
    total_input <- input_dim + coord_dim + cov_dim

    self$enc1 <- torch::nn_linear(total_input, hidden_dim)
    self$enc2 <- torch::nn_linear(hidden_dim,  hidden_dim)

    # Graph message passing: alpha gates the neighbourhood signal
    self$msg   <- torch::nn_linear(hidden_dim, hidden_dim)
    self$alpha <- torch::nn_parameter(torch::torch_tensor(0.05))

    self$dec1 <- torch::nn_linear(hidden_dim,  hidden_dim)
    self$dec2 <- torch::nn_linear(hidden_dim,  input_dim)   # outputs delta

    # Learnable residual scale: sigmoid(-4) ~ 0.018 initial; cap at 0.5
    self$res_raw <- torch::nn_parameter(torch::torch_tensor(-4.0))

    self$act  <- torch::nn_relu()
    self$drop <- torch::nn_dropout(0.10)

    # Learned mask token replaces corrupted trait values in the input
    self$mask_token <- torch::nn_parameter(torch::torch_zeros(c(1L, input_dim)))
  },

  forward = function(x, coords, covs, adj) {
    combined <- torch::torch_cat(list(x, coords, covs), dim = 2L)

    h <- self$enc1(combined)
    h <- self$act(h)
    h <- self$drop(h)
    h <- self$enc2(h)
    h <- self$act(h)

    # One-step message passing: h <- h + alpha * msg(A h)
    m <- torch::torch_matmul(adj, h)
    h <- h + self$alpha * self$msg(m)

    h     <- self$dec1(h)
    h     <- self$act(h)
    delta <- self$dec2(h)

    # Bounded residual scale in (0, 0.5)
    rs <- torch::torch_sigmoid(self$res_raw) * 0.5

    list(delta = delta, rs = rs)
  }
)
