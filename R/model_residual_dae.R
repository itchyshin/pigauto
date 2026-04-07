# ResidualPhyloDAE — internal nn_module, not exported.
#
# Architecture:
#   Input  : x (n x p), coords (n x k), covs (n x cov_dim)
#   Encoder: two linear layers with ReLU + dropout
#   Message: one-step graph message passing with learnable strength alpha
#   Decoder: two linear layers -> delta (n x p)
#   Output : (1-rs) * BM_baseline + rs * delta,  rs = sigmoid(res_raw) * 0.5
#
# res_raw is a per-column vector (length p), so each latent column has
# its own learnable gate.  Initialised at -1.0 (sigmoid(-1)*0.5 ≈ 0.135)
# to allow gradient flow through rs*delta; the shrinkage penalty on delta
# drives unused columns toward zero correction.

ResidualPhyloDAE <- torch::nn_module(
  "ResidualPhyloDAE",
  initialize = function(input_dim, hidden_dim, coord_dim, cov_dim,
                        per_column_rs = TRUE) {
    total_input <- input_dim + coord_dim + cov_dim

    self$enc1 <- torch::nn_linear(total_input, hidden_dim)
    self$enc2 <- torch::nn_linear(hidden_dim,  hidden_dim)

    # Graph message passing: alpha gates the neighbourhood signal
    self$msg   <- torch::nn_linear(hidden_dim, hidden_dim)
    self$alpha <- torch::nn_parameter(torch::torch_tensor(0.05))

    self$dec1 <- torch::nn_linear(hidden_dim,  hidden_dim)
    self$dec2 <- torch::nn_linear(hidden_dim,  input_dim)   # outputs delta

    # Learnable residual scale: sigmoid(res_raw) * 0.5 caps correction at 50%.
    # Per-column vector allows each latent column to have its own gate.
    # Init at -1.0 → sigmoid(-1)*0.5 ≈ 0.135 — enough for reconstruction
    # gradient to flow through rs*delta, while shrinkage on delta still
    # drives unused columns toward zero correction.
    self$per_column_rs <- per_column_rs
    if (per_column_rs) {
      self$res_raw <- torch::nn_parameter(
        torch::torch_full(input_dim, -1.0, dtype = torch::torch_float())
      )
    } else {
      self$res_raw <- torch::nn_parameter(torch::torch_tensor(-1.0))
    }

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
    if (self$per_column_rs) {
      rs <- torch::torch_sigmoid(self$res_raw)$unsqueeze(1L) * 0.5
    } else {
      rs <- torch::torch_sigmoid(self$res_raw) * 0.5
    }

    list(delta = delta, rs = rs)
  }
)
