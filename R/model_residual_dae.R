# ResidualPhyloDAE — internal nn_module, not exported.
#
# Architecture:
#   Input  : x (n x p), coords (n x k), covs (n x cov_dim)
#   Encoder: two linear layers with ReLU + dropout
#   Message: n_gnn_layers steps of graph message passing, each with
#            layer norm, ReLU, dropout, and per-layer learnable alpha
#   Decoder: two linear layers -> delta (n x p)
#   Output : (1-rs) * BM_baseline + rs * delta,  rs = sigmoid(res_raw) * cap
#
# res_raw is a per-column vector (length p), so each latent column has
# its own learnable gate.  Continuous columns init at ~0.135 effective
# gate; discrete columns init near 0.  Safety comes from lambda_gate
# regularisation, not the architectural cap.

ResidualPhyloDAE <- torch::nn_module(
  "ResidualPhyloDAE",
  initialize = function(input_dim, hidden_dim, coord_dim, cov_dim,
                        per_column_rs = TRUE, n_gnn_layers = 2L,
                        gate_cap = 0.8) {
    total_input <- input_dim + coord_dim + cov_dim

    self$enc1 <- torch::nn_linear(total_input, hidden_dim)
    self$enc2 <- torch::nn_linear(hidden_dim,  hidden_dim)

    # Multi-layer graph message passing: each layer has its own
    # linear transform, layer norm, and learnable alpha gate.
    self$n_gnn_layers <- as.integer(n_gnn_layers)
    self$msg_layers  <- torch::nn_module_list(lapply(
      seq_len(n_gnn_layers), function(i) torch::nn_linear(hidden_dim, hidden_dim)
    ))
    self$layer_norms <- torch::nn_module_list(lapply(
      seq_len(n_gnn_layers), function(i) torch::nn_layer_norm(hidden_dim)
    ))
    # Per-layer learnable alpha gates (init 0.05)
    self$alphas <- torch::nn_module_list(lapply(
      seq_len(n_gnn_layers), function(i) {
        m <- torch::nn_module("AlphaHolder",
          initialize = function() {
            self$val <- torch::nn_parameter(torch::torch_tensor(0.05))
          }
        )
        m()
      }
    ))

    self$dec1 <- torch::nn_linear(hidden_dim,  hidden_dim)
    self$dec2 <- torch::nn_linear(hidden_dim,  input_dim)   # outputs delta

    # Learnable residual scale: sigmoid(res_raw) * gate_cap.
    # Per-column vector allows each latent column to have its own gate.
    self$per_column_rs <- per_column_rs
    self$gate_cap      <- gate_cap
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

    # Multi-step message passing with gated residual connections
    for (l in seq_len(self$n_gnn_layers)) {
      m <- torch::torch_matmul(adj, h)
      m <- self$msg_layers[[l]](m)
      m <- self$layer_norms[[l]](m)
      m <- self$act(m)
      m <- self$drop(m)
      h <- h + self$alphas[[l]]$val * m
    }

    h     <- self$dec1(h)
    h     <- self$act(h)
    delta <- self$dec2(h)

    # Bounded residual scale in (0, gate_cap)
    if (self$per_column_rs) {
      rs <- torch::torch_sigmoid(self$res_raw)$unsqueeze(1L) * self$gate_cap
    } else {
      rs <- torch::torch_sigmoid(self$res_raw) * self$gate_cap
    }

    list(delta = delta, rs = rs)
  }
)
