# ResidualPhyloDAE -- internal nn_module, not exported.
#
# NAME CLARIFICATION. "Residual" in the class name refers to the
# ResNet-style residual skip connections used INSIDE each GNN layer
# (h_next = norm(h + alpha * msg)). It does NOT mean that the module
# predicts a statistical residual (y - baseline). The output `delta`
# is a full per-cell prediction; it is blended with the baseline
# externally via `(1 - rs) * baseline + rs * delta`, and the loss
# (MSE / BCE / cross-entropy, dispatched by trait type in
# compute_mixed_loss) is computed on the BLEND, not on (y - baseline).
# The shrinkage penalty lambda_shrink * MSE(delta - baseline) is a
# regularisation choice that keeps `delta` close to the baseline when
# the GNN provides no useful correction -- it is not a statistical
# residual model.
#
# Architecture:
#   Input  : x (n_obs x p), coords (n_species x k), covs (n_obs x cov_dim)
#   Encoder: two linear layers with ReLU + dropout
#   Message: n_gnn_layers steps of graph message passing, each with
#            layer norm, ReLU, dropout, and per-layer learnable alpha
#            (ResNet-style skip: h_next = norm(h + alpha * msg)).
#            Two modes: simple adjacency multiplication (default) or
#            attention-based message passing (use_attention = TRUE).
#            Attention uses scaled dot-product with a learnable
#            log-adjacency bias so the model is initialised close to
#            the phylogenetic prior but can learn to deviate.
#            When multi_obs=TRUE, observations are aggregated to species
#            before message passing, then broadcast back.
#   Decoder: two linear layers -> delta (n_obs x p)
#   Output : (1-rs) * baseline + rs * delta,  rs = sigmoid(res_raw) * cap
#
# res_raw is a per-column vector (length p), so each latent column has
# its own learnable blend gate.  Continuous columns init at ~0.135
# effective gate; discrete columns init near 0.  Safety comes from
# lambda_gate regularisation, not the architectural cap.

ResidualPhyloDAE <- torch::nn_module(
  "ResidualPhyloDAE",
  initialize = function(input_dim, hidden_dim, coord_dim, cov_dim,
                        per_column_rs = TRUE, n_gnn_layers = 2L,
                        gate_cap = 0.8, use_attention = FALSE) {
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

    # Attention-based message passing (optional) --------------------------------
    self$use_attention <- use_attention
    if (use_attention) {
      attn_dim <- as.integer(max(hidden_dim %/% 4L, 16L))
      self$attn_dim <- attn_dim
      self$query_layers <- torch::nn_module_list(lapply(
        seq_len(n_gnn_layers), function(i) torch::nn_linear(hidden_dim, attn_dim)
      ))
      self$key_layers <- torch::nn_module_list(lapply(
        seq_len(n_gnn_layers), function(i) torch::nn_linear(hidden_dim, attn_dim)
      ))
      self$value_layers <- torch::nn_module_list(lapply(
        seq_len(n_gnn_layers), function(i) torch::nn_linear(hidden_dim, hidden_dim)
      ))
      # Learnable scale for adjacency bias (initialised to 1.0)
      self$attn_bias_scale <- torch::nn_module_list(lapply(
        seq_len(n_gnn_layers), function(i) {
          m <- torch::nn_module("BiasScale",
            initialize = function() {
              self$val <- torch::nn_parameter(torch::torch_tensor(1.0))
            }
          )
          m()
        }
      ))
      self$attn_drop <- torch::nn_dropout(0.10)
    }

    self$dec1 <- torch::nn_linear(hidden_dim,  hidden_dim)
    self$dec2 <- torch::nn_linear(hidden_dim,  input_dim)   # outputs delta

    # Learnable blend gate: sigmoid(res_raw) * gate_cap. Historically
    # named "residual scale" because it weights the GNN branch of the
    # blend; it is NOT a statistical residual coefficient.
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

  forward = function(x, coords, covs, adj, obs_to_species = NULL) {
    # x:     (n_obs x p)      -- input traits (possibly with mask tokens)
    # coords: (n_species x k) -- spectral coordinates (species-level)
    # covs:  (n_obs x cov_dim) -- covariates (baseline + mask indicator)
    # adj:   (n_species x n_species) -- phylogenetic adjacency
    # obs_to_species: long tensor (n_obs,) -- maps obs to species index.
    #   NULL when single-obs (n_obs = n_species).

    multi_obs <- !is.null(obs_to_species)

    if (multi_obs) {
      # Expand species-level coords to observation level
      coords_obs <- coords$index_select(1L, obs_to_species)
    } else {
      coords_obs <- coords
    }

    combined <- torch::torch_cat(list(x, coords_obs, covs), dim = 2L)

    h <- self$enc1(combined)
    h <- self$act(h)
    h <- self$drop(h)
    h <- self$enc2(h)
    h <- self$act(h)

    # ---- Message passing (species level) ------------------------------------
    if (multi_obs) {
      # Aggregate observations to species: scatter_mean
      n_species <- as.integer(adj$size(1))
      hidden    <- as.integer(h$size(2))
      h_species <- scatter_mean(h, obs_to_species, n_species)
    } else {
      h_species <- h
    }

    # Precompute log-adjacency for attention bias (only when using attention)
    if (self$use_attention) {
      log_adj <- torch::torch_log(adj$clamp(min = 1e-8))
    }

    for (l in seq_len(self$n_gnn_layers)) {
      if (self$use_attention) {
        # Attention-based message passing with phylogenetic prior
        Q <- self$query_layers[[l]](h_species)
        K <- self$key_layers[[l]](h_species)
        V <- self$value_layers[[l]](h_species)

        # Scaled dot-product attention + adjacency bias
        attn_scores <- torch::torch_matmul(Q, K$t()) / sqrt(self$attn_dim)
        attn_scores <- attn_scores + self$attn_bias_scale[[l]]$val * log_adj

        attn_weights <- torch::nnf_softmax(attn_scores, dim = 2L)
        attn_weights <- self$attn_drop(attn_weights)

        m <- torch::torch_matmul(attn_weights, V)
      } else {
        m <- torch::torch_matmul(adj, h_species)
        m <- self$msg_layers[[l]](m)
      }
      m <- self$layer_norms[[l]](m)
      m <- self$act(m)
      m <- self$drop(m)
      h_species <- h_species + self$alphas[[l]]$val * m
    }

    if (multi_obs) {
      # Expand back from species to observation level
      h <- h_species$index_select(1L, obs_to_species)
    } else {
      h <- h_species
    }

    h     <- self$dec1(h)
    h     <- self$act(h)
    delta <- self$dec2(h)

    # Bounded blend gate in (0, gate_cap)
    if (self$per_column_rs) {
      rs <- torch::torch_sigmoid(self$res_raw)$unsqueeze(1L) * self$gate_cap
    } else {
      rs <- torch::torch_sigmoid(self$res_raw) * self$gate_cap
    }

    list(delta = delta, rs = rs)
  }
)


# ---- Internal: scatter_mean for aggregation ----------------------------------
# Computes the mean of src rows for each group defined by index.
# src:   (n_obs x dim) tensor
# index: (n_obs,) long tensor with values in 1..n_groups (1-indexed, R convention)
# n_groups: integer, number of output groups (species)
# Returns: (n_groups x dim) tensor

scatter_mean <- function(src, index, n_groups) {
  dim <- as.integer(src$size(2))
  device <- src$device

  # Sum via scatter_add
  out <- torch::torch_zeros(c(n_groups, dim), dtype = src$dtype, device = device)
  idx_expand <- index$unsqueeze(2L)$expand(c(-1L, dim))
  out$scatter_add_(1L, idx_expand, src)

  # Count per group
  ones <- torch::torch_ones(c(src$size(1), 1L), dtype = src$dtype, device = device)
  counts <- torch::torch_zeros(c(n_groups, 1L), dtype = src$dtype, device = device)
  idx_count <- index$unsqueeze(2L)
  counts$scatter_add_(1L, idx_count, ones)
  counts <- counts$clamp(min = 1)  # avoid division by zero

  out / counts
}
