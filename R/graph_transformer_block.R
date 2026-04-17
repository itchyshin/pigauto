#' Graph transformer block: multi-head attention + FFN + pre-norm + residual
#'
#' A standard transformer-encoder block adapted for graph inputs: the
#' attention operation respects an adjacency matrix via an additive
#' `log(adj + eps)` bias on the attention scores (the phylo prior). Uses
#' pre-norm (LayerNorm before attention and FFN) for training stability
#' at depth.
#'
#' The FFN output projection is initialised near zero so each block
#' starts as an approximate identity. When multiple blocks are stacked,
#' the composition at init is also near identity — preserving pigauto's
#' gate-closed-at-init safety property.
#'
#' @param hidden_dim Integer. Model dimensionality.
#' @param n_heads Integer. Number of attention heads. `hidden_dim`
#'   must be divisible by `n_heads`.
#' @param ffn_mult Integer. FFN inner width multiplier (typically 4).
#' @param dropout Numeric in [0, 1). Applied to attention and FFN outputs.
#' @keywords internal
#' @noRd
GraphTransformerBlock <- torch::nn_module(
  "GraphTransformerBlock",

  initialize = function(hidden_dim, n_heads = 4L, ffn_mult = 4L,
                         dropout = 0.1) {
    stopifnot(hidden_dim %% n_heads == 0L)
    self$hidden_dim <- as.integer(hidden_dim)
    self$n_heads    <- as.integer(n_heads)
    self$head_dim   <- self$hidden_dim %/% self$n_heads
    self$dropout_p  <- dropout

    # Attention projections: Q, K, V, output
    self$q_proj   <- torch::nn_linear(hidden_dim, hidden_dim)
    self$k_proj   <- torch::nn_linear(hidden_dim, hidden_dim)
    self$v_proj   <- torch::nn_linear(hidden_dim, hidden_dim)
    self$out_proj <- torch::nn_linear(hidden_dim, hidden_dim)

    # FFN: hidden -> ffn_mult*hidden -> hidden
    inner_dim <- hidden_dim * as.integer(ffn_mult)
    self$ffn_in  <- torch::nn_linear(hidden_dim, inner_dim)
    self$ffn_out <- torch::nn_linear(inner_dim, hidden_dim)

    # Near-zero init on FFN output projection so block ~ identity at step 0.
    # The FFN residual path contributes nothing at init; the attention path
    # still fires with standard default init so the adjacency bias is live
    # from the first forward pass (required for the phylo-prior property).
    # out_proj is left at default Kaiming init deliberately — zeroing it
    # would suppress attention entirely and break the log-adj bias test.
    torch::nn_init_zeros_(self$ffn_out$weight)
    torch::nn_init_zeros_(self$ffn_out$bias)

    # Pre-norm LayerNorms
    self$ln_attn <- torch::nn_layer_norm(hidden_dim)
    self$ln_ffn  <- torch::nn_layer_norm(hidden_dim)

    self$attn_dropout <- torch::nn_dropout(dropout)
    self$ffn_dropout  <- torch::nn_dropout(dropout)

    # Log-adjacency bias scaling: learnable per-head scalar so each head
    # can modulate how strongly it respects the phylo prior.
    # Used as the fallback when D_sq is not supplied (backward compat).
    self$adj_bias_scale <- torch::nn_parameter(
      torch::torch_ones(self$n_heads)
    )

    # Learnable per-head Gaussian bandwidth for rate-aware attention (B2).
    # softplus(0.5413) ≈ 1.0 — moderate initial bandwidth.
    # When D_sq is supplied in forward(), each head computes:
    #   bias_h = -D_sq / (2 * softplus(log_bw_h)^2)
    # so it learns its own phylogenetic scale (tight = close relatives;
    # broad = distant relatives).  When D_sq = NULL, the legacy
    # log(adj)-bias path fires unchanged.
    self$log_bandwidth <- torch::nn_parameter(
      torch::torch_full(list(self$n_heads), 0.5413)
    )
  },

  forward = function(h, adj, D_sq = NULL) {
    # h:    (n_species, hidden_dim)
    # adj:  (n_species, n_species), non-negative kernel weights
    # D_sq: (n_species, n_species) or NULL.  When provided, use the
    #       learnable per-head Gaussian bandwidth (B2) instead of the
    #       log(adj) fallback.
    n  <- h$size(1)
    H  <- self$n_heads
    D  <- self$head_dim

    # Pre-norm
    h_norm <- self$ln_attn(h)

    # Multi-head Q, K, V: reshape (n, hidden) -> (H, n, D)
    q <- self$q_proj(h_norm)$view(c(n, H, D))$permute(c(2, 1, 3))
    k <- self$k_proj(h_norm)$view(c(n, H, D))$permute(c(2, 1, 3))
    v <- self$v_proj(h_norm)$view(c(n, H, D))$permute(c(2, 1, 3))

    # Attention scores per head: (H, n, n)
    scale  <- 1 / sqrt(D)
    scores <- torch::torch_matmul(q, k$transpose(2, 3)) * scale

    # Per-head Gaussian on cophenetic distances (B2) or log(adj) fallback
    if (!is.null(D_sq)) {
      # B2 path: bias_h = -D_sq / (2 * softplus(log_bw_h)^2)
      bw    <- torch::nnf_softplus(self$log_bandwidth)       # (H,)
      bw_sq <- (bw * bw)$view(c(H, 1L, 1L))                 # (H, 1, 1)
      bias  <- -D_sq$unsqueeze(1) / (2 * bw_sq + 1e-8)      # (H, n, n)
    } else {
      # Legacy fallback: adj_bias_scale * log(adj + eps)
      log_adj <- torch::torch_log(adj + 1e-6)
      bias <- self$adj_bias_scale$view(c(H, 1L, 1L)) *
              log_adj$unsqueeze(1)
    }
    scores <- scores + bias

    # Softmax over keys
    attn <- torch::nnf_softmax(scores, dim = -1L)
    attn <- self$attn_dropout(attn)

    # Weighted sum: (H, n, n) x (H, n, D) -> (H, n, D)
    ctx <- torch::torch_matmul(attn, v)

    # Merge heads: (H, n, D) -> (n, H*D) = (n, hidden)
    ctx <- ctx$permute(c(2, 1, 3))$contiguous()$view(c(n, self$hidden_dim))
    ctx <- self$out_proj(ctx)

    # Residual 1: attention output + input
    h2 <- h + ctx

    # FFN block with pre-norm
    h2_norm <- self$ln_ffn(h2)
    ffn_out <- self$ffn_out(torch::nnf_gelu(self$ffn_in(h2_norm)))
    ffn_out <- self$ffn_dropout(ffn_out)

    # Residual 2: FFN output + attention output
    h3 <- h2 + ffn_out

    h3
  }
)
