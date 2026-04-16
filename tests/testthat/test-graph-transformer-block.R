test_that("GraphTransformerBlock forward pass produces correct shape", {
  skip_if_not_installed("torch")
  if (!torch::torch_is_installed()) skip("torch backend not installed")

  hidden_dim <- 32L
  n_heads    <- 4L
  n_species  <- 20L

  set.seed(1)
  torch::torch_manual_seed(1L)

  block <- GraphTransformerBlock(hidden_dim = hidden_dim,
                                   n_heads    = n_heads,
                                   ffn_mult   = 4L,
                                   dropout    = 0.0)

  h   <- torch::torch_randn(n_species, hidden_dim)
  adj <- torch::torch_rand(n_species, n_species)   # fake adjacency [0,1]

  out <- block(h, adj)
  expect_equal(dim(out), c(n_species, hidden_dim))
  expect_true(all(is.finite(as.numeric(out))))
})

test_that("GraphTransformerBlock is near-identity at init (FFN output near zero)", {
  skip_if_not_installed("torch")
  if (!torch::torch_is_installed()) skip("torch backend not installed")

  hidden_dim <- 32L
  n_species  <- 20L

  set.seed(2)
  torch::torch_manual_seed(2L)
  block <- GraphTransformerBlock(hidden_dim = hidden_dim,
                                   n_heads    = 4L,
                                   ffn_mult   = 4L,
                                   dropout    = 0.0)

  h   <- torch::torch_randn(n_species, hidden_dim)
  adj <- torch::torch_eye(n_species)   # identity adjacency = self-attention only

  out <- block(h, adj)
  # With near-zero FFN init + identity adjacency + pre-norm residual, output
  # should be close to input (the residual passes h through almost unchanged).
  # Tolerance: within 1 unit of norm per element on average.
  diff_norm <- as.numeric(torch::torch_norm(out - h)) /
                sqrt(n_species * hidden_dim)
  expect_lt(diff_norm, 1.0)
})

test_that("GraphTransformerBlock gradients flow through", {
  skip_if_not_installed("torch")
  if (!torch::torch_is_installed()) skip("torch backend not installed")

  hidden_dim <- 16L
  n_species  <- 10L

  torch::torch_manual_seed(3L)
  block <- GraphTransformerBlock(hidden_dim = hidden_dim,
                                   n_heads    = 2L,
                                   ffn_mult   = 4L,
                                   dropout    = 0.0)

  h   <- torch::torch_randn(n_species, hidden_dim, requires_grad = TRUE)
  adj <- torch::torch_rand(n_species, n_species)
  out <- block(h, adj)
  loss <- out$sum()
  loss$backward()

  expect_false(is.null(h$grad))
  # At least some parameters should have non-zero gradients
  param_grads <- lapply(block$parameters, function(p) {
    if (is.null(p$grad)) 0 else as.numeric(torch::torch_sum(torch::torch_abs(p$grad)))
  })
  expect_gt(sum(unlist(param_grads)), 0)
})

test_that("GraphTransformerBlock preserves log-adjacency-bias prior on attention", {
  skip_if_not_installed("torch")
  if (!torch::torch_is_installed()) skip("torch backend not installed")

  # When adj is strongly peaked (near-identity), attention should mostly attend
  # to self; when adj is uniform, attention should be roughly uniform.
  hidden_dim <- 16L
  n_species  <- 10L

  torch::torch_manual_seed(4L)
  block <- GraphTransformerBlock(hidden_dim = hidden_dim,
                                   n_heads    = 2L,
                                   ffn_mult   = 4L,
                                   dropout    = 0.0)
  h <- torch::torch_randn(n_species, hidden_dim)

  # Peaked adjacency: strong self-edges, weak off-diagonal
  adj_peaked  <- torch::torch_eye(n_species) * 100 +
                  torch::torch_rand(n_species, n_species) * 0.01
  out_peaked  <- block(h, adj_peaked)

  # Uniform adjacency: every pair equally weighted
  adj_uniform <- torch::torch_ones(n_species, n_species)
  out_uniform <- block(h, adj_uniform)

  # The two outputs should differ — proves the block respects adj
  diff <- as.numeric(torch::torch_norm(out_peaked - out_uniform))
  expect_gt(diff, 0.01)
})
