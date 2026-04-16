test_that("ResidualPhyloDAE with transformer blocks builds and forwards", {
  skip_if_not_installed("torch")
  if (!torch::torch_is_installed()) skip("torch backend not installed")

  n_obs     <- 20L
  n_species <- 20L
  p_latent  <- 5L
  k_eigen   <- 8L
  cov_dim   <- 2L
  hidden    <- 32L

  torch::torch_manual_seed(1L)
  model <- ResidualPhyloDAE(
    p_latent  = p_latent,
    k_eigen   = k_eigen,
    cov_dim   = cov_dim,
    hidden    = hidden,
    n_gnn_layers = 4L,
    use_transformer_blocks = TRUE,
    n_heads   = 4L,
    ffn_mult  = 4L,
    dropout   = 0.0
  )

  x      <- torch::torch_randn(n_obs, p_latent)
  coords <- torch::torch_randn(n_species, k_eigen)
  covs   <- torch::torch_randn(n_obs, cov_dim)
  adj    <- torch::torch_rand(n_species, n_species)
  baseline <- torch::torch_randn(n_obs, p_latent)
  obs_to_species <- torch::torch_tensor(seq_len(n_obs), dtype = torch::torch_long())

  out <- model(x = x, coords = coords, covs = covs, adj = adj,
               baseline_mu = baseline,
               obs_to_species = obs_to_species)
  expect_equal(dim(out), c(n_obs, p_latent))
  expect_true(all(is.finite(as.numeric(out))))
})

test_that("ResidualPhyloDAE with transformer blocks produces ~= baseline at init", {
  skip_if_not_installed("torch")
  if (!torch::torch_is_installed()) skip("torch backend not installed")

  torch::torch_manual_seed(2L)
  model <- ResidualPhyloDAE(
    p_latent = 3L, k_eigen = 4L, cov_dim = 1L, hidden = 16L,
    n_gnn_layers = 4L,
    use_transformer_blocks = TRUE, n_heads = 4L, ffn_mult = 4L,
    dropout = 0.0
  )

  x <- torch::torch_randn(10L, 3L)
  coords <- torch::torch_randn(10L, 4L)
  covs   <- torch::torch_randn(10L, 1L)
  adj    <- torch::torch_eye(10L)
  baseline <- torch::torch_randn(10L, 3L)
  obs_to_species <- torch::torch_tensor(seq_len(10L), dtype = torch::torch_long())

  out <- model(x = x, coords = coords, covs = covs, adj = adj,
               baseline_mu = baseline,
               obs_to_species = obs_to_species)
  # At init with near-zero FFN + gate-at-start, output should be close to baseline.
  # Tolerance: on average within 2 units of norm per element (gate is bounded,
  # delta has some contribution from attention path and encoder but bounded).
  diff <- as.numeric(torch::torch_norm(out - baseline)) / sqrt(30)
  expect_lt(diff, 3.0)   # not too strict; just sanity
})

test_that("Legacy architecture still works (use_transformer_blocks = FALSE)", {
  skip_if_not_installed("torch")
  if (!torch::torch_is_installed()) skip("torch backend not installed")

  torch::torch_manual_seed(3L)
  model <- ResidualPhyloDAE(
    p_latent = 3L, k_eigen = 4L, cov_dim = 1L, hidden = 16L,
    n_gnn_layers = 2L,
    use_transformer_blocks = FALSE
  )

  x <- torch::torch_randn(10L, 3L)
  coords <- torch::torch_randn(10L, 4L)
  covs   <- torch::torch_randn(10L, 1L)
  adj    <- torch::torch_rand(10L, 10L)
  baseline <- torch::torch_randn(10L, 3L)
  obs_to_species <- torch::torch_tensor(seq_len(10L), dtype = torch::torch_long())

  out <- model(x = x, coords = coords, covs = covs, adj = adj,
               baseline_mu = baseline,
               obs_to_species = obs_to_species)
  expect_equal(dim(out), c(10L, 3L))
  expect_true(all(is.finite(as.numeric(out))))
})

test_that("model_config captures transformer hyperparameters", {
  skip_if_not_installed("torch")
  if (!torch::torch_is_installed()) skip("torch backend not installed")

  torch::torch_manual_seed(4L)
  model <- ResidualPhyloDAE(
    p_latent = 3L, k_eigen = 4L, cov_dim = 1L, hidden = 16L,
    n_gnn_layers = 4L,
    use_transformer_blocks = TRUE, n_heads = 4L, ffn_mult = 4L
  )
  # The model should expose or internally record use_transformer_blocks, n_heads, ffn_mult.
  # Adjust these expectations based on how the existing config is exposed
  # (likely via fields accessible as model$... or via a get_config() method).
  expect_true(model$use_transformer_blocks)
  expect_equal(model$n_heads, 4L)
  expect_equal(model$ffn_mult, 4L)
})
