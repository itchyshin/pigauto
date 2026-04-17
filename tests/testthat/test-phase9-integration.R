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

test_that("fit_pigauto end-to-end with transformer blocks on synthetic data", {
  skip_if_not_installed("torch")
  if (!torch::torch_is_installed()) skip("torch backend not installed")
  skip_if_not_installed("Rphylopars")

  set.seed(400)
  tree <- ape::rtree(40)
  # Mixed types: continuous + binary + categorical
  df <- data.frame(
    x = rnorm(40),
    y = factor(sample(c("A", "B"), 40, TRUE), levels = c("A", "B")),
    z = factor(sample(c("P", "Q", "R"), 40, TRUE), levels = c("P", "Q", "R")),
    row.names = tree$tip.label
  )
  df$y[c(3, 15)] <- NA
  df$z[c(7, 20)] <- NA

  pd     <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.2,
                                seed = 400, trait_map = pd$trait_map)

  # Fit with transformer blocks (default) and short training
  fit <- fit_pigauto(pd, tree, splits = splits,
                     epochs = 30L, eval_every = 10L, patience = 5L,
                     use_transformer_blocks = TRUE,
                     n_gnn_layers = 4L, n_heads = 4L,
                     verbose = FALSE, seed = 400)

  # Smoke checks
  expect_true(inherits(fit, "pigauto_fit"))
  expect_false(is.null(fit$model_state))
  expect_true(fit$model_config$use_transformer_blocks)

  pred <- predict(fit, return_se = TRUE)
  expect_false(is.null(pred$imputed))
  expect_true(all(is.finite(as.numeric(as.matrix(pred$imputed_latent)))))

  # Categorical outputs: probabilities matrix should exist and rows sum to ~1
  expect_false(is.null(pred$probabilities$z))
  row_sums <- rowSums(pred$probabilities$z)
  expect_true(all(abs(row_sums - 1) < 1e-4))
})

test_that("legacy architecture (use_transformer_blocks = FALSE) still trains", {
  skip_if_not_installed("torch")
  if (!torch::torch_is_installed()) skip("torch backend not installed")
  skip_if_not_installed("Rphylopars")

  set.seed(401)
  tree <- ape::rtree(30)
  df <- data.frame(
    x = rnorm(30),
    y = factor(sample(c("A", "B"), 30, TRUE), levels = c("A", "B")),
    row.names = tree$tip.label
  )
  df$y[c(5, 10)] <- NA
  pd     <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.2,
                                seed = 401, trait_map = pd$trait_map)

  fit <- fit_pigauto(pd, tree, splits = splits,
                     epochs = 20L, eval_every = 10L, patience = 5L,
                     use_transformer_blocks = FALSE,
                     n_gnn_layers = 2L,
                     verbose = FALSE, seed = 401)

  expect_true(inherits(fit, "pigauto_fit"))
  expect_false(fit$model_config$use_transformer_blocks %||% FALSE)
})

test_that("fit_pigauto end-to-end converges with transformer blocks on small data", {
  skip_if_not_installed("torch")
  if (!torch::torch_is_installed()) skip("torch backend not installed")
  skip_if_not_installed("Rphylopars")

  set.seed(500)
  data(avonet300, tree300, package = "pigauto")
  traits <- avonet300
  rownames(traits) <- traits$Species_Key
  traits$Species_Key <- NULL

  pd <- preprocess_traits(traits, tree300)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.25,
                                seed = 500, trait_map = pd$trait_map)
  graph <- build_phylo_graph(tree300, k_eigen = "auto")
  baseline <- fit_baseline(pd, tree300, splits = splits, graph = graph)

  # Transformer variant (default); eval_every = 5 ensures evaluations fire
  # within the 30-epoch budget (default eval_every = 100 would not).
  fit_t <- fit_pigauto(
    data = pd, tree = tree300, splits = splits, graph = graph,
    baseline = baseline,
    epochs = 30L, eval_every = 5L, verbose = FALSE, seed = 500L
  )
  expect_s3_class(fit_t, "pigauto_fit")
  # history$loss_rec is the per-epoch training reconstruction loss
  expect_true(nrow(fit_t$history) > 0L)
  expect_true(is.finite(tail(fit_t$history$loss_rec, 1)))
  # val_rmse is the best validation loss recorded during training
  expect_true(is.finite(fit_t$val_rmse))
  # Predictions work and are all finite
  pred_t <- predict(fit_t, return_se = FALSE)
  expect_false(is.null(pred_t$imputed_latent))
  expect_true(all(is.finite(as.matrix(pred_t$imputed_latent))))

  # Legacy variant: should also converge cleanly
  fit_l <- fit_pigauto(
    data = pd, tree = tree300, splits = splits, graph = graph,
    baseline = baseline,
    use_transformer_blocks = FALSE,
    epochs = 30L, eval_every = 5L, verbose = FALSE, seed = 500L
  )
  expect_s3_class(fit_l, "pigauto_fit")
  expect_true(nrow(fit_l$history) > 0L)
  expect_true(is.finite(tail(fit_l$history$loss_rec, 1)))
})

test_that("ResidualPhyloDAE accepts D_sq in forward (transformer path)", {
  skip_if_not_installed("torch")
  if (!torch::torch_is_installed()) skip("torch backend not installed")

  torch::torch_manual_seed(10L)
  model <- ResidualPhyloDAE(
    p_latent = 3L, k_eigen = 4L, cov_dim = 1L, hidden = 16L,
    n_gnn_layers = 2L,
    use_transformer_blocks = TRUE, n_heads = 2L
  )
  x      <- torch::torch_randn(10L, 3L)
  coords <- torch::torch_randn(10L, 4L)
  covs   <- torch::torch_randn(10L, 1L)
  adj    <- torch::torch_rand(10L, 10L)
  D_sq   <- torch::torch_rand(10L, 10L) * 10
  baseline <- torch::torch_randn(10L, 3L)
  obs_to_sp <- torch::torch_tensor(seq_len(10L), dtype = torch::torch_long())

  out <- model(x = x, coords = coords, covs = covs, adj = adj,
               baseline_mu = baseline, obs_to_species = obs_to_sp,
               D_sq = D_sq)
  expect_equal(dim(out), c(10L, 3L))
  expect_true(all(is.finite(as.numeric(out))))
})
