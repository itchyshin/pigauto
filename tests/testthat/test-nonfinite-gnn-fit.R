test_that("non-finite final GNN outputs fail during fitting with learning rate guidance", {
  skip_if_no_libtorch()
  testthat::local_mocked_bindings(cuda_is_available = function() FALSE,
    backends_mps_is_available = function() FALSE, .package = "torch")
  expect_identical(get_device()$type, "cpu")
  set.seed(41)
  tree <- ape::rtree(12)
  traits <- data.frame(x = seq_len(12) / 10, row.names = tree$tip.label)
  data <- preprocess_traits(traits, tree)
  graph <- build_phylo_graph(tree)
  run_fit <- function(lr, splits = NULL, epochs = 3L) fit_pigauto(data, tree, graph = graph,
    splits = splits, gnn = TRUE, epochs = epochs, eval_every = 1L, patience = 2L,
    lr = lr, warmup_epochs = 1L, hidden_dim = 8L, n_gnn_layers = 1L,
    use_transformer_blocks = FALSE, use_attention = TRUE,
    corruption_ramp = 0L, edge_dropout = 0, corruption_rate = 0.8,
    dropout = 0, verbose = FALSE, seed = 41, predict_method = "per_column",
    refine_steps = 1L)
  expect_error(run_fit(1e10), "non-finite.*learning rate")
  good <- run_fit(0.003)
  expect_s3_class(good, "pigauto_fit")
  expect_true(all(is.finite(predict(good, return_se = FALSE)$imputed_latent)))
  untrained <- run_fit(0.003, epochs = 0L)
  expect_equal(nrow(untrained$history), 0L)
  expect_true(all(is.finite(predict(untrained, return_se = FALSE)$imputed_latent)))
  splits <- make_missing_splits(data$X_scaled, missing_frac = 0.25,
                                val_frac = 0.5, seed = 41)
  expect_error(run_fit(1e10, splits), "non-finite.*learning rate")
  expect_warning(restored <- run_fit(1e4, splits), "Small validation set")
  expect_true(is.finite(restored$val_rmse))
  expect_true(all(is.finite(predict(restored, return_se = FALSE)$imputed_latent)))
})
