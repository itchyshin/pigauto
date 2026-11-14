test_that("early stopping restores the validation-best model state", {
  skip_if_no_libtorch()
  testthat::local_mocked_bindings(cuda_is_available = function() FALSE,
    backends_mps_is_available = function() FALSE, .package = "torch")
  expect_identical(get_device()$type, "cpu")
  set.seed(41)
  tree <- ape::rtree(12)
  traits <- data.frame(x = seq_len(12) / 10, row.names = tree$tip.label)
  data <- preprocess_traits(traits, tree)
  graph <- build_phylo_graph(tree)
  splits <- make_missing_splits(data$X_scaled, missing_frac = 0.25,
                                val_frac = 0.5, seed = 41)
  run_fit <- function(epochs) fit_pigauto(data, tree, graph = graph, splits = splits,
    gnn = TRUE, epochs = epochs, eval_every = 1L, patience = 2L, lr = 10,
    warmup_epochs = 1L, hidden_dim = 8L, n_gnn_layers = 1L,
    use_transformer_blocks = FALSE, use_attention = TRUE,
    corruption_ramp = 0L, edge_dropout = 0, corruption_rate = 0.8,
    dropout = 0, verbose = FALSE, seed = 41, predict_method = "per_column",
    refine_steps = 1L)
  expect_warning(one <- run_fit(1L), "Small validation set")
  expect_warning(three <- run_fit(3L), "Small validation set")
  expect_equal(one$history$lr[1], three$history$lr[1])
  expect_equal(one$history$loss_rec[1], three$history$loss_rec[1], tolerance = 1e-6)
  expect_equal(one$history$val_loss[1], three$history$val_loss[1], tolerance = 1e-6)
  expect_equal(nrow(three$history), 3L)
  expect_true(all(three$history$val_loss[-1] > three$history$val_loss[1]))
  expect_equal(three$val_rmse, one$val_rmse, tolerance = 1e-6)
  expect_identical(names(one$model_state), names(three$model_state))
  difference <- vapply(names(one$model_state), function(nm) {
    max(abs(as.array(one$model_state[[nm]]) - as.array(three$model_state[[nm]])))
  }, numeric(1))
  expect_lt(max(difference), 1e-6)
})
