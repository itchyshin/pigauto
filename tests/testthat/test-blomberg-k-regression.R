test_that("Blomberg K is returned through the helper and the public fit", {
  skip_if_not_installed("phytools")
  set.seed(6401)
  tree <- ape::rtree(32)
  traits <- data.frame(x = sin(seq_len(32) / 3) + seq_len(32) / 20,
                       row.names = tree$tip.label)
  traits$x[1] <- NA_real_
  data <- preprocess_traits(traits, tree)
  observed <- !is.na(data$X_scaled[, 1])
  observed_tree <- ape::keep.tip(tree, tree$tip.label[observed])
  values <- stats::setNames(data$X_scaled[observed, 1], tree$tip.label[observed])
  values <- values[observed_tree$tip.label]
  tm <- data$trait_map[[1]]
  reference_k <- phytools::phylosig(observed_tree, values,
    method = "K", test = FALSE)
  reference_k <- if (is.list(reference_k)) reference_k$K else reference_k
  actual <- tryCatch(.phylo_signal_lambda_continuous(tm, data, tree, "blomberg_k", 20L),
                     error = function(e) NA_real_)
  expect_true(is.finite(actual))
  expect_equal(actual, as.numeric(reference_k), tolerance = 1e-10)
  lambda <- .phylo_signal_lambda_continuous(tm, data, tree, "lambda", 20L)
  reference_lambda <- phytools::phylosig(observed_tree, values,
    method = "lambda", test = FALSE)$lambda
  expect_equal(lambda, as.numeric(reference_lambda), tolerance = 1e-10)

  result <- impute(traits, tree, gnn = FALSE, verbose = FALSE, seed = 18L,
                   missing_frac = 0, phylo_signal_gate = TRUE,
                   phylo_signal_method = "blomberg_k")
  expect_equal(result$fit$phylo_signal_method, "blomberg_k")
  expect_true(is.finite(result$fit$phylo_signal_per_trait[["x"]]))
  expect_equal(unname(result$fit$phylo_signal_per_trait[["x"]]),
               as.numeric(reference_k), tolerance = 1e-10)
})
