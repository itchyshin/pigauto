test_that("phylogenetic signal calculation can exclude internal holdouts", {
  data <- list(
    X_scaled = matrix(as.numeric(1:12), nrow = 6L, ncol = 2L),
    trait_map = list(list(name = "y", latent_cols = 1:2))
  )
  original <- data$X_scaled
  splits <- list(val_idx = c(2L, 7L), test_idx = c(5L, 12L))

  signal_data <- .phylo_signal_training_data(data, splits)

  expect_true(all(is.na(signal_data$X_scaled[c(2L, 5L, 7L, 12L)])))
  expect_identical(signal_data$X_scaled[-c(2L, 5L, 7L, 12L)],
                   original[-c(2L, 5L, 7L, 12L)])
  expect_identical(data$X_scaled, original)
})

test_that("phylogenetic signal training mask accepts absent split partitions", {
  data <- list(X_scaled = matrix(1:6, ncol = 1L))
  signal_data <- .phylo_signal_training_data(data, list(val_idx = 2L))
  expect_true(is.na(signal_data$X_scaled[2L]))
  expect_identical(signal_data$X_scaled[-2L], data$X_scaled[-2L])
})

test_that("signal estimates equal an explicitly training-masked fit", {
  skip_if_not_installed("phytools")
  tree <- ape::rtree(50L)
  set.seed(20914L)
  values <- as.numeric(ape::rTraitCont(tree)[tree$tip.label])
  names(values) <- tree$tip.label
  values[47:50] <- c(-20, 20, -30, 30)
  data <- list(
    X_scaled = matrix(values, ncol = 1L,
      dimnames = list(tree$tip.label, "y")),
    species_names = tree$tip.label,
    trait_map = list(list(name = "y", type = "continuous", latent_cols = 1L))
  )
  splits <- list(val_idx = 49:50, test_idx = 47:48)
  manually_masked <- .phylo_signal_training_data(data, splits)

  signal_with_splits <- compute_phylo_signal_per_trait(data, tree,
    method = "lambda", min_tips = 20L, splits = splits)
  signal_pre_masked <- compute_phylo_signal_per_trait(data = manually_masked,
    tree = tree, method = "lambda", min_tips = 20L)

  # Independent oracle: explicitly subset the tree and named training values,
  # then call phytools directly without pigauto's masking/estimation helpers.
  training_names <- tree$tip.label[-c(47:50)]
  training_tree <- ape::keep.tip(tree, training_names)
  training_values <- values[training_tree$tip.label]
  direct_lambda <- phytools::phylosig(training_tree, training_values,
    method = "lambda", test = FALSE)$lambda

  expect_equal(signal_with_splits, signal_pre_masked)
  expect_equal(unname(signal_with_splits), unname(direct_lambda))
  expect_identical(data$X_scaled[, 1L], values)
})

test_that("public phylogenetic-signal gate defaults match the documented policy", {
  expect_identical(formals(fit_pigauto)$phylo_signal_gate, FALSE)
  expect_identical(formals(impute)$phylo_signal_gate, FALSE)
})
