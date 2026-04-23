# tests/testthat/test-phylo-signal-gate.R
# Smoke canary for phylogenetic-signal gate.
# See specs/2026-04-23-phylo-signal-gate-design.md.

test_that("compute_phylo_signal_per_trait returns lambda for continuous traits", {
  skip_if_not_installed("phytools")
  set.seed(2026L)
  n <- 200L
  tree <- ape::rcoal(n)
  # Strong-signal trait: simulate under BM
  t_strong <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  # Weak-signal trait: white noise
  t_weak   <- stats::rnorm(n)
  traits <- data.frame(strong = t_strong, weak = t_weak,
                        row.names = tree$tip.label)
  data_obj <- pigauto::preprocess_traits(traits, tree)
  lambdas <- pigauto:::compute_phylo_signal_per_trait(
    data = data_obj, tree = tree, method = "lambda")
  expect_type(lambdas, "double")
  expect_named(lambdas, c("strong", "weak"))
  expect_gt(as.numeric(lambdas["strong"]), 0.7)   # strong BM signal
  expect_lt(as.numeric(lambdas["weak"]),   0.3)   # weak/no signal
})

test_that("compute_phylo_signal_per_trait handles constant column as NA", {
  skip_if_not_installed("phytools")
  set.seed(1L)
  n <- 50L
  tree <- ape::rcoal(n)
  traits <- data.frame(const = rep(3.14, n),
                        row.names = tree$tip.label)
  data_obj <- pigauto::preprocess_traits(traits, tree)
  lambdas <- pigauto:::compute_phylo_signal_per_trait(
    data = data_obj, tree = tree, method = "lambda")
  expect_true(is.na(as.numeric(lambdas["const"])))
})

test_that("compute_phylo_signal_per_trait returns NA below min_tips", {
  skip_if_not_installed("phytools")
  set.seed(1L)
  tree <- ape::rcoal(15L)  # below default min_tips = 20
  traits <- data.frame(x = stats::rnorm(15L),
                        row.names = tree$tip.label)
  data_obj <- pigauto::preprocess_traits(traits, tree)
  lambdas <- pigauto:::compute_phylo_signal_per_trait(
    data = data_obj, tree = tree, method = "lambda", min_tips = 20L)
  expect_true(is.na(as.numeric(lambdas["x"])))
})
