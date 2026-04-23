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

test_that("fit_pigauto(phylo_signal_gate = TRUE) stores phylo_signal slots", {
  skip_if_not_installed("phytools")
  data("avonet300", package = "pigauto")
  data("tree300",   package = "pigauto")
  df <- avonet300
  if ("Species_Key" %in% colnames(df)) {
    rownames(df) <- df$Species_Key; df$Species_Key <- NULL
  }
  set.seed(2026L)
  df$Mass[sample(300, 30)] <- NA_real_
  res <- pigauto::impute(df, tree300,
                           phylo_signal_gate = TRUE,
                           phylo_signal_threshold = 0.2,
                           safety_floor = TRUE,
                           epochs = 50L, n_imputations = 1L,
                           verbose = FALSE, seed = 2026L)
  fit <- res$fit
  expect_true(!is.null(fit$phylo_signal_per_trait))
  expect_true(!is.null(fit$phylo_gate_triggered))
  expect_equal(fit$phylo_signal_method, "lambda")
  expect_equal(fit$phylo_signal_threshold, 0.2)
  # AVONET traits all have strong phylogenetic signal; none should be gated
  expect_false(any(fit$phylo_gate_triggered, na.rm = TRUE))
})

test_that("fit_pigauto gates weak-signal traits to (0, 0, 1) exactly", {
  skip_if_not_installed("phytools")
  set.seed(2026L)
  n <- 200L
  tree <- ape::rcoal(n)
  # White-noise trait — lambda ~ 0
  traits <- data.frame(noise = stats::rnorm(n),
                        row.names = tree$tip.label)
  # Mask some cells
  traits$noise[sample(n, 30)] <- NA_real_
  res <- pigauto::impute(traits, tree,
                           phylo_signal_gate = TRUE,
                           phylo_signal_threshold = 0.2,
                           safety_floor = TRUE,
                           epochs = 30L, n_imputations = 1L,
                           verbose = FALSE, seed = 2026L)
  fit <- res$fit
  expect_true(isTRUE(as.logical(fit$phylo_gate_triggered["noise"])))
  # r_cal_bm = 0, r_cal_gnn = 0, r_cal_mean = 1 on the noise latent col
  expect_equal(as.numeric(fit$r_cal_bm["noise"]),   0, tolerance = 1e-10)
  expect_equal(as.numeric(fit$r_cal_gnn["noise"]),  0, tolerance = 1e-10)
  expect_equal(as.numeric(fit$r_cal_mean["noise"]), 1, tolerance = 1e-10)
})

test_that("impute(phylo_signal_gate = TRUE) is the default and threads to fit", {
  skip_if_not_installed("phytools")
  data("avonet300", package = "pigauto")
  data("tree300",   package = "pigauto")
  df <- avonet300
  if ("Species_Key" %in% colnames(df)) {
    rownames(df) <- df$Species_Key; df$Species_Key <- NULL
  }
  set.seed(2026L)
  df$Mass[sample(300, 30)] <- NA_real_
  r1 <- pigauto::impute(df, tree300, epochs = 30L, n_imputations = 1L,
                          verbose = FALSE, seed = 2026L)
  expect_equal(r1$fit$phylo_signal_method, "lambda")
  r2 <- pigauto::impute(df, tree300, phylo_signal_gate = FALSE,
                          epochs = 30L, n_imputations = 1L,
                          verbose = FALSE, seed = 2026L)
  expect_true(all(is.na(r2$fit$phylo_signal_per_trait)))
  expect_false(any(r2$fit$phylo_gate_triggered, na.rm = TRUE))
})
