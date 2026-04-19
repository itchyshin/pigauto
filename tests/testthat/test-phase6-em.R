# Phase 6 EM tests.
#
# Core invariants being guarded:
#   1. em_iterations = 0L produces v0.9.1-identical baseline output.
#   2. em_iterations = 1L through the EM wrapper produces the same numerical
#      baseline as the single-pass path (degenerate, but verifies the
#      plumbing).
#   3. em_iterations = 5L runs to completion on mixed-type data without
#      error and attaches em_state to the baseline / OVR output.
#   4. em_iterations = 5L works on a K = 3 OVR categorical fit.

test_that("em_iterations = 0L: fit_baseline byte-identical to default call", {
  skip_if_not_installed("Rphylopars")
  skip_if_not(pigauto:::joint_mvn_available())

  # n = 50 and median-split binary so phylopars gets a balanced, non-singular
  # fit. At n = 30 with `> 0` split, trait2 is too imbalanced and Rphylopars'
  # internal solve() hits a near-singular matrix. This is a preexisting edge
  # case in the single-pass path (not a Phase 6 regression).
  set.seed(7L)
  tree <- ape::rcoal(50L)
  df   <- pigauto:::simulate_bm_traits(tree, n_traits = 2L, seed = 7L)
  df$trait2 <- factor(ifelse(df$trait2 > median(df$trait2), "1", "0"),
                       levels = c("0", "1"))
  df$trait1[1:10]  <- NA
  df$trait2[11:20] <- NA

  pd <- pigauto::preprocess_traits(df, tree, log_transform = FALSE)
  sp <- pigauto::make_missing_splits(pd$X_scaled, missing_frac = 0.2,
                                       seed = 1L, trait_map = pd$trait_map)
  gr <- pigauto::build_phylo_graph(tree, k_eigen = 4L)

  bl_default <- pigauto::fit_baseline(pd, tree, splits = sp, graph = gr)
  bl_em0     <- pigauto::fit_baseline(pd, tree, splits = sp, graph = gr,
                                       em_iterations = 0L)

  expect_equal(bl_em0$mu, bl_default$mu)
  expect_equal(bl_em0$se, bl_default$se)
})

test_that("em_iterations = 1L: threshold-joint EM wrapper matches single-pass", {
  skip_if_not_installed("Rphylopars")
  skip_if_not(pigauto:::joint_mvn_available())

  set.seed(2L)
  tree <- ape::rcoal(40L)
  df   <- pigauto:::simulate_bm_traits(tree, n_traits = 2L, seed = 2L)
  df$trait2 <- factor(ifelse(df$trait2 > 0, "1", "0"), levels = c("0", "1"))
  df$trait1[1:4] <- NA
  df$trait2[5:8] <- NA

  pd <- pigauto::preprocess_traits(df, tree, log_transform = FALSE)
  sp <- pigauto::make_missing_splits(pd$X_scaled, missing_frac = 0.2,
                                       seed = 1L, trait_map = pd$trait_map)

  out0    <- pigauto:::fit_joint_threshold_baseline(pd, tree, splits = sp)
  out_em1 <- pigauto:::fit_joint_threshold_baseline_em(pd, tree, splits = sp,
                                                        em_iterations = 1L)

  expect_equal(out_em1$mu_liab, out0$mu_liab)
  expect_equal(out_em1$se_liab, out0$se_liab)
  # em_state is attached only in the EM path
  expect_false(is.null(out_em1$em_state))
  expect_true(is.null(out0$em_state))
})

test_that("em_iterations = 5L runs to completion on mixed-type data", {
  skip_if_not_installed("Rphylopars")
  skip_if_not(pigauto:::joint_mvn_available())

  set.seed(7L)
  tree <- ape::rcoal(50L)
  df   <- pigauto:::simulate_bm_traits(tree, n_traits = 2L, seed = 7L)
  df$trait2 <- factor(ifelse(df$trait2 > median(df$trait2), "1", "0"),
                       levels = c("0", "1"))
  df$trait1[1:10]  <- NA
  df$trait2[11:20] <- NA

  pd <- pigauto::preprocess_traits(df, tree, log_transform = FALSE)
  sp <- pigauto::make_missing_splits(pd$X_scaled, missing_frac = 0.2,
                                       seed = 1L, trait_map = pd$trait_map)
  gr <- pigauto::build_phylo_graph(tree, k_eigen = 4L)

  bl <- pigauto::fit_baseline(pd, tree, splits = sp, graph = gr,
                                em_iterations = 5L)
  expect_false(is.null(bl))
  expect_equal(dim(bl$mu), c(nrow(pd$X_scaled), ncol(pd$X_scaled)))
  # mu is finite wherever X_scaled has observations
  obs_cells <- !is.na(pd$X_scaled)
  expect_true(all(is.finite(bl$mu[obs_cells])))
})

test_that("em_iterations = 5L works on K = 3 OVR categorical", {
  skip_if_not_installed("Rphylopars")
  skip_if_not(pigauto:::joint_mvn_available())

  set.seed(13L)
  tree <- ape::rcoal(60L)
  df   <- pigauto:::simulate_bm_traits(tree, n_traits = 1L, seed = 13L)
  df$cat <- factor(sample(c("A", "B", "C"), 60, replace = TRUE),
                    levels = c("A", "B", "C"))
  df$cat[1:10] <- NA

  pd <- pigauto::preprocess_traits(df, tree, log_transform = FALSE)
  sp <- pigauto::make_missing_splits(pd$X_scaled, missing_frac = 0.2,
                                       seed = 1L, trait_map = pd$trait_map)
  gr <- pigauto::build_phylo_graph(tree, k_eigen = 4L)

  bl <- pigauto::fit_baseline(pd, tree, splits = sp, graph = gr,
                                em_iterations = 5L)
  expect_false(is.null(bl))
})

test_that("em_iterations = 0L via impute() matches default impute() baseline", {
  skip_if_not_installed("Rphylopars")
  skip_if_not(pigauto:::joint_mvn_available())
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "libtorch not installed")

  set.seed(7L)
  tree <- ape::rcoal(50L)
  df   <- pigauto:::simulate_bm_traits(tree, n_traits = 2L, seed = 7L)
  df$trait2 <- factor(ifelse(df$trait2 > median(df$trait2), "1", "0"),
                       levels = c("0", "1"))
  df$trait1[1:10] <- NA

  r_default <- pigauto::impute(df, tree, log_transform = FALSE,
                                 epochs = 20L, verbose = FALSE, seed = 1L)
  r_em0     <- pigauto::impute(df, tree, log_transform = FALSE,
                                 epochs = 20L, verbose = FALSE, seed = 1L,
                                 em_iterations = 0L)

  expect_equal(r_em0$baseline$mu, r_default$baseline$mu)
  expect_equal(r_em0$baseline$se, r_default$baseline$se)
})

test_that("build_liability_matrix sd_prior_vec = rep(1, n) matches default", {
  skip_if_not_installed("Rphylopars")

  set.seed(4L)
  tree <- ape::rcoal(30L)
  df   <- pigauto:::simulate_bm_traits(tree, n_traits = 2L, seed = 4L)
  df$trait2 <- factor(ifelse(df$trait2 > 0, "1", "0"), levels = c("0", "1"))
  df$trait1[1:5] <- NA
  pd <- pigauto::preprocess_traits(df, tree, log_transform = FALSE)

  m1 <- pigauto:::build_liability_matrix(pd)
  m2 <- pigauto:::build_liability_matrix(pd,
    sd_prior_vec = rep(1, length(m1$liab_cols)))
  expect_identical(m1$X_liab, m2$X_liab)

  # sd_prior != 1 should shift the liability
  m3 <- pigauto:::build_liability_matrix(pd,
    sd_prior_vec = rep(0.5, length(m1$liab_cols)))
  # Binary column (trait2) should differ; continuous column should not
  # (continuous cols don't go through estep_liability).
  bin_col <- which(m1$liab_cols %in% {
    trait_map_entries <- pd$trait_map[vapply(pd$trait_map, function(tm) tm$type, character(1)) == "binary"]
    if (length(trait_map_entries) > 0L) trait_map_entries[[1]]$latent_cols else integer(0)
  })
  if (length(bin_col) == 1L) {
    expect_false(identical(m1$X_liab[, bin_col], m3$X_liab[, bin_col]))
  }
})

test_that("rel_frobenius + extract_liability_variances sanity", {
  expect_equal(pigauto:::rel_frobenius(c(1, 1), c(1, 1)), 0)
  expect_true(is.infinite(pigauto:::rel_frobenius(c(1, 1), c(0, 0))))

  fake_mat <- list(pars = list(phylocov = matrix(c(2, 0.5, 0.5, 3), 2, 2)))
  expect_equal(pigauto:::extract_liability_variances(fake_mat), c(2, 3))

  fake_sig <- list(pars = list(sigma2 = 4.0))
  expect_equal(pigauto:::extract_liability_variances(fake_sig), 4)
})
