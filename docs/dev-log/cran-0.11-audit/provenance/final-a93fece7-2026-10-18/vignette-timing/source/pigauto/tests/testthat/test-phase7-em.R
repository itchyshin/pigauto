# Phase 7 EM tests.
#
# Invariants guarded:
#   1. em_offdiag = FALSE preserves Phase 6 behaviour at every
#      em_iterations value (byte-identical regression).
#   2. em_offdiag = TRUE is silently clamped to FALSE when em_iterations
#      < 2L (no previous Σ to condition on).
#   3. em_offdiag = TRUE, em_iterations >= 2L on correlated-binary data
#      runs to completion, em_state$em_offdiag == TRUE, and produces
#      different liability posteriors vs Phase 6.
#   4. build_conditional_prior handles missing-pattern fallback.

test_that("em_offdiag = FALSE at em_iterations = 0 matches default fit_baseline", {
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

  bl_default <- pigauto::fit_baseline(pd, tree, splits = sp, graph = gr)
  bl_off0    <- pigauto::fit_baseline(pd, tree, splits = sp, graph = gr,
                                        em_iterations = 0L, em_offdiag = FALSE)
  bl_off0T   <- pigauto::fit_baseline(pd, tree, splits = sp, graph = gr,
                                        em_iterations = 0L, em_offdiag = TRUE)

  expect_equal(bl_off0$mu,  bl_default$mu)
  expect_equal(bl_off0T$mu, bl_default$mu)  # offdiag clamped at em=0
})

test_that("em_offdiag = TRUE at em_iterations = 1 silently clamps to FALSE", {
  skip_if_not_installed("Rphylopars")
  skip_if_not(pigauto:::joint_mvn_available())

  set.seed(7L)
  tree <- ape::rcoal(50L)
  df   <- pigauto:::simulate_bm_traits(tree, n_traits = 2L, seed = 7L)
  df$trait2 <- factor(ifelse(df$trait2 > median(df$trait2), "1", "0"),
                       levels = c("0", "1"))
  df$trait1[1:10] <- NA

  pd <- pigauto::preprocess_traits(df, tree, log_transform = FALSE)
  sp <- pigauto::make_missing_splits(pd$X_scaled, missing_frac = 0.2,
                                       seed = 1L, trait_map = pd$trait_map)
  gr <- pigauto::build_phylo_graph(tree, k_eigen = 4L)

  bl_off1F <- pigauto::fit_baseline(pd, tree, splits = sp, graph = gr,
                                      em_iterations = 1L, em_offdiag = FALSE)
  bl_off1T <- pigauto::fit_baseline(pd, tree, splits = sp, graph = gr,
                                      em_iterations = 1L, em_offdiag = TRUE)

  # Both should produce identical mu since iter 1 is plug-in regardless.
  expect_equal(bl_off1F$mu, bl_off1T$mu)
})

test_that("em_offdiag = TRUE, em_iterations >= 2L runs on correlated binary", {
  skip_if_not_installed("Rphylopars")
  skip_if_not(pigauto:::joint_mvn_available())

  set.seed(42L)
  tree <- ape::rcoal(100L)
  z <- ape::rTraitCont(tree, model = "BM")
  df <- data.frame(
    b1 = factor(ifelse(z + stats::rnorm(100, sd = 0.2) > 0, "1", "0"),
                 levels = c("0", "1")),
    b2 = factor(ifelse(z + stats::rnorm(100, sd = 0.2) > 0, "1", "0"),
                 levels = c("0", "1")),
    c1 = z + stats::rnorm(100, sd = 0.5),
    row.names = tree$tip.label
  )
  df$b1[1:15]  <- NA
  df$b2[16:30] <- NA

  pd <- pigauto::preprocess_traits(df, tree, log_transform = FALSE)
  sp <- pigauto::make_missing_splits(pd$X_scaled, missing_frac = 0.2,
                                       seed = 1L, trait_map = pd$trait_map)

  out_p6 <- pigauto:::fit_joint_threshold_baseline_em(pd, tree, splits = sp,
                                                      em_iterations = 5L,
                                                      em_offdiag = FALSE)
  out_p7 <- pigauto:::fit_joint_threshold_baseline_em(pd, tree, splits = sp,
                                                      em_iterations = 5L,
                                                      em_offdiag = TRUE)

  expect_true(out_p7$em_state$em_offdiag)
  expect_false(out_p6$em_state$em_offdiag)
  expect_true(out_p7$em_state$iterations_run >= 1L)
  # Phase 7 should produce different posteriors than Phase 6 on correlated
  # data. (At iter 1 they agree; difference shows after iter 2+.)
  expect_false(identical(out_p6$mu_liab, out_p7$mu_liab))
})

test_that("build_conditional_prior sanity + missing-pattern fallback", {
  # Uncorrelated Σ: conditional prior = unconditional.
  Sigma <- diag(c(2, 3, 5))
  L     <- matrix(stats::rnorm(9), 3, 3)
  cp <- pigauto:::build_conditional_prior(Sigma, L)
  expect_equal(cp$mu_prior, matrix(0, 3, 3))
  expect_equal(cp$sd_prior[1, ], sqrt(diag(Sigma)))

  # All-NA row falls back to unconditional (mu=0, sd=sqrt(Σ_jj)).
  Sigma2 <- matrix(c(1, 0.6, 0.3,
                     0.6, 1, 0.2,
                     0.3, 0.2, 1), 3, 3)
  L_na   <- matrix(rep(NA_real_, 3), nrow = 1)
  cp_na  <- pigauto:::build_conditional_prior(Sigma2, L_na)
  expect_equal(as.numeric(cp_na$mu_prior), rep(0, 3))
  expect_equal(as.numeric(cp_na$sd_prior), sqrt(diag(Sigma2)))

  # Partial-NA: observed subset drives the conditioning.
  L_partial <- matrix(c(NA, 0.8, 0.2), nrow = 1)
  cp_part   <- pigauto:::build_conditional_prior(Sigma2, L_partial)
  # mu_prior[1, 1] uses obs cols 2, 3. Value depends on Σ; just verify
  # it's finite and non-zero.
  expect_true(is.finite(cp_part$mu_prior[1, 1]))
  expect_false(cp_part$mu_prior[1, 1] == 0)

  # All-observed: standard conditional-MVN formula produces finite output.
  L_obs <- matrix(c(0.5, 0.5, 0.5), nrow = 1)
  cp_obs <- pigauto:::build_conditional_prior(Sigma2, L_obs)
  expect_true(all(is.finite(cp_obs$mu_prior)))
  expect_true(all(cp_obs$sd_prior > 0))
})
