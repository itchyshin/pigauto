# Tests for S5 of the joint-lambda-default lane (docs/dev-log/2026-09-22-joint-lambda-alignment.md):
# the lambda_mode = "estimate" default flip and its plumbing through
# impute() / fit_pigauto() / multi_impute() / multi_impute_trees().
#
# These tests deliberately favour gnn = FALSE where possible (see
# test-gnn-off.R's header comment): gnn = FALSE makes zero torch:: calls, so
# this file does not need skip_if_no_libtorch() at the top level.
#
# fit_baseline() returns `lambda_per_trait` / `lambda_block` and accepts
# `lambda_fixed` for predict-time rebuilds. Tests that read those fields keep a
# skip_if(..., "baseline has no lambda fields") guard so a missing field shows
# up as a visible skip rather than an unrelated error.

# ---- (i) defaults ------------------------------------------------------

test_that("[lambda-default] impute/fit_pigauto/multi_impute default lambda_mode to 'estimate'", {
  expect_identical(eval(formals(impute)$lambda_mode)[1], "estimate")
  expect_identical(eval(formals(fit_pigauto)$lambda_mode)[1], "estimate")
  expect_identical(eval(formals(multi_impute)$lambda_mode)[1], "estimate")
})

# ---- (ii) impute(gnn = FALSE) model_config ------------------------------

.lambda_default_fixture <- function(seed = 50L, n = 30L) {
  set.seed(seed)
  tree <- ape::rcoal(n)
  df <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                    row.names = tree$tip.label)
  df$x1[c(2L, 7L, 15L)] <- NA
  df$x2[c(3L, 11L)] <- NA
  list(tree = tree, df = df)
}

test_that("[lambda-default] impute(gnn = FALSE) stores model_config$lambda_mode == 'estimate'", {
  fx <- .lambda_default_fixture()
  res <- suppressWarnings(impute(fx$df, fx$tree, gnn = FALSE,
                                  missing_frac = 0.2, verbose = FALSE,
                                  seed = 50L))

  expect_identical(res$fit$model_config$lambda_mode, "estimate")

  lpt <- res$fit$model_config$lambda_per_trait
  skip_if(is.null(lpt), "baseline has no lambda fields")
  expect_true(is.numeric(lpt))
  expect_length(lpt, ncol(res$data$X_scaled))
})

# ---- (iii) predict rebuild -----------------------------------------------

# gnn = FALSE blend (frozen contract,
# docs/dev-log/arc/2026-09-18-gnn-off-contract.md): pred = r_bm * mu +
# r_mean * mean_baseline_per_col. Named distinctly from test-gnn-off.R's
# own `.gnn_off_blend()` helper (both files are sourced into one
# environment by testthat::test_dir()) but computes the same thing.
.lambda_gnn_off_blend <- function(fit, mu) {
  latent_names <- names(fit$r_cal_bm)
  p <- length(latent_names)
  r_mean <- fit$r_cal_mean %||% stats::setNames(rep(0, p), latent_names)
  r_mean <- r_mean[latent_names]
  mean_bl <- fit$mean_baseline_per_col %||% stats::setNames(rep(0, p), latent_names)
  mean_bl <- mean_bl[latent_names]
  sweep(mu, 2, fit$r_cal_bm[latent_names], `*`) +
    matrix(r_mean * mean_bl, nrow = nrow(mu), ncol = p, byrow = TRUE)
}

test_that("[lambda-default] predict rebuild: predict() reproduces the fit-time baseline at estimated lambda", {
  fx  <- .lambda_default_fixture(seed = 60L, n = 60L)
  pd  <- preprocess_traits(fx$df, fx$tree)
  spl <- make_missing_splits(pd$X_scaled, missing_frac = 0.2, seed = 60L,
                              trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, fx$tree, splits = spl, gnn = FALSE,
                      lambda_mode = "estimate", verbose = FALSE, seed = 60L)

  skip_if(is.null(fit$model_config$lambda_per_trait), "baseline has no lambda fields")

  # A genuine predict() call: the previous version of this test only
  # replayed `lambda_fixed` through fit_baseline() and never called
  # predict() at all, so a mutant that estimates lambda correctly but
  # never applies it to a live prediction would have passed silently
  # (Rose review, 2026-09-23, "TESTS THAT CAN FAIL").
  pred <- predict(fit, return_se = FALSE)

  miss <- is.na(pd$X_scaled)
  expected <- .lambda_gnn_off_blend(fit, fit$baseline_full$mu)
  expect_equal(pred$imputed_latent[miss], expected[miss], tolerance = 1e-6)

  # S4 contract: fit_baseline(..., lambda_fixed = <named numeric>) rebuilds
  # the SAME baseline at previously-estimated per-trait lambda values
  # without re-running ML estimation for each column -- this is the
  # mechanism predict() uses to rebuild a baseline consistent with fit time
  # (e.g. after save/load, or for a fresh prediction call) rather than
  # re-estimating lambda from scratch.
  bl_rebuilt <- fit_baseline(pd, fx$tree,
                              lambda_fixed = fit$model_config$lambda_per_trait)
  expect_equal(bl_rebuilt$mu, fit$baseline_full$mu, tolerance = 1e-8)
})

# ---- (iv) per-tree lambda in multi_impute_trees --------------------------

test_that("[lambda-default] multi_impute_trees stores different lambda_per_trait per tree", {
  # Moderate, non-boundary phylogenetic signal under tree1 (lambda_true =
  # 0.5, mixed phylo/iid): a pure-noise or pure-BM trait pushes the ML
  # estimate to a boundary (0 or 1) under BOTH trees regardless of tree
  # identity, which defeats the point of this test (checked: iid traits
  # give 0.005/0.005 on both trees; lambda_true = 1 traits give 0.995/0.995
  # on both). A moderate mixed signal lands in the interior for tree1 and
  # is pulled toward the boundary for tree2 (whose OWN correlation matrix
  # is a shrunk version of tree1's), which is where per-tree lambda
  # actually differs.
  set.seed(70L)
  n <- 40L
  tree1 <- ape::rcoal(n)
  R1 <- stats::cov2cor(ape::vcv.phylo(tree1))
  R1 <- R1[tree1$tip.label, tree1$tip.label]
  lambda_true <- 0.5
  phylo1 <- as.numeric(t(chol(R1)) %*% rnorm(n))
  phylo2 <- as.numeric(t(chol(R1)) %*% rnorm(n))
  x1 <- sqrt(lambda_true) * phylo1 + sqrt(1 - lambda_true) * rnorm(n)
  x2 <- sqrt(lambda_true) * phylo2 + sqrt(1 - lambda_true) * rnorm(n)
  df <- data.frame(x1 = x1, x2 = x2, row.names = tree1$tip.label)
  df$x1[c(2L, 7L, 15L, 22L)] <- NA
  df$x2[c(3L, 11L, 28L)] <- NA

  # A second, phylogenetically DIFFERENT tree (correlation shrunk toward
  # I by transform_tree_pagel), same tip set -- so a lambda estimated per
  # tree should differ between the two.
  tree2 <- pigauto:::transform_tree_pagel(tree1, 0.2)
  trees <- list(tree1, tree2)
  class(trees) <- "multiPhylo"

  mi <- suppressWarnings(multi_impute_trees(
    df, trees, m_per_tree = 1L, gnn = FALSE,
    missing_frac = 0.2, verbose = FALSE, seed = 70L
  ))

  lpt <- mi$lambda_per_trait_by_tree
  skip_if(is.null(lpt) || any(vapply(lpt, is.null, logical(1))),
          "baseline has no lambda fields")
  expect_length(lpt, 2L)
  expect_false(isTRUE(all.equal(lpt[[1L]], lpt[[2L]])))
})

# ---- (v) estimated lambda is applied (Rose review, 2026-09-23) -----------

# Required change 8 from docs/dev-log/lambda-default/rose-final-review.md:
# a dispatcher-level test that fails if lambda is ESTIMATED but not
# actually APPLIED to the prediction. Section 4 of that review found no
# such test: every existing lambda-dispatch/lambda-default test either
# checked the reported `lambda_per_trait` value (not that it reaches mu)
# or replayed `lambda_fixed` against itself (trivially self-consistent
# even when lambda is silently 1 everywhere).
test_that("[lambda-default] estimated lambda is applied", {
  skip_if_not_installed("Matrix")
  set.seed(80L)
  n <- 200L
  tree <- ape::rcoal(n)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  Lc <- chol(R)
  # lambda_true = 0.3 is safely clear of phylo_signal_threshold's default
  # 0.2 boundary -- but phylo_signal_gate/safety_floor are still disabled
  # below (the "pure traditional-stats arm", per test-gnn-off.R) so this
  # test isolates lambda_mode's effect on the baseline from that separate
  # low-signal safety mechanism.
  lambda_true <- 0.3
  x1 <- sqrt(lambda_true) * as.numeric(t(Lc) %*% rnorm(n)) +
    sqrt(1 - lambda_true) * rnorm(n)
  x2 <- sqrt(lambda_true) * as.numeric(t(Lc) %*% rnorm(n)) +
    sqrt(1 - lambda_true) * rnorm(n)
  df <- data.frame(x1 = x1, x2 = x2, row.names = tree$tip.label)
  df$x1[sample(n, round(0.3 * n))] <- NA
  df$x2[sample(n, round(0.3 * n))] <- NA

  res_est <- suppressWarnings(impute(df, tree, gnn = FALSE,
                                      missing_frac = 0.2, verbose = FALSE,
                                      seed = 80L, lambda_mode = "estimate",
                                      safety_floor = FALSE,
                                      phylo_signal_gate = FALSE))
  res_f1  <- suppressWarnings(impute(df, tree, gnn = FALSE,
                                      missing_frac = 0.2, verbose = FALSE,
                                      seed = 80L, lambda_mode = "fixed_1",
                                      safety_floor = FALSE,
                                      phylo_signal_gate = FALSE))

  # (a) default ("estimate") and fixed_1 must give different completed
  # values on data simulated under lambda_true = 0.3 (a non-trivial
  # distance from the fixed_1 = 1 assumption).
  expect_false(isTRUE(all.equal(res_est$completed$x1, res_f1$completed$x1)))

  lpt <- res_est$fit$model_config$lambda_per_trait
  skip_if(is.null(lpt), "baseline has no lambda fields")

  # (b) the estimate-mode prediction itself must equal an INDEPENDENT
  # bm_impute_col() call at the fit's own reported lambda_hat -- not
  # merely a self-consistency check against another call that could share
  # the same bug (see the "predict rebuild" test above for that mutant
  # class).
  lambda_hat_x1 <- unname(lpt["x1"])
  spp <- rownames(res_est$data$X_scaled)
  ref <- bm_impute_col(res_est$data$X_scaled[spp, "x1"], R[spp, spp],
                        lambda = lambda_hat_x1)

  pred_x1 <- res_est$prediction$imputed_latent[spp, "x1"]
  expect_equal(unname(pred_x1), unname(ref$mu), tolerance = 1e-3)
})

# ---- (vi) NEWS mentions lambda early --------------------------------------

test_that("[lambda-default] NEWS.md mentions lambda in its first 60 lines", {
  news_path <- system.file("..", "NEWS.md", package = "pigauto")
  if (!nzchar(news_path) || !file.exists(news_path)) {
    # Not installed with NEWS.md (e.g. devtools::load_all from source):
    # fall back to the package source tree.
    news_path <- testthat::test_path("..", "..", "NEWS.md")
  }
  skip_if_not(file.exists(news_path), "NEWS.md not found")
  head_lines <- readLines(news_path, n = 60L, warn = FALSE)
  expect_true(any(grepl("lambda", head_lines, ignore.case = TRUE)))
})
