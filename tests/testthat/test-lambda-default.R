# Tests for S5 of the joint-lambda-default lane (docs/dev-log/2026-09-22-joint-lambda-alignment.md):
# the lambda_mode = "estimate" default flip and its plumbing through
# impute() / fit_pigauto() / multi_impute() / multi_impute_trees().
#
# These tests deliberately favour gnn = FALSE where possible (see
# test-gnn-off.R's header comment): gnn = FALSE makes zero torch:: calls, so
# this file does not need skip_if_no_libtorch() at the top level.
#
# S4 (R/fit_baseline.R, R/joint_mvn_solver.R, R/joint_*_baseline.R,
# R/ovr_categorical.R) is what populates `lambda_per_trait` / `lambda_block`
# on fit_baseline()'s return value and adds `lambda_fixed` for predict-time
# rebuilds. Written before S4 landed in this shared worktree; verified
# passing (not skipping) once it did. Tests that depend on those fields
# still guard with skip_if(..., "S4 not landed") rather than assume landing
# order, so this file stays green regardless of merge order with S4.

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
  skip_if(is.null(lpt), "S4 not landed")
  expect_true(is.numeric(lpt))
  expect_length(lpt, ncol(res$data$X_scaled))
})

# ---- (iii) predict rebuild -----------------------------------------------

test_that("[lambda-default] predict rebuild: baseline at predict time matches fit time to 1e-8", {
  fx <- .lambda_default_fixture(seed = 60L)
  pd <- preprocess_traits(fx$df, fx$tree)
  bl <- fit_baseline(pd, fx$tree, lambda_mode = "estimate")

  skip_if(is.null(bl$lambda_per_trait), "S4 not landed")

  # S4 contract: fit_baseline(..., lambda_fixed = <named numeric>) rebuilds
  # the SAME baseline at previously-estimated per-trait lambda values
  # without re-running ML estimation for each column -- this is the
  # mechanism predict() uses to rebuild a baseline consistent with fit time
  # (e.g. after save/load, or for a fresh prediction call) rather than
  # re-estimating lambda from scratch.
  bl_rebuilt <- fit_baseline(pd, fx$tree, lambda_fixed = bl$lambda_per_trait)
  expect_equal(bl_rebuilt$mu, bl$mu, tolerance = 1e-8)
  expect_equal(bl_rebuilt$se, bl$se, tolerance = 1e-8)
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
          "S4 not landed")
  expect_length(lpt, 2L)
  expect_false(isTRUE(all.equal(lpt[[1L]], lpt[[2L]])))
})

# ---- (v) NEWS mentions lambda early --------------------------------------

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
