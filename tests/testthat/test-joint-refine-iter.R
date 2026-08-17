# Tests for `joint_refine_iter` (R/joint_mvn_solver.R::fit_joint_solver()'s
# `max_iter` cross-trait EM cell-refinement, threaded through fit_baseline() /
# fit_pigauto() / impute() the same way as `joint_solver`). See
# docs/dev-log/2026-08-17-sigma-recovery-results.md.

make_two_trait_bm_data <- function(seed = 5, n = 40) {
  set.seed(seed)
  tree <- ape::rcoal(n)
  C <- ape::vcv(tree)
  L <- chol(C)
  x1 <- as.numeric(t(L) %*% stats::rnorm(n))
  x2 <- 0.8 * x1 + as.numeric(t(L) %*% stats::rnorm(n)) * 0.3
  df <- data.frame(t1 = x1, t2 = x2)
  rownames(df) <- tree$tip.label
  df$t1[sample(n, 4)] <- NA
  df$t2[sample(n, 5)] <- NA

  pd <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.1,
                                val_frac = 0.5, seed = 1,
                                trait_map = pd$trait_map)
  graph <- build_phylo_graph(tree, k_eigen = 4L)
  list(tree = tree, pd = pd, splits = splits, graph = graph, df = df)
}

test_that("joint_refine_iter default (0L) is byte-identical whether passed explicitly or not", {
  skip_on_cran()
  d <- make_two_trait_bm_data(seed = 5)

  bl_default <- fit_baseline(d$pd, d$tree, splits = d$splits, graph = d$graph)
  bl_explicit <- fit_baseline(d$pd, d$tree, splits = d$splits, graph = d$graph,
                               joint_refine_iter = 0L)

  expect_identical(bl_default$mu, bl_explicit$mu)
  expect_identical(bl_default$se, bl_explicit$se)
})

test_that("joint_refine_iter = 3L runs, returns finite mu/se of identical shape, and differs from 0L", {
  skip_on_cran()
  d <- make_two_trait_bm_data(seed = 5)

  bl0 <- fit_baseline(d$pd, d$tree, splits = d$splits, graph = d$graph,
                       joint_refine_iter = 0L)
  bl3 <- fit_baseline(d$pd, d$tree, splits = d$splits, graph = d$graph,
                       joint_refine_iter = 3L)

  expect_identical(dim(bl3$mu), dim(bl0$mu))
  expect_identical(dim(bl3$se), dim(bl0$se))
  expect_true(all(is.finite(bl3$mu)))
  expect_true(all(is.finite(bl3$se)))

  # Cross-trait EM refinement actually moves the fit on correlated
  # multi-trait data -- not a coincidental byte-match with 0L.
  expect_gt(max(abs(bl3$mu - bl0$mu)), 1e-4)
})

test_that("joint_refine_iter validation rejects a negative value", {
  skip_on_cran()
  d <- make_two_trait_bm_data(seed = 5)

  expect_error(
    fit_baseline(d$pd, d$tree, splits = d$splits, graph = d$graph,
                 joint_refine_iter = -1L),
    "non-negative integer"
  )
})

test_that("joint_refine_iter validation rejects a non-numeric value", {
  skip_on_cran()
  d <- make_two_trait_bm_data(seed = 5)

  expect_error(
    fit_baseline(d$pd, d$tree, splits = d$splits, graph = d$graph,
                 joint_refine_iter = "3"),
    "non-negative integer"
  )
})

test_that("fit_pigauto() validates joint_refine_iter", {
  skip_on_cran()
  d <- make_two_trait_bm_data(seed = 5)

  expect_error(
    fit_pigauto(d$pd, d$tree, splits = d$splits, graph = d$graph,
                joint_refine_iter = -1L, epochs = 2L, verbose = FALSE),
    "non-negative integer"
  )
})

test_that("impute() validates joint_refine_iter", {
  skip_on_cran()
  d <- make_two_trait_bm_data(seed = 5)

  expect_error(
    impute(d$df, d$tree, joint_refine_iter = -1L, epochs = 2L, verbose = FALSE),
    "non-negative integer"
  )
})

test_that("impute() threads joint_refine_iter = 2L through to a working fit", {
  skip_if_no_libtorch()
  skip_on_cran()

  d <- make_two_trait_bm_data(seed = 5)

  result <- impute(d$df, d$tree, joint_refine_iter = 2L,
                   epochs = 5L, verbose = FALSE, seed = 5)

  expect_s3_class(result, "pigauto_result")
})
