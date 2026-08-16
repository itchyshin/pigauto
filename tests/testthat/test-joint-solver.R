# Tests for the joint_solver = c("inhouse", "rphylopars") dispatcher
# (R/joint_mvn_solver.R::fit_joint_solver()) and its threading through
# fit_baseline() / fit_pigauto() / impute(). Diagnosis:
# docs/dev-log/2026-08-16-continuous-gap-diagnosis.md

make_two_trait_bm_data <- function(seed = 42, n = 30) {
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
  list(tree = tree, pd = pd, splits = splits, graph = graph)
}

test_that("joint_solver default ('inhouse') is unchanged whether passed explicitly or not", {
  skip_on_cran()
  d <- make_two_trait_bm_data(seed = 42)

  bl_default <- fit_baseline(d$pd, d$tree, splits = d$splits, graph = d$graph)
  bl_explicit <- fit_baseline(d$pd, d$tree, splits = d$splits, graph = d$graph,
                               joint_solver = "inhouse")

  expect_identical(bl_default$mu, bl_explicit$mu)
  expect_identical(bl_default$se, bl_explicit$se)
})

test_that("joint_solver = 'rphylopars' runs and returns finite mu/se of the right shape", {
  skip_on_cran()
  skip_if_not_installed("Rphylopars")
  d <- make_two_trait_bm_data(seed = 42)

  bl_rp <- fit_baseline(d$pd, d$tree, splits = d$splits, graph = d$graph,
                         joint_solver = "rphylopars")
  bl_ih <- fit_baseline(d$pd, d$tree, splits = d$splits, graph = d$graph,
                         joint_solver = "inhouse")

  expect_equal(dim(bl_rp$mu), dim(bl_ih$mu))
  expect_true(all(is.finite(bl_rp$mu)))
  expect_true(all(is.finite(bl_rp$se)))

  # Different estimators (single-pass in-house vs converged REML) -- they
  # should differ measurably on this data, not coincide to machine precision.
  expect_gt(max(abs(bl_rp$mu - bl_ih$mu)), 1e-4)
})

test_that("joint_solver = 'rphylopars' output shape/colnames match the in-house contract", {
  skip_on_cran()
  skip_if_not_installed("Rphylopars")
  d <- make_two_trait_bm_data(seed = 43)

  bl_rp <- fit_baseline(d$pd, d$tree, splits = d$splits, graph = d$graph,
                         joint_solver = "rphylopars")
  bl_ih <- fit_baseline(d$pd, d$tree, splits = d$splits, graph = d$graph,
                         joint_solver = "inhouse")

  expect_identical(dim(bl_rp$mu), dim(bl_ih$mu))
  expect_identical(dim(bl_rp$se), dim(bl_ih$se))
  expect_identical(dimnames(bl_rp$mu), dimnames(bl_ih$mu))
  expect_identical(dimnames(bl_rp$se), dimnames(bl_ih$se))
})

test_that("joint_solver = 'rphylopars' falls back to in-house with a warning on failure", {
  skip_on_cran()
  skip_if_not_installed("Rphylopars")
  testthat::local_mocked_bindings(
    phylopars = function(...) stop("boom"),
    .package = "Rphylopars"
  )
  d <- make_two_trait_bm_data(seed = 44)

  bl_ih <- fit_baseline(d$pd, d$tree, splits = d$splits, graph = d$graph,
                         joint_solver = "inhouse")
  expect_warning(
    bl_fallback <- fit_baseline(d$pd, d$tree, splits = d$splits, graph = d$graph,
                                 joint_solver = "rphylopars"),
    "falling back to the in-house solver"
  )
  expect_identical(bl_fallback$mu, bl_ih$mu)
  expect_identical(bl_fallback$se, bl_ih$se)
})

test_that("impute() threads joint_solver = 'rphylopars' through to a working fit", {
  skip_if_no_libtorch()
  skip_on_cran()
  skip_if_not_installed("Rphylopars")

  set.seed(45)
  n <- 20L
  tree <- ape::rcoal(n)
  C <- ape::vcv(tree)
  L <- chol(C)
  x1 <- as.numeric(t(L) %*% stats::rnorm(n))
  x2 <- 0.8 * x1 + as.numeric(t(L) %*% stats::rnorm(n)) * 0.3
  df <- data.frame(t1 = x1, t2 = x2)
  rownames(df) <- tree$tip.label

  result <- impute(df, tree, joint_solver = "rphylopars",
                   epochs = 5L, verbose = FALSE, seed = 45)

  expect_s3_class(result, "pigauto_result")
})
