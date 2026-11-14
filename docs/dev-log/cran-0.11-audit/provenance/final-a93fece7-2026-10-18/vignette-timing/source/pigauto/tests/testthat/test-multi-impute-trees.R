skip_if_no_libtorch()

# P1-11: multi_impute_trees() now accepts draws_method = c("mc_dropout",
# "conformal") (default "mc_dropout", back compat). Tests use tiny trees,
# epochs = 5L, m_per_tree = 2L, T = 2 to keep runtime low.

# ---- Helpers ----------------------------------------------------------------

make_trees_test_data <- function(n = 40L, miss_frac = 0.15, seed = 321L) {
  set.seed(seed)
  tree1 <- ape::rtree(n)
  tree2 <- ape::rtree(n)
  sp    <- tree1$tip.label
  df <- data.frame(
    row.names = sp,
    tr1 = abs(stats::rnorm(n)) + 0.5,
    tr2 = abs(stats::rnorm(n)) + 0.5
  )
  n_cells   <- n * 2L
  n_missing <- max(1L, round(n_cells * miss_frac))
  idx <- sample.int(n_cells, n_missing)
  m   <- as.matrix(df)
  m[idx] <- NA
  df <- as.data.frame(m)
  rownames(df) <- sp

  trees <- list(tree1, tree2)
  class(trees) <- "multiPhylo"
  list(trees = trees, df = df)
}

quick_mit <- function(draws_method = "mc_dropout", share_gnn = TRUE,
                      m_per_tree = 2L, seed = 5L) {
  td <- make_trees_test_data(seed = seed)
  mi <- multi_impute_trees(
    traits       = td$df,
    trees        = td$trees,
    m_per_tree   = m_per_tree,
    draws_method = draws_method,
    share_gnn    = share_gnn,
    epochs       = 5L,
    missing_frac = 0.25,
    verbose      = FALSE,
    seed         = seed
  )
  list(mi = mi, td = td)
}

# ---- 1. Default draws_method is "mc_dropout" (back compat) -----------------

test_that("multi_impute_trees defaults to draws_method = 'mc_dropout'", {
  setup <- quick_mit(seed = 11L)
  mi <- setup$mi

  expect_s3_class(mi, "pigauto_mi_trees")
  expect_identical(mi$draws_method, "mc_dropout")
  expect_equal(mi$m, 4L)
  expect_equal(mi$n_trees, 2L)
  expect_equal(mi$m_per_tree, 2L)
  expect_equal(length(mi$datasets), 4L)
  for (d in mi$datasets) {
    expect_equal(dim(d), dim(setup$td$df))
  }
})

# ---- 2. draws_method = "conformal" (share_gnn = TRUE, the default path) ----

test_that("draws_method = 'conformal' returns m datasets per tree, correct shape, non-identical draws (share_gnn = TRUE)", {
  setup <- quick_mit(draws_method = "conformal", share_gnn = TRUE, seed = 12L)
  mi <- setup$mi
  td <- setup$td

  expect_s3_class(mi, "pigauto_mi_trees")
  expect_identical(mi$draws_method, "conformal")
  expect_equal(mi$m, 4L)
  expect_equal(length(mi$datasets), 4L)
  for (d in mi$datasets) {
    expect_equal(dim(d), dim(td$df))
    expect_named(d, names(td$df))
  }
  expect_equal(mi$tree_index, c(1L, 1L, 2L, 2L))

  # Non-identical draws within the same tree (k=1 vs k=2 for tree 1).
  d1 <- mi$datasets[[1]]
  d2 <- mi$datasets[[2]]
  imputed_rows <- which(is.na(td$df$tr1))
  expect_gt(length(imputed_rows), 0L)
  expect_false(isTRUE(all.equal(
    d1$tr1[imputed_rows], d2$tr1[imputed_rows]
  )))

  # Observed cells preserved identically across all draws.
  obs_rows <- which(!is.na(td$df$tr1))
  expect_equal(d1$tr1[obs_rows], td$df$tr1[obs_rows])
  expect_equal(d2$tr1[obs_rows], td$df$tr1[obs_rows])
})

# ---- 3. draws_method = "conformal" (share_gnn = FALSE, per-tree impute()) --

test_that("draws_method = 'conformal' returns m datasets per tree, correct shape, non-identical draws (share_gnn = FALSE)", {
  setup <- quick_mit(draws_method = "conformal", share_gnn = FALSE, seed = 13L)
  mi <- setup$mi
  td <- setup$td

  expect_s3_class(mi, "pigauto_mi_trees")
  expect_identical(mi$draws_method, "conformal")
  expect_equal(mi$m, 4L)
  expect_equal(length(mi$datasets), 4L)
  for (d in mi$datasets) {
    expect_equal(dim(d), dim(td$df))
  }

  d1 <- mi$datasets[[1]]
  d2 <- mi$datasets[[2]]
  imputed_rows <- which(is.na(td$df$tr1))
  expect_gt(length(imputed_rows), 0L)
  expect_false(isTRUE(all.equal(
    d1$tr1[imputed_rows], d2$tr1[imputed_rows]
  )))
})

# ---- 4. print() reports the correct draws_method label ----------------------

test_that("print.pigauto_mi_trees reports the draws_method actually used", {
  setup <- quick_mit(draws_method = "conformal", share_gnn = TRUE, seed = 14L)
  out <- capture.output(print(setup$mi))
  expect_true(any(grepl("conformal-width sampling", out, fixed = TRUE)))
})
