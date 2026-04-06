# Tests for TabPFN baseline integration
# These tests are skipped when reticulate or tabimpute are unavailable.

tabpfn_available <- function() {
  requireNamespace("reticulate", quietly = TRUE) &&
    tryCatch({
      reticulate::import("tabimpute.interface")
      TRUE
    }, error = function(e) FALSE)
}

make_test_data <- function(n = 40, p = 2, seed = 42) {
  set.seed(seed)
  tree <- ape::rtree(n)
  sp   <- tree$tip.label
  df   <- data.frame(
    row.names = sp,
    tr1 = abs(stats::rnorm(n)) + 0.5,
    tr2 = abs(stats::rnorm(n)) + 0.5
  )
  list(tree = tree, df = df)
}

test_that("fit_baseline_tabpfn errors without reticulate", {
  td <- make_test_data()
  pd <- preprocess_traits(td$df, td$tree)
  # This test verifies the input validation path
  expect_error(
    fit_baseline_tabpfn(td$df),
    "pigauto_data"
  )
})

test_that("fit_baseline_tabpfn returns correct structure", {
  skip_if_not(tabpfn_available(), "tabimpute Python package not available")

  td  <- make_test_data()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 1)

  bl <- fit_baseline_tabpfn(pd, splits = spl, envname = NULL)

  expect_type(bl, "list")
  expect_named(bl, c("mu", "se"))
  expect_true(is.matrix(bl$mu))
  expect_equal(dim(bl$mu), c(40L, 2L))
  expect_equal(dim(bl$se), c(40L, 2L))
  expect_equal(rownames(bl$mu), pd$species_names)
  expect_equal(colnames(bl$mu), pd$trait_names)
})

test_that("fit_baseline_tabpfn predictions are finite", {
  skip_if_not(tabpfn_available(), "tabimpute Python package not available")

  td  <- make_test_data()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 2)

  bl <- fit_baseline_tabpfn(pd, splits = spl, envname = NULL)
  expect_true(all(is.finite(bl$mu)))
})

test_that("fit_baseline_tabpfn works with evaluate_imputation", {
  skip_if_not(tabpfn_available(), "tabimpute Python package not available")

  td  <- make_test_data()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 3)

  bl <- fit_baseline_tabpfn(pd, splits = spl, envname = NULL)
  ev <- evaluate_imputation(bl$mu, pd$X_scaled, spl)

  expect_s3_class(ev, "data.frame")
  expect_true(all(c("split", "trait", "n", "rmse", "pearson_r") %in% names(ev)))
  expect_true(all(ev$rmse > 0))
})

test_that("fit_baseline_tabpfn can serve as baseline for fit_pigauto", {
  skip_if_not(tabpfn_available(), "tabimpute Python package not available")

  td  <- make_test_data()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 4)

  bl_tab <- fit_baseline_tabpfn(pd, splits = spl, envname = NULL)
  graph  <- build_phylo_graph(td$tree)

  # TabPFN baseline fed into the GNN
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     graph = graph, baseline = bl_tab,
                     epochs = 20L, eval_every = 10L, patience = 5L,
                     verbose = FALSE, seed = 1)

  expect_s3_class(fit, "pigauto_fit")
  pred <- predict(fit, return_se = FALSE)
  expect_true(all(is.finite(pred)))
})

test_that("se_method='none' returns NA matrix", {
  skip_if_not(tabpfn_available(), "tabimpute Python package not available")

  td  <- make_test_data()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 5)

  bl <- fit_baseline_tabpfn(pd, splits = spl, se_method = "none",
                             envname = NULL)
  expect_true(all(is.na(bl$se)))
})
