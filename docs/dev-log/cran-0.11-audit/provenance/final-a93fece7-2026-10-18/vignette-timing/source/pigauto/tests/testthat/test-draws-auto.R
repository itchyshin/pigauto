# multi_impute(draws_method = "auto"): posterior when the data fit the
# posterior route's contract, otherwise conformal with a message.

da_tree_df <- function(n = 30L, factor_col = FALSE) {
  set.seed(31L)
  tree <- ape::rtree(n)
  df <- data.frame(mass = ape::rTraitCont(tree) + 10,
                   wing = ape::rTraitCont(tree) + 10,
                   row.names = tree$tip.label)
  if (factor_col) {
    df$diet <- factor(sample(c("a", "b"), n, replace = TRUE))
  }
  df$mass[c(2, 5, 9)] <- NA
  df$wing[c(4, 7)] <- NA
  list(tree = tree, df = df)
}

test_that("draws_method defaults to \"auto\"", {
  expect_identical(eval(formals(multi_impute)$draws_method)[1], "auto")
})

test_that("auto resolves to posterior for continuous-only, one-row-per-species data", {
  x <- da_tree_df()
  expect_message(
    out <- pigauto:::.mi_resolve_draws_auto(
      x$df, x$tree, species_col = NULL, trait_types = NULL,
      multi_proportion_groups = NULL, log_transform = TRUE,
      covariates = NULL, verbose = TRUE),
    "using \"posterior\"")
  expect_identical(out, "posterior")
  expect_silent(pigauto:::.mi_resolve_draws_auto(
    x$df, x$tree, species_col = NULL, trait_types = NULL,
    multi_proportion_groups = NULL, log_transform = TRUE,
    covariates = NULL, verbose = FALSE))
})

test_that("auto falls back to conformal, with a reason, when the posterior route cannot apply", {
  x <- da_tree_df(factor_col = TRUE)
  expect_message(
    out <- pigauto:::.mi_resolve_draws_auto(
      x$df, x$tree, species_col = NULL, trait_types = NULL,
      multi_proportion_groups = NULL, log_transform = TRUE,
      covariates = NULL, verbose = FALSE),
    "not continuous \\(diet: binary\\).*refuse them")
  expect_identical(out, "conformal")

  y <- da_tree_df()
  covs <- data.frame(t = rnorm(nrow(y$df)), row.names = rownames(y$df))
  expect_message(
    out <- pigauto:::.mi_resolve_draws_auto(
      y$df, y$tree, species_col = NULL, trait_types = NULL,
      multi_proportion_groups = NULL, log_transform = TRUE,
      covariates = covs, verbose = FALSE),
    "`covariates` were supplied")
  expect_identical(out, "conformal")

  expect_message(
    out <- pigauto:::.mi_resolve_draws_auto(
      y$df, y$tree, species_col = "sp", trait_types = NULL,
      multi_proportion_groups = NULL, log_transform = TRUE,
      covariates = NULL, verbose = FALSE),
    "`species_col` was supplied")
  expect_identical(out, "conformal")
})

test_that("default multi_impute() on continuous data gives poolable posterior draws", {
  skip_on_cran()
  x <- da_tree_df()
  ctl <- list(n_chains = 2L, burnin = 60L, n_iter = 120L, keep_draws = 40L,
              auto_extend = FALSE)
  mi <- suppressWarnings(suppressMessages(
    multi_impute(x$df, x$tree, m = 3L, posterior_control = ctl, seed = 1,
                 verbose = FALSE)))
  expect_identical(mi$draws_method, "posterior")
  expect_s3_class(mi, "pigauto_posterior_mi")
  fits <- with_imputations(mi, function(d) stats::lm(mass ~ wing, data = d))
  pooled <- pool_mi(fits)
  expect_true(all(is.finite(pooled$estimate)))
})

test_that("default multi_impute() on mixed-type data falls back to conformal", {
  skip_on_cran()
  x <- da_tree_df(factor_col = TRUE)
  expect_message(
    mi <- suppressWarnings(multi_impute(x$df, x$tree, m = 2L, gnn = FALSE,
                                        verbose = FALSE, seed = 1)),
    "using \"conformal\"")
  expect_identical(mi$draws_method, "conformal")
  expect_error(with_imputations(mi, function(d) stats::lm(mass ~ wing, data = d)))
})
