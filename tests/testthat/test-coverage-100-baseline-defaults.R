test_that("baseline helper defaults agree with public and explicit routes", {
  skip_if_not_installed("Matrix")
  tree <- ape::read.tree(text = paste0(
    "(((a:1,b:1):1,(c:1,d:1):1):1,",
    "((e:1,f:1):1,(g:1,h:1):1):1);"
  ))
  traits <- data.frame(x = c(1, NA, 4, 3, 7, 5, 8, 6),
                       y = c(2, 5, 3, 6, 4, 8, 7, 9),
                       row.names = tree$tip.label)
  data <- preprocess_traits(traits, tree, log_transform = FALSE)

  dispatch <- .fit_baseline_dispatch(data, tree, lambda_mode = "fixed_1")
  public <- fit_baseline(data, tree, lambda_mode = "fixed_1")
  default_core <- .fit_baseline_core(data, tree, lambda_mode = "fixed_1")
  explicit_core <- .fit_baseline_core(data, tree, lambda_mode = "fixed_1",
                                      predict_method = "exact")

  expect_identical(dispatch, public)
  expect_identical(default_core$predict_method_used, "exact")
  expect_equal(default_core$mu, explicit_core$mu)
  expect_equal(default_core$se, explicit_core$se)
  expect_true(all(is.finite(default_core$mu)))
  expect_true(all(is.finite(default_core$se)))
})

test_that("implicit exact fallback is quiet while explicit exact warns", {
  tree <- ape::read.tree(text = paste0(
    "(((a:1,b:1):1,(c:1,d:1):1):1,",
    "((e:1,f:1):1,(g:1,h:1):1):1);"
  ))
  traits <- data.frame(x = c(1, NA, 4, 3, 7, 5, 8, 6),
                       y = c(2, 5, 3, 6, 4, 8, 7, 9),
                       row.names = tree$tip.label)
  data <- preprocess_traits(traits, tree, log_transform = FALSE)
  real_require_namespace <- base::requireNamespace
  local_mocked_bindings(
    requireNamespace = function(package, ...) {
      if (package == "Matrix") return(FALSE)
      real_require_namespace(package, ...)
    },
    .package = "base"
  )

  expect_no_warning(default <- suppressMessages(
    .fit_baseline_core(data, tree, lambda_mode = "fixed_1")
  ))
  expect_warning(
    explicit <- .fit_baseline_core(data, tree, lambda_mode = "fixed_1",
                                   predict_method = "exact"),
    'predict_method = "exact" was not usable'
  )
  expect_equal(default$mu, explicit$mu)
  expect_equal(default$se, explicit$se)
  expect_true(all(is.finite(default$mu)))
  expect_true(all(is.finite(default$se)))
})

test_that("mixed baseline routes retain ordinal metadata and per-route results", {
  skip_if_not_installed("Matrix")
  withr::local_seed(28)
  tree <- ape::rtree(20)
  i <- seq_len(20)
  traits <- data.frame(
    grade = ordered(rep(c("low", "mid", "high", "mid"), 5),
                    levels = c("low", "mid", "high")),
    x = sin(i / 2), row.names = tree$tip.label
  )
  data <- preprocess_traits(traits, tree,
                            trait_types = c(grade = "ordinal", x = "continuous"),
                            log_transform = FALSE)
  splits <- list(val_idx = 22L, test_idx = 23L)
  mixed <- fit_baseline(data, tree, splits = splits,
                        predict_route = c(grade = "exact", x = "per_column"),
                        lambda_mode = "fixed_1", seed = 42)
  exact <- fit_baseline(data, tree, splits = splits, predict_method = "exact",
                        lambda_mode = "fixed_1")
  per_column <- fit_baseline(data, tree, splits = splits,
                             predict_method = "per_column", lambda_mode = "fixed_1")

  expect_identical(mixed$ordinal_path_chosen, c(`1` = "threshold_joint"))
  expect_identical(mixed$predict_method_by_trait,
                   c(grade = "exact", x = "per_column"))
  expect_equal(mixed$mu[, 1], exact$mu[, 1])
  expect_equal(mixed$se[, 1], exact$se[, 1])
  expect_equal(mixed$mu[, 2], per_column$mu[, 2])
  expect_equal(mixed$se[, 2], per_column$se[, 2])
  expect_true(all(is.finite(mixed$mu)))
  expect_true(all(is.finite(mixed$se)))
})
