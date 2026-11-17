test_that("constant ordinal correlations are NA without warnings", {
  truth <- matrix(c(1, 2, 3, 4), ncol = 1)
  constant <- matrix(rep(1, 4), ncol = 1)
  map <- list(list(name = "trait", type = "ordinal", latent_cols = 1L))
  split <- list(val_idx = 1:4, test_idx = integer())

  for (reverse in c(FALSE, TRUE)) {
    observed <- if (reverse) constant else truth
    predicted <- if (reverse) truth else constant
    expect_no_warning(out <- evaluate_imputation(predicted, observed, split,
                                                 trait_map = map))
    expect_true(is.na(out$spearman_rho))
    expect_no_warning(scored <- pigauto:::eval_test_cells(
      predicted, NULL, observed, 1:4, rep(1L, 4), map, 4L, "pigauto"
    ))
    expect_true(is.na(scored$value[scored$metric == "spearman_rho"]))
    expect_equal(scored$value[scored$metric == "rmse"], sqrt(3.5))
  }
})

test_that("ordinal scores retain defined rank correlations", {
  truth <- matrix(c(1, 2, 3, 4), ncol = 1)
  predicted <- matrix(c(4, 3, 1, 2), ncol = 1)
  map <- list(list(name = "trait", type = "ordinal", latent_cols = 1L))
  out <- evaluate_imputation(predicted, truth,
                             list(val_idx = 1:4, test_idx = integer()),
                             trait_map = map)
  expected <- cor(as.vector(truth), as.vector(predicted), method = "spearman")
  expect_equal(out$spearman_rho, expected)
  scored <- pigauto:::eval_test_cells(
    predicted, NULL, truth, 1:4, rep(1L, 4), map, 4L, "pigauto"
  )
  expect_equal(scored$value[scored$metric == "spearman_rho"], expected)
})

test_that("ordinal constant detection uses finite pairs and preserves ties", {
  map <- list(list(name = "trait", type = "ordinal", latent_cols = 1L))
  truth <- matrix(c(1, 2, 3, Inf, NA), ncol = 1)
  predicted <- matrix(c(2, 2, 2, 3, 4), ncol = 1)
  split <- list(val_idx = 1:5, test_idx = integer())
  expect_no_warning(out <- evaluate_imputation(predicted, truth, split,
                                               trait_map = map))
  expect_true(is.na(out$spearman_rho))
  expect_equal(out$n, 3L)

  truth <- matrix(c(1, 1, 2, 3), ncol = 1)
  predicted <- matrix(c(3, 2, 2, 1), ncol = 1)
  out <- evaluate_imputation(predicted, truth,
                             list(val_idx = 1:4, test_idx = integer()),
                             trait_map = map)
  expect_equal(out$spearman_rho, cor(as.vector(truth), as.vector(predicted),
                                     method = "spearman"))
})

test_that("constant numeric correlations are NA without warnings", {
  truth <- matrix(c(1, 2, 3, 4), ncol = 1,
                  dimnames = list(NULL, "trait"))
  constant <- truth * 0 + 1
  split <- list(val_idx = 1:4, test_idx = integer())
  for (type in c("continuous", "count", "proportion")) {
    map <- list(list(name = "trait", type = type, latent_cols = 1L))
    for (reverse in c(FALSE, TRUE)) {
      observed <- if (reverse) constant else truth
      predicted <- if (reverse) truth else constant
      expect_no_warning(out <- evaluate_imputation(predicted, observed, split,
                                                   trait_map = map))
      expect_true(is.na(out$pearson_r))
      expect_no_warning(scored <- pigauto:::eval_test_cells(
        predicted, NULL, observed, 1:4, rep(1L, 4), map, 4L, "pigauto"
      ))
      expect_true(is.na(scored$value[scored$metric == "pearson_r"]))
    }
  }
  expect_no_warning(out <- evaluate_imputation(constant, truth, split))
  expect_true(is.na(out$pearson_r))
  expect_no_warning(scored <- pigauto:::eval_test_cells(
    constant, NULL, truth, 1:4, rep(1L, 4), NULL, 4L, "pigauto"
  ))
  expect_true(is.na(scored$value[scored$metric == "pearson_r"]))
})

test_that("zero-inflated count correlations handle constant predictions", {
  truth <- cbind(gate = rep(1, 4), magnitude = log1p(1:4))
  predicted <- cbind(gate = rep(0, 4), magnitude = rep(log(2), 4))
  map <- list(list(name = "trait", type = "zi_count", latent_cols = 1:2,
                   mean = 0, sd = 1))
  split <- list(val_idx = 1:4, test_idx = integer())
  expect_no_warning(out <- evaluate_imputation(predicted, truth, split,
                                               trait_map = map))
  expect_true(is.na(out$pearson_r))
  expect_no_warning(scored <- pigauto:::eval_test_cells(
    predicted, NULL, truth, 1:4, rep(1L, 4), map, 4L, "pigauto"
  ))
  expect_true(is.na(scored$value[scored$metric == "pearson_r"]))
})

test_that("nonconstant numeric correlations remain unchanged", {
  truth <- matrix(c(1, 2, 3, 4), ncol = 1)
  predicted <- matrix(c(1, 3, 2, 5), ncol = 1)
  map <- list(list(name = "trait", type = "continuous", latent_cols = 1L))
  out <- evaluate_imputation(predicted, truth,
                             list(val_idx = 1:4, test_idx = integer()),
                             trait_map = map)
  expect_equal(out$pearson_r, cor(as.vector(truth), as.vector(predicted)))
})

test_that("nonfinite scoring correlations preserve base R semantics", {
  x <- c(Inf, Inf, Inf, Inf)
  y <- c(1, 2, 3, 4)
  expect_identical(suppressWarnings(pigauto:::cor_for_scoring(x, y)),
                   suppressWarnings(stats::cor(x, y)))
  x <- c(1, NA, 2, 3)
  expect_identical(pigauto:::cor_for_scoring(x, y), stats::cor(x, y))
})

test_that("zero-inflated constant truth also has undefined correlation", {
  truth <- cbind(gate = rep(1, 4), magnitude = rep(log(2), 4))
  predicted <- cbind(gate = rep(0, 4), magnitude = log1p(1:4))
  map <- list(list(name = "trait", type = "zi_count", latent_cols = 1:2,
                   mean = 0, sd = 1))
  split <- list(val_idx = 1:4, test_idx = integer())
  expect_no_warning(out <- evaluate_imputation(predicted, truth, split,
                                               trait_map = map))
  expect_true(is.na(out$pearson_r))
  expect_no_warning(scored <- pigauto:::eval_test_cells(
    predicted, NULL, truth, 1:4, rep(1L, 4), map, 4L, "pigauto"
  ))
  expect_true(is.na(scored$value[scored$metric == "pearson_r"]))
})
