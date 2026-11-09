.analysis_mi_guard_fixture <- function() {
  t <- seq_len(16) / 16
  data.frame(x = replace(sin(3 * t), c(3, 8), NA_real_),
    y = cos(4 * t), z = t, a = t^2,
    group = factor(rep(c("a", "b"), each = 8)))
}

test_that("analysis MI rejects incomplete non-target roles before backend loading", {
  for (column in c("y", "z", "a", "group")) {
    data <- .analysis_mi_guard_fixture()
    data[[column]][2] <- NA
    formula <- if (column == "group") y ~ x + z + (1 | group) else y ~ x + z
    model <- if (column == "group") "lmer" else "lm"
    expect_error(multi_impute_analysis(data, formula, missing = "x",
      auxiliary = "a", model = model, m = 2L, seed = 17L),
      paste0("Only `x` may contain missing values; incomplete columns: ", column))
  }
  data <- .analysis_mi_guard_fixture()
  data$a[2] <- NaN
  expect_error(multi_impute_analysis(data, y ~ x + z, missing = "x",
    auxiliary = "a", m = 2L), "Only `x` may contain missing values")
  data$a[2] <- Inf
  expect_error(multi_impute_analysis(data, y ~ x + z, missing = "x",
    auxiliary = "a", m = 2L), "must contain only finite values")
})

test_that("analysis MI preserves complete inputs and fills its single target", {
  data <- .analysis_mi_guard_fixture()
  original <- data
  result <- multi_impute_analysis(data, y ~ x + z, missing = "x",
    auxiliary = "a", m = 2L, seed = 17L)
  expect_s3_class(result, "pigauto_analysis_mi")
  expect_length(result$datasets, 2L)
  expect_identical(data, original)
  observed <- !is.na(data$x)
  for (completed in result$datasets) {
    expect_true(all(is.finite(completed$x)))
    expect_identical(completed$x[observed], data$x[observed])
    expect_identical(completed[setdiff(names(data), "x")],
                     data[setdiff(names(data), "x")])
  }
})
