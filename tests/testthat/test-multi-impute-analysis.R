make_analysis_mi_data <- function(n = 180L, seed = 1201L) {
  set.seed(seed)
  z <- stats::rnorm(n)
  x <- 0.4 * z + stats::rnorm(n, sd = 0.8)
  y <- 0.5 + 0.7 * x - 0.35 * z + stats::rnorm(n, sd = 0.9)
  data <- data.frame(y = y, x = x, z = z)
  data$x[seq(3L, n, by = 4L)] <- NA_real_
  data
}

test_that("Bayesian Normal analysis-aware MI preserves data and varies draws", {
  data <- make_analysis_mi_data()
  mi <- multi_impute_analysis(
    data, y ~ x + z, missing = "x", model = "lm", m = 8L, seed = 44L
  )

  expect_s3_class(mi, "pigauto_analysis_mi")
  expect_s3_class(mi, "pigauto_mi")
  expect_identical(mi$engine, "bayes_norm")
  expect_named(mi$dependency_versions, "R")
  expect_true(nzchar(mi$package_version))
  expect_identical(mi$m, 8L)
  expect_length(mi$datasets, 8L)
  expect_identical(mi$control, list())

  missing <- is.na(data$x)
  for (completed in mi$datasets) {
    expect_identical(dim(completed), dim(data))
    expect_identical(names(completed), names(data))
    expect_identical(rownames(completed), rownames(data))
    expect_identical(completed$x[!missing], data$x[!missing])
    expect_identical(completed[c("y", "z")], data[c("y", "z")])
    expect_true(all(is.finite(completed$x)))
  }
  draws <- vapply(mi$datasets, function(x) x$x[which(missing)[1]], numeric(1))
  expect_gt(stats::sd(draws), 0)
})

test_that("Bayesian Normal analysis-aware MI is seed reproducible", {
  data <- make_analysis_mi_data(seed = 1202L)
  a <- multi_impute_analysis(
    data, y ~ x + z, "x", model = "lm", m = 3L, seed = 9L
  )
  b <- multi_impute_analysis(
    data, y ~ x + z, "x", model = "lm", m = 3L, seed = 9L
  )
  c <- multi_impute_analysis(
    data, y ~ x + z, "x", model = "lm", m = 3L, seed = 10L
  )
  expect_identical(a$datasets, b$datasets)
  expect_false(identical(a$datasets, c$datasets))
})

test_that("analysis-aware MI composes with fit and Rubin pooling", {
  data <- make_analysis_mi_data(n = 400L, seed = 1203L)
  mi <- multi_impute_analysis(
    data, y ~ x + z, "x", model = "lm", m = 50L, seed = 81L
  )
  fits <- with_imputations(
    mi, function(d) stats::lm(y ~ x + z, data = d), .progress = FALSE
  )
  pooled <- pool_mi(fits)
  x_row <- pooled[pooled$term == "x", , drop = FALSE]

  expect_s3_class(fits, "pigauto_mi_fits")
  expect_s3_class(pooled, "pigauto_pooled")
  expect_equal(x_row$estimate, 0.7, tolerance = 0.15)
  expect_lt(x_row$conf.low, 0.7)
  expect_gt(x_row$conf.high, 0.7)
})

test_that("analysis-aware MI prints its narrow provenance", {
  data <- make_analysis_mi_data(n = 80L, seed = 1204L)
  mi <- multi_impute_analysis(
    data, y ~ x + z, "x", model = "lm", m = 2L, seed = 2L
  )
  out <- capture.output(print(mi))
  expect_match(paste(out, collapse = "\n"), "experimental analysis-aware")
  expect_match(paste(out, collapse = "\n"), "bayes_norm")
  expect_match(paste(out, collapse = "\n"), "y ~ x \\+ z")
})

test_that("analysis-aware MI enforces the narrow public contract", {
  data <- make_analysis_mi_data(n = 80L, seed = 1205L)

  expect_error(
    multi_impute_analysis(data, y ~ log(x) + z, "x", model = "lm", m = 2L),
    "additive, untransformed"
  )
  expect_error(
    multi_impute_analysis(data, y ~ x * z, "x", model = "lm", m = 2L),
    "additive, untransformed"
  )
  expect_error(
    multi_impute_analysis(data, y ~ z, "x", model = "lm", m = 2L),
    "main-effect predictor"
  )
  expect_error(
    multi_impute_analysis(data, y ~ 0 + x + z, "x", model = "lm", m = 2L),
    "with an intercept"
  )
  integer_data <- data
  integer_data$x <- as.integer(round(integer_data$x))
  expect_error(
    multi_impute_analysis(
      integer_data, y ~ x + z, "x", model = "lm", m = 2L
    ),
    "double-precision"
  )
  data$z[1] <- NA_real_
  expect_error(
    multi_impute_analysis(data, y ~ x + z, "x", model = "lm", m = 2L),
    "Only `x` may contain"
  )
  data$z[1] <- Inf
  expect_error(
    multi_impute_analysis(data, y ~ x + z, "x", model = "lm", m = 2L),
    "only finite"
  )
})

test_that("analysis-aware MI validates controls without silent truncation", {
  data <- make_analysis_mi_data(n = 80L, seed = 1206L)
  data$binary <- as.integer(data$y > stats::median(data$y))

  expect_identical(
    pigauto:::.analysis_mi_control(list(), "glm"),
    list(numit = 20L, rjlimit = 100000L)
  )
  expect_identical(
    pigauto:::.analysis_mi_control(list(), "lmer"),
    list(nburn = 1000L, nbetween = 100L)
  )

  expect_error(
    multi_impute_analysis(
      data, binary ~ x + z, "x", model = "glm", m = 2L,
      control = list(numit = 2.5)
    ),
    "whole number"
  )
  expect_error(
    multi_impute_analysis(
      data, binary ~ x + z, "x", model = "glm", m = 2L,
      control = list(unknown = 1L)
    ),
    "Unsupported `control`"
  )
  expect_error(
    multi_impute_analysis(data, y ~ x + z, "x", model = "lm", m = 2.5),
    "whole number"
  )
})

test_that("Bayesian Normal route rejects a rank-deficient imputation model", {
  data <- make_analysis_mi_data(n = 80L, seed = 1207L)
  data$duplicate_y <- data$y
  expect_error(
    multi_impute_analysis(
      data, y ~ x + duplicate_y, "x", model = "lm", m = 2L
    ),
    "rank deficient"
  )
})

test_that("Bayesian Normal predictive moments match the analytic posterior", {
  data <- make_analysis_mi_data(n = 120L, seed = 1211L)
  missing_row <- which(is.na(data$x))[1]
  datasets <- pigauto:::.analysis_mi_bayes_norm(
    data, y ~ x + z, "x", character(), m = 5000L, seed = 91L
  )

  observed <- !is.na(data$x)
  X <- stats::model.matrix(~ y + z, data = data)
  Xo <- X[observed, , drop = FALSE]
  fit <- stats::lm.fit(Xo, data$x[observed])
  df <- sum(observed) - ncol(Xo)
  sse <- sum(fit$residuals^2)
  xm <- X[missing_row, , drop = FALSE]
  expected_mean <- as.numeric(xm %*% fit$coefficients)
  leverage <- as.numeric(xm %*% solve(crossprod(Xo), t(xm)))
  expected_variance <- sse / (df - 2) * (1 + leverage)
  draws <- vapply(datasets, function(x) x$x[missing_row], numeric(1))

  expect_equal(mean(draws), expected_mean, tolerance = 0.04)
  expect_equal(stats::var(draws), expected_variance,
               tolerance = 0.08 * expected_variance)
})

test_that("Bayesian Normal route rejects ill-conditioned designs", {
  data <- make_analysis_mi_data(n = 80L, seed = 1212L)
  data$almost_y <- data$y + seq_along(data$y) * 1e-14
  expect_error(
    multi_impute_analysis(
      data, y ~ x + almost_y, "x", model = "lm", m = 2L
    ),
    "ill-conditioned|rank deficient"
  )
})

test_that("engine warnings are captured for provenance", {
  captured <- pigauto:::.analysis_mi_capture_warnings({
    warning("engine warning", call. = FALSE)
    7L
  })
  expect_identical(captured$value, 7L)
  expect_identical(captured$warnings, "engine warning")
  expect_identical(captured$messages, character())
  expect_identical(captured$output, character())
  expect_error(
    pigauto:::.analysis_mi_abort_on_warnings(
      captured$warnings, captured$messages, "glm"
    ),
    "no inference-ready object was returned"
  )
})

test_that("analysis-aware MI rejects unsupported types and auxiliary overlap", {
  data <- make_analysis_mi_data(n = 80L, seed = 1213L)
  data$character_z <- as.character(round(data$z, 2))
  expect_error(
    multi_impute_analysis(
      data, y ~ x + character_z, "x", model = "lm", m = 2L
    ),
    "must be numeric"
  )
  expect_error(
    multi_impute_analysis(
      data, y ~ x + z, "x", model = "lm", m = 2L, auxiliary = "z"
    ),
    "must not duplicate formula variables"
  )
})

test_that("glm route has an informative optional-package guard", {
  skip_if(requireNamespace("smcfcs", quietly = TRUE),
          "smcfcs is installed; integration test covers this route")
  data <- make_analysis_mi_data(n = 80L, seed = 1208L)
  data$binary <- as.integer(data$y > stats::median(data$y))
  expect_error(
    multi_impute_analysis(
      data, binary ~ x + z, "x", model = "glm", m = 2L
    ),
    "requires optional package 'smcfcs'"
  )
})

test_that("smcfcs route returns valid completed datasets", {
  skip_if_not_installed("smcfcs")
  data <- make_analysis_mi_data(n = 100L, seed = 1209L)
  data$binary <- as.integer(data$y > stats::median(data$y))
  mi <- multi_impute_analysis(
    data, binary ~ x + z, "x", model = "glm", m = 2L, seed = 5L,
    control = list(numit = 2L, rjlimit = 10000L)
  )
  expect_identical(mi$engine, "smcfcs")
  expect_length(mi$datasets, 2L)
  expect_true(all(vapply(mi$datasets, function(x) all(is.finite(x$x)),
                         logical(1))))
})

test_that("lmer route enforces one random intercept and returns valid data", {
  skip_if_not_installed("lme4")
  data <- make_analysis_mi_data(n = 80L, seed = 1210L)
  data$group <- factor(rep(seq_len(20L), each = 4L))

  expect_error(
    multi_impute_analysis(
      data, y ~ x + z + (x | group), "x", model = "lmer", m = 2L
    ),
    "requires exactly one"
  )

  if (!requireNamespace("jomo", quietly = TRUE)) {
    expect_error(
      multi_impute_analysis(
        data, y ~ x + z + (1 | group), "x", model = "lmer", m = 2L
      ),
      "requires optional package 'jomo'"
    )
    return(invisible())
  }

  mi <- multi_impute_analysis(
    data, y ~ x + z + (1 | group), "x", model = "lmer", m = 2L, seed = 6L,
    control = list(nburn = 20L, nbetween = 10L)
  )
  expect_identical(mi$engine, "jomo")
  expect_length(mi$datasets, 2L)
  expect_true(all(vapply(mi$datasets, function(x) all(is.finite(x$x)),
                         logical(1))))
  expect_true(all(vapply(mi$datasets, function(x) identical(x$group, data$group),
                         logical(1))))
})
