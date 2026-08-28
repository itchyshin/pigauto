make_analysis_mi_provenance_fixture <- function(m = 2L) {
  set.seed(20260825)
  dat <- data.frame(
    y = stats::rnorm(30L),
    x = stats::rnorm(30L),
    z = stats::rnorm(30L)
  )
  dat$x[seq(3L, 30L, by = 5L)] <- NA_real_
  multi_impute_analysis(
    dat, y ~ x + z, missing = "x", model = "lm", m = m, seed = 17L
  )
}

make_mi_provenance_datasets <- function() {
  x <- seq_len(8L)
  list(
    data.frame(y = 1 + 0.4 * x + c(0.2, -0.1, 0.3, -0.2, 0.1, -0.3, 0.2, -0.1),
               x = x),
    data.frame(y = 1.1 + 0.4 * x + c(-0.1, 0.2, -0.2, 0.3, -0.3, 0.1, -0.1, 0.2),
               x = x + 0.05)
  )
}

test_that("analysis-aware MI carries durable fitting provenance", {
  mi <- make_analysis_mi_provenance_fixture()

  expect_s3_class(mi, "pigauto_analysis_mi")
  expect_identical(mi$mi_workflow, "pigauto_analysis_mi_v1")

  mi_path <- tempfile(fileext = ".rds")
  fits_path <- tempfile(fileext = ".rds")
  on.exit(unlink(c(mi_path, fits_path)), add = TRUE)
  saveRDS(mi, mi_path)
  restored_mi <- readRDS(mi_path)
  expect_s3_class(restored_mi, "pigauto_analysis_mi")
  expect_identical(restored_mi$mi_workflow, "pigauto_analysis_mi_v1")

  fits <- with_imputations(
    restored_mi, function(d) stats::lm(y ~ x + z, data = d), .progress = FALSE
  )
  expect_s3_class(fits, "pigauto_mi_fits")
  expect_identical(attr(fits, "mi_workflow"), "pigauto_analysis_mi_v1")
  saveRDS(fits, fits_path)
  restored_fits <- readRDS(fits_path)
  expect_identical(attr(restored_fits, "mi_workflow"), "pigauto_analysis_mi_v1")
  expect_s3_class(pool_mi(restored_fits), "pigauto_pooled")
})

test_that("with_imputations rejects mixed diagnostic classes before admission", {
  diagnostic <- make_analysis_mi_provenance_fixture()
  class(diagnostic) <- c("pigauto_diagnostic_mi", class(diagnostic))
  tree <- make_analysis_mi_provenance_fixture()
  class(tree) <- c("pigauto_tree_sensitivity_diagnostic", class(tree))
  called <- FALSE
  f <- function(d) {
    called <<- TRUE
    stats::lm(y ~ x + z, data = d)
  }

  expect_error(
    with_imputations(diagnostic, f, .progress = FALSE),
    "prediction-diagnostic"
  )
  expect_false(called)
  expect_error(
    with_imputations(tree, f, .progress = FALSE),
    "tree-sensitivity"
  )
  expect_false(called)
})

test_that("with_imputations refuses diagnostic, tree, and legacy MI inputs", {
  datasets <- make_mi_provenance_datasets()
  diagnostic <- structure(
    list(datasets = datasets, mi_workflow = "pigauto_diagnostic_mi"),
    class = c("pigauto_diagnostic_mi", "pigauto_mi", "list")
  )
  tree <- structure(
    list(datasets = datasets, mi_workflow = "pigauto_tree_sensitivity_diagnostic"),
    class = c("pigauto_tree_sensitivity_diagnostic", "pigauto_mi_trees",
              "pigauto_mi", "list")
  )
  legacy <- structure(list(datasets = datasets), class = c("pigauto_mi", "list"))

  expect_error(
    with_imputations(diagnostic, stats::lm, .progress = FALSE),
    "prediction-diagnostic.*multi_impute_analysis"
  )
  expect_error(
    with_imputations(tree, stats::lm, .progress = FALSE),
    "tree-sensitivity.*multi_impute_analysis"
  )
  expect_error(
    with_imputations(legacy, stats::lm, .progress = FALSE),
    "(?i)legacy.*multi_impute_analysis"
  )
  expect_error(
    with_imputations(datasets, stats::lm, .progress = FALSE),
    "multi_impute_analysis"
  )
})

test_that("pool_mi checks known workflow provenance before extraction", {
  datasets <- make_mi_provenance_datasets()
  fitted <- lapply(datasets, function(d) stats::lm(y ~ x, data = d))

  diagnostic <- structure(
    fitted, mi_workflow = "pigauto_diagnostic_mi",
    class = c("pigauto_diagnostic_mi", "list")
  )
  tree <- structure(
    fitted, mi_workflow = "pigauto_tree_sensitivity_diagnostic",
    class = c("pigauto_tree_sensitivity_diagnostic", "list")
  )
  legacy <- structure(fitted, class = c("pigauto_mi_fits", "list"))

  expect_error(pool_mi(diagnostic), "(?i)prediction-diagnostic")
  expect_error(pool_mi(tree), "(?i)tree-sensitivity")
  expect_error(pool_mi(legacy), "(?i)legacy.*multi_impute_analysis")

  warnings <- character()
  pooled <- withCallingHandlers(
    pool_mi(fitted),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(pooled, "pigauto_pooled")
  expect_length(warnings, 1L)
  expect_match(warnings, "unverified provenance")
})

test_that("pool_mi rejects forged and legacy provenance before extractors", {
  fitted <- lapply(make_mi_provenance_datasets(), function(d) stats::lm(y ~ x, data = d))
  called <- c(coef = FALSE, vcov = FALSE, df = FALSE)
  extractors <- list(
    coef_fun = function(f) {
      called[["coef"]] <<- TRUE
      stats::coef(f)
    },
    vcov_fun = function(f) {
      called[["vcov"]] <<- TRUE
      stats::vcov(f)
    },
    df_fun = function(f) {
      called[["df"]] <<- TRUE
      stats::df.residual(f)
    }
  )
  expect_invalid <- function(x, pattern) {
    called[] <<- FALSE
    expect_error(do.call(pool_mi, c(list(x), extractors)), pattern)
    expect_false(any(called))
  }

  forged <- structure(
    fitted, mi_workflow = "pigauto_analysis_mi_v1",
    class = "list"
  )
  legacy_mi <- structure(fitted, class = c("pigauto_mi", "list"))
  legacy_tree <- structure(fitted, class = c("pigauto_mi_trees", "pigauto_mi", "list"))
  unknown <- structure(fitted, mi_workflow = "unknown_workflow")

  expect_invalid(forged, "inconsistent provenance")
  expect_invalid(legacy_mi, "(?i)legacy")
  expect_invalid(legacy_tree, "(?i)legacy|tree")
  expect_invalid(unknown, "unknown")
})

test_that("bare lists warn only after successful pooling", {
  fitted <- lapply(make_mi_provenance_datasets(), function(d) stats::lm(y ~ x, data = d))
  warnings <- character()
  expect_error(
    withCallingHandlers(
      pool_mi(fitted, df_fun = function(f) stop("bad complete-data df")),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    "bad complete-data df"
  )
  expect_length(warnings, 0L)
})

test_that("analysis-aware fitting retains captured failures for pooling", {
  mi <- make_analysis_mi_provenance_fixture(m = 3L)
  calls <- 0L
  fit_warnings <- character()
  fits <- withCallingHandlers(
    with_imputations(
      mi,
      function(d) {
        calls <<- calls + 1L
        if (calls == 2L) stop("synthetic failure")
        stats::lm(y ~ x + z, data = d)
      },
      .progress = FALSE,
      .on_error = "continue"
    ),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(fit_warnings, 1L)
  expect_match(fit_warnings, "1 of 3 fits failed")
  expect_identical(attr(fits, "mi_workflow"), "pigauto_analysis_mi_v1")
  expect_equal(attr(fits, "failed"), 2L)
  expect_s3_class(fits[[2]], "pigauto_mi_error")
  pool_warnings <- character()
  pooled <- withCallingHandlers(
    pool_mi(fits),
    warning = function(w) {
      pool_warnings <<- c(pool_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(pool_warnings, 1L)
  expect_match(pool_warnings, "Dropping 1 fit")
  expect_s3_class(pooled, "pigauto_pooled")
})

test_that("diagnostic and tree provenance survive RDS refusal", {
  datasets <- make_mi_provenance_datasets()
  objects <- list(
    diagnostic = structure(
      list(datasets = datasets, mi_workflow = "pigauto_diagnostic_mi"),
      class = c("pigauto_diagnostic_mi", "pigauto_mi", "list")
    ),
    tree = structure(
      list(datasets = datasets, mi_workflow = "pigauto_tree_sensitivity_diagnostic"),
      class = c("pigauto_tree_sensitivity_diagnostic", "pigauto_mi_trees",
                "pigauto_mi", "list")
    )
  )
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  for (name in names(objects)) {
    saveRDS(objects[[name]], path)
    restored <- readRDS(path)
    expect_identical(restored$mi_workflow, objects[[name]]$mi_workflow)
    expect_s3_class(restored, class(objects[[name]])[[1]])
    expect_error(with_imputations(restored, stats::lm, .progress = FALSE),
                 if (name == "diagnostic") "prediction-diagnostic" else "tree-sensitivity")
  }
})
