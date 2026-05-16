# Unit tests for script/gha/ helpers.
# These exercise the helpers via testthat so they're covered by the
# regular `devtools::test()` run, even though they live under script/.

.gha_path <- function(...) {
  # tests/testthat run from tests/testthat/; the helpers are at
  # ../../script/gha/<name>.R relative to the test file.
  file.path("..", "..", "script", "gha", ...)
}

test_that("[gha-config] PIGAUTO_CI_CONFIG loads with documented defaults", {
  e <- new.env()
  source(.gha_path("_ci_config.R"), local = e)
  cfg <- e$PIGAUTO_CI_CONFIG
  expect_type(cfg, "list")
  expect_equal(cfg$subset_n,          2000L)
  expect_equal(cfg$n_imputations,     10L)
  expect_equal(cfg$missing_frac,      0.30)
  expect_equal(cfg$seed,              2026L)
  expect_equal(cfg$pool_method,       "median")
  expect_false(cfg$clamp_outliers)
  expect_true(cfg$phylo_signal_gate)
  expect_equal(cfg$mc_cores,          4L)
})

test_that("[gha-config] .normalize_eval() projects to the canonical schema", {
  e <- new.env()
  source(.gha_path("_ci_config.R"), local = e)
  norm <- e$.normalize_eval

  raw <- data.frame(
    trait          = c("mass", "diet"),
    type           = c("continuous", "categorical"),
    imputation_idx = c(1L, 1L),
    rmse           = c(0.41, NA_real_),
    mae            = c(0.30, NA_real_),
    pearson_r      = c(0.85, NA_real_),
    accuracy       = c(NA_real_, 0.72),
    brier          = c(NA_real_, 0.18),
    time_sec       = c(120, 120),
    stringsAsFactors = FALSE
  )

  out <- norm(raw, dataset = "avonet", method = "pigauto_ci")
  expected_cols <- c("dataset", "trait", "type", "method",
                     "imputation_idx", "rmse", "mae", "pearson_r",
                     "accuracy", "brier", "time_sec")
  expect_equal(colnames(out), expected_cols)
  expect_equal(out$dataset, c("avonet", "avonet"))
  expect_equal(out$method,  c("pigauto_ci", "pigauto_ci"))
  expect_equal(nrow(out), 2L)
})
