# Silent-fallback honesty warnings (arc/silent-fallback-honesty).
#
# 1. P1-8: user covariates never reach the joint MVN / threshold-joint
#    baselines -- warn instead of silently ignoring them.
# 2. Small-validation conformal ceiling: at n_val < 19 the 95% split-
#    conformal target is arithmetically unreachable (ceiling is
#    n_val / (n_val + 1)); the warning must say so, not just "noisy".

test_that("joint MVN baseline warning describes covariate omission", {
  set.seed(42L)
  n <- 24L
  tree <- ape::rcoal(n)
  # >= 2 continuous traits => use_continuous_joint fires (single-obs,
  # no multi_proportion, default lambda_mode).
  df <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                   row.names = tree$tip.label)
  df$x1[c(2L, 9L)] <- NA
  df$x2[c(4L, 11L)] <- NA
  covs <- data.frame(env = stats::rnorm(n), row.names = tree$tip.label)
  pd <- preprocess_traits(df, tree, covariates = covs)

  warnings <- character()
  fit <- withCallingHandlers(
    fit_baseline(pd, tree, joint_solver = "inhouse"),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("omitted from the selected joint MVN baseline route", warnings)))
  expect_true(all(fit$path == "joint_mvn"))
})

test_that("joint MVN baseline is silent without covariates", {
  set.seed(43L)
  n <- 24L
  tree <- ape::rcoal(n)
  df <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                   row.names = tree$tip.label)
  df$x1[c(2L, 9L)] <- NA
  df$x2[c(4L, 11L)] <- NA
  pd <- preprocess_traits(df, tree)

  expect_no_warning(fit_baseline(pd, tree))
})

test_that("threshold-joint baseline warning names its route", {
  set.seed(451L)
  n <- 36L
  tree <- ape::rcoal(n)
  df <- data.frame(x = stats::rnorm(n),
                   b = factor(ifelse(stats::rnorm(n) > 0, "yes", "no")),
                   row.names = tree$tip.label)
  df$x[c(3L, 10L)] <- NA_real_
  df$b[c(5L, 12L)] <- NA
  covs <- data.frame(env = stats::rnorm(n), row.names = tree$tip.label)
  pd <- preprocess_traits(df, tree, covariates = covs)
  warnings <- character()
  fit <- withCallingHandlers(
    fit_baseline(pd, tree, joint_solver = "inhouse"),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("omitted from the selected joint threshold-joint baseline route", warnings)))
  expect_true(any(fit$path == "threshold_joint"))
})

test_that("OVR-only categorical baseline omits covariates without the joint warning", {
  set.seed(452L)
  n <- 36L
  tree <- ape::rcoal(n)
  df <- data.frame(cat = factor(sample(c("a", "b", "c"), n, replace = TRUE)),
                   row.names = tree$tip.label)
  df$cat[c(4L, 13L)] <- NA
  covs <- data.frame(env = stats::rnorm(n), row.names = tree$tip.label)
  pd <- preprocess_traits(df, tree, covariates = covs)
  warnings <- character()
  fit <- withCallingHandlers(
    fit_baseline(pd, tree, joint_solver = "inhouse"),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(all(fit$path == "ovr_categorical"))
  expect_false(any(grepl("omitted from the selected joint", warnings)))
})

test_that("mixed numeric and threshold routes retain covariates on the numeric fallback", {
  set.seed(453L)
  n <- 36L
  tree <- ape::rcoal(n)
  df <- data.frame(x = stats::rnorm(n),
                   b = factor(ifelse(stats::rnorm(n) > 0, "yes", "no")),
                   row.names = tree$tip.label)
  df$x[c(3L, 10L)] <- NA_real_
  df$b[c(5L, 12L)] <- NA
  covs <- data.frame(env = stats::rnorm(n), row.names = tree$tip.label)
  pd <- preprocess_traits(df, tree, covariates = covs)
  cov_fit_calls <- 0L
  trace("bm_impute_col_with_cov",
        tracer = function() cov_fit_calls <<- cov_fit_calls + 1L,
        print = FALSE, where = asNamespace("pigauto"))
  on.exit(try(untrace("bm_impute_col_with_cov", where = asNamespace("pigauto")),
              silent = TRUE), add = TRUE)
  warnings <- character()
  fit <- withCallingHandlers(
    fit_baseline(pd, tree, lambda_mode = "cv", lambda_fixed = c(x = 1),
                 discrete_lambda = "fixed_1", joint_solver = "inhouse",
                 predict_method = "per_column"),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(unname(fit$path[["x"]]), "per_column_bm")
  expect_identical(unname(fit$path[["b"]]), "threshold_joint")
  expect_gt(cov_fit_calls, 0L)
  expect_true(any(grepl("omitted from the selected joint threshold-joint baseline route", warnings)))
})

test_that("[conformal-ceiling] calibrate_gates names the n_val < 19 ceiling", {
  set.seed(2033L)
  n <- 30L
  tm <- list(list(name = "x", type = "continuous", latent_cols = 1L,
                  mean = 0, sd = 1))
  mu    <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "x"))
  delta <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "x"))
  truth <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "x"))
  val   <- matrix(FALSE, n, 1L)
  val[1:12, 1L] <- TRUE   # 12 val cells: above min_val_cells (10), below 19

  # 12 >= min_val_cells so the "Small validation set" branch does NOT
  # fire; the unreachable-ceiling branch must fire on its own.
  expect_warning(
    pigauto:::calibrate_gates(
      trait_map = tm, mu_cal = mu, delta_cal = delta,
      X_truth_r = truth, val_mask_mat = val,
      gate_grid = seq(0, 1, 0.1), gate_cap = 1,
      safety_floor = TRUE,
      mean_baseline_per_col = c(x = 0),
      simplex_step = 0.1,
      latent_names = "x", verbose = FALSE, seed = 2033L),
    regexp = "NOT achievable.*n_val / \\(n_val \\+ 1\\)"
  )
})

test_that("[conformal-ceiling] no ceiling warning at n_val >= 19", {
  set.seed(2034L)
  n <- 60L
  tm <- list(list(name = "x", type = "continuous", latent_cols = 1L,
                  mean = 0, sd = 1))
  mu    <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "x"))
  delta <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "x"))
  truth <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "x"))
  val   <- matrix(FALSE, n, 1L)
  val[1:20, 1L] <- TRUE   # 20 val cells: ceiling 20/21 > 0.95

  expect_no_warning(
    pigauto:::calibrate_gates(
      trait_map = tm, mu_cal = mu, delta_cal = delta,
      X_truth_r = truth, val_mask_mat = val,
      gate_grid = seq(0, 1, 0.1), gate_cap = 1,
      safety_floor = TRUE,
      mean_baseline_per_col = c(x = 0),
      simplex_step = 0.1,
      latent_names = "x", verbose = FALSE, seed = 2034L)
  )
})

test_that("[conformal-ceiling] binary traits do not trigger the ceiling warning", {
  # Conformal scores are not computed for binary/categorical, so the
  # ceiling message would be misleading there.
  set.seed(2035L)
  n <- 30L
  tm <- list(list(name = "b", type = "binary", latent_cols = 1L,
                  levels = c("no", "yes"), mean = 0, sd = 1))
  mu    <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "b"))
  delta <- matrix(rnorm(n), n, 1L, dimnames = list(NULL, "b"))
  truth <- matrix(rbinom(n, 1L, 0.5), n, 1L, dimnames = list(NULL, "b"))
  val   <- matrix(FALSE, n, 1L)
  val[1:12, 1L] <- TRUE   # 12 val cells, binary

  expect_no_warning(
    pigauto:::calibrate_gates(
      trait_map = tm, mu_cal = mu, delta_cal = delta,
      X_truth_r = truth, val_mask_mat = val,
      gate_grid = seq(0, 1, 0.1), gate_cap = 1,
      safety_floor = TRUE,
      mean_baseline_per_col = c(b = 0),
      simplex_step = 0.1,
      latent_names = "b", verbose = FALSE, seed = 2035L)
  )
})

test_that("[split-val] impute() exposes conformal_split_val explicitly", {
  expect_true("conformal_split_val" %in% names(formals(impute)))
  set.seed(2036L)
  n <- 30L
  tree <- ape::rcoal(n)
  df <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n),
                   row.names = tree$tip.label)
  df$x1[sample(n, 6L)] <- NA
  df$x2[sample(n, 6L)] <- NA
  res <- suppressWarnings(
    impute(df, tree, epochs = 5L, verbose = FALSE, seed = 2036L,
           conformal_split_val = TRUE)
  )
  expect_s3_class(res, "pigauto_result")
  expect_false(anyNA(res$completed))
})
