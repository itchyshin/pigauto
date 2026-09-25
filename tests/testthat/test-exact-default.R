# predict_method = "exact" default flip (S3, docs/dev-log/exact-default/
# S3-default-report.md). Covers: defaults resolving to "exact" at the three
# user-facing entry points, $model_config$predict_method_used, the quiet
# once-per-session fallback message vs the explicit warning, and a discrete
# (binary + categorical) accuracy check under the new default.

mixed_fixture <- function(n = 60, seed) {
  set.seed(seed)
  tree <- ape::rtree(n)
  sp <- tree$tip.label
  df <- data.frame(
    row.names = sp,
    mass = abs(stats::rnorm(n)) + 0.5,
    wing = abs(stats::rnorm(n)) * 2 + 1,
    migr = factor(sample(c("no", "yes"), n, replace = TRUE)),
    diet = factor(sample(c("herb", "carn", "omni"), n, replace = TRUE))
  )
  list(tree = tree, df = df)
}

# Masked-cell accuracy for the binary (migr) and categorical (diet) traits
# combined. Cells belonging to one trait are held out together (mask_missing.R),
# so checking the first latent column of each trait for NA-after-masking
# identifies every masked row for that trait.
discrete_accuracy <- function(bl, pd, spl) {
  X_masked <- pd$X_scaled
  X_masked[spl$val_idx]  <- NA
  X_masked[spl$test_idx] <- NA

  migr_col <- pd$trait_map$migr$latent_cols
  rows_migr  <- which(is.na(X_masked[, migr_col]) & !is.na(pd$X_scaled[, migr_col]))
  truth_migr <- pd$X_scaled[rows_migr, migr_col]
  pred_migr  <- as.numeric(stats::plogis(bl$mu[rows_migr, migr_col]) > 0.5)

  diet_cols <- pd$trait_map$diet$latent_cols
  rows_diet  <- which(is.na(X_masked[, diet_cols[1]]) & !is.na(pd$X_scaled[, diet_cols[1]]))
  truth_diet <- apply(pd$X_scaled[rows_diet, diet_cols, drop = FALSE], 1, which.max)
  pred_diet  <- apply(bl$mu[rows_diet, diet_cols, drop = FALSE], 1, which.max)

  correct <- sum(pred_migr == truth_migr) + sum(pred_diet == truth_diet)
  total   <- length(rows_migr) + length(rows_diet)
  list(acc = correct / total, n = total)
}

test_that("[exact-default] impute()/fit_pigauto()/fit_baseline() default to predict_method = 'auto'", {
  # S5b (docs/dev-log/exact-default/S5b-route-choice-report.md): the
  # dispatcher default moved from "exact" to "auto" (per-trait choice on
  # the validation split); "exact" is still the SECOND enum entry (the
  # concrete route .fit_baseline_core() runs when the auto path falls
  # through with no validation evidence).
  expect_identical(eval(formals(fit_baseline)$predict_method)[1], "auto")
  expect_identical(eval(formals(fit_pigauto)$predict_method)[1], "auto")
  expect_identical(eval(formals(impute)$predict_method)[1], "auto")
})

test_that("[exact-default] model_config$predict_method_used is 'exact' on a small fit", {
  skip_if_not_installed("Matrix")
  set.seed(1)
  n <- 30L
  tree <- ape::rcoal(n)
  df <- data.frame(a = stats::rnorm(n), b = stats::rnorm(n), row.names = tree$tip.label)
  df$a[1:5] <- NA; df$b[6:10] <- NA
  # Pinned to predict_method = "exact" (S5b): this test is about the
  # EXACT route's own $model_config$predict_method_used bookkeeping, not
  # about what the "auto" default happens to pick on this particular
  # random draw. See tests/testthat/test-route-choice.R for "auto" itself.
  res <- suppressWarnings(impute(df, tree, epochs = 5L, verbose = FALSE, seed = 1,
                                  predict_method = "exact"))
  expect_identical(res$fit$model_config$predict_method_used, "exact")
})

test_that("[exact-default] oversized/mocked exact under the default: message, per_column, predict_method_used", {
  skip_if_not_installed("Matrix")
  pigauto:::.pigauto_exact_fallback_reset()
  testthat::local_mocked_bindings(exact_conditional_mvn = function(...) NULL,
                                   .package = "pigauto")
  set.seed(2)
  n <- 25L; K <- 3L
  tree <- ape::rcoal(n)
  df <- data.frame(a = stats::rnorm(n), b = stats::rnorm(n), c = stats::rnorm(n),
                   row.names = tree$tip.label)
  df$a[1:5] <- NA; df$b[6:10] <- NA; df$c[11:15] <- NA
  pd <- preprocess_traits(df, tree)

  # Left at the true default (S5b: "auto"), no predict_method argument: with
  # splits = NULL there are no validation cells for "auto" to compare
  # routes with, so it runs the exact route internally with
  # predict_method_explicit = FALSE (docs/dev-log/exact-default/
  # S5b-route-choice-report.md, "No splits ... use exact for every trait")
  # -- same message-not-warning behaviour as the pre-S5b "exact" default.
  # $predict_method_used is now the constant "auto" marker;
  # $predict_method_by_trait carries the real per-trait outcome, which
  # still shows "per_column" here because of the mocked fallback.
  bl <- expect_no_warning(fit_baseline(pd, tree))
  expect_identical(bl$predict_method_used, "auto")
  expect_true(all(bl$predict_method_by_trait == "per_column"))

  pigauto:::.pigauto_exact_fallback_reset()
  msgs <- character(0)
  withCallingHandlers(
    fit_baseline(pd, tree),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(length(msgs) >= 1L)
  expect_match(msgs[1], "predict_method")
})

test_that("[exact-default] explicit predict_method = 'exact' on the same mock still warns", {
  skip_if_not_installed("Matrix")
  pigauto:::.pigauto_exact_fallback_reset()
  testthat::local_mocked_bindings(exact_conditional_mvn = function(...) NULL,
                                   .package = "pigauto")
  set.seed(3)
  n <- 25L
  tree <- ape::rcoal(n)
  df <- data.frame(a = stats::rnorm(n), b = stats::rnorm(n), row.names = tree$tip.label)
  df$a[1:5] <- NA; df$b[6:10] <- NA
  pd <- preprocess_traits(df, tree)

  expect_warning(
    bl <- fit_baseline(pd, tree, predict_method = "exact"),
    "falling back to the per-column path"
  )
  expect_identical(bl$predict_method_used, "per_column")
})

test_that("[exact-default] discrete accuracy under exact is not worse than per_column", {
  skip_if_not_installed("Matrix")
  n_seeds <- 10L
  acc_exact <- numeric(n_seeds)
  acc_pc    <- numeric(n_seeds)
  for (s in seq_len(n_seeds)) {
    fx  <- mixed_fixture(n = 60, seed = 200 + s)
    pd  <- preprocess_traits(fx$df, fx$tree)
    spl <- make_missing_splits(pd$X_scaled, seed = s, trait_map = pd$trait_map)
    # Pinned to predict_method = "exact" (S5b): this test compares the
    # concrete exact and per_column routes against each other, not "auto"
    # (the new default) against per_column.
    bl_exact <- fit_baseline(pd, fx$tree, splits = spl, predict_method = "exact")
    bl_pc    <- fit_baseline(pd, fx$tree, splits = spl, predict_method = "per_column")
    acc_exact[s] <- discrete_accuracy(bl_exact, pd, spl)$acc
    acc_pc[s]    <- discrete_accuracy(bl_pc, pd, spl)$acc
  }
  # Measured (2026-09-25): mean exact accuracy 0.411, mean per_column 0.396,
  # across 10 seeds on a random (no real phylogenetic signal) mixed fixture.
  expect_gte(mean(acc_exact), mean(acc_pc) - 0.02)
})
