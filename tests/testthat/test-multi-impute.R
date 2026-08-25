skip_if_no_libtorch()

# Tests for multi_impute(), with_imputations(), and pool_mi().
# Uses tiny synthetic data (small tree, short training) to keep runtime low.

# ---- Helpers ----------------------------------------------------------------

make_mi_test_data <- function(n = 40, p = 3, miss_frac = 0.15, seed = 123) {
  set.seed(seed)
  tree <- ape::rtree(n)
  sp   <- tree$tip.label
  df <- data.frame(
    row.names = sp,
    tr1 = abs(stats::rnorm(n)) + 0.5,
    tr2 = abs(stats::rnorm(n)) + 0.5,
    tr3 = abs(stats::rnorm(n)) + 0.5
  )
  # Punch holes in the matrix so we have actual missing cells to impute.
  n_cells   <- n * p
  n_missing <- max(1L, round(n_cells * miss_frac))
  idx <- sample.int(n_cells, n_missing)
  m   <- as.matrix(df)
  m[idx] <- NA
  df <- as.data.frame(m)
  rownames(df) <- sp
  list(tree = tree, df = df, n_missing = n_missing)
}

quick_mi <- function(m = 3L, seed = 123, draws_method = "conformal") {
  td <- make_mi_test_data(seed = seed)
  mi <- multi_impute(
    traits        = td$df,
    tree          = td$tree,
    m             = m,
    draws_method  = draws_method,
    epochs        = 20L,
    missing_frac  = 0.25,
    verbose       = FALSE,
    seed          = seed,
    eval_every    = 10L,
    patience      = 5L
  )
  # Include traits_with_na so tests can identify originally-missing cells
  list(mi = mi, td = td, traits_with_na = td$df)
}

pool_unverified <- function(...) {
  warnings <- character()
  result <- withCallingHandlers(
    pool_mi(...),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warnings, 1L)
  expect_identical(
    warnings,
    "pool_mi() received a bare fit list with unverified provenance; Rubin arithmetic will run, but the imputation workflow was not validated by pigauto."
  )
  result
}


# ---- 1. multi_impute() structural shape -------------------------------------

test_that("multi_impute returns pigauto_mi with M datasets of the right shape", {
  setup <- quick_mi(m = 3L, seed = 11)
  mi <- setup$mi

  expect_s3_class(mi, "pigauto_mi")
  expect_s3_class(mi, "pigauto_diagnostic_mi")
  expect_equal(mi$m, 3L)
  expect_true(is.list(mi$datasets))
  expect_equal(length(mi$datasets), 3L)

  # Each dataset should have the same shape as the input.
  for (d in mi$datasets) {
    expect_s3_class(d, "data.frame")
    expect_equal(nrow(d), nrow(setup$td$df))
    expect_equal(ncol(d), ncol(setup$td$df))
    expect_equal(names(d), names(setup$td$df))
    expect_equal(rownames(d), rownames(setup$td$df))
  }

  # Slots expected by downstream helpers.
  expect_true("imputed_mask" %in% names(mi))
  expect_true("fit"          %in% names(mi))
  expect_true("tree"         %in% names(mi))
  expect_identical(mi$mi_workflow, "pigauto_diagnostic_mi")
})


# ---- 2. Dataset consistency across imputations ------------------------------

test_that("all M datasets have matching rownames, colnames, and column classes", {
  setup <- quick_mi(m = 3L, seed = 12)
  mi <- setup$mi

  ref <- mi$datasets[[1]]
  ref_rn <- rownames(ref)
  ref_cn <- colnames(ref)
  ref_cls <- vapply(ref, function(x) class(x)[1], character(1))

  for (d in mi$datasets) {
    expect_equal(rownames(d), ref_rn)
    expect_equal(colnames(d), ref_cn)
    expect_equal(vapply(d, function(x) class(x)[1], character(1)), ref_cls)
  }
})


# ---- 3. Observed cells are identical across imputations ---------------------

test_that("cells that were observed in the input are identical across datasets", {
  setup <- quick_mi(m = 3L, seed = 13)
  mi <- setup$mi
  orig <- setup$td$df

  # For every originally observed cell, all M datasets should return
  # exactly the input value.
  obs_mask <- !is.na(orig)
  for (d in mi$datasets) {
    for (j in seq_len(ncol(orig))) {
      col_obs <- obs_mask[, j]
      expect_equal(
        d[col_obs, j],
        orig[col_obs, j],
        tolerance = 0,  # exact match expected for observed cells
        info = paste("column", colnames(orig)[j])
      )
    }
  }
})


# ---- 4. Imputed cells vary across draws -------------------------------------

test_that("cells that were originally missing vary across at least some datasets", {
  setup <- quick_mi(m = 3L, seed = 14)
  mi <- setup$mi
  miss <- is.na(setup$td$df)

  # Collect the values at originally-missing positions for each draw.
  if (sum(miss) == 0L) skip("no missing cells in test data")

  mat_stack <- do.call(rbind, lapply(mi$datasets, function(d) {
    as.numeric(as.matrix(d)[miss])
  }))
  # Should have M rows, n_missing cols
  expect_equal(nrow(mat_stack), 3L)

  # At least one missing cell should have non-zero variance across M draws.
  # (BM posterior draws guarantee this even when calibrated gates are 0.)
  col_sds <- apply(mat_stack, 2, stats::sd, na.rm = TRUE)
  expect_true(any(col_sds > 1e-10),
              info = "draws should produce varying imputations")
})


# ---- 5. Numerical correctness of Rubin's rules against hand calculation ----

test_that("pool_mi() matches a hand-computed Rubin's rules reference (M=3, p=2)", {
  # Three synthetic fits with known coef and vcov.
  make_fake <- function(beta, var_diag) {
    V <- diag(var_diag)
    dimnames(V) <- list(names(beta), names(beta))
    list(beta = beta, V = V)
  }
  fits <- list(
    make_fake(c(a = 1.00, b = 2.00), c(0.010, 0.040)),
    make_fake(c(a = 1.20, b = 2.10), c(0.020, 0.050)),
    make_fake(c(a = 0.90, b = 1.80), c(0.015, 0.030))
  )

  pooled <- pool_unverified(
    fits,
    coef_fun = function(f) f$beta,
    vcov_fun = function(f) f$V
  )

  M  <- 3L
  beta_a <- c(1.00, 1.20, 0.90)
  beta_b <- c(2.00, 2.10, 1.80)
  theta_a <- mean(beta_a)
  theta_b <- mean(beta_b)
  W_a <- mean(c(0.010, 0.020, 0.015))
  W_b <- mean(c(0.040, 0.050, 0.030))
  B_a <- stats::var(beta_a)
  B_b <- stats::var(beta_b)
  T_a <- W_a + (1 + 1 / M) * B_a
  T_b <- W_b + (1 + 1 / M) * B_b
  r_a <- (1 + 1 / M) * B_a / W_a
  r_b <- (1 + 1 / M) * B_b / W_b
  df_a <- (M - 1) * (1 + 1 / r_a)^2
  df_b <- (M - 1) * (1 + 1 / r_b)^2

  # Order of rows in pooled should match input coef names.
  row_a <- which(pooled$term == "a")
  row_b <- which(pooled$term == "b")

  expect_equal(pooled$estimate[row_a],  theta_a,    tolerance = 1e-10)
  expect_equal(pooled$estimate[row_b],  theta_b,    tolerance = 1e-10)
  expect_equal(pooled$std.error[row_a], sqrt(T_a),  tolerance = 1e-10)
  expect_equal(pooled$std.error[row_b], sqrt(T_b),  tolerance = 1e-10)
  expect_equal(pooled$riv[row_a],       r_a,        tolerance = 1e-10)
  expect_equal(pooled$riv[row_b],       r_b,        tolerance = 1e-10)
  expect_equal(pooled$df[row_a],        df_a,       tolerance = 1e-10)
  expect_equal(pooled$df[row_b],        df_b,       tolerance = 1e-10)

  # fmi bounded in [0, 1] and positive because B > 0.
  expect_true(all(pooled$fmi > 0 & pooled$fmi < 1))
})


test_that("pool_mi() handles Rubin variance boundaries exactly", {
  make_fake <- function(beta, V) list(beta = beta, V = V)
  coef_fun <- function(f) f$beta
  vcov_fun <- function(f) f$V
  named_vcov <- function(x) {
    V <- diag(x)
    dimnames(V) <- list(c("a", "b"), c("a", "b"))
    V
  }

  same <- replicate(
    3L,
    make_fake(c(a = 1, b = 2), named_vcov(c(0.25, 0.5))),
    simplify = FALSE
  )
  classical <- pool_unverified(same, coef_fun = coef_fun, vcov_fun = vcov_fun)
  expect_equal(classical$riv, c(0, 0))
  expect_equal(classical$fmi, c(0, 0))
  expect_true(all(is.infinite(classical$df)))

  adjusted <- pool_unverified(
    same,
    coef_fun = coef_fun,
    vcov_fun = vcov_fun,
    df_fun = function(f) 10
  )
  expected_df <- 10 * 11 / 13
  expect_equal(adjusted$df, rep(expected_df, 2), tolerance = 1e-12)
  expect_equal(adjusted$fmi, rep(2 / (expected_df + 3), 2),
               tolerance = 1e-12)
  expect_true(all(adjusted$fmi > 0))

  deterministic <- replicate(
    3L,
    make_fake(c(a = 1, b = 2), named_vcov(c(0, 0))),
    simplify = FALSE
  )
  zero <- pool_unverified(
    deterministic,
    coef_fun = coef_fun,
    vcov_fun = vcov_fun,
    df_fun = function(f) 10
  )
  expect_equal(zero$std.error, c(0, 0))
  expect_equal(zero$riv, c(0, 0))
  expect_equal(zero$fmi, c(0, 0))
  expect_true(all(is.infinite(zero$df)))
  expect_equal(zero$conf.low, zero$estimate)
  expect_equal(zero$conf.high, zero$estimate)
  expect_true(all(is.na(zero$statistic)))
  expect_true(all(is.na(zero$p.value)))

  between_only <- list(
    make_fake(c(a = 0), matrix(0, 1, 1, dimnames = list("a", "a"))),
    make_fake(c(a = 1), matrix(0, 1, 1, dimnames = list("a", "a"))),
    make_fake(c(a = 2), matrix(0, 1, 1, dimnames = list("a", "a")))
  )
  pure_between <- pool_unverified(
    between_only,
    coef_fun = coef_fun,
    vcov_fun = vcov_fun,
    df_fun = function(f) 10
  )
  expect_true(is.infinite(pure_between$riv))
  expect_equal(pure_between$fmi, 1)
  expect_equal(pure_between$df, 2)
})


test_that("pool_mi() accepts Matrix covariance objects and validates them", {
  make_fake <- function(V) list(beta = c(a = 1, b = 2), V = V)
  V <- matrix(c(0.2, 0.03, 0.03, 0.4), 2, 2,
              dimnames = list(c("a", "b"), c("a", "b")))
  base <- pool_unverified(
    list(make_fake(V), make_fake(V)),
    coef_fun = function(f) f$beta,
    vcov_fun = function(f) f$V
  )
  sparse_V <- Matrix::Matrix(V, sparse = TRUE)
  sparse <- pool_unverified(
    list(make_fake(sparse_V), make_fake(sparse_V)),
    coef_fun = function(f) f$beta,
    vcov_fun = function(f) f$V
  )
  expect_equal(sparse, base)

  expect_error(
    pool_mi(
      list(make_fake(matrix(1, 2, 3)), make_fake(V)),
      coef_fun = function(f) f$beta, vcov_fun = function(f) f$V
    ),
    "square"
  )
  V_na <- V
  V_na[1, 1] <- NA_real_
  expect_error(
    pool_mi(
      list(make_fake(V_na), make_fake(V_na)),
      coef_fun = function(f) f$beta, vcov_fun = function(f) f$V
    ),
    "non-finite"
  )
  V_negative <- V
  V_negative[1, 1] <- -0.1
  expect_error(
    pool_mi(
      list(make_fake(V_negative), make_fake(V_negative)),
      coef_fun = function(f) f$beta, vcov_fun = function(f) f$V
    ),
    "negative"
  )
  V_asymmetric <- V
  V_asymmetric[1, 2] <- 0.2
  expect_error(
    pool_mi(
      list(make_fake(V_asymmetric), make_fake(V_asymmetric)),
      coef_fun = function(f) f$beta, vcov_fun = function(f) f$V
    ),
    "symmetric"
  )
  V_names <- V
  dimnames(V_names) <- list(c("a", "c"), c("a", "c"))
  expect_error(
    pool_mi(
      list(make_fake(V_names), make_fake(V_names)),
      coef_fun = function(f) f$beta, vcov_fun = function(f) f$V
    ),
    "align"
  )
})


test_that("pool_mi() supports a strict tidy extractor contract", {
  fits <- list(list(i = 1), list(i = 2))
  tidy_fun <- function(f) {
    data.frame(
      term = c("a", "b"),
      estimate = c(f$i, 2 * f$i),
      std.error = c(0.2, 0.4)
    )
  }
  pooled <- pool_unverified(fits, tidy_fun = tidy_fun)
  expect_equal(pooled$term, c("a", "b"))
  expect_equal(pooled$estimate, c(1.5, 3))

  expect_error(
    pool_mi(fits, tidy_fun = tidy_fun, coef_fun = function(f) f$i),
    "cannot be combined"
  )
  expect_error(
    pool_mi(fits, tidy_fun = function(f) {
      data.frame(term = c("a", "a"), estimate = c(1, 2),
                 std.error = c(0.1, 0.1))
    }),
    "duplicate"
  )
  expect_error(
    pool_mi(fits, tidy_fun = function(f) {
      data.frame(term = "a", estimate = 1, std.error = -0.1)
    }),
    "non-negative"
  )
})


test_that("automatic fixed-effect adapters cover lm and glm", {
  set.seed(1901)
  dat <- data.frame(
    y = rnorm(60),
    x = rnorm(60),
    z = rnorm(60),
    g = factor(rep(seq_len(12), each = 5))
  )
  dat$binary <- stats::rbinom(60, 1, stats::plogis(0.2 + 0.5 * dat$x))

  lm_fit <- stats::lm(y ~ x + z, data = dat)
  glm_fit <- stats::glm(binary ~ x + z, data = dat,
                        family = stats::binomial())
  expect_equal(pool_unverified(list(lm_fit, lm_fit))$term,
               names(stats::coef(lm_fit)))
  expect_equal(pool_unverified(list(glm_fit, glm_fit))$term,
               names(stats::coef(glm_fit)))
})


test_that("glmmTMB adapter selects conditional fixed effects only", {
  skip_if_not_installed("glmmTMB")
  set.seed(1904)
  n <- 240L
  dat <- data.frame(x = rnorm(n), z = rnorm(n))
  mu <- exp(0.2 + 0.4 * dat$x)
  dat$y <- stats::rnbinom(n, mu = mu, size = exp(0.3 - 0.2 * dat$x))
  dat$y[stats::rbinom(n, 1, stats::plogis(-1 + 0.5 * dat$z)) == 1L] <- 0L

  fit <- suppressWarnings(glmmTMB::glmmTMB(
    y ~ x,
    ziformula = ~ z,
    dispformula = ~ x,
    family = glmmTMB::nbinom2,
    data = dat
  ))
  fixed <- glmmTMB::fixef(fit)
  expect_gt(length(fixed$zi), 0L)
  expect_gt(length(fixed$disp), 0L)

  pooled <- pool_unverified(list(fit, fit))
  expect_equal(pooled$term, names(fixed$cond))
  expect_equal(pooled$estimate, unname(fixed$cond))
  expect_equal(nrow(pooled), length(fixed$cond))
  expect_lt(nrow(pooled), sum(lengths(fixed)))
})


test_that("automatic fixed-effect adapters cover nlme gls and lme", {
  skip_if_not_installed("nlme")
  set.seed(1901)
  dat <- data.frame(
    y = rnorm(60), x = rnorm(60), z = rnorm(60),
    g = factor(rep(seq_len(12), each = 5))
  )
  gls_fit <- nlme::gls(y ~ x + z, data = dat)
  lme_fit <- nlme::lme(y ~ x + z, random = ~1 | g, data = dat)
  expect_equal(pool_unverified(list(gls_fit, gls_fit))$term,
               names(stats::coef(gls_fit)))
  expect_equal(pool_unverified(list(lme_fit, lme_fit))$term,
               names(nlme::fixef(lme_fit)))
})


test_that("automatic fixed-effect adapter covers lme4 Matrix covariance", {
  skip_if_not_installed("lme4")
  set.seed(1901)
  dat <- data.frame(
    y = rnorm(60), x = rnorm(60), z = rnorm(60),
    g = factor(rep(seq_len(12), each = 5))
  )
  mer_fit <- suppressMessages(lme4::lmer(y ~ x + z + (1 | g), data = dat))
  mer_pooled <- pool_unverified(list(mer_fit, mer_fit))
  expect_equal(mer_pooled$term, names(lme4::fixef(mer_fit)))
  expect_true(all(is.finite(mer_pooled$std.error)))

  custom_coef <- pool_unverified(
    list(mer_fit, mer_fit),
    coef_fun = function(f) lme4::fixef(f)
  )
  custom_vcov <- pool_unverified(
    list(mer_fit, mer_fit),
    vcov_fun = function(f) stats::vcov(f)
  )
  expect_equal(custom_coef, mer_pooled)
  expect_equal(custom_vcov, mer_pooled)
})


test_that("automatic drmTMB adapter pools all distributional fixed effects", {
  skip_if_not_installed("drmTMB")
  drm_fit <- getExportedValue("drmTMB", "drmTMB")
  bf <- getExportedValue("drmTMB", "bf")
  set.seed(1902)
  dat <- data.frame(x = seq(-1, 1, length.out = 40))
  dat$y <- 0.3 + 0.7 * dat$x + rnorm(40, sd = exp(-0.2 + 0.1 * dat$x))
  fit1 <- drm_fit(bf(y ~ x, sigma ~ x), data = dat)
  dat$y <- dat$y + rnorm(40, sd = 0.01)
  fit2 <- drm_fit(bf(y ~ x, sigma ~ x), data = dat)

  pooled <- pool_unverified(list(fit1, fit2))
  expect_equal(
    pooled$term,
    c("mu:(Intercept)", "mu:x", "sigma:(Intercept)", "sigma:x")
  )
  expect_true(all(is.finite(pooled$estimate)))
  expect_true(all(is.finite(pooled$std.error)))
})


test_that("automatic gllvmTMB adapter uses fixed-effect tidy rows only", {
  skip_if_not_installed("gllvmTMB")
  simulate_site_trait <- getExportedValue("gllvmTMB", "simulate_site_trait")
  gllvm_fit <- getExportedValue("gllvmTMB", "gllvmTMB")
  set.seed(1903)
  sim <- simulate_site_trait(
    n_sites = 30, n_species = 8, n_traits = 2,
    mean_species_per_site = 5,
    Lambda_B = matrix(c(0.5, -0.2), nrow = 2),
    psi_B = c(0.2, 0.2), sigma2_eps = 0.1, seed = 1903
  )
  fit <- suppressMessages(suppressWarnings(gllvm_fit(
    value ~ 0 + trait + latent(0 + trait | site, d = 1),
    data = sim$data, silent = TRUE
  )))

  tidy_method <- utils::getS3method(
    "tidy", "gllvmTMB_multi", envir = asNamespace("gllvmTMB")
  )
  fixed <- tidy_method(fit, effects = "fixed")

  # Saved fits must remain poolable in a fresh session where gllvmTMB has not
  # already been attached or loaded. The adapter is responsible for loading
  # the namespace that owns the class-specific tidy method.
  unloadNamespace("gllvmTMB")
  on.exit(loadNamespace("gllvmTMB"), add = TRUE)
  expect_false("gllvmTMB" %in% loadedNamespaces())
  pooled <- pool_unverified(list(fit, fit))
  expect_true("gllvmTMB" %in% loadedNamespaces())
  expect_equal(pooled$term, fixed$term)
  expect_equal(pooled$estimate, fixed$estimate)
  expect_equal(pooled$std.error, fixed$std.error)
  expect_false(any(grepl("^sd_|loglambda|cutpoint", pooled$term)))
})


# ---- 10. pool_mi() rejects MCMCglmm fits ------------------------------------

test_that("pool_mi() errors cleanly when passed MCMCglmm fits", {
  # Fake MCMCglmm-classed objects -- no need to actually run MCMCglmm.
  fake <- structure(list(), class = "MCMCglmm")
  expect_error(
    pool_mi(list(fake, fake)),
    regexp = "MCMCglmm|posterior"
  )
})


# ---- 11. pool_mi() errors on inconsistent coefficient names -----------------

test_that("pool_mi() errors when fits have inconsistent coefficient names", {
  make_fake <- function(beta, var_diag) {
    V <- diag(var_diag)
    dimnames(V) <- list(names(beta), names(beta))
    list(beta = beta, V = V)
  }
  fits <- list(
    make_fake(c(a = 1.0, b = 2.0), c(0.01, 0.02)),
    make_fake(c(a = 1.1, c = 2.1), c(0.01, 0.02))   # 'c' instead of 'b'
  )
  expect_error(
    pool_mi(
      fits,
      coef_fun = function(f) f$beta,
      vcov_fun = function(f) f$V
    ),
    regexp = "names differ|Rubin"
  )
})


# ---- 12. multi_impute_trees() structural shape --------------------------------

test_that("multi_impute_trees returns pigauto_mi_trees with T*m datasets", {
  # Build 3 small random trees sharing the same tip labels
  n <- 40
  set.seed(200)
  tree1 <- ape::rtree(n)
  sp <- tree1$tip.label

  # Make 2 more trees by randomly perturbing edge lengths
  tree2 <- tree1; tree2$edge.length <- tree1$edge.length * stats::runif(length(tree1$edge.length), 0.5, 1.5)
  tree3 <- tree1; tree3$edge.length <- tree1$edge.length * stats::runif(length(tree1$edge.length), 0.5, 1.5)
  trees <- list(tree1, tree2, tree3)
  class(trees) <- "multiPhylo"

  df <- data.frame(
    row.names = sp,
    tr1 = abs(stats::rnorm(n)) + 0.5,
    tr2 = abs(stats::rnorm(n)) + 0.5,
    tr3 = abs(stats::rnorm(n)) + 0.5
  )
  # Punch some holes
  m <- as.matrix(df)
  idx <- sample.int(n * 3, 15)
  m[idx] <- NA
  df <- as.data.frame(m); rownames(df) <- sp

  # share_gnn = FALSE exercises the legacy per-tree-fit path, which this
  # test asserts on via $fits. suppressWarnings() silences the expected
  # low-M runtime warning triggered by the tiny T*m used for test speed.
  mi_t <- suppressWarnings(multi_impute_trees(
    traits       = df,
    trees        = trees,
    m_per_tree   = 2L,
    epochs       = 20L,
    missing_frac = 0.25,
    verbose      = FALSE,
    seed         = 200L,
    share_gnn    = FALSE,
    eval_every   = 10L,
    patience     = 5L
  ))

  # Class hierarchy

  expect_s3_class(mi_t, "pigauto_mi_trees")
  expect_s3_class(mi_t, "pigauto_tree_sensitivity_diagnostic")
  expect_identical(mi_t$mi_workflow, "pigauto_tree_sensitivity_diagnostic")
  expect_s3_class(mi_t, "pigauto_mi")

  # Structure
  expect_equal(mi_t$m, 6L)           # 3 trees * 2 imputations
  expect_equal(mi_t$n_trees, 3L)
  expect_equal(mi_t$m_per_tree, 2L)
  expect_equal(length(mi_t$datasets), 6L)
  expect_equal(length(mi_t$tree_index), 6L)
  expect_equal(mi_t$tree_index, c(1L, 1L, 2L, 2L, 3L, 3L))

  # Each dataset has correct shape
  for (d in mi_t$datasets) {
    expect_s3_class(d, "data.frame")
    expect_equal(nrow(d), n)
    expect_equal(ncol(d), 3L)
    expect_equal(names(d), c("tr1", "tr2", "tr3"))
  }

  # Fits list has one per tree (legacy behaviour, share_gnn = FALSE)
  expect_equal(length(mi_t$fits), 3L)
  for (f in mi_t$fits) {
    expect_s3_class(f, "pigauto_fit")
  }
})


# ---- 13. draws_method = "conformal" explicitly ----------------------------

test_that("draws_method='conformal' stores method and produces non-zero between-dataset variance", {
  setup <- quick_mi(m = 5L, seed = 77L, draws_method = "conformal")
  mi    <- setup$mi

  expect_equal(mi$draws_method, "conformal")

  # Collect values at originally-missing positions
  miss <- is.na(setup$traits_with_na[["tr1"]])
  if (sum(miss) == 0L) skip("no missing cells in test data")

  mat_stack <- do.call(rbind, lapply(mi$datasets, function(d) {
    as.numeric(d[["tr1"]][miss])
  }))
  col_sds <- apply(mat_stack, 2, stats::sd, na.rm = TRUE)
  expect_true(any(col_sds > 1e-10),
              info = "conformal draws should produce non-zero between-dataset variance")

  # Observed cells must be identical across all M datasets
  obs_vals <- lapply(mi$datasets, function(d) d[["tr1"]][!miss])
  for (i in seq_along(obs_vals)[-1])
    expect_identical(obs_vals[[i]], obs_vals[[1]],
                     info = paste("dataset", i, "should not alter observed cells"))
})

test_that("multi_impute() multi-obs aligns datasets when input is shuffled", {
  # Regression test for the multi_impute path of the row-alignment bug
  # found 2026-04-26 (commit a3e6d39).  preprocess_traits() reorders rows
  # internally to tree-tip order, and build_completed() needs the inverse
  # mapping to write predictions back to the user's input row order.  This
  # test exercises both draws_method paths (default conformal + mc_dropout)
  # to make sure both code paths in multi_impute.R were patched.
  set.seed(20260426)
  tree <- ape::rtree(15)
  # Two obs per species, species column SHUFFLED.
  sp_shuffled <- sample(rep(tree$tip.label, each = 2L))
  sp_truth <- setNames(rnorm(15, mean = 5, sd = 3), tree$tip.label)
  trait_full <- sp_truth[sp_shuffled] + rnorm(length(sp_shuffled), sd = 0.05)
  df <- data.frame(species = sp_shuffled, trait = trait_full,
                    stringsAsFactors = FALSE)

  # Mask one of each species's two observations.
  mask_idx <- vapply(unique(df$species),
                      function(sp) sample(which(df$species == sp), 1L),
                      integer(1))
  df_obs <- df
  df_obs$trait[mask_idx] <- NA

  for (method in c("conformal", "mc_dropout")) {
    mi <- multi_impute(
      traits        = df_obs,
      tree          = tree,
      species_col   = "species",
      m             = 2L,
      draws_method  = method,
      epochs        = 30L,
      missing_frac  = 0.0,
      verbose       = FALSE,
      seed          = 1L
    )

    for (k in seq_along(mi$datasets)) {
      d <- mi$datasets[[k]]
      observed_idx <- setdiff(seq_len(nrow(df)), mask_idx)
      # Observed rows preserved exactly (sanity check).
      expect_equal(d$trait[observed_idx], df$trait[observed_idx],
                    tolerance = 1e-8,
                    label = sprintf("draws_method=%s, dataset %d, observed", method, k))

      # Imputed rows must align with their own species.  We allow
      # generous tolerance because mc_dropout adds Gaussian noise.
      imp_pred <- d$trait[mask_idx]
      imp_truth <- sp_truth[as.character(df$species[mask_idx])]
      r <- cor(imp_pred, imp_truth)
      expect_gt(r, 0.85,
                 label = sprintf("cor(pred, truth) for draws_method=%s dataset %d (got %.3f)",
                                  method, k, r))
    }
  }
})

test_that("multi_impute(draws_method='conformal') produces non-zero between-draw variance on shuffled input", {
  # Regression test for the row-alignment bug in .sample_conformal_draw
  # found by adversarial review (commit fixing it: pending).
  # The bug: imputed_mask is in user-input row order, but pred$imputed is
  # in internal (tree-tip) order.  .sample_conformal_draw was indexing
  # imp[[nm]][rows] with rows from imputed_mask, so on shuffled-input
  # data it perturbed the WRONG cells, and build_completed re-aligned
  # the unperturbed (deterministic) values back into the user's missing
  # cells — producing identical datasets across all M draws (zero
  # between-imputation variance).
  set.seed(20260428)
  tree <- ape::rtree(15)
  # Single-obs but ROWNAMES are SHUFFLED relative to tree$tip.label so
  # preprocess_traits has to reorder.
  shuffled_sp <- sample(tree$tip.label)
  df <- data.frame(
    row.names = shuffled_sp,
    tr = abs(rnorm(15)) + 0.5
  )
  # Mask 5 of 15 cells.
  miss_idx <- sample(15, 5)
  df_obs <- df
  df_obs$tr[miss_idx] <- NA

  mi <- multi_impute(
    traits        = df_obs,
    tree          = tree,
    m             = 5L,
    draws_method  = "conformal",
    epochs        = 30L,
    missing_frac  = 0.0,
    verbose       = FALSE,
    seed          = 1L
  )

  # Pull values at the originally-missing cells from each of the 5 datasets.
  draws_at_miss <- sapply(mi$datasets, function(d) d$tr[miss_idx])
  # Compute SD across the 5 draws at each masked cell.
  per_cell_sd <- apply(draws_at_miss, 1L, stats::sd)

  # Pre-fix bug behaviour: per_cell_sd ≡ 0 for all cells (zero variance
  # because the wrong cells were perturbed).  Post-fix: the SD should be
  # strictly positive at every masked cell because conformal-width
  # Gaussian sampling was correctly applied.
  expect_true(all(per_cell_sd > 1e-8),
              info = sprintf(
                "All 5 masked cells should have non-zero SD across draws. Got: %s",
                paste(round(per_cell_sd, 6), collapse = ", ")))
})

test_that("conformal draws perturb multi_proportion groups on the simplex", {
  comp_cols <- c("a", "b", "c")
  tm <- list(
    name = "comp", type = "multi_proportion", latent_cols = 1:3,
    levels = comp_cols, input_cols = comp_cols, n_latent = 3L,
    mean = c(0, 0, 0), sd = c(1, 1, 1)
  )
  pred <- list(
    imputed = data.frame(
      a = c(1 / 3, 0.2), b = c(1 / 3, 0.3), c = c(1 / 3, 0.5)
    ),
    imputed_latent = matrix(c(0, 0, 0, -0.2, 0, 0.2),
                            nrow = 2, byrow = TRUE),
    se_latent = matrix(c(0.3, 0.3, 0.3, 0.3, 0.3, 0.3),
                       nrow = 2, byrow = TRUE),
    probabilities = list(),
    conformal_scores = c(comp = NA_real_)
  )
  imputed_mask <- matrix(
    c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
    nrow = 2, byrow = TRUE,
    dimnames = list(NULL, comp_cols)
  )

  se_lat <- compute_latent_se(list(pred$imputed_latent), list(tm),
                              pred$se_latent, c("s1", "s2"))
  expect_equal(unname(se_lat[, tm$latent_cols]), unname(pred$se_latent),
               tolerance = 0)

  draw1 <- .sample_conformal_draw(pred, imputed_mask, list(tm), seed_i = 1L)
  draw2 <- .sample_conformal_draw(pred, imputed_mask, list(tm), seed_i = 2L)

  first1 <- unlist(draw1[1, comp_cols], use.names = FALSE)
  first2 <- unlist(draw2[1, comp_cols], use.names = FALSE)
  second <- unlist(draw1[2, comp_cols], use.names = FALSE)

  expect_gt(max(abs(first1 - first2)), 1e-10)
  expect_equal(sum(first1), 1, tolerance = 1e-12)
  expect_equal(sum(first2), 1, tolerance = 1e-12)
  expect_equal(second,
               unlist(pred$imputed[2, comp_cols], use.names = FALSE),
               tolerance = 0)
})
