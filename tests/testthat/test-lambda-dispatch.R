# Tests for the S4 dispatcher slice (feat/joint-lambda-default):
# fit_baseline()'s lambda_mode = "estimate" default now keeps
# continuous-family columns on the joint MVN / threshold-joint path
# (lambda-aware via R/joint_mvn_solver.R's lambda_cols), instead of
# discarding that fit's continuous output for a separate per-column
# re-fit. Design: docs/dev-log/2026-09-22-joint-lambda-alignment.md
# section 7; dispatch changes: R/fit_baseline.R, R/joint_mvn_baseline.R,
# R/joint_threshold_baseline.R.

# ---- shared fixture --------------------------------------------------

# 3 continuous BM(lambda = 0.3) traits + 1 binary + 1 three-level
# categorical, all correlated with the SAME tree, 30% NA per column.
make_lambda_dispatch_fixture <- function(seed = 900L, n = 80L,
                                          lambda_true = 0.3,
                                          miss_frac = 0.3) {
  set.seed(seed)
  tree <- ape::rcoal(n)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  Lc <- chol(R)
  gen_bm_lambda <- function(lam) {
    phylo_part <- as.numeric(t(Lc) %*% stats::rnorm(n))
    iid_part <- stats::rnorm(n)
    sqrt(lam) * phylo_part + sqrt(1 - lam) * iid_part
  }
  c1 <- gen_bm_lambda(lambda_true)
  c2 <- gen_bm_lambda(lambda_true)
  c3 <- gen_bm_lambda(lambda_true)

  # Binary + categorical: threshold a phylo-correlated latent for some
  # signal. lambda has no meaning for these types (B iii cut); this is
  # for realism only.
  lat_b <- as.numeric(t(Lc) %*% stats::rnorm(n))
  b1 <- factor(ifelse(lat_b > stats::median(lat_b), "yes", "no"))
  lat_k <- as.numeric(t(Lc) %*% stats::rnorm(n))
  brks <- stats::quantile(lat_k, c(0, 1 / 3, 2 / 3, 1))
  k1 <- factor(cut(lat_k, brks, labels = c("A", "B", "C"),
                   include.lowest = TRUE))

  df <- data.frame(c1 = c1, c2 = c2, c3 = c3, b1 = b1, k1 = k1,
                    row.names = tree$tip.label)
  for (col in c("c1", "c2", "c3", "b1", "k1")) {
    df[sample(n, round(miss_frac * n)), col] <- NA
  }
  list(tree = tree, df = df)
}

fit_dispatch_baseline <- function(fx, ..., splits_seed = 1L) {
  pd <- preprocess_traits(fx$df, fx$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = splits_seed,
                              trait_map = pd$trait_map)
  list(pd = pd, splits = spl,
       bl = fit_baseline(pd, fx$tree, splits = spl, ...))
}

# ---- joint path under estimate --------------------------------------

test_that("[lambda-dispatch] joint path under estimate", {
  fx <- make_lambda_dispatch_fixture()
  fit <- NULL
  expect_no_warning(
    fit <- fit_dispatch_baseline(fx, lambda_mode = "estimate")
  )
  bl <- fit$bl
  cont_names <- c("c1", "c2", "c3")
  expect_true(all(bl$path[cont_names] %in% c("joint_mvn", "threshold_joint")))
})

# ---- lambda_per_trait populated ---------------------------------------

test_that("[lambda-dispatch] lambda_per_trait populated", {
  fx <- make_lambda_dispatch_fixture()
  fit <- fit_dispatch_baseline(fx, lambda_mode = "estimate")
  pd <- fit$pd; bl <- fit$bl

  expect_true(is.numeric(bl$lambda_per_trait))
  expect_length(bl$lambda_per_trait, ncol(pd$X_scaled))
  expect_identical(names(bl$lambda_per_trait), colnames(pd$X_scaled))

  cont_names <- colnames(pd$X_scaled)[c(pd$trait_map$c1$latent_cols,
                                         pd$trait_map$c2$latent_cols,
                                         pd$trait_map$c3$latent_cols)]
  # Per-trait lambda comes from ml_lambda_for_col() (R/bm_internal.R),
  # whose grid spans [0.005, 0.995] -- not [0.01, 0.99] (that narrower
  # interval is the SOLVER's own block-lambda optimize() range, a
  # different estimator; see R/joint_mvn_solver.R's .mvn_resolve_lambda()).
  lam_cont <- bl$lambda_per_trait[cont_names]
  expect_true(all(lam_cont >= 0.005 & lam_cont <= 0.995))
  expect_false(all(lam_cont == 1))

  disc_names <- setdiff(colnames(pd$X_scaled), cont_names)
  expect_true(all(bl$lambda_per_trait[disc_names] == 1))
})

# ---- covariates keep lambda --------------------------------------------

test_that("[lambda-dispatch] covariates keep lambda", {
  # Isolate the per-column covariate path (bm_impute_col_with_cov()): a
  # single continuous trait keeps bm_cols == 1, which is below
  # use_continuous_joint's >= 2 threshold and has no binary/ordinal for
  # use_threshold_joint, so no joint dispatch fires and no P1-8
  # "covariates ignored by the joint baseline" warning competes with the
  # assertion below.
  set.seed(902L)
  n <- 80L
  tree <- ape::rcoal(n)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  Lc <- chol(R)
  phylo_part <- as.numeric(t(Lc) %*% stats::rnorm(n))
  iid_part <- stats::rnorm(n)
  y <- sqrt(0.3) * phylo_part + sqrt(0.7) * iid_part
  df <- data.frame(c1 = y, row.names = tree$tip.label)
  df$c1[sample(n, round(0.3 * n))] <- NA
  covs <- data.frame(env1 = stats::rnorm(n), env2 = stats::rnorm(n),
                      row.names = tree$tip.label)
  pd <- preprocess_traits(df, tree, covariates = covs)
  spl <- make_missing_splits(pd$X_scaled, seed = 1, trait_map = pd$trait_map)

  bl <- NULL
  expect_no_warning(
    bl <- fit_baseline(pd, tree, splits = spl, lambda_mode = "estimate")
  )
  cont_name <- colnames(pd$X_scaled)[pd$trait_map$c1$latent_cols]
  expect_false(isTRUE(unname(bl$lambda_per_trait[cont_name]) == 1))
})

# ---- lambda_fixed rebuild -----------------------------------------------

test_that("[lambda-dispatch] lambda_fixed rebuild", {
  fx <- make_lambda_dispatch_fixture()
  fit <- fit_dispatch_baseline(fx, lambda_mode = "estimate")
  pd <- fit$pd; spl <- fit$splits; bl <- fit$bl

  bl_rebuilt <- fit_baseline(pd, fx$tree, splits = spl,
                              lambda_mode = "estimate",
                              lambda_fixed = bl$lambda_per_trait)
  expect_equal(bl_rebuilt$mu, bl$mu, tolerance = 1e-8)
  expect_equal(bl_rebuilt$se, bl$se, tolerance = 1e-8)
})

# ---- cv and bayes still per-column --------------------------------------

test_that("[lambda-dispatch] cv and bayes still per-column", {
  skip_on_cran()
  fx <- make_lambda_dispatch_fixture()
  fit_cv <- fit_dispatch_baseline(fx, lambda_mode = "cv")
  cont_names <- c("c1", "c2", "c3")
  expect_true(all(fit_cv$bl$path[cont_names] == "per_column_bm"))
})

# ---- fixed_1 dispatcher reference ---------------------------------------

test_that("[lambda-dispatch] fixed_1 dispatcher reference", {
  ref <- readRDS(testthat::test_path("fixtures",
                                      "lambda_fixed1_reference_ab02e31.rds"))
  pd <- preprocess_traits(ref$df, ref$tree)
  bl <- fit_baseline(pd, ref$tree, splits = ref$splits, lambda_mode = "fixed_1")

  expect_equal(bl$mu, ref$baseline$mu, tolerance = 1e-12)
  expect_equal(bl$se, ref$baseline$se, tolerance = 1e-12)
  expect_identical(bl$path, ref$baseline$path)
})
