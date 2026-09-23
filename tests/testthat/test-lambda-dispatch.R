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

# ---- partial lambda_fixed defaults missing columns to 1 (Rose review) -----

# Required change 2 from docs/dev-log/lambda-default/rose-final-review.md:
# `fit_baseline()`'s roxygen for `lambda_fixed` says "Columns not present
# in `lambda_fixed` keep their lambda = 1 default", but a PARTIAL vector
# used to error ("numeric 'lambda' vector must have all entries in [0,
# 1]") because R/joint_mvn_baseline.R and R/joint_threshold_baseline.R
# both subset `lambda_fixed` by column name without filling in the
# missing names first, leaving NA. Threshold-joint variant (this fixture
# has a binary + categorical column, so `use_threshold_joint` fires):
test_that("[lambda-dispatch] partial lambda_fixed defaults missing columns to 1 (threshold-joint)", {
  fx <- make_lambda_dispatch_fixture()
  fit <- fit_dispatch_baseline(fx, lambda_mode = "estimate")
  pd <- fit$pd; spl <- fit$splits; bl <- fit$bl

  partial <- bl$lambda_per_trait["c1"]  # names only ONE of the three continuous columns

  bl_partial <- NULL
  expect_no_error(
    bl_partial <- fit_baseline(pd, fx$tree, splits = spl,
                                lambda_fixed = partial)
  )
  bl_fixed1 <- fit_baseline(pd, fx$tree, splits = spl, lambda_mode = "fixed_1")

  other_cols <- setdiff(colnames(pd$X_scaled), names(partial))
  expect_equal(bl_partial$mu[, other_cols], bl_fixed1$mu[, other_cols],
               tolerance = 1e-8)
  expect_equal(bl_partial$se[, other_cols], bl_fixed1$se[, other_cols],
               tolerance = 1e-8)
})

# Continuous-only variant (>= 2 continuous columns, no binary/ordinal),
# which dispatches to `use_continuous_joint` / `fit_joint_mvn_baseline()`
# instead -- the OTHER call site the review flagged
# (R/joint_mvn_baseline.R:~129).
test_that("[lambda-dispatch] partial lambda_fixed defaults missing columns to 1 (joint MVN)", {
  set.seed(903L)
  n <- 80L
  tree <- ape::rcoal(n)
  R <- stats::cov2cor(ape::vcv.phylo(tree))
  R <- R[tree$tip.label, tree$tip.label]
  Lc <- chol(R)
  gen_bm_lambda <- function(lam) {
    phylo_part <- as.numeric(t(Lc) %*% stats::rnorm(n))
    iid_part <- stats::rnorm(n)
    sqrt(lam) * phylo_part + sqrt(1 - lam) * iid_part
  }
  df <- data.frame(c1 = gen_bm_lambda(0.3), c2 = gen_bm_lambda(0.3),
                    row.names = tree$tip.label)
  df$c1[sample(n, round(0.3 * n))] <- NA
  df$c2[sample(n, round(0.3 * n))] <- NA

  pd  <- preprocess_traits(df, tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 1, trait_map = pd$trait_map)

  bl <- fit_baseline(pd, tree, splits = spl, lambda_mode = "estimate")
  expect_identical(unname(bl$path[c("c1", "c2")]), c("joint_mvn", "joint_mvn"))

  partial <- bl$lambda_per_trait["c1"]

  bl_partial <- NULL
  expect_no_error(
    bl_partial <- fit_baseline(pd, tree, splits = spl, lambda_fixed = partial)
  )
  bl_fixed1 <- fit_baseline(pd, tree, splits = spl, lambda_mode = "fixed_1")

  expect_equal(unname(bl_partial$mu[, "c2"]), unname(bl_fixed1$mu[, "c2"]),
               tolerance = 1e-8)
  expect_equal(unname(bl_partial$se[, "c2"]), unname(bl_fixed1$se[, "c2"]),
               tolerance = 1e-8)
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

# ---- ordinal stays at lambda = 1 under estimate (Rose review, 2026-09-23) --

# 2 continuous BM(lambda = 0.3) traits + 1 three-level ordinal trait, all
# correlated with the SAME tree, 30% NA per column. Small enough (K = 3)
# that R/fit_baseline.R's "Per-trait ordinal path selection" block
# (~line 606) exercises all three candidates (threshold_joint, the
# BM-via-MVN alternative, and the K-class-OVR LP alternative) across
# seeds, matching the regime Rose's review reproduced the leak in
# (docs/dev-log/lambda-default/rose-final-review.md, section 2, "ord.R").
make_ordinal_lambda_fixture <- function(seed, n = 100L, lambda_true = 0.3,
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

  lat_o <- as.numeric(t(Lc) %*% stats::rnorm(n))
  brks  <- stats::quantile(lat_o, c(0, 1 / 3, 2 / 3, 1))
  o1 <- factor(cut(lat_o, brks, labels = c("1", "2", "3"),
                    include.lowest = TRUE), ordered = TRUE)

  df <- data.frame(c1 = c1, c2 = c2, o1 = o1, row.names = tree$tip.label)
  for (col in c("c1", "c2", "o1")) {
    df[sample(n, round(miss_frac * n)), col] <- NA
  }
  list(tree = tree, df = df)
}

test_that("[lambda-dispatch] ordinal stays at lambda = 1 under estimate", {
  for (seed in 900:905) {
    fx <- make_ordinal_lambda_fixture(seed)
    pd <- preprocess_traits(fx$df, fx$tree)
    spl <- make_missing_splits(pd$X_scaled, seed = seed,
                                trait_map = pd$trait_map)

    bl_fixed <- fit_baseline(pd, fx$tree, splits = spl, lambda_mode = "fixed_1")
    bl_est   <- fit_baseline(pd, fx$tree, splits = spl, lambda_mode = "estimate")

    o1_col <- colnames(pd$X_scaled)[pd$trait_map$o1$latent_cols]

    # Documented contract (R/fit_baseline.R roxygen, "Per-type lambda
    # dispatch"): ordinal columns stay at lambda = 1 in EVERY path under
    # "estimate", so their mu/se must be bit-identical to the fixed_1 fit
    # -- not just close, and `lambda_per_trait["o1"]` must actually equal
    # the lambda that was used (1), not just report 1 while a different
    # lambda leaked into mu (the bug this test guards against).
    expect_identical(unname(bl_fixed$mu[, o1_col]), unname(bl_est$mu[, o1_col]),
                      info = paste("seed", seed))
    expect_identical(unname(bl_fixed$se[, o1_col]), unname(bl_est$se[, o1_col]),
                      info = paste("seed", seed))
    expect_true(unname(bl_est$lambda_per_trait[o1_col]) == 1,
                info = paste("seed", seed))
  }
})

test_that("[lambda-dispatch] a fully observed column reports its estimated lambda, not 1", {
  set.seed(3)
  tr <- ape::rcoal(60)
  P <- matrix(rexp(180), 60, 3); P <- P / rowSums(P)
  df <- data.frame(a = P[, 1], b = P[, 2], c = P[, 3], x = rnorm(60),
                   row.names = tr$tip.label)
  df[1:8, c("a", "b", "c")] <- NA
  pd <- preprocess_traits(df, tr, multi_proportion_groups = list(comp = c("a", "b", "c")))
  bl <- fit_baseline(pd, tr, lambda_mode = "estimate")
  # x is iid noise with no missing cells: its lambda estimate sits at the lower bound
  expect_lt(unname(bl$lambda_per_trait["x"]), 0.5)
  bl1 <- fit_baseline(pd, tr, lambda_mode = "fixed_1")
  expect_equal(unname(bl1$lambda_per_trait["x"]), 1)
})
