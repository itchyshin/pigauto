# Tests for arc/lambda-per-type: per-trait-type lambda_mode dispatch in
# fit_baseline().
#
# Context: before this change, `lambda_mode != "fixed_1"` set
# `force_per_column <- TRUE`, which disabled the joint MVN, threshold-joint,
# AND OVR-categorical baselines for ALL traits -- not just the BM-eligible
# continuous-family columns the lambda machinery actually applies to.
# Measured cost on AVONET (script/bench_avonet_lambda_modes.md on branch
# arc/bace-comparators; docs/dev-log/2026-08-16-external-comparison-results.md
# on the handover branch): lambda_mode = "bayes" improves every continuous
# trait but drops Trophic.Level (categorical) accuracy from 0.789 to 0.600
# because categorical falls from the joint/OVR baseline to plain label
# propagation. Discrete traits have no lambda concept and were always fit
# at lambda = 1 regardless of `lambda_mode`.
#
# These tests confirm: (1) `lambda_mode = "fixed_1"` is byte-identical to
# pre-change behaviour (numerical regression pin, verified independently
# against origin/main's R/fit_baseline.R during implementation); (2) for
# `lambda_mode != "fixed_1"`, binary/ordinal/categorical columns still come
# from the threshold-joint / OVR-categorical baseline (identical to the
# fixed_1 run, and distinguishable from a pure label-propagation reference);
# (3) continuous-family columns route through the per-column
# `bm_impute_col(..., lambda = bm_lambda)` path instead.
#
# `joint_mvn_available()` (R/joint_mvn_baseline.R) returns TRUE
# unconditionally -- the joint MVN / threshold-joint / OVR-categorical
# baselines use the in-house solver (R/joint_mvn_solver.R) by default
# (`joint_solver = "inhouse"`, the default used throughout this file) and
# have no Rphylopars dependency. Rphylopars is only needed when a caller
# explicitly passes `joint_solver = "rphylopars"`, which none of these
# tests do, so no `skip_if_not_installed("Rphylopars")` guard is required.

make_mixed_data_lpt <- function(n = 30, seed = 100) {
  set.seed(seed)
  tree <- ape::rtree(n)
  sp   <- tree$tip.label
  df <- data.frame(
    row.names = sp,
    mass = abs(stats::rnorm(n)) + 0.5,               # continuous
    wing = abs(stats::rnorm(n)) * 2 + 1,              # continuous
    migr = factor(sample(c("no", "yes"), n, replace = TRUE)),          # binary
    diet = factor(sample(c("herb", "carn", "omni"), n, replace = TRUE)) # categorical
  )
  list(tree = tree, df = df)
}

# Replicates fit_baseline()'s phylogenetic label-propagation formula
# (R/fit_baseline.R, binary/categorical LP blocks) independently, as a
# "pure LP" oracle to confirm the joint/OVR baseline -- not LP -- produced
# the discrete-trait output under lambda_mode != "fixed_1".
lp_reference_binary <- function(tree, X_masked, col) {
  D <- ape::cophenetic.phylo(tree)
  D <- D[rownames(X_masked), rownames(X_masked)]
  sigma_lp <- stats::median(D) * 0.5
  sim <- exp(-(D^2) / (2 * sigma_lp^2))
  diag(sim) <- 0
  vals <- X_masked[, col]
  observed <- !is.na(vals)
  sim_obs <- sim[, observed, drop = FALSE]
  rw <- rowSums(sim_obs)
  rw[rw < 1e-10] <- 1e-10
  probs <- as.numeric(sim_obs %*% vals[observed]) / rw
  probs <- pmin(pmax(probs, 0.01), 0.99)
  pigauto:::logit(probs)
}

# ---- fixed_1 regression pin -----------------------------------------------

test_that("[lambda-per-type] lambda_mode = 'fixed_1' is unchanged (numerical regression)", {
  td  <- make_mixed_data_lpt()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 1, trait_map = pd$trait_map)

  bl <- fit_baseline(pd, td$tree, splits = spl, lambda_mode = "fixed_1")

  expect_true(is.matrix(bl$mu))
  expect_true(all(is.finite(bl$mu)))
  expect_equal(dim(bl$mu), dim(pd$X_scaled))

  mass_col <- pd$trait_map$mass$latent_cols
  migr_col <- pd$trait_map$migr$latent_cols
  diet_cols <- pd$trait_map$diet$latent_cols

  # Pinned values computed from this same call; a change here means
  # lambda_mode = "fixed_1" no longer reproduces pre-change fit_baseline()
  # output. Verified independently against origin/main's R/fit_baseline.R
  # (byte-identical mu/se) during implementation of arc/lambda-per-type.
  expect_equal(
    unname(bl$mu[1:5, mass_col]),
    c(-0.22414418, 0.56787943, -2.29468316, 1.85821583, 1.42645745),
    tolerance = 1e-6
  )
  expect_equal(
    unname(bl$mu[1:5, migr_col]),
    c(-1.31010818, 1.31010818, -1.31010818, 1.31010818, 1.31010818),
    tolerance = 1e-6
  )
  expect_equal(
    unname(bl$mu[1:3, diet_cols]),
    matrix(c(-1.74161896, -0.43151078, -0.43151078,
             -1.74161896, -1.74161896, -1.74161896,
             -0.43151078, -1.74161896, -1.74161896), nrow = 3),
    tolerance = 1e-6
  )
})

# ---- lambda_mode != "fixed_1": discrete traits keep joint/OVR -------------

test_that("[lambda-per-type] lambda_mode = 'estimate' leaves binary/categorical baseline identical to fixed_1", {
  td  <- make_mixed_data_lpt()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 1, trait_map = pd$trait_map)

  bl_fixed <- fit_baseline(pd, td$tree, splits = spl, lambda_mode = "fixed_1")
  bl_est   <- fit_baseline(pd, td$tree, splits = spl, lambda_mode = "estimate")

  migr_col  <- pd$trait_map$migr$latent_cols
  diet_cols <- pd$trait_map$diet$latent_cols

  # The threshold-joint / OVR-categorical fits never see bm_lambda (the
  # underlying solvers have no lambda argument), so their binary/
  # categorical output must be bit-identical regardless of lambda_mode.
  expect_equal(bl_est$mu[, migr_col], bl_fixed$mu[, migr_col])
  expect_equal(bl_est$mu[, diet_cols], bl_fixed$mu[, diet_cols])
  expect_equal(bl_est$se[, migr_col], bl_fixed$se[, migr_col])
  expect_equal(bl_est$se[, diet_cols], bl_fixed$se[, diet_cols])

  # And it is NOT the plain label-propagation result -- confirms the
  # joint/OVR baseline actually fired rather than silently falling back to
  # LP (the pre-fix behaviour for lambda_mode != "fixed_1").
  X <- pd$X_scaled
  X[spl$val_idx]  <- NA
  X[spl$test_idx] <- NA
  lp_migr <- lp_reference_binary(td$tree, X, migr_col)
  expect_false(isTRUE(all.equal(unname(bl_est$mu[, migr_col]), unname(lp_migr),
                                 tolerance = 1e-4)))
})

test_that("[lambda-per-type] lambda_mode = 'bayes' also leaves categorical baseline identical to fixed_1", {
  td  <- make_mixed_data_lpt()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 1, trait_map = pd$trait_map)

  bl_fixed <- fit_baseline(pd, td$tree, splits = spl, lambda_mode = "fixed_1")
  bl_bayes <- fit_baseline(pd, td$tree, splits = spl, lambda_mode = "bayes")

  diet_cols <- pd$trait_map$diet$latent_cols
  migr_col  <- pd$trait_map$migr$latent_cols
  expect_equal(bl_bayes$mu[, diet_cols], bl_fixed$mu[, diet_cols])
  expect_equal(bl_bayes$mu[, migr_col], bl_fixed$mu[, migr_col])
})

# ---- lambda_mode != "fixed_1": continuous-family takes the lambda path ----

test_that("[lambda-per-type] lambda_mode = 'estimate' routes continuous columns through the per-column lambda path", {
  td  <- make_mixed_data_lpt()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 1, trait_map = pd$trait_map)

  bl_fixed <- fit_baseline(pd, td$tree, splits = spl, lambda_mode = "fixed_1")
  bl_est   <- fit_baseline(pd, td$tree, splits = spl, lambda_mode = "estimate")

  mass_col <- pd$trait_map$mass$latent_cols
  wing_col <- pd$trait_map$wing$latent_cols

  # lambda_mode actually changed the continuous-family baseline (it did
  # NOT silently keep the fixed_1 joint-fit output).
  expect_false(isTRUE(all.equal(bl_est$mu[, mass_col], bl_fixed$mu[, mass_col])))
  expect_false(isTRUE(all.equal(bl_est$mu[, wing_col], bl_fixed$mu[, wing_col])))

  # And it matches a direct per-column bm_impute_col(lambda = "estimate")
  # call on the same masked column + phylogenetic correlation matrix that
  # fit_baseline() uses internally -- i.e. continuous columns went through
  # the lambda-aware per-column BM path, not the (lambda-blind) joint fit.
  X <- pd$X_scaled
  X[spl$val_idx]  <- NA
  X[spl$test_idx] <- NA
  R_phy <- pigauto:::phylo_cor_matrix(td$tree)[rownames(X), rownames(X)]

  direct_mass <- pigauto:::bm_impute_col(X[, mass_col], R_phy, lambda = "estimate")
  direct_wing <- pigauto:::bm_impute_col(X[, wing_col], R_phy, lambda = "estimate")

  expect_equal(unname(bl_est$mu[, mass_col]), unname(direct_mass$mu), tolerance = 1e-8)
  expect_equal(unname(bl_est$se[, mass_col]), unname(direct_mass$se), tolerance = 1e-8)
  expect_equal(unname(bl_est$mu[, wing_col]), unname(direct_wing$mu), tolerance = 1e-8)
})
