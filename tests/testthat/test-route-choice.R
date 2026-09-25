# predict_method = "auto" per-trait route choice (S5b, docs/dev-log/
# exact-default/S5b-route-choice-report.md). "auto" fits the baseline with
# both the "exact" and "per_column" routes and picks, per trait, whichever
# has the lower loss on that trait's validation cells.

# ---- shared DGP helper ------------------------------------------------
# Simulate K columns from vec(L) ~ MVN(0, Sigma %x% R(lambda)), the exact
# route's own generative assumption: R(lambda)[i, j] = lambda * R[i, j] for
# i != j, 1 on the diagonal (Pagel's lambda on the correlation scale); Sigma
# is the K x K cross-trait covariance (compound symmetry via `rho`).
sim_joint_mvn <- function(tree, R, lambda, Sigma, seed) {
  set.seed(seed)
  n <- length(tree$tip.label)
  K <- nrow(Sigma)
  R_lambda <- lambda * R
  diag(R_lambda) <- 1
  Lr <- chol(R_lambda)   # Lr^T Lr = R_lambda
  Lc <- chol(Sigma)      # Lc^T Lc = Sigma
  Z <- matrix(stats::rnorm(n * K), n, K)
  L <- t(Lr) %*% Z %*% Lc
  dimnames(L) <- list(tree$tip.label, paste0("t", seq_len(K)))
  L
}

mask_frac <- function(df, frac, seed) {
  set.seed(seed)
  n <- nrow(df)
  for (j in seq_along(df)) df[sample(n, floor(frac * n)), j] <- NA
  df
}

# ---- (a) defaults ------------------------------------------------------

test_that("[route] impute()/fit_pigauto()/fit_baseline() default to predict_method = 'auto'", {
  expect_identical(eval(formals(fit_baseline)$predict_method)[1], "auto")
  expect_identical(eval(formals(fit_pigauto)$predict_method)[1], "auto")
  expect_identical(eval(formals(impute)$predict_method)[1], "auto")
})

# ---- (b) exact clearly wins: correlated continuous traits --------------

test_that("[route] auto chooses exact for all continuous traits when exact clearly wins", {
  skip_if_not_installed("Matrix")
  set.seed(42)
  n <- 200L
  tree <- ape::rtree(n)
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  R <- pigauto:::phylo_cor_matrix(tree)[tree$tip.label, tree$tip.label]
  K <- 3L
  Sigma <- matrix(0.6, K, K); diag(Sigma) <- 1
  L <- sim_joint_mvn(tree, R, lambda = 0.7, Sigma = Sigma, seed = 100)
  df <- mask_frac(as.data.frame(L), frac = 0.2, seed = 101)

  pd   <- preprocess_traits(df, tree)
  spl  <- make_missing_splits(pd$X_scaled, trait_map = pd$trait_map, seed = 1)
  bl   <- fit_baseline(pd, tree, spl)

  expect_identical(bl$predict_method_used, "auto")
  expect_true(all(bl$predict_method_by_trait == "exact"))
})

# ---- (c) independent traits at true lambda = 1 + low-signal count ------

test_that("[route] auto chooses per_column for at least one continuous trait when a low-signal count pulls lambda_block below 1, without hurting its validation MSE", {
  skip_if_not_installed("Matrix")
  set.seed(7)
  n <- 300L
  tree <- ape::rtree(n)
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  R <- pigauto:::phylo_cor_matrix(tree)[tree$tip.label, tree$tip.label]

  # c1, c2: independent (Sigma = I), true lambda = 1 (full phylogenetic signal).
  L_cont <- sim_joint_mvn(tree, R, lambda = 1.0, Sigma = diag(2), seed = 200)
  # cnt: near-white-noise phylogenetic signal (lambda = 0.05), exponentiated
  # and rounded to non-negative counts so preprocess_traits() auto-detects
  # it as a "count" trait. Low signal here is what pulls the joint fit's
  # single shared lambda_block below 1 (the mechanism reported in
  # docs/dev-log/exact-default/S5-benchmark-round1.md).
  cnt_L <- sim_joint_mvn(tree, R, lambda = 0.05, Sigma = matrix(1, 1, 1), seed = 201)
  cnt_vals <- pmax(0L, round(exp(0.3 * cnt_L[, 1] + 1)))

  df <- data.frame(c1 = L_cont[, 1], c2 = L_cont[, 2], cnt = as.integer(cnt_vals),
                    row.names = tree$tip.label)
  df <- mask_frac(df, frac = 0.2, seed = 202)

  pd  <- preprocess_traits(df, tree)
  spl <- make_missing_splits(pd$X_scaled, trait_map = pd$trait_map, seed = 2)

  bl_auto  <- fit_baseline(pd, tree, spl)
  bl_exact <- fit_baseline(pd, tree, spl, predict_method = "exact")
  bl_pc    <- fit_baseline(pd, tree, spl, predict_method = "per_column")

  expect_identical(bl_exact$lambda_block < 1, TRUE)  # the low-signal count pulls it down

  cont_names <- c("c1", "c2")
  chosen <- bl_auto$predict_method_by_trait[cont_names]
  expect_true(any(chosen == "per_column"))

  # Per-trait validation MSE on the z-scored latent scale, using the same
  # linear-index decode fit_baseline() uses internally for splits$val_idx.
  X_truth <- pd$X_scaled
  n_obs   <- nrow(X_truth)
  val_idx <- spl$val_idx
  val_col <- ((val_idx - 1L) %/% n_obs) + 1L
  val_row <- ((val_idx - 1L) %% n_obs) + 1L

  mse_for <- function(bl, tm) {
    col <- tm$latent_cols[1]
    keep <- val_col == col
    vr <- val_row[keep]
    truth <- X_truth[vr, col]
    pred  <- bl$mu[vr, col]
    ok <- is.finite(truth) & is.finite(pred)
    mean((pred[ok] - truth[ok])^2)
  }

  results <- character(0)
  for (nm in cont_names) {
    tm <- pd$trait_map[[nm]]
    mse_exact <- mse_for(bl_exact, tm)
    mse_pc    <- mse_for(bl_pc, tm)
    mse_auto  <- mse_for(bl_auto, tm)
    results <- c(results, sprintf("%s: chosen=%s exact=%.4f per_column=%.4f",
                                   nm, chosen[[nm]], mse_exact, mse_pc))
    # auto's MSE for this trait must not be worse than the route it did NOT
    # pick, i.e. it must equal whichever of the two is lower.
    expect_lte(mse_auto, max(mse_exact, mse_pc) + 1e-8)
    if (identical(chosen[[nm]], "per_column")) {
      expect_lte(mse_pc, mse_exact + 1e-8)
    }
  }
  # Report which traits chose what and their MSEs (visible with
  # testthat::test_file(..., reporter = "summary") / -v).
  message(paste(results, collapse = "; "))
})

# ---- (d) predict_route forces the route -------------------------------

test_that("[route] predict_route forces the route and reproduces a forced-route fit within 1e-10", {
  skip_if_not_installed("Matrix")
  set.seed(11)
  n <- 60L
  tree <- ape::rtree(n)
  df <- data.frame(
    row.names = tree$tip.label,
    c1 = stats::rnorm(n), c2 = stats::rnorm(n),
    b1 = factor(sample(c("no", "yes"), n, replace = TRUE)),
    k1 = factor(sample(letters[1:3], n, replace = TRUE))
  )
  df <- mask_frac(df, frac = 0.2, seed = 12)
  pd  <- preprocess_traits(df, tree)
  spl <- make_missing_splits(pd$X_scaled, trait_map = pd$trait_map, seed = 3)

  trait_names <- vapply(pd$trait_map, function(tm) tm$name, character(1))

  route_all_exact <- stats::setNames(rep("exact", length(trait_names)), trait_names)
  bl_forced_exact <- fit_baseline(pd, tree, spl, predict_route = route_all_exact)
  bl_exact        <- fit_baseline(pd, tree, spl, predict_method = "exact")
  expect_equal(bl_forced_exact$mu, bl_exact$mu, tolerance = 1e-10)
  expect_equal(bl_forced_exact$se, bl_exact$se, tolerance = 1e-10)
  expect_identical(bl_forced_exact$predict_method_by_trait, route_all_exact)

  route_all_pc <- stats::setNames(rep("per_column", length(trait_names)), trait_names)
  bl_forced_pc <- fit_baseline(pd, tree, spl, predict_route = route_all_pc)
  bl_pc        <- fit_baseline(pd, tree, spl, predict_method = "per_column")
  expect_equal(bl_forced_pc$mu, bl_pc$mu, tolerance = 1e-10)
  expect_equal(bl_forced_pc$se, bl_pc$se, tolerance = 1e-10)

  # Mixed route: c1/b1/k1 -> exact, c2 -> per_column. Each trait's columns
  # must come from the matching single-route fit exactly.
  route_mixed <- route_all_exact
  route_mixed[["c2"]] <- "per_column"
  bl_mixed <- fit_baseline(pd, tree, spl, predict_route = route_mixed)
  c2_col <- pd$trait_map$c2$latent_cols
  c1_col <- pd$trait_map$c1$latent_cols
  expect_equal(bl_mixed$mu[, c2_col], bl_pc$mu[, c2_col], tolerance = 1e-10)
  expect_equal(bl_mixed$mu[, c1_col], bl_exact$mu[, c1_col], tolerance = 1e-10)
  expect_identical(bl_mixed$predict_method_by_trait, route_mixed)
})

# ---- (e) production refit reuses the validation-split choice -----------

test_that("[route] the production refit in fit_pigauto() reuses the validation-split route choice", {
  skip_if_not_installed("Matrix")
  set.seed(21)
  n <- 200L
  tree <- ape::rtree(n)
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  R <- pigauto:::phylo_cor_matrix(tree)[tree$tip.label, tree$tip.label]
  L_cont <- sim_joint_mvn(tree, R, lambda = 1.0, Sigma = diag(2), seed = 220)
  cnt_L <- sim_joint_mvn(tree, R, lambda = 0.05, Sigma = matrix(1, 1, 1), seed = 221)
  cnt_vals <- pmax(0L, round(exp(0.3 * cnt_L[, 1] + 1)))
  df <- data.frame(c1 = L_cont[, 1], c2 = L_cont[, 2], cnt = as.integer(cnt_vals),
                    row.names = tree$tip.label)
  df <- mask_frac(df, frac = 0.2, seed = 222)

  res <- suppressWarnings(impute(df, tree, gnn = FALSE, missing_frac = 0.25,
                                  seed = 23, verbose = FALSE))
  fit <- res$fit

  expect_identical(fit$model_config$predict_method_used, "auto")
  expect_identical(fit$model_config$predict_method_by_trait,
                    fit$baseline$predict_method_by_trait)
  # baseline_full (splits = NULL, no validation cells of its own) must have
  # REUSED baseline's (splits-based, real validation cells) per-trait choice
  # rather than re-deciding with zero evidence (which would default every
  # trait to "exact").
  expect_identical(fit$baseline_full$predict_method_by_trait,
                    fit$baseline$predict_method_by_trait)
})

# ---- (f) no splits -> every trait "exact" -------------------------------

test_that("[route] splits = NULL gives every trait 'exact'", {
  skip_if_not_installed("Matrix")
  set.seed(31)
  n <- 60L
  tree <- ape::rtree(n)
  df <- data.frame(
    row.names = tree$tip.label,
    c1 = stats::rnorm(n), c2 = stats::rnorm(n),
    b1 = factor(sample(c("no", "yes"), n, replace = TRUE))
  )
  df <- mask_frac(df, frac = 0.2, seed = 32)
  pd <- preprocess_traits(df, tree)

  bl <- fit_baseline(pd, tree, splits = NULL)
  expect_identical(bl$predict_method_used, "auto")
  expect_true(all(bl$predict_method_by_trait == "exact"))

  # Same for an empty-val-cells splits object (val_idx length 0).
  spl_empty <- make_missing_splits(pd$X_scaled, trait_map = pd$trait_map, seed = 4)
  spl_empty$val_idx <- integer(0)
  bl2 <- fit_baseline(pd, tree, spl_empty)
  expect_identical(bl2$predict_method_used, "auto")
  expect_true(all(bl2$predict_method_by_trait == "exact"))
})
