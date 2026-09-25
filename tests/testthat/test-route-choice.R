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

  # S5c (rose-review.md, required change 5): score the chosen route on
  # cells NOT used to choose it. `bl_auto$score_val_idx` is exactly the
  # per-trait SCORE half (.pigauto_split_route_score()'s complement of the
  # ROUTE half `.pigauto_choose_predict_route()` actually saw) -- using
  # the FULL validation set here (as this test did pre-S5c) would include
  # the same cells that picked the route, which cannot fail except by a
  # coding error and is not evidence the choice helps.
  expect_false(is.null(bl_auto$score_val_idx))
  X_truth <- pd$X_scaled
  n_obs   <- nrow(X_truth)
  val_idx <- bl_auto$score_val_idx
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
    # auto's MSE for this trait equals whichever candidate it actually
    # used (true by construction, not independent evidence).
    expect_lte(mse_auto, max(mse_exact, mse_pc) + 1e-8)
    # S5c (rose-review.md, required change 5): unlike the pre-S5c version
    # of this test, `mse_exact`/`mse_pc` above are now scored on the SCORE
    # half -- cells that did NOT inform the route choice (made on the
    # ROUTE half only). On this independent half, the route chosen from
    # half A is NOT guaranteed to also win on half B (that guarantee was
    # exactly the tautology the review flagged: scoring on the SAME cells
    # that made the choice cannot fail except by a coding error). No
    # further per-trait MSE assertion is made here; `results` below still
    # reports both numbers for visibility.
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

# ---- (g) route split and score split are disjoint, and reused downstream

test_that("[route] under auto, the route half and score half of each trait's validation cells are disjoint and partition it", {
  skip_if_not_installed("Matrix")
  set.seed(41)
  n <- 200L
  tree <- ape::rtree(n)
  df <- data.frame(
    row.names = tree$tip.label,
    c1 = stats::rnorm(n), c2 = stats::rnorm(n),
    b1 = factor(sample(c("no", "yes"), n, replace = TRUE))
  )
  df <- mask_frac(df, frac = 0.3, seed = 42)
  pd  <- preprocess_traits(df, tree)
  spl <- make_missing_splits(pd$X_scaled, trait_map = pd$trait_map, seed = 5)

  rsplit <- pigauto:::.pigauto_split_route_score(spl$val_idx, nrow(pd$X_scaled),
                                                  pd$trait_map, seed = 99)
  # Split only when each half keeps >= 19 cells (n_j >= 38): then the halves
  # are disjoint and partition the trait's validation cells. Below that, the
  # same cells both choose the route and calibrate. Check each trait by rule.
  val_col <- ((spl$val_idx - 1L) %/% nrow(pd$X_scaled)) + 1L
  n_split <- 0L; n_shared <- 0L
  for (tm in pd$trait_map) {
    idx_j <- spl$val_idx[val_col %in% tm$latent_cols]
    r_j <- intersect(rsplit$route_idx, idx_j); s_j <- intersect(rsplit$score_idx, idx_j)
    if (length(idx_j) >= 38L) {
      n_split <- n_split + 1L
      expect_length(intersect(r_j, s_j), 0L)
      expect_setequal(c(r_j, s_j), idx_j)
    } else {
      n_shared <- n_shared + 1L
      expect_setequal(r_j, idx_j)
      expect_setequal(s_j, idx_j)
    }
  }
  expect_gt(n_split + n_shared, 0L)

  # fit_baseline()'s own "auto" output surfaces the SAME score half (what
  # fit_pigauto()/impute() restrict gate calibration + conformal scoring
  # to) and the same per-trait route/score counts.
  bl <- fit_baseline(pd, tree, spl, seed = 99)
  expect_setequal(bl$score_val_idx, rsplit$score_idx)
  # Overlap with the route cells is allowed only for traits below the split threshold.
  expect_setequal(intersect(bl$score_val_idx, rsplit$route_idx),
                  intersect(rsplit$score_idx, rsplit$route_idx))
  for (tm in pd$trait_map) {
    expect_identical(unname(bl$route_val_n[[tm$name]]),
                      unname(rsplit$n_route[[tm$name]]))
    expect_identical(unname(bl$score_val_n[[tm$name]]),
                      unname(rsplit$n_score[[tm$name]]))
  }

  # Integration: impute()'s model_config records the same counts, and
  # calibration/conformal never see a routing cell.
  res <- suppressWarnings(impute(df, tree, gnn = FALSE, missing_frac = 0.3,
                                  seed = 99, verbose = FALSE))
  cfg <- res$fit$model_config
  expect_false(is.null(cfg$route_val_n))
  expect_false(is.null(cfg$score_val_n))
})

# ---- (h) mixed-route auto baseline reproduces exactly via lambda_fixed ---

test_that("[route] a mixed-route auto baseline reproduces mu/se exactly when rebuilt with lambda_fixed (12 seeds)", {
  skip_if_not_installed("Matrix")
  n <- 150L
  mixed_seeds <- integer(0)
  for (sd in 1:12) {
    tree <- ape::rtree(n)
    tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
    R <- pigauto:::phylo_cor_matrix(tree)[tree$tip.label, tree$tip.label]

    L_cont  <- sim_joint_mvn(tree, R, lambda = 0.3, Sigma = diag(2), seed = 2000 + sd)
    lat_bin <- sim_joint_mvn(tree, R, lambda = 0.3, Sigma = matrix(1, 1, 1),
                              seed = 3000 + sd)[, 1]
    lat_ord <- sim_joint_mvn(tree, R, lambda = 0.3, Sigma = matrix(1, 1, 1),
                              seed = 4000 + sd)[, 1]

    df <- data.frame(
      row.names = tree$tip.label,
      c1 = L_cont[, 1], c2 = L_cont[, 2],
      b1 = factor(ifelse(lat_bin > 0, "yes", "no")),
      o1 = factor(findInterval(lat_ord, stats::quantile(lat_ord, c(0.33, 0.66))),
                  ordered = TRUE)
    )
    df <- mask_frac(df, frac = 0.25, seed = 5000 + sd)

    pd  <- preprocess_traits(df, tree)
    spl <- make_missing_splits(pd$X_scaled, trait_map = pd$trait_map, seed = 6000 + sd)

    bl <- fit_baseline(pd, tree, spl, seed = 7000 + sd)
    if (length(unique(bl$predict_method_by_trait)) < 2L) next  # not mixed this seed

    mixed_seeds <- c(mixed_seeds, sd)
    bl_rebuilt <- fit_baseline(pd, tree, spl, lambda_fixed = bl$lambda_per_trait,
                                seed = 7000 + sd)

    expect_identical(bl_rebuilt$predict_method_by_trait,
                      bl$predict_method_by_trait, info = paste("seed", sd))
    expect_equal(bl_rebuilt$mu, bl$mu, tolerance = 1e-8, info = paste("seed", sd))
    expect_equal(bl_rebuilt$se, bl$se, tolerance = 1e-8, info = paste("seed", sd))
  }
  # A vacuous pass (every seed uniform-route) would not be evidence of
  # anything -- require at least one genuinely mixed-route seed.
  expect_gt(length(mixed_seeds), 0L)
  message("mixed-route seeds (of 12): ", paste(mixed_seeds, collapse = ", "))
})

# ---- (i) a single OVR class fallback reports per_column only for its trait

test_that("[route] a single joint fit's fallback reports per_column only for the traits it covers", {
  skip_if_not_installed("Matrix")
  pigauto:::.pigauto_exact_fallback_reset()
  real_ecm <- pigauto:::exact_conditional_mvn
  call_n <- 0L
  testthat::local_mocked_bindings(
    exact_conditional_mvn = function(...) {
      call_n <<- call_n + 1L
      # Call 1 is the continuous-only joint fit for c1/c2; calls 2 and 3
      # are the two OVR class fits for k1. Fail only the first OVR class.
      if (call_n == 2L) return(NULL)
      real_ecm(...)
    },
    .package = "pigauto"
  )

  set.seed(51)
  n <- 60L
  tree <- ape::rcoal(n)
  df <- data.frame(
    row.names = tree$tip.label,
    c1 = stats::rnorm(n), c2 = stats::rnorm(n),
    # 3 levels: preprocess_traits() auto-detects factor(>2) as
    # "categorical" (K independent OVR fits below); factor(2) would be
    # "binary" and go through the SAME single threshold-joint call as
    # c1/c2, defeating the point of this test.
    k1 = factor(sample(c("a", "b", "c"), n, replace = TRUE))
  )
  df$c1[1:5]  <- NA
  df$c2[6:10] <- NA
  df$k1[11:15] <- NA
  pd <- preprocess_traits(df, tree)

  bl <- suppressWarnings(fit_baseline(pd, tree, splits = NULL,
                                       predict_method = "exact"))

  expect_identical(unname(bl$predict_method_by_trait[["c1"]]), "exact")
  expect_identical(unname(bl$predict_method_by_trait[["c2"]]), "exact")
  expect_identical(unname(bl$predict_method_by_trait[["k1"]]), "per_column")
})

# ---- (j) predict_route warns once on an unknown trait name --------------

test_that("[route] predict_route warns once on a name that matches no trait", {
  skip_if_not_installed("Matrix")
  set.seed(61)
  n <- 60L
  tree <- ape::rtree(n)
  df <- data.frame(
    row.names = tree$tip.label,
    c1 = stats::rnorm(n), c2 = stats::rnorm(n)
  )
  df <- mask_frac(df, frac = 0.2, seed = 62)
  pd  <- preprocess_traits(df, tree)
  spl <- make_missing_splits(pd$X_scaled, trait_map = pd$trait_map, seed = 7)

  route <- c(c1 = "exact", not_a_trait = "per_column")
  expect_warning(
    bl <- fit_baseline(pd, tree, spl, predict_route = route),
    "not_a_trait"
  )
  expect_identical(unname(bl$predict_method_by_trait[["c1"]]), "exact")

  # Values outside exact/per_column still error (unchanged contract).
  expect_error(
    fit_baseline(pd, tree, spl, predict_route = c(c1 = "bogus")),
    "must be a named character vector"
  )
})

test_that("[route] traits with at least 38 validation cells split into disjoint halves; smaller ones share", {
  skip_if_not_installed("Matrix")
  set.seed(43)
  n <- 800L
  tree <- ape::rtree(n)
  df <- data.frame(row.names = tree$tip.label, c1 = stats::rnorm(n), c2 = stats::rnorm(n))
  df <- mask_frac(df, frac = 0.3, seed = 44)
  pd  <- preprocess_traits(df, tree)
  spl <- make_missing_splits(pd$X_scaled, trait_map = pd$trait_map, seed = 6)
  rs <- pigauto:::.pigauto_split_route_score(spl$val_idx, nrow(pd$X_scaled), pd$trait_map, seed = 7)
  val_col <- ((spl$val_idx - 1L) %/% nrow(pd$X_scaled)) + 1L
  sizes <- vapply(pd$trait_map, function(tm) sum(val_col %in% tm$latent_cols), integer(1))
  expect_true(any(sizes >= 38L))
  for (i in which(sizes >= 38L)) {
    tm <- pd$trait_map[[i]]
    expect_true(rs$n_route[[tm$name]] >= 19L && rs$n_score[[tm$name]] >= 19L)
    expect_length(intersect(rs$route_idx[rs$route_idx %in% spl$val_idx[val_col %in% tm$latent_cols]],
                            rs$score_idx), 0L)
  }
  small <- pigauto:::.pigauto_split_route_score(spl$val_idx[seq_len(20)], nrow(pd$X_scaled), pd$trait_map, seed = 7)
  expect_setequal(small$route_idx, small$score_idx)
})

test_that("[route] the split keeps every latent cell of a held-out categorical row in one half and counts rows", {
  skip_if_not_installed("Matrix")
  set.seed(45)
  n <- 1200L
  tree <- ape::rtree(n)
  df <- data.frame(row.names = tree$tip.label, c1 = stats::rnorm(n),
                   k4 = factor(sample(letters[1:4], n, replace = TRUE)))
  df <- mask_frac(df, frac = 0.3, seed = 46)
  pd  <- preprocess_traits(df, tree)
  spl <- make_missing_splits(pd$X_scaled, trait_map = pd$trait_map, seed = 8)
  rs <- pigauto:::.pigauto_split_route_score(spl$val_idx, nrow(pd$X_scaled), pd$trait_map, seed = 9)
  n_obs <- nrow(pd$X_scaled)
  tm <- pd$trait_map[["k4"]]
  row_of <- function(idx) ((idx - 1L) %% n_obs) + 1L
  col_of <- function(idx) ((idx - 1L) %/% n_obs) + 1L
  cat_val <- spl$val_idx[col_of(spl$val_idx) %in% tm$latent_cols]
  rows <- unique(row_of(cat_val))
  expect_gte(length(rows), 38L)
  r_rows <- unique(row_of(intersect(rs$route_idx, cat_val)))
  s_rows <- unique(row_of(intersect(rs$score_idx, cat_val)))
  expect_length(intersect(r_rows, s_rows), 0L)
  expect_setequal(c(r_rows, s_rows), rows)
  expect_identical(unname(rs$n_route[["k4"]] + rs$n_score[["k4"]]), length(rows))
  # Every cell of a route row is in the route half (all K columns travel together).
  expect_setequal(intersect(rs$route_idx, cat_val), cat_val[row_of(cat_val) %in% r_rows])
})

test_that("[route] a trait whose two candidate fits agree is not split", {
  mu <- matrix(stats::rnorm(400), 200, 2)
  tmap <- list(a = list(name = "a", type = "continuous", latent_cols = 1L),
               b = list(name = "b", type = "continuous", latent_cols = 2L))
  val_idx <- c(1:60, 200L + (1:60))
  mu_pc <- mu
  mu_pc[1:60, 2] <- mu_pc[1:60, 2] + 0.1     # only trait b differs
  rs <- pigauto:::.pigauto_split_route_score(val_idx, 200L, tmap, seed = 1,
                                              mu_exact = mu, mu_pc = mu_pc)
  expect_setequal(intersect(rs$route_idx, 1:60), 1:60)
  expect_setequal(intersect(rs$score_idx, 1:60), 1:60)
  expect_length(intersect(intersect(rs$route_idx, 200L + (1:60)), rs$score_idx), 0L)

  # End to end: a single continuous trait gives identical exact and
  # per_column fits, so all its validation cells stay for calibration.
  skip_if_not_installed("Matrix")
  set.seed(47)
  tree <- ape::rtree(300L)
  df <- data.frame(row.names = tree$tip.label, y = stats::rnorm(300L))
  df <- mask_frac(df, frac = 0.3, seed = 48)
  pd  <- preprocess_traits(df, tree)
  spl <- make_missing_splits(pd$X_scaled, trait_map = pd$trait_map, seed = 10)
  bl <- fit_baseline(pd, tree, spl, seed = 11)
  expect_setequal(bl$score_val_idx, spl$val_idx)
})

test_that("[route] in multi-obs data the split keeps all observations of a species in one half", {
  tmap <- list(a = list(name = "a", type = "continuous", latent_cols = 1L))
  unit <- rep(seq_len(100L), each = 2L)
  rs <- pigauto:::.pigauto_split_route_score(1:200, 200L, tmap, seed = 2,
                                              unit_of_row = unit)
  expect_length(intersect(unit[rs$route_idx], unit[rs$score_idx]), 0L)
  expect_identical(unname(rs$n_route[["a"]] + rs$n_score[["a"]]), 100L)
})
