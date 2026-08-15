test_that("joint_mvn_available() reflects the in-house solver (always TRUE; no Rphylopars)", {
  # Default joint MVN path is R/joint_mvn_solver.R — always available.
  # Do not gate this on requireNamespace("Rphylopars").
  # 3c31f1d contract: identical TRUE (not requireNamespace("Rphylopars")).
  res <- joint_mvn_available()
  expect_type(res, "logical")
  expect_length(res, 1L)
  expect_identical(res, TRUE)
})

test_that("fit_joint_mvn_baseline recovers cross-trait structure on correlated BM", {
  skip_on_cran()
  set.seed(42)
  tree <- ape::rcoal(50)
  # Simulate two correlated BM traits: x1 ~ BM, x2 = 0.8 * x1 + noise (also BM)
  C <- ape::vcv(tree)
  L <- chol(C)
  x1 <- as.numeric(t(L) %*% stats::rnorm(50))
  x2 <- 0.8 * x1 + as.numeric(t(L) %*% stats::rnorm(50)) * 0.3
  df <- data.frame(t1 = x1, t2 = x2)
  rownames(df) <- tree$tip.label

  # Punch holes: 20% of t2 values missing, ALL t1 observed
  miss <- sample(50, 10)
  df$t2[miss] <- NA

  pd <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.1,
                                val_frac = 0.5, seed = 1,
                                trait_map = pd$trait_map)
  graph <- build_phylo_graph(tree, k_eigen = 4L)

  res <- fit_joint_mvn_baseline(pd, tree, splits = splits, graph = graph)

  # Shape
  expect_equal(dim(res$mu), c(50L, 2L))
  expect_equal(dim(res$se), c(50L, 2L))
  expect_false(any(is.na(res$mu[miss, "t2"])))

  # At observed, non-masked cells mu should equal the observed latent value.
  # Exclude cells that splits also masked (val/test mask can overlap observed rows).
  n <- nrow(pd$X_scaled)
  hold_idx   <- c(splits$val_idx, splits$test_idx)
  hold_rows  <- ((hold_idx - 1L) %% n) + 1L
  hold_cols  <- ((hold_idx - 1L) %/% n) + 1L
  masked_t2_rows <- hold_rows[hold_cols == which(colnames(pd$X_scaled) == "t2")]
  obs_unmasked   <- setdiff(seq_len(50L)[-miss], masked_t2_rows)
  expect_equal(res$mu[obs_unmasked, "t2"],
               pd$X_scaled[obs_unmasked, "t2"],
               tolerance = 1e-6)

  # Joint baseline should beat naive column-mean imputation on the held-out
  # cells because it can exploit the 0.8 correlation with t1.
  # Ground truth for miss cells: scale x2[miss] using the same mean/sd as pd.
  tm2 <- pd$trait_map[["t2"]]
  x2_true_latent <- (x2[miss] - tm2$mean) / tm2$sd
  col_mean_pred <- mean(pd$X_scaled[-miss, "t2"], na.rm = TRUE)
  joint_rmse <- sqrt(mean((res$mu[miss, "t2"] - x2_true_latent)^2))
  naive_rmse <- sqrt(mean((col_mean_pred       - x2_true_latent)^2))
  expect_lt(joint_rmse, naive_rmse)
})

test_that("joint baseline matches per-column BM on single-trait data", {
  skip_on_cran()
  set.seed(7)
  tree <- ape::rcoal(30)
  # Simulate a single BM trait directly via Cholesky of vcv(tree)
  C  <- ape::vcv(tree)
  L  <- chol(C)
  x1 <- as.numeric(t(L) %*% stats::rnorm(30))
  df <- data.frame(trait = x1)
  rownames(df) <- tree$tip.label
  df$trait[sample(30, 6)] <- NA

  pd     <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.1,
                                val_frac = 0.5, seed = 1,
                                trait_map = pd$trait_map)
  graph  <- build_phylo_graph(tree, k_eigen = 4L)

  # Path A: current per-column BM
  bl_old <- fit_baseline(pd, tree, splits = splits, graph = graph)

  # Path B: joint MVN (single trait should collapse to the same math)
  bl_new <- fit_joint_mvn_baseline(pd, tree, splits = splits, graph = graph)

  # Report max absolute diff for diagnostic purposes even when within tolerance.
  max_mu_diff <- max(abs(bl_new$mu - bl_old$mu), na.rm = TRUE)
  max_se_diff <- max(abs(bl_new$se - bl_old$se), na.rm = TRUE)
  message(sprintf("single-trait back-compat: max |mu diff| = %.2e, max |se diff| = %.2e",
                  max_mu_diff, max_se_diff))

  # Tolerance reflects in-house joint MVN vs per-column BM (Cholesky nugget);
  # they should agree within ~1e-3 in practice for n=30.  SE may diverge
  # slightly more than mu because the two paths use different variance
  # estimators, so SE tolerance is kept at 1e-2.
  expect_equal(bl_new$mu, bl_old$mu, tolerance = 1e-3)
  expect_equal(bl_new$se, bl_old$se, tolerance = 1e-2)
})

test_that("joint baseline does not leak val/test cells into BM fit", {
  skip_on_cran()
  set.seed(11)
  tree <- ape::rcoal(30)
  # Two independent BM traits (no cross-correlation — so the joint is just
  # two independent BMs, making the leak detection clean)
  C <- ape::vcv(tree)
  L <- chol(C)
  x1 <- as.numeric(t(L) %*% stats::rnorm(30))
  x2 <- as.numeric(t(L) %*% stats::rnorm(30))
  df <- data.frame(t1 = x1, t2 = x2)
  rownames(df) <- tree$tip.label

  pd <- preprocess_traits(df, tree)
  # Two splits that differ only in which cells are held out; the BM
  # prediction at cells held out in ONE but observed in the OTHER must
  # differ between the two runs.
  splits_a <- make_missing_splits(pd$X_scaled, missing_frac = 0.2,
                                  val_frac = 0.5, seed = 1,
                                  trait_map = pd$trait_map)
  splits_b <- make_missing_splits(pd$X_scaled, missing_frac = 0.2,
                                  val_frac = 0.5, seed = 99,
                                  trait_map = pd$trait_map)
  graph <- build_phylo_graph(tree, k_eigen = 4L)

  bl_a <- fit_joint_mvn_baseline(pd, tree, splits = splits_a, graph = graph)
  bl_b <- fit_joint_mvn_baseline(pd, tree, splits = splits_b, graph = graph)

  # Cells held out in splits_a but NOT in splits_b — at these cells
  # bl_a had to predict (observed cell was masked), but bl_b saw the
  # observed value, so predictions should differ (observed != smoothed).
  held_a <- c(splits_a$val_idx, splits_a$test_idx)
  held_b <- c(splits_b$val_idx, splits_b$test_idx)
  diff_cells <- setdiff(held_a, held_b)
  if (length(diff_cells) >= 3) {
    vals_a <- bl_a$mu[diff_cells]
    vals_b <- bl_b$mu[diff_cells]
    expect_gt(max(abs(vals_a - vals_b)), 1e-6)
  } else {
    # Rare — splits aligned too closely. Test is still valid but useless;
    # emit a skip with a note.
    skip("not enough differentiating cells between splits_a and splits_b")
  }
})

test_that("fit_baseline dispatches to joint MVN when trait count >= 2", {
  skip_on_cran()
  set.seed(17)
  tree <- ape::rcoal(30)
  # Two independent BM traits
  C <- ape::vcv(tree)
  L <- chol(C)
  x1 <- as.numeric(t(L) %*% stats::rnorm(30))
  x2 <- as.numeric(t(L) %*% stats::rnorm(30))
  df <- data.frame(t1 = x1, t2 = x2)
  rownames(df) <- tree$tip.label
  # Punch holes to make imputation meaningful
  df$t1[sample(30, 4)] <- NA
  df$t2[sample(30, 5)] <- NA

  pd <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.1,
                                val_frac = 0.5, seed = 1,
                                trait_map = pd$trait_map)
  graph <- build_phylo_graph(tree, k_eigen = 4L)

  # Call fit_baseline() — should internally use joint MVN path
  bl <- fit_baseline(pd, tree, splits = splits, graph = graph)
  expect_equal(dim(bl$mu), c(30L, 2L))
  expect_equal(dim(bl$se), c(30L, 2L))
  expect_false(any(is.na(bl$mu)))

  # Compare to the direct joint call — should be identical (same path)
  bl_direct <- fit_joint_mvn_baseline(pd, tree, splits = splits, graph = graph)
  expect_equal(bl$mu, bl_direct$mu, tolerance = 1e-8)
  expect_equal(bl$se, bl_direct$se, tolerance = 1e-8)
})

test_that("fit_baseline falls back cleanly when joint_mvn_available() is FALSE", {
  skip_on_cran()
  # Use testthat::local_mocked_bindings to shadow joint_mvn_available()
  # for this test only. testthat 3.1+ syntax. Production always returns
  # TRUE (in-house solver); this keeps the per-column fallback path honest.
  testthat::local_mocked_bindings(
    joint_mvn_available = function() FALSE,
    .package = "pigauto"
  )
  # Sanity check the mock took effect
  expect_false(joint_mvn_available())

  set.seed(19)
  tree <- ape::rcoal(30)
  C <- ape::vcv(tree)
  L <- chol(C)
  x1 <- as.numeric(t(L) %*% stats::rnorm(30))
  x2 <- as.numeric(t(L) %*% stats::rnorm(30))
  df <- data.frame(t1 = x1, t2 = x2)
  rownames(df) <- tree$tip.label
  df$t1[sample(30, 5)] <- NA

  pd <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.1,
                                val_frac = 0.5, seed = 1,
                                trait_map = pd$trait_map)
  graph <- build_phylo_graph(tree, k_eigen = 4L)

  # Should succeed without the joint path by falling back to per-column BM
  bl <- fit_baseline(pd, tree, splits = splits, graph = graph)
  expect_equal(dim(bl$mu), c(30L, 2L))
  expect_false(any(is.na(bl$mu)))
  expect_false(any(is.na(bl$se)))
})

# P0 B-Blk1: max_iter=0 Henderson path must return BM-style per-tip SE,
# not one empirical SD recycled for every missing tip.
test_that("K>=2 max_iter=0 missing-cell SE is per-tip BM-style, not one empirical SD", {
  skip_if_not_installed("Matrix")
  tree <- ape::read.tree(text = paste0(
    "(((close_obs:0.01,close_miss:0.01):0.04,(a:0.03,b:0.03):0.02):0.95,",
    "((far_miss:0.90,c:0.05):0.05,(d:0.4,e:0.4):0.1):0.05);"
  ))
  tips <- tree$tip.label
  n <- length(tips)
  # Two weakly correlated BM-ish columns; only t1 has the designed holes.
  set.seed(20260808L)
  t1 <- setNames(as.numeric(tips == "close_obs") * 2 + rnorm(n, sd = 0.15), tips)
  t2 <- setNames(rnorm(n), tips)
  L <- cbind(t1 = t1, t2 = t2)
  L[match(c("close_miss", "far_miss"), rownames(L)), "t1"] <- NA_real_

  fit <- fit_mvn_bm_inhouse(L, tree = tree, max_iter = 0L)
  se1 <- sqrt(fit$anc_var[, "t1"])
  se_close <- unname(se1["close_miss"])
  se_far   <- unname(se1["far_miss"])

  expect_true(is.finite(se_close) && se_close > 0)
  expect_true(is.finite(se_far) && se_far > 0)
  # Observed cells keep SE = 0.
  expect_equal(unname(se1[setdiff(tips, c("close_miss", "far_miss"))]),
               rep(0, n - 2L))
  # Heteroscedastic: closer-to-observed tip has strictly smaller SE.
  expect_lt(se_close, se_far)
  # Not a single recycled column SD (the pre-fix Henderson contract).
  expect_false(isTRUE(all.equal(se_close, se_far)))

  # Same ordering as dense BM σ√(1−h_i) on this column.
  R <- phylo_cor_matrix(tree)
  dense <- bm_impute_col(L[, "t1"], R)
  expect_lt(unname(dense$se[tips == "close_miss"]),
            unname(dense$se[tips == "far_miss"]))
})
