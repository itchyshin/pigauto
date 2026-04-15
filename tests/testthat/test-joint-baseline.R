test_that("joint_mvn_available() detects Rphylopars correctly", {
  # We expect TRUE in dev environments where Rphylopars is installed
  res <- joint_mvn_available()
  expect_type(res, "logical")
  expect_length(res, 1L)
  expect_equal(res, requireNamespace("Rphylopars", quietly = TRUE))
})

test_that("fit_joint_mvn_baseline recovers cross-trait structure on correlated BM", {
  skip_on_cran()
  skip_if_not(joint_mvn_available())
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
  skip_if_not(joint_mvn_available())
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

  # Tolerance reflects Rphylopars' REML solver vs our Cholesky nugget; they
  # should agree within ~1e-3 in practice for n=30.  SE may diverge slightly
  # more than mu because the two paths use different variance estimators
  # (profile REML vs plug-in σ²), so SE tolerance is kept at 1e-2.
  expect_equal(bl_new$mu, bl_old$mu, tolerance = 1e-3)
  expect_equal(bl_new$se, bl_old$se, tolerance = 1e-2)
})
