test_that("pooling requires two fits and preserves Rubin arithmetic at that boundary", {
  first <- stats::lm(y ~ 1, data = data.frame(y = c(1, 2, 4, 5)))
  second <- stats::lm(y ~ 1, data = data.frame(y = c(2, 4, 6, 8)))
  expect_error(pool_mi(list(first)), "Need at least 2 fits to pool; got 1")
  expect_warning(pooled <- pool_mi(list(first, second)), "unverified provenance")
  # Means 3 and 5, within variances 5/6 and 5/3, between variance 2.
  expect_equal(pooled$estimate, 4)
  expect_equal(pooled$std.error, sqrt(5 / 4 + 3))
  expect_identical(attr(pooled, "m"), 2L)
})

test_that("GNN-off conformal metadata stays NULL for its supported methods", {
  set.seed(883)
  tree <- ape::rtree(60)
  traits <- data.frame(x = sin(seq_len(60)), row.names = tree$tip.label)
  data <- preprocess_traits(traits, tree)
  splits <- make_missing_splits(data$X_scaled, missing_frac = .5,
    val_frac = .8, trait_map = data$trait_map, seed = 22)
  expect_error(fit_pigauto(data, tree, splits = splits, gnn = FALSE,
    conformal_method = "mondrian", verbose = FALSE),
    "not supported when gnn = FALSE")
  for (method in c("split", "bootstrap")) {
    fit <- fit_pigauto(data, tree, splits = splits, gnn = FALSE,
      conformal_method = method, conformal_split_val = FALSE,
      conformal_bootstrap_B = 20L, verbose = FALSE, seed = 31)
    expect_identical(fit$conformal_method, method)
    expect_null(fit$conformal_mondrian)
    expect_length(fit$conformal_scores, 1L)
    expect_true(all(is.finite(unlist(fit$conformal_scores))))
    expect_length(fit$r_cal_gnn, ncol(data$X_scaled))
    expect_true(all(fit$r_cal_gnn == 0))
  }
})
