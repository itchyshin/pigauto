test_that("fit_ovr_categorical_fits returns K probability vectors per tip", {
  skip_if_not_installed("Rphylopars")
  set.seed(100)
  tree <- ape::rtree(30)
  df <- data.frame(
    x = rnorm(30),
    z = factor(sample(c("P", "Q", "R"), 30, TRUE), levels = c("P", "Q", "R")),
    row.names = tree$tip.label
  )
  df$z[c(5, 15, 25)] <- NA
  pd <- preprocess_traits(df, tree)

  probs <- fit_ovr_categorical_fits(pd, tree, trait_name = "z", splits = NULL)

  expect_equal(dim(probs), c(30, 3))
  expect_equal(colnames(probs), c("z=P", "z=Q", "z=R"))
  expect_true(all(is.finite(probs) | is.na(probs)))
  finite_vals <- probs[!is.na(probs) & is.finite(probs)]
  expect_true(all(finite_vals >= 0 & finite_vals <= 1))
})

test_that("fit_baseline routes categorical through OVR K-fits when Rphylopars available", {
  skip_if_not_installed("Rphylopars")
  set.seed(200)
  tree <- ape::rtree(40)
  df <- data.frame(
    x = rnorm(40),
    y = factor(sample(c("A","B"), 40, TRUE), levels = c("A","B")),
    z = factor(sample(c("P","Q","R","S"), 40, TRUE),
               levels = c("P","Q","R","S")),
    row.names = tree$tip.label
  )
  df$z[c(5, 15, 25)] <- NA
  pd <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.2,
                                 seed = 200, trait_map = pd$trait_map)
  bl <- fit_baseline(pd, tree, splits = splits)

  cat_col_names <- grep("^z=", colnames(pd$X_scaled), value = TRUE)
  expect_length(cat_col_names, 4L)
  row_probs_sums <- rowSums(exp(bl$mu[, cat_col_names]))
  expect_true(all(abs(row_probs_sums - 1) < 1e-6))
  expect_true(all(is.finite(bl$mu[, cat_col_names])))
})

test_that("decode_ovr_categorical normalises K probs to valid log-prob distribution", {
  probs <- matrix(c(0.8, 0.3, 0.1,
                     0.2, 0.9, 0.4,
                     0.1, 0.2, 0.7), nrow = 3, byrow = TRUE)
  log_probs <- decode_ovr_categorical(probs)
  expect_true(all(abs(rowSums(exp(log_probs)) - 1) < 1e-6))
  expect_equal(which.max(log_probs[1, ]), 1L)

  # All-NA row gets uniform
  probs_na <- probs; probs_na[1, ] <- NA
  lp <- decode_ovr_categorical(probs_na)
  expect_equal(lp[1, ], rep(log(1/3), 3), tolerance = 1e-9)
  expect_true(all(abs(rowSums(exp(lp)) - 1) < 1e-6))

  # Single-NA col gets residual
  probs_na2 <- probs; probs_na2[2, 2] <- NA
  lp2 <- decode_ovr_categorical(probs_na2)
  expect_true(all(is.finite(lp2)))
  expect_true(all(abs(rowSums(exp(lp2)) - 1) < 1e-6))
})
