test_that("build_liability_matrix fills observed binary cells via truncated-Gaussian E-step", {
  skip_if_not_installed("Rphylopars")
  set.seed(42)
  tree <- ape::rtree(10)
  df <- data.frame(
    x1 = rnorm(10),
    y  = factor(c("A","B","A","B","A","B","A","B","A","B"),
                levels = c("A","B")),
    row.names = tree$tip.label
  )
  df$y[c(2, 5)] <- NA  # two missing binary cells
  pd <- preprocess_traits(df, tree)

  out  <- build_liability_matrix(pd, splits = NULL)
  liab <- out$X_liab

  # Shape: n_species x n_liab_cols (1 continuous + 1 binary = 2)
  expect_equal(dim(liab), c(10, 2))
  expect_equal(colnames(liab), c("x1", "y"))

  # Continuous col: observed values pass through (z-scored already)
  expect_equal(unname(liab[, "x1"]), unname(pd$X_scaled[, "x1"]))

  # Binary col: observed cells get posterior mean from estep_liability_binary
  # A = 0 -> liability < 0 -> posterior mean is negative
  # B = 1 -> liability > 0 -> posterior mean is positive
  y_obs_rows <- which(!is.na(df$y))
  for (i in y_obs_rows) {
    if (df$y[i] == "A") expect_lt(liab[i, "y"], 0)
    if (df$y[i] == "B") expect_gt(liab[i, "y"], 0)
  }
  # Missing cells: NA
  expect_true(all(is.na(liab[c(2, 5), "y"])))

  expect_named(out, c("X_liab", "liab_cols", "liab_types"))
  expect_equal(out$liab_types, c("continuous", "binary"))
})

test_that("fit_joint_threshold_baseline returns liability-scale posterior for every tip", {
  skip_if_not_installed("Rphylopars")
  set.seed(7)
  tree <- ape::rtree(30)
  df <- data.frame(
    x1 = rnorm(30),
    x2 = rnorm(30),
    row.names = tree$tip.label
  )
  df$x1[c(3, 17)] <- NA
  pd <- preprocess_traits(df, tree)

  res <- fit_joint_threshold_baseline(pd, tree, splits = NULL)

  # Expected shape: list(mu_liab, se_liab) each n_species x n_liab_cols
  expect_named(res, c("mu_liab", "se_liab", "liab_cols", "liab_types"))
  expect_equal(dim(res$mu_liab), c(30, 2))
  expect_equal(dim(res$se_liab), c(30, 2))
  expect_true(all(is.finite(res$mu_liab)))
  expect_true(all(res$se_liab >= 0))

  # Observed x1 cells should have posterior close to the truth
  obs_idx <- which(!is.na(df$x1))
  expect_lt(max(abs(res$mu_liab[obs_idx, "x1"] -
                      pd$X_scaled[obs_idx, "x1"])), 0.5)
})

test_that("fit_joint_threshold_baseline handles mixed continuous + binary", {
  skip_if_not_installed("Rphylopars")
  set.seed(21)
  tree <- ape::rtree(40)
  df <- data.frame(
    x = rnorm(40),
    y = factor(sample(c("A", "B"), 40, TRUE), levels = c("A", "B")),
    row.names = tree$tip.label
  )
  df$y[c(5, 10, 15)] <- NA
  pd <- preprocess_traits(df, tree)

  res <- fit_joint_threshold_baseline(pd, tree, splits = NULL)

  expect_equal(dim(res$mu_liab), c(40, 2))
  expect_equal(res$liab_types, c("continuous", "binary"))
  # All-finite posterior (no NA columns since observed cells exist in both)
  expect_true(all(is.finite(res$mu_liab)))
  expect_true(all(res$se_liab >= 0))
  # Posterior for missing binary cells is between the A and B class means
  # on the liability scale (cross-trait correlation keeps them in range).
  expect_true(all(abs(res$mu_liab[c(5, 10, 15), "y"]) < 5))
})

test_that("fit_joint_threshold_baseline handles all-NA liability columns gracefully", {
  skip_if_not_installed("Rphylopars")
  set.seed(33)
  tree <- ape::rtree(15)
  df <- data.frame(
    x = rnorm(15),
    y = factor(sample(c("A", "B"), 15, TRUE), levels = c("A", "B")),
    row.names = tree$tip.label
  )
  # Make ALL y cells missing so the binary liability column is all-NA
  df$y <- factor(rep(NA, 15), levels = c("A", "B"))
  pd <- preprocess_traits(df, tree)

  res <- fit_joint_threshold_baseline(pd, tree, splits = NULL)

  # Continuous column fits normally; binary column stays all-NA in output
  expect_equal(dim(res$mu_liab), c(15, 2))
  expect_true(all(is.finite(res$mu_liab[, "x"])))
  expect_true(all(is.na(res$mu_liab[, "y"])))
  expect_true(all(is.na(res$se_liab[, "y"])))
})

test_that("build_liability_matrix masks val/test split cells to NA", {
  skip_if_not_installed("Rphylopars")
  set.seed(11)
  tree <- ape::rtree(12)
  df <- data.frame(
    x = rnorm(12),
    y = factor(sample(c("A", "B"), 12, replace = TRUE), levels = c("A", "B")),
    row.names = tree$tip.label
  )
  pd     <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.25,
                                 seed = 11, trait_map = pd$trait_map)
  out <- build_liability_matrix(pd, splits = splits)

  # Every cell in val_idx or test_idx should be NA in X_liab
  hold_idx <- c(splits$val_idx, splits$test_idx)
  n <- nrow(pd$X_scaled); p <- ncol(pd$X_scaled)
  row_i <- ((hold_idx - 1L) %% n) + 1L
  col_j <- ((hold_idx - 1L) %/% n) + 1L
  # Only cells whose column is in out$liab_cols survive in X_liab
  keep  <- col_j %in% out$liab_cols
  row_i <- row_i[keep]
  col_j <- col_j[keep]
  local_col <- match(col_j, out$liab_cols)
  for (k in seq_along(row_i)) {
    expect_true(is.na(out$X_liab[row_i[k], local_col[k]]))
  }
})

test_that("fit_joint_threshold_baseline treats 1-non-NA columns as unfit", {
  skip_if_not_installed("Rphylopars")
  set.seed(44)
  tree <- ape::rtree(20)
  df <- data.frame(
    x = rnorm(20),
    y = factor(sample(c("A", "B"), 20, TRUE), levels = c("A", "B")),
    row.names = tree$tip.label
  )
  # Leave exactly 1 observed y; NA the rest.
  keep_row <- 7L
  y_keep   <- df$y[keep_row]
  df$y     <- factor(rep(NA, 20), levels = c("A", "B"))
  df$y[keep_row] <- y_keep
  pd <- preprocess_traits(df, tree)

  res <- fit_joint_threshold_baseline(pd, tree, splits = NULL)

  # Continuous column should fit normally
  expect_true(all(is.finite(res$mu_liab[, "x"])))
  # Binary column has only 1 non-NA -> treated as unfit -> all-NA output
  expect_true(all(is.na(res$mu_liab[, "y"])))
  expect_true(all(is.na(res$se_liab[, "y"])))
})

test_that("decode_binary_liability returns logit(pnorm(mu/sqrt(1+se^2)))", {
  # Large positive liability -> P near 1 -> large positive logit
  res <- decode_binary_liability(mu_liab = 2, se_liab = 0.1)
  expect_lt(abs(res$p - pnorm(2 / sqrt(1 + 0.01))), 1e-6)
  expect_lt(abs(res$mu_logit - qlogis(res$p)), 1e-6)
  expect_gt(res$mu_logit, 2)

  # Large negative liability -> P near 0 -> large negative logit
  res <- decode_binary_liability(mu_liab = -2, se_liab = 0.1)
  expect_lt(abs(res$p - pnorm(-2 / sqrt(1 + 0.01))), 1e-6)
  expect_lt(res$mu_logit, -2)

  # Zero liability -> P = 0.5 -> logit = 0
  res <- decode_binary_liability(mu_liab = 0, se_liab = 1)
  expect_equal(res$p, 0.5, tolerance = 1e-9)
  expect_equal(res$mu_logit, 0, tolerance = 1e-9)

  # Clipping guard: huge liability should not produce Inf logit
  res <- decode_binary_liability(mu_liab = 50, se_liab = 0)
  expect_equal(res$mu_logit, qlogis(0.99))

  # Vectorised: accepts numeric vectors
  res <- decode_binary_liability(mu_liab = c(-3, 0, 3),
                                  se_liab = c(0.1, 1, 0.1))
  expect_length(res$mu_logit, 3)
  expect_length(res$p, 3)
})
