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
