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

  # Expected shape: list with at least (mu_liab, se_liab, liab_cols,
  # liab_types); Phase 6 additively attaches phylopars_fit + fit_cols_idx so
  # downstream EM code can extract Σ without re-running phylopars.
  expect_true(all(c("mu_liab", "se_liab", "liab_cols", "liab_types")
                   %in% names(res)))
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

test_that("fit_baseline uses threshold-joint path for mixed continuous+binary", {
  skip_if_not_installed("Rphylopars")
  set.seed(1)
  tree <- ape::rtree(40)
  df <- data.frame(
    x  = rnorm(40),
    y  = factor(sample(c("A", "B"), 40, replace = TRUE),
                levels = c("A", "B")),
    row.names = tree$tip.label
  )
  df$x[c(5, 20)] <- NA
  df$y[c(10, 30)] <- NA
  pd     <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.1,
                                 seed = 1, trait_map = pd$trait_map)
  bl <- fit_baseline(pd, tree, splits = splits)

  # Continuous col: mu on z-score scale, finite
  expect_true(all(is.finite(bl$mu[, "x"])))
  # Binary col: mu on logit scale, within clipped range
  expect_true(all(is.finite(bl$mu[, "y"])))
  expect_true(all(bl$mu[, "y"] >= qlogis(0.01) - 1e-6))
  expect_true(all(bl$mu[, "y"] <= qlogis(0.99) + 1e-6))
})

test_that("fit_baseline falls back to LP when Rphylopars not installed", {
  skip_if_not(joint_mvn_available(),
              "This test checks the fallback path; run with Rphylopars removed to exercise it.")
  # When Rphylopars IS installed, simulate absence by stubbing joint_mvn_available
  set.seed(2)
  tree <- ape::rtree(30)
  df <- data.frame(
    x = rnorm(30),
    y = factor(sample(c("A","B"), 30, TRUE), levels = c("A","B")),
    row.names = tree$tip.label
  )
  pd <- preprocess_traits(df, tree)

  # Stub joint_mvn_available to FALSE for this test
  orig <- pigauto:::joint_mvn_available
  assignInNamespace("joint_mvn_available", function() FALSE, ns = "pigauto")
  on.exit(assignInNamespace("joint_mvn_available", orig, ns = "pigauto"))

  bl <- fit_baseline(pd, tree, splits = NULL)

  # LP path returns logit(probs) for binary; within clip range
  expect_true(all(is.finite(bl$mu[, "y"])))
  expect_true(all(bl$mu[, "y"] >= qlogis(0.01) - 1e-6))
  expect_true(all(bl$mu[, "y"] <= qlogis(0.99) + 1e-6))
  # Continuous col: still finite (per-column BM fallback)
  expect_true(all(is.finite(bl$mu[, "x"])))
})

test_that("build_liability_matrix excludes categorical cols from joint liability", {
  skip_if_not_installed("Rphylopars")
  set.seed(50)
  tree <- ape::rtree(15)
  df <- data.frame(
    x = rnorm(15),
    z = factor(sample(c("P", "Q", "R"), 15, replace = TRUE),
               levels = c("P", "Q", "R")),
    row.names = tree$tip.label
  )
  pd <- preprocess_traits(df, tree)
  out <- build_liability_matrix(pd, splits = NULL)
  # Only the 1 continuous col should be in the liability matrix; the K=3
  # categorical cols are handled separately by fit_ovr_categorical_fits().
  expect_equal(ncol(out$X_liab), 1L)
  expect_equal(out$liab_types, "continuous")
})

# ---- B3: ordinal threshold baseline tests ----

test_that("build_liability_matrix applies ordinal E-step (not passthrough)", {
  skip_if_not_installed("Rphylopars")
  set.seed(99)
  tree <- ape::rtree(15)
  df <- data.frame(
    x = rnorm(15),
    o = ordered(sample(c("low", "med", "high"), 15, TRUE),
                levels = c("low", "med", "high")),
    row.names = tree$tip.label
  )
  df$o[c(3, 7)] <- NA
  pd <- preprocess_traits(df, tree)
  out <- build_liability_matrix(pd, splits = NULL)

  # Should have 2 cols (x continuous + o ordinal)
  expect_equal(ncol(out$X_liab), 2)
  expect_true("ordinal" %in% out$liab_types)
  # Ordinal col: non-NA observed cells should differ from the z-scored
  # integer passthrough (the E-step produces different values)
  o_col <- which(out$liab_types == "ordinal")
  obs_rows <- which(!is.na(df$o))
  # The E-step posterior should be finite and in a reasonable range
  expect_true(all(is.finite(out$X_liab[obs_rows, o_col])))
  # Values should NOT exactly equal the z-scored integers from X_scaled
  o_latent_col <- pd$trait_map[[2]]$latent_cols
  z_passthrough <- pd$X_scaled[obs_rows, o_latent_col]
  expect_false(all(abs(out$X_liab[obs_rows, o_col] - z_passthrough) < 1e-10))
  # Missing cells should be NA
  expect_true(all(is.na(out$X_liab[c(3, 7), o_col])))
})

test_that("fit_baseline uses threshold-joint for ordinal + continuous", {
  skip_if_not_installed("Rphylopars")
  set.seed(100)
  tree <- ape::rtree(30)
  df <- data.frame(
    x = rnorm(30),
    o = ordered(sample(c("low", "med", "high"), 30, TRUE),
                levels = c("low", "med", "high")),
    row.names = tree$tip.label
  )
  df$o[c(5, 15, 25)] <- NA
  pd <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.2,
                                 seed = 100, trait_map = pd$trait_map)
  bl <- fit_baseline(pd, tree, splits = splits)

  # Output should be finite for all species
  o_col <- pd$trait_map[[2]]$latent_cols
  expect_true(all(is.finite(bl$mu[, o_col])))
  # Ordinal mu should be within plausible z-scored integer range
  expect_true(all(abs(bl$mu[, o_col]) < 5))
  # Continuous col should also be finite
  expect_true(all(is.finite(bl$mu[, "x"])))
})

test_that("decode_ordinal_liability returns z-scored integer class", {
  skip_if_not_installed("Rphylopars")
  set.seed(101)
  tree <- ape::rtree(10)
  df <- data.frame(
    x = rnorm(10),
    o = ordered(rep(c("A", "B", "C", "D", "E"), 2),
                levels = c("A", "B", "C", "D", "E")),
    row.names = tree$tip.label
  )
  pd <- preprocess_traits(df, tree)
  tm_o <- NULL
  for (tm in pd$trait_map) if (tm$type == "ordinal") { tm_o <- tm; break }

  # liability_info for K=5: thresholds = 0.5, 1.5, 2.5, 3.5
  # Liability value 2.0 falls in interval [1.5, 2.5] = class 3 (0-indexed = 2)
  dec <- decode_ordinal_liability(mu_liab = 2.0, se_liab = 0.1, tm = tm_o)
  expected_z <- (2 - tm_o$mean) / tm_o$sd
  expect_equal(dec$mu_z, expected_z, tolerance = 1e-6)
  expect_equal(dec$se_z, 0)

  # Extreme high liability -> highest class (K=5, 0-indexed = 4)
  dec_high <- decode_ordinal_liability(mu_liab = 100, se_liab = 0, tm = tm_o)
  expected_high <- (4 - tm_o$mean) / tm_o$sd
  expect_equal(dec_high$mu_z, expected_high, tolerance = 1e-6)

  # Extreme low liability -> lowest class (0-indexed = 0)
  dec_low <- decode_ordinal_liability(mu_liab = -100, se_liab = 0, tm = tm_o)
  expected_low <- (0 - tm_o$mean) / tm_o$sd
  expect_equal(dec_low$mu_z, expected_low, tolerance = 1e-6)
})

test_that("ordinal falls back to per-column BM when Rphylopars unavailable", {
  skip_if_not(joint_mvn_available(),
              "This test checks the fallback path")
  set.seed(102)
  tree <- ape::rtree(20)
  df <- data.frame(
    x = rnorm(20),
    o = ordered(sample(c("A", "B", "C"), 20, TRUE),
                levels = c("A", "B", "C")),
    row.names = tree$tip.label
  )
  df$o[c(3, 10)] <- NA
  pd <- preprocess_traits(df, tree)

  # Stub joint_mvn_available to FALSE
  orig <- pigauto:::joint_mvn_available
  assignInNamespace("joint_mvn_available", function() FALSE, ns = "pigauto")
  on.exit(assignInNamespace("joint_mvn_available", orig, ns = "pigauto"))

  bl <- fit_baseline(pd, tree, splits = NULL)

  # Should still produce finite output (per-column BM fallback)
  o_col <- pd$trait_map[[2]]$latent_cols
  expect_true(all(is.finite(bl$mu[, o_col])))
  expect_true(all(is.finite(bl$mu[, "x"])))
})
