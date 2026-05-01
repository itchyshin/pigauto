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

# ---- Phase F: per-trait ordinal baseline path selection ----
# Opus adversarial review #6 (2026-04-30).  The threshold-joint baseline
# under-performs per-column BM-via-MVN at K=3 ordinals (see the AVONET
# Migration regression bisected to commit a541dbd).  fit_baseline()
# now computes both paths and picks the lower-val-MSE one per ordinal
# trait, exposing the choice via $ordinal_path_chosen.

test_that("fit_baseline picks per-trait ordinal path and reports it", {
  skip_if_not_installed("Rphylopars")
  set.seed(2026)
  tree <- ape::rtree(40)
  df <- data.frame(
    x = rnorm(40),
    o = ordered(sample(c("low", "med", "high"), 40, TRUE),
                levels = c("low", "med", "high")),
    row.names = tree$tip.label
  )
  pd <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.30,
                                 seed = 2026, trait_map = pd$trait_map)
  bl <- fit_baseline(pd, tree, splits = splits)

  o_col <- pd$trait_map[[2]]$latent_cols
  # Output finite for all species
  expect_true(all(is.finite(bl$mu[, o_col])))
  # Selection field present and recorded for the ordinal column
  expect_true("ordinal_path_chosen" %in% names(bl))
  expect_true(as.character(o_col) %in% names(bl$ordinal_path_chosen))
  expect_true(bl$ordinal_path_chosen[[as.character(o_col)]] %in%
              c("threshold_joint", "bm_mvn"))
})

test_that("fit_baseline ordinal selection picks BM when threshold-joint loses on val", {
  skip_if_not_installed("Rphylopars")
  set.seed(2027)
  tree <- ape::rtree(40)
  df <- data.frame(
    x = rnorm(40),
    o = ordered(sample(c("low", "med", "high"), 40, TRUE),
                levels = c("low", "med", "high")),
    row.names = tree$tip.label
  )
  pd <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.30,
                                 seed = 2027, trait_map = pd$trait_map)

  # Run fit_baseline as usual (selection enabled)
  bl <- fit_baseline(pd, tree, splits = splits)
  o_col <- pd$trait_map[[2]]$latent_cols
  chosen <- bl$ordinal_path_chosen[[as.character(o_col)]]

  # Compute the BOTH-paths predictions independently and verify the
  # selection picked the lower-val-MSE path (no monkey-patching: the
  # selection logic is a pure function of the predictions and val cells).
  truth_full <- pd$X_scaled
  n_obs <- nrow(truth_full)
  val_idx  <- splits$val_idx
  val_col  <- ((val_idx - 1L) %/% n_obs) + 1L
  val_row  <- ((val_idx - 1L) %% n_obs) + 1L
  val_rows <- val_row[val_col == o_col]
  expect_gt(length(val_rows), 0L)
  # Compute alternative prediction: per-column BM on z-scored col with
  # val cells masked.
  X_masked <- pd$X_scaled
  X_masked[splits$val_idx]  <- NA
  X_masked[splits$test_idx] <- NA
  R_phy <- pigauto:::phylo_cor_matrix(tree)
  R_phy <- R_phy[tree$tip.label, tree$tip.label]
  bm_alt <- pigauto:::bm_impute_col(X_masked[, o_col], R_phy)
  truth_o <- truth_full[val_rows, o_col]
  finite_t <- is.finite(truth_o)
  bm_mse <- mean((bm_alt$mu[val_rows[finite_t]] - truth_o[finite_t])^2)
  tj_mse <- mean((bl$mu[val_rows[finite_t], o_col] - truth_o[finite_t])^2)

  if (chosen == "bm_mvn") {
    # bl$mu = bm_alt; tj_mse here is BM MSE; nothing extra to check.
    expect_equal(unname(bl$mu[, o_col]), unname(bm_alt$mu),
                 tolerance = 1e-9,
                 info = "bm_mvn was chosen so bl$mu must match bm_alt")
  } else {
    # threshold-joint kept; verify it was not strictly worse than BM
    # by recomputing both MSEs using bl$mu (which is the threshold
    # path) and the independent bm_alt.
    expect_lte(tj_mse, bm_mse + 1e-9,
               label = "threshold-joint kept => its val MSE <= BM val MSE")
  }
})

# ---- C1: AVONET-300 Migration K=3 regression bottled as fixture ----------
# This test locks in the Phase F win against future refactors.  Pre-Phase-F
# (commit a541dbd .. f43edb6^), AVONET-300 Migration's threshold-joint
# baseline regressed RMSE from 0.879 to 0.975 vs the LP-equivalent path
# (see useful/MEMO_2026-04-29_phase6_migration_bisect.md).  Phase F's
# per-trait selection now picks the BM-via-MVN path on Migration whenever
# threshold-joint loses on val.
#
# This test asserts the SELECTION OUTCOME on AVONET-300, not the absolute
# accuracy.  If a future refactor reintroduces the regression silently
# (e.g. by changing the val-MSE comparison logic, or swapping the BM path
# implementation), this test will fail at the AVONET regression bench
# script's expense before the user-facing benchmark drift surfaces.

test_that("[C1] AVONET-300 Migration ordinal: BM-via-MVN beats threshold-joint on val (Phase F regression bottle)", {
  skip_if_not_installed("Rphylopars")
  data("avonet300", package = "pigauto")
  data("tree300",   package = "pigauto")
  set.seed(42)
  traits <- avonet300
  rownames(traits) <- traits$Species_Key
  traits$Species_Key <- NULL

  pd     <- preprocess_traits(traits, tree300)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.25, seed = 42,
                                 trait_map = pd$trait_map)
  bl <- fit_baseline(pd, tree300, splits = splits)

  # Locate the Migration latent column (it's an ordered factor)
  mig_col <- NULL
  for (tm in pd$trait_map) {
    if (identical(tm$name, "Migration") && tm$type == "ordinal") {
      mig_col <- tm$latent_cols
      break
    }
  }
  expect_false(is.null(mig_col),
               info = "AVONET-300 must have a Migration ordinal trait")

  expect_true("ordinal_path_chosen" %in% names(bl),
              info = "fit_baseline must expose $ordinal_path_chosen on AVONET-300 (single-obs)")
  chosen <- bl$ordinal_path_chosen[[as.character(mig_col)]]
  expect_identical(chosen, "bm_mvn",
                   info = paste0(
                     "[C1] On AVONET-300 K=3 Migration ordinal, the per-trait ",
                     "selection MUST pick BM-via-MVN over threshold-joint. ",
                     "If this fails, Phase F has regressed -- check ",
                     "useful/MEMO_2026-04-29_phase6_migration_bisect.md and ",
                     "useful/MEMO_2026-04-30_ordinal_path_selection.md."))

  # Sanity: the absolute lift on Migration RMSE vs LP-only should be
  # near-zero (Phase F restores parity with the LP path), not the
  # pre-Phase-F regression of -0.085.  Check on test-set RMSE so this
  # asserts the user-facing outcome.
  n   <- nrow(pd$X_scaled)
  ri  <- ((splits$test_idx - 1L) %% n) + 1L
  cj  <- ((splits$test_idx - 1L) %/% n) + 1L
  test_mig_rows <- ri[cj == mig_col]
  truth <- pd$X_scaled[test_mig_rows, mig_col]
  pred  <- bl$mu[test_mig_rows, mig_col]
  rmse_phaseF <- sqrt(mean((truth - pred)^2))
  # Pre-Phase-F regression value was 0.975; Phase F restores 0.890 (and
  # the LP path gives 0.890 too).  Assert RMSE is comfortably below the
  # regression threshold of 0.93 (halfway between 0.890 and 0.975).
  expect_lt(rmse_phaseF, 0.93,
            label = sprintf("[C1] Migration test RMSE = %.3f", rmse_phaseF))
})

# Phase F (2026-05-01): LP added as a third option in the per-trait
# ordinal path selection.  Tests below verify mechanical correctness:
# the LP path is reachable, produces finite z-scaled predictions, and
# does not cause the existing AVONET-300 regression bottle to fail.

test_that("[Phase F] fit_baseline ordinal selection lists 'lp' as a valid choice", {
  skip_if_not_installed("Rphylopars")
  set.seed(2031L)
  tree <- ape::rtree(40L)
  df <- data.frame(
    x = stats::rnorm(40L),
    o = ordered(sample(c("low", "med", "high"), 40L, TRUE),
                levels = c("low", "med", "high")),
    row.names = tree$tip.label
  )
  pd     <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.30,
                                 seed = 2031L,
                                 trait_map = pd$trait_map)
  bl <- fit_baseline(pd, tree, splits = splits)

  o_col <- pd$trait_map[[2L]]$latent_cols
  expect_true("ordinal_path_chosen" %in% names(bl))
  chosen <- bl$ordinal_path_chosen[[as.character(o_col)]]
  # Phase F adds "lp" to the legal set of selection outcomes.
  expect_true(chosen %in% c("threshold_joint", "bm_mvn", "lp"),
              info = "Phase F: chosen path must be one of the three options")
  # Whichever path was selected, predictions must be finite and on the
  # same z-scale (sanity).
  expect_true(all(is.finite(bl$mu[, o_col])))
  expect_true(all(is.finite(bl$se[, o_col])))
  expect_true(all(abs(bl$mu[, o_col]) < 5),
              info = "Phase F: ordinal latent predictions must stay within ~5 SD")
})

test_that("[Phase F] LP option produces a sensible z-scale prediction when chosen", {
  skip_if_not_installed("Rphylopars")
  # K=3 ordinal where the threshold-joint baseline is misspecified.
  # If LP wins, verify its predictions are bounded by E[class] in [1, K]
  # which on the integer-z scale is bounded approximately in [-z_max, +z_max]
  # with z_max = (K - mean(1:K)) / sd(1:K).  For K=3 uniform: mean=2,
  # sd=1, so z_max = 1.0 (predictions bounded in [-1, +1]).
  set.seed(2032L)
  n_tip <- 60L
  tree  <- ape::rtree(n_tip)
  cls <- sample(c("low", "med", "high"), n_tip, replace = TRUE,
                prob = c(0.4, 0.2, 0.4))
  df <- data.frame(
    x = stats::rnorm(n_tip),
    o = ordered(cls, levels = c("low", "med", "high")),
    row.names = tree$tip.label
  )
  pd     <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.30,
                                 seed = 2032L,
                                 trait_map = pd$trait_map)
  bl <- fit_baseline(pd, tree, splits = splits)

  o_col <- NULL
  for (tm in pd$trait_map) {
    if (tm$type == "ordinal") { o_col <- tm$latent_cols; break }
  }
  expect_false(is.null(o_col),
               info = "K=3 ordered factor must preprocess as ordinal")
  chosen <- bl$ordinal_path_chosen[[as.character(o_col)]]
  expect_true(chosen %in% c("threshold_joint", "bm_mvn", "lp"))
  expect_true(all(is.finite(bl$mu[, o_col])))
  if (chosen == "lp") {
    # LP-via-OVR: E[class] bounded in [1, K]; z-scaled bounded in
    # roughly [(1 - mean) / sd, (K - mean) / sd]; for K=3 uniform
    # ~[-1, +1] but with non-uniform priors it can be slightly outside.
    # Use a generous bound (+/- 1.5) to allow for class-frequency shifts.
    rng <- range(bl$mu[, o_col])
    expect_gte(rng[1], -1.5)
    expect_lte(rng[2],  1.5)
  }
})
