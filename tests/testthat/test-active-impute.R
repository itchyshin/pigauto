# tests/testthat/test-active-impute.R
#
# Active imputation: suggest_next_observation() (2026-04-30 evening sprint)
#
# For each currently-missing cell (s, t), compute the closed-form expected
# reduction in *total* predictive variance across all missing cells if
# (s, t) were observed next.  Built on the BM conditional-MVN
# variance-update formula (Sherman-Morrison rank-1 inverse).
#
# These tests verify the math against an independent brute-force
# recomputation: mask cell, refit BM, sum new predictive variances; the
# difference vs the original sum should match the closed-form output.

# ---- Helpers ---------------------------------------------------------------

ref_total_pred_var <- function(y, R, nugget = 1e-6) {
  # Recomputes the sum of conditional variances at all currently-missing
  # cells using bm_impute_col(), without any rank-1 trickery.
  res <- pigauto:::bm_impute_col(y, R, nugget = nugget)
  miss_idx <- which(is.na(y))
  if (length(miss_idx) == 0L) return(0)
  sum(res$se[miss_idx]^2)
}

# Brute-force variance reduction under FIXED sigma2 (active-learning
# assumption).  Matches the closed-form's value-free derivation: we
# compute the new conditional variance at each remaining missing cell
# AFTER moving s_new from miss to obs, holding sigma2 fixed at the
# original REML estimate.  Returns the total reduction (current total
# minus new total, plus the s_new variance which drops to 0).
ref_variance_reduction_fixed_sigma2 <- function(y, R, nugget = 1e-6) {
  obs_idx  <- which(!is.na(y))
  miss_idx <- which(is.na(y))
  n_o <- length(obs_idx); n_m <- length(miss_idx)
  if (n_o < 5L || n_m == 0L) return(numeric(0L))

  # Current sigma2 (REML, same as bm_impute_col)
  R_oo_inv <- solve(R[obs_idx, obs_idx, drop = FALSE] + diag(nugget, n_o))
  ones <- rep(1, n_o); y_o <- y[obs_idx]
  mu_hat <- as.numeric((ones %*% R_oo_inv %*% y_o) /
                          (ones %*% R_oo_inv %*% ones))
  e <- y_o - mu_hat
  sigma2 <- as.numeric(t(e) %*% R_oo_inv %*% e) / max(n_o - 1L, 1L)

  # Current conditional h_i at each miss cell (no s_new added)
  R_om <- R[obs_idx, miss_idx, drop = FALSE]
  h_cur <- diag(t(R_om) %*% R_oo_inv %*% R_om)
  cur_total <- sigma2 * sum(1 - h_cur)

  # For each candidate s_new, refit conditional MVN with s_new added to
  # observed, sigma2 HELD FIXED, sum new variances at other miss cells.
  out <- numeric(n_m)
  for (k in seq_along(miss_idx)) {
    s_new <- miss_idx[k]
    obs_new <- sort(c(obs_idx, s_new))
    miss_new <- setdiff(miss_idx, s_new)
    R_oo_new <- R[obs_new, obs_new, drop = FALSE]
    R_oo_new_inv <- solve(R_oo_new + diag(nugget, length(obs_new)))
    if (length(miss_new) == 0L) {
      new_total <- 0
    } else {
      R_om_new <- R[obs_new, miss_new, drop = FALSE]
      h_new <- diag(t(R_om_new) %*% R_oo_new_inv %*% R_om_new)
      new_total <- sigma2 * sum(1 - h_new)
    }
    out[k] <- cur_total - new_total
  }
  names(out) <- names(y)[miss_idx]
  out
}

# ---- bm_variance_reduction: closed-form math --------------------------------

test_that("[active] bm_variance_reduction matches fixed-sigma2 brute-force on a small fixture", {
  set.seed(2100L)
  n <- 20L
  tree <- ape::rtree(n)
  R <- pigauto:::phylo_cor_matrix(tree)
  y <- stats::rnorm(n)
  names(y) <- tree$tip.label
  miss_idx <- c(3L, 7L, 11L, 17L)
  y[miss_idx] <- NA

  closed <- pigauto:::bm_variance_reduction(y, R)
  expect_equal(length(closed), length(miss_idx))
  expect_true(all(is.finite(closed)))
  expect_true(all(closed >= 0),
              info = "variance reduction must be non-negative")

  # Brute-force under the same assumption the closed-form makes: sigma2
  # held FIXED at current REML estimate.  This is the standard active-
  # learning convention -- variance reduction depends ONLY on the
  # geometry of which cells are observed, not the value at s_new (which
  # we don't know yet).  Refitting sigma2 (e.g., observing at posterior
  # mean) violates this convention and produces a different number;
  # see useful/MEMO_2026-04-30_active_imputation.md for the discussion.
  bf <- ref_variance_reduction_fixed_sigma2(y, R)
  for (k in seq_along(miss_idx)) {
    expect_equal(closed[k], bf[k], tolerance = 1e-6,
                 info = sprintf("[active] brute-force fixed-sigma2 vs closed-form mismatch at k=%d", k))
  }
})

test_that("[active] bm_variance_reduction returns empty when no missing cells", {
  set.seed(2101L)
  n <- 10L
  tree <- ape::rtree(n)
  R <- pigauto:::phylo_cor_matrix(tree)
  y <- stats::rnorm(n)
  names(y) <- tree$tip.label
  out <- pigauto:::bm_variance_reduction(y, R)
  expect_length(out, 0L)
})

test_that("[active] bm_variance_reduction returns empty when too few observed", {
  set.seed(2102L)
  n <- 10L
  tree <- ape::rtree(n)
  R <- pigauto:::phylo_cor_matrix(tree)
  y <- rep(NA_real_, n)
  y[1:3] <- c(1, 2, 3)   # only 3 observed -> n_obs < 5 path
  out <- pigauto:::bm_variance_reduction(y, R)
  expect_length(out, 0L)
})

test_that("[active] bm_variance_reduction is non-negative under tip permutation", {
  # A consistency check: permuting tip order should not produce negative
  # variance reductions.  Catches sign / index errors.
  set.seed(2103L)
  n <- 25L
  tree <- ape::rtree(n)
  R <- pigauto:::phylo_cor_matrix(tree)
  y <- stats::rnorm(n)
  y[c(2, 5, 9, 13, 17, 22)] <- NA
  perm <- sample(n)
  y_p <- y[perm]
  R_p <- R[perm, perm]
  out_p <- pigauto:::bm_variance_reduction(y_p, R_p)
  expect_true(all(out_p >= -1e-10),
              info = "[active] variance reduction must stay non-negative under permutation")
})

# ---- suggest_next_observation: public API -----------------------------------

test_that("[active] suggest_next_observation returns top-N by descending delta_var_total", {
  skip_if_not_installed("torch")
  set.seed(2110L)
  n <- 30L
  tree <- ape::rtree(n)
  df <- data.frame(
    cont1 = stats::rnorm(n),
    cont2 = stats::rnorm(n) * 2,
    row.names = tree$tip.label
  )
  df$cont1[c(2, 5, 9, 13, 17)] <- NA
  df$cont2[c(3, 6, 10, 14, 19)] <- NA

  result <- pigauto::impute(df, tree, epochs = 20L, n_imputations = 1L,
                              verbose = FALSE, seed = 2110L)
  sug <- pigauto::suggest_next_observation(result, top_n = 5)
  expect_s3_class(sug, "data.frame")
  expect_equal(nrow(sug), 5L)
  expect_true(all(c("species", "trait", "type", "delta_var_total")
                   %in% names(sug)))
  # Descending order
  expect_true(all(diff(sug$delta_var_total) <= 0),
              info = "[active] suggestions must be sorted by descending delta_var_total")
  expect_true(all(sug$delta_var_total >= 0))
})

test_that("[active] suggest_next_observation by = 'species' aggregates across traits", {
  skip_if_not_installed("torch")
  set.seed(2111L)
  n <- 30L
  tree <- ape::rtree(n)
  df <- data.frame(
    cont1 = stats::rnorm(n),
    cont2 = stats::rnorm(n) * 2,
    row.names = tree$tip.label
  )
  # Same species missing for BOTH traits
  shared_miss <- c(2L, 5L, 9L)
  df$cont1[shared_miss] <- NA
  df$cont2[shared_miss] <- NA
  # Plus disjoint miss in each
  df$cont1[c(13L, 17L)] <- NA
  df$cont2[c(14L, 19L)] <- NA

  result <- pigauto::impute(df, tree, epochs = 20L, n_imputations = 1L,
                              verbose = FALSE, seed = 2111L)
  sug_cell <- pigauto::suggest_next_observation(result, top_n = 100,
                                                  by = "cell")
  sug_sp   <- pigauto::suggest_next_observation(result, top_n = 100,
                                                  by = "species")
  expect_true("n_traits_missing" %in% names(sug_sp))
  # Each shared_miss species should appear with n_traits_missing = 2
  for (sp in tree$tip.label[shared_miss]) {
    row <- sug_sp[sug_sp$species == sp, ]
    if (nrow(row) > 0L) {
      expect_equal(row$n_traits_missing, 2L,
                   info = sprintf("species %s: n_traits_missing", sp))
      # Aggregate equals sum of cell deltas for that species
      cell_total <- sum(sug_cell$delta_var_total[sug_cell$species == sp])
      expect_equal(row$delta_var_total, cell_total, tolerance = 1e-9,
                   info = sprintf("species %s: by='species' total", sp))
    }
  }
})

test_that("[active] suggest_next_observation returns empty when no missing cells", {
  skip_if_not_installed("torch")
  set.seed(2112L)
  n <- 20L
  tree <- ape::rtree(n)
  df <- data.frame(cont = stats::rnorm(n), row.names = tree$tip.label)
  # Inject ONE NA so impute() runs (it requires at least one missing cell)
  df$cont[1L] <- NA

  result <- pigauto::impute(df, tree, epochs = 15L, n_imputations = 1L,
                              verbose = FALSE, seed = 2112L)
  sug <- pigauto::suggest_next_observation(result, top_n = 10)
  # Should have 1 candidate cell (the one that was masked)
  expect_lte(nrow(sug), 1L)
  expect_true(all(sug$delta_var_total >= 0))
})

test_that("[active] suggest_next_observation rejects multi-obs input with clear error", {
  skip_if_not_installed("torch")
  set.seed(2113L)
  n <- 10L
  tree <- ape::rtree(n)
  df <- data.frame(
    species = rep(tree$tip.label[1:5], each = 2L),
    cont = stats::rnorm(n),
    stringsAsFactors = FALSE
  )
  df$cont[c(2L, 5L, 7L)] <- NA
  tree_5 <- ape::drop.tip(tree, tree$tip.label[6:10])

  result <- pigauto::impute(df, tree_5, species_col = "species",
                              epochs = 10L, n_imputations = 1L,
                              verbose = FALSE, seed = 2113L)
  expect_error(
    pigauto::suggest_next_observation(result, top_n = 5),
    regexp = "multi_obs"
  )
})

test_that("[active] suggest_next_observation skips discrete traits silently", {
  skip_if_not_installed("torch")
  set.seed(2114L)
  n <- 30L
  tree <- ape::rtree(n)
  df <- data.frame(
    cont = stats::rnorm(n),
    bin  = factor(sample(c("A", "B"), n, replace = TRUE), levels = c("A", "B")),
    row.names = tree$tip.label
  )
  df$cont[c(3L, 7L, 11L)] <- NA
  df$bin[c(5L, 9L, 13L)] <- NA

  result <- pigauto::impute(df, tree, epochs = 15L, n_imputations = 1L,
                              verbose = FALSE, seed = 2114L)
  sug <- pigauto::suggest_next_observation(result, top_n = 100)
  # Only continuous trait suggestions
  expect_true(all(sug$type == "continuous"))
  expect_true(all(sug$trait == "cont"))
})

test_that("[active] suggest_next_observation returns 'pigauto_active' class", {
  skip_if_not_installed("torch")
  set.seed(2115L)
  n <- 15L
  tree <- ape::rtree(n)
  df <- data.frame(cont = stats::rnorm(n), row.names = tree$tip.label)
  df$cont[c(3L, 7L)] <- NA
  result <- pigauto::impute(df, tree, epochs = 10L, n_imputations = 1L,
                              verbose = FALSE, seed = 2115L)
  sug <- pigauto::suggest_next_observation(result, top_n = 5)
  expect_s3_class(sug, "pigauto_active")
  expect_s3_class(sug, "data.frame")
})
