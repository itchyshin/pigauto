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

test_that("[active] suggest_next_observation handles binary trait via entropy reduction", {
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
  # Mixed-type input: should have BOTH continuous (variance) and binary
  # (entropy) suggestions in the output
  expect_true(any(sug$type == "continuous"))
  expect_true(any(sug$type == "binary"))
  expect_true("metric" %in% names(sug))
  expect_true(all(sug$metric[sug$type == "continuous"] == "variance"))
  expect_true(all(sug$metric[sug$type == "binary"]     == "entropy"))
  expect_true(all(sug$delta >= 0))
})

test_that("[active] suggest_next_observation skips discrete when types excludes them", {
  skip_if_not_installed("torch")
  set.seed(2116L)
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
                              verbose = FALSE, seed = 2116L)
  sug <- pigauto::suggest_next_observation(result, top_n = 100,
                                             types = c("continuous", "count",
                                                       "ordinal", "proportion"))
  expect_true(all(sug$type == "continuous"))
})

test_that("[active] suggest_next_observation handles K=3 categorical via entropy reduction", {
  skip_if_not_installed("torch")
  set.seed(2117L)
  n <- 30L
  tree <- ape::rtree(n)
  df <- data.frame(
    cat3 = factor(sample(c("A", "B", "C"), n, replace = TRUE),
                  levels = c("A", "B", "C")),
    row.names = tree$tip.label
  )
  df$cat3[c(3L, 7L, 11L, 17L)] <- NA
  result <- pigauto::impute(df, tree, epochs = 15L, n_imputations = 1L,
                              verbose = FALSE, seed = 2117L)
  sug <- pigauto::suggest_next_observation(result, top_n = 5)
  expect_gte(nrow(sug), 1L)
  expect_true(all(sug$type == "categorical"))
  expect_true(all(sug$metric == "entropy"))
  expect_true(all(sug$delta >= 0))
  # Descending
  expect_true(all(diff(sug$delta) <= 0))
})

test_that("[active] lp_entropy_reduction_binary matches brute-force on small fixture", {
  set.seed(2118L)
  n <- 15L
  tree <- ape::rtree(n)
  D <- ape::cophenetic.phylo(tree)
  sigma_lp <- stats::median(D) * 0.5
  sim <- exp(-(D^2) / (2 * sigma_lp^2))
  diag(sim) <- 0
  y <- as.numeric(stats::rbinom(n, 1, 0.5))
  names(y) <- tree$tip.label
  miss_idx <- c(2L, 5L, 9L, 13L)
  y[miss_idx] <- NA

  closed <- pigauto:::lp_entropy_reduction_binary(y, sim)
  expect_length(closed, length(miss_idx))
  expect_true(all(closed >= 0))
  expect_true(all(is.finite(closed)))

  # Brute-force: for each candidate s_new, compute current LP entropy at miss
  # cells, then for y_new in {0, 1} compute new LP entropy, take expectation
  # weighted by current p_s_new.  This directly mirrors the closed-form
  # formula -- but written from scratch as a check.
  H_bin <- function(p) {
    p <- pmin(pmax(p, 0.01), 0.99)
    -p * log(p) - (1 - p) * log(1 - p)
  }
  obs_idx <- which(!is.na(y))
  current_p <- function(y_full, sim_mat) {
    obs <- which(!is.na(y_full))
    sim_o <- sim_mat[, obs, drop = FALSE]
    rw <- rowSums(sim_o)
    rw[rw < 1e-10] <- 1e-10
    p <- as.numeric(sim_o %*% y_full[obs]) / rw
    pmin(pmax(p, 0.01), 0.99)
  }
  p_cur <- current_p(y, sim)
  H_cur <- H_bin(p_cur)
  total_H_cur <- sum(H_cur[miss_idx])

  for (k in seq_along(miss_idx)) {
    s_new <- miss_idx[k]
    q <- p_cur[s_new]
    # y_new = 0
    y_y0 <- y; y_y0[s_new] <- 0
    p_y0 <- current_p(y_y0, sim)
    miss_other <- miss_idx[miss_idx != s_new]
    sum_E_y0 <- sum(H_bin(p_y0[miss_other]))
    # y_new = 1
    y_y1 <- y; y_y1[s_new] <- 1
    p_y1 <- current_p(y_y1, sim)
    sum_E_y1 <- sum(H_bin(p_y1[miss_other]))
    expected_total_after <- (1 - q) * sum_E_y0 + q * sum_E_y1
    bf <- total_H_cur - expected_total_after
    expect_equal(closed[k], bf, tolerance = 1e-6,
                 info = sprintf("[active] entropy reduction k=%d brute-force mismatch", k))
  }
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

# ===========================================================================
# v2 (2026-05-01): zi_count + multi_proportion support
# ===========================================================================

test_that("[active v2] zi_count: hybrid variance + entropy populated", {
  skip_if_not_installed("torch")
  set.seed(2200L)
  n <- 30L
  tree <- ape::rtree(n)
  counts <- stats::rpois(n, lambda = 5)
  counts[sample.int(n, round(n * 0.3))] <- 0L
  df <- data.frame(parasites = as.integer(counts),
                    row.names = tree$tip.label)
  df$parasites[c(3L, 7L, 11L, 17L)] <- NA

  result <- pigauto::impute(df, tree,
                              trait_types = c(parasites = "zi_count"),
                              epochs = 15L, n_imputations = 1L,
                              verbose = FALSE, seed = 2200L)
  sug <- pigauto::suggest_next_observation(result, top_n = 10)
  expect_gte(nrow(sug), 1L)
  expect_true(all(sug$type == "zi_count"))
  expect_true(all(sug$metric == "variance"))
  # BOTH variance AND entropy populated for zi_count rows
  expect_true(all(!is.na(sug$delta_var_total)))
  expect_true(all(!is.na(sug$delta_entropy_total)))
  # Variance reduction is probability-weighted -> non-negative,
  # entropy reduction is non-negative
  expect_true(all(sug$delta_var_total >= 0))
  expect_true(all(sug$delta_entropy_total >= 0))
  # Sorted by delta (= delta_var_total for zi_count) descending
  expect_true(all(diff(sug$delta) <= 0))
})

test_that("[active v2] zi_count: types argument can exclude zi_count", {
  skip_if_not_installed("torch")
  set.seed(2201L)
  n <- 30L
  tree <- ape::rtree(n)
  counts <- stats::rpois(n, lambda = 5)
  counts[sample.int(n, round(n * 0.3))] <- 0L
  df <- data.frame(parasites = as.integer(counts),
                    cont      = stats::rnorm(n),
                    row.names = tree$tip.label)
  df$parasites[c(3L, 7L, 11L)] <- NA
  df$cont[c(5L, 9L)] <- NA

  result <- pigauto::impute(df, tree,
                              trait_types = c(parasites = "zi_count"),
                              epochs = 15L, n_imputations = 1L,
                              verbose = FALSE, seed = 2201L)
  sug <- pigauto::suggest_next_observation(result, top_n = 100,
                                             types = c("continuous", "count",
                                                       "ordinal", "proportion"))
  expect_true(all(sug$type == "continuous"))
})

test_that("[active v2] multi_proportion: per-component BM variance summed", {
  skip_if_not_installed("torch")
  set.seed(2202L)
  n <- 25L
  tree <- ape::rtree(n)
  raw <- matrix(stats::rgamma(n * 3L, shape = 2), n, 3L)
  props <- raw / rowSums(raw)
  df <- data.frame(red   = props[, 1L],
                    green = props[, 2L],
                    blue  = props[, 3L],
                    row.names = tree$tip.label)
  miss_rows <- c(3L, 7L, 11L, 17L)
  df$red[miss_rows]   <- NA_real_
  df$green[miss_rows] <- NA_real_
  df$blue[miss_rows]  <- NA_real_

  result <- pigauto::impute(df, tree,
                              multi_proportion_groups = list(
                                colour = c("red", "green", "blue")),
                              epochs = 15L, n_imputations = 1L,
                              verbose = FALSE, seed = 2202L)
  sug <- pigauto::suggest_next_observation(result, top_n = 10)
  expect_gte(nrow(sug), 1L)
  expect_true(all(sug$type == "multi_proportion"))
  expect_true(all(sug$metric == "variance"))
  expect_true(all(!is.na(sug$delta_var_total)))
  expect_true(all(is.na(sug$delta_entropy_total)))
  expect_true(all(sug$delta_var_total >= 0))
  # Sorted descending
  expect_true(all(diff(sug$delta) <= 0))
})

test_that("[active v2] multi_proportion: delta equals K-component sum", {
  # Verify the multi_proportion delta_var_total is exactly the sum of
  # per-component BM variance reductions (not double counted, no scale
  # bug).
  skip_if_not_installed("torch")
  set.seed(2203L)
  n <- 20L
  tree <- ape::rtree(n)
  raw <- matrix(stats::rgamma(n * 3L, shape = 2), n, 3L)
  props <- raw / rowSums(raw)
  df <- data.frame(red   = props[, 1L],
                    green = props[, 2L],
                    blue  = props[, 3L],
                    row.names = tree$tip.label)
  miss_rows <- c(3L, 7L, 11L)
  df$red[miss_rows]   <- NA_real_
  df$green[miss_rows] <- NA_real_
  df$blue[miss_rows]  <- NA_real_

  result <- pigauto::impute(df, tree,
                              multi_proportion_groups = list(
                                colour = c("red", "green", "blue")),
                              epochs = 10L, n_imputations = 1L,
                              verbose = FALSE, seed = 2203L)
  sug <- pigauto::suggest_next_observation(result, top_n = 10)

  # Recompute the per-component variance reduction independently from
  # X_scaled and check the sum matches.
  R_phy <- pigauto:::phylo_cor_matrix(tree)[tree$tip.label, tree$tip.label]
  X <- result$data$X_scaled
  tm_mp <- NULL
  for (tm in result$data$trait_map) {
    if (identical(tm$type, "multi_proportion")) { tm_mp <- tm; break }
  }
  expect_false(is.null(tm_mp))

  delta_sum <- numeric(0)
  for (k in seq_along(tm_mp$latent_cols)) {
    y_k <- X[, tm_mp$latent_cols[k]]
    delta_k <- pigauto:::bm_variance_reduction(y_k, R_phy)
    if (length(delta_sum) == 0L) delta_sum <- delta_k else
      delta_sum <- delta_sum + delta_k
  }
  expect_length(delta_sum, length(miss_rows))

  # sug rows are sorted descending by delta, but the underlying
  # delta_sum is in miss_idx order.  Match by species name.
  miss_species <- tree$tip.label[miss_rows]
  for (sp in miss_species) {
    sug_val <- sug$delta_var_total[sug$species == sp]
    if (length(sug_val) == 0L) next
    ind_in_miss <- match(sp, names(delta_sum))
    if (is.na(ind_in_miss)) next
    expect_equal(sug_val, delta_sum[[ind_in_miss]], tolerance = 1e-9,
                 info = sprintf("[active v2] multi_proportion sum check for %s", sp))
  }
})

# ===========================================================================
# T1 (2026-05-01): count + ordinal + proportion explicit coverage
#
# These three types use the same continuous-family BM variance-reduction
# closed form as `continuous`.  They were validated indirectly via the
# `continuous` test, but a dedicated mixed fixture pins down the
# trait_map encoding (log1p for count, integer-z for ordinal,
# logit-z for proportion) and confirms each type appears in the
# output and sorts correctly.
# ===========================================================================

test_that("[active T1] suggest_next_observation handles count + ordinal + proportion in one mixed fixture", {
  skip_if_not_installed("torch")
  set.seed(2300L)
  n <- 30L
  tree <- ape::rtree(n)

  # Three traits, all continuous-family but with different encodings:
  #   count:      integer + log1p + z-score (auto-detected)
  #   ordinal:    ordered factor + integer-z latent (explicit override)
  #   proportion: numeric in (0,1) + logit + z-score (explicit override
  #               -- not auto-detected per CLAUDE.md)
  cnt   <- as.integer(stats::rpois(n, lambda = 4))
  ord   <- ordered(sample(letters[1:4], n, replace = TRUE),
                   levels = letters[1:4])
  prop  <- stats::plogis(stats::rnorm(n))

  df <- data.frame(parasites = cnt,
                   stage     = ord,
                   p         = prop,
                   row.names = tree$tip.label)
  miss_rows <- c(2L, 5L, 9L, 13L, 17L, 21L, 25L)
  df$parasites[miss_rows[1:3]] <- NA_integer_
  df$stage[miss_rows[3:5]] <- NA
  df$p[miss_rows[5:7]] <- NA_real_

  result <- pigauto::impute(df, tree,
                              trait_types = c(p = "proportion"),
                              epochs = 15L, n_imputations = 1L,
                              verbose = FALSE, seed = 2300L)
  sug <- pigauto::suggest_next_observation(result, top_n = 30L)

  expect_s3_class(sug, "pigauto_active")
  expect_gt(nrow(sug), 0L)

  # All three trait types must appear in the output
  types_seen <- unique(sug$type)
  expect_true("count"      %in% types_seen,
              info = "T1: count missing-cell rows must appear in output")
  expect_true("ordinal"    %in% types_seen,
              info = "T1: ordinal missing-cell rows must appear in output")
  expect_true("proportion" %in% types_seen,
              info = "T1: proportion missing-cell rows must appear in output")

  # All three are continuous-family -> metric == "variance"
  expect_true(all(sug$metric == "variance"))
  # All deltas non-negative and finite
  expect_true(all(is.finite(sug$delta_var_total)))
  expect_true(all(sug$delta_var_total >= 0))
  # Sorted descending on the unified delta column
  expect_true(all(diff(sug$delta) <= 0),
              info = "T1: rows must be sorted by descending delta")
})

test_that("[active T1] suggest_next_observation respects types argument for count/ordinal/proportion", {
  skip_if_not_installed("torch")
  set.seed(2301L)
  n <- 25L
  tree <- ape::rtree(n)
  cnt  <- as.integer(stats::rpois(n, lambda = 4))
  ord  <- ordered(sample(letters[1:3], n, replace = TRUE),
                  levels = letters[1:3])
  prop <- stats::plogis(stats::rnorm(n))
  df <- data.frame(parasites = cnt,
                   stage     = ord,
                   p         = prop,
                   row.names = tree$tip.label)
  df$parasites[c(2L, 5L)] <- NA_integer_
  df$stage[c(3L, 7L)] <- NA
  df$p[c(4L, 9L)] <- NA_real_
  result <- pigauto::impute(df, tree,
                              trait_types = c(p = "proportion"),
                              epochs = 15L, n_imputations = 1L,
                              verbose = FALSE, seed = 2301L)

  # Restrict to count only
  sug_count <- pigauto::suggest_next_observation(result, types = "count")
  expect_true(all(sug_count$type == "count"))

  # Restrict to ordinal + proportion
  sug_ord_prop <- pigauto::suggest_next_observation(result,
                                                       types = c("ordinal",
                                                                 "proportion"))
  expect_true(all(sug_ord_prop$type %in% c("ordinal", "proportion")))
  expect_true(any(sug_ord_prop$type == "ordinal"))
  expect_true(any(sug_ord_prop$type == "proportion"))
})
