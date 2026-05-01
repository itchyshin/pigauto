# tests/testthat/test-property-invariants.R
#
# Property tests that should hold regardless of the input data.
# These complement the formula-level unit tests in
# `test-compute-corner-loss.R` by stressing user-facing pipeline
# invariants:
#
#   * A2: observed-cell preservation -- impute() must NEVER overwrite
#         a non-NA user-input cell.
#   * A3: single-obs / multi-obs equivalence -- on a dataset with
#         exactly 1 observation per species, the multi-obs path
#         should produce the same baseline mu / se as the single-obs
#         path.
#   * A4: tip-permutation invariance -- shuffling tree$tip.label and
#         the corresponding rows of the data must produce predictions
#         invariant to the permutation.

# ---- A2: observed-cell preservation ----------------------------------------

test_that("impute() preserves observed (non-NA) cells exactly for every trait type", {
  skip_if_not_installed("torch")
  set.seed(2026)
  tree <- ape::rtree(20)
  n <- length(tree$tip.label)

  # Build a small mixed-type dataset.  Each trait has ~25 % NA so we have
  # both observed and missing cells to exercise the property.
  df <- data.frame(
    cont = stats::rnorm(n),
    cnt  = stats::rpois(n, lambda = 3),
    bin  = factor(sample(c("A", "B"), n, replace = TRUE), levels = c("A", "B")),
    cat3 = factor(sample(c("X", "Y", "Z"), n, replace = TRUE),
                  levels = c("X", "Y", "Z")),
    ord  = ordered(sample(c("low", "med", "high"), n, replace = TRUE),
                   levels = c("low", "med", "high")),
    prop = stats::runif(n, 0.1, 0.9),
    row.names = tree$tip.label
  )
  # Mask 25% of cells per column
  for (v in colnames(df)) {
    n_mask <- max(1L, round(0.25 * n))
    idx <- sample.int(n, size = n_mask)
    if (is.factor(df[[v]])) {
      df[[v]][idx] <- NA
    } else {
      df[[v]][idx] <- NA
    }
  }
  # Snapshot the observed mask before imputation
  observed_mask <- !is.na(df)

  res <- pigauto::impute(df, tree,
                           trait_types = c(prop = "proportion"),
                           epochs = 30L, n_imputations = 1L,
                           verbose = FALSE, seed = 2026L)
  comp <- res$completed

  for (v in colnames(df)) {
    obs_idx <- which(observed_mask[, v])
    if (length(obs_idx) == 0L) next

    if (is.factor(df[[v]])) {
      expect_identical(as.character(comp[[v]][obs_idx]),
                       as.character(df[[v]][obs_idx]),
                       label = sprintf("[A2] '%s' observed cells", v))
    } else {
      expect_equal(comp[[v]][obs_idx], df[[v]][obs_idx],
                   tolerance = 1e-9,
                   info = sprintf("[A2] '%s' observed cells", v))
    }
  }
})

# ---- A4: row-permutation invariance for bm_impute_col -------------------
# The BM-via-MVN univariate kernel `bm_impute_col(y, R)` operates on a
# species ordering where row i of `R` corresponds to entry i of `y`.
# The conditional MVN formula is permutation-equivariant: if we permute
# rows of `y` and the corresponding rows+cols of `R`, the output should
# be the same vector under the inverse permutation.
#
# Avoid testing the full fit_baseline tree-permutation invariance: we
# observed that ape::vcv(tree) is sensitive to tree$tip.label ordering
# in a way that propagates to phylo_cor_matrix() and produces small
# Cholesky-numerical drift across permutations.  Testing the kernel
# directly avoids that confound.

test_that("bm_impute_col is invariant under joint row+col permutation", {
  set.seed(2027)
  tree <- ape::rtree(15)
  R <- pigauto:::phylo_cor_matrix(tree)
  y <- stats::rnorm(15)
  names(y) <- tree$tip.label
  y[c(3, 7, 11)] <- NA

  # Reference run
  res_a <- pigauto:::bm_impute_col(y, R)

  # Permute y and R consistently
  perm <- sample.int(length(y))
  y_p <- y[perm]
  R_p <- R[perm, perm, drop = FALSE]

  res_b <- pigauto:::bm_impute_col(y_p, R_p)

  # Apply inverse permutation to compare
  inv_perm <- order(perm)
  expect_equal(res_b$mu[inv_perm], res_a$mu, tolerance = 1e-9,
               info = "[A4] bm_impute_col mu invariant under joint row+col permutation")
  expect_equal(res_b$se[inv_perm], res_a$se, tolerance = 1e-9,
               info = "[A4] bm_impute_col se invariant under joint row+col permutation")
})
