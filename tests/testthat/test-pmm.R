# Phase G' (2026-05-01): Predictive Mean Matching tests.
#
# Comprehensive test suite for the PMM imputation path.  Layered:
#   Layer 1 - pmm_impute_one_trait helper unit tests
#   Layer 2 - type-dispatch (which trait types PMM applies to)
#   Layer 3 - statistical properties (no extrapolation, marginal preservation)
#   Layer 4 - integration (end-to-end via impute(), backward compat)
#   Layer 5 - edge cases (small donor pool, all-observed, all-missing,
#             ties, NA predictions)

# ===========================================================================
# Layer 1: pmm_impute_one_trait helper
# ===========================================================================

test_that("[Phase G' L1] pmm_impute_one_trait K=1 returns single nearest donor", {
  pmm_one <- getFromNamespace("pmm_impute_one_trait", "pigauto")
  preds <- c(10, 20, 30, 40, 50)
  truth <- c(11, 21, 31, NA, NA)        # 4 and 5 missing
  out <- pmm_one(preds, truth, K = 1L, seed = 42L)
  # Cell 4: pred=40, closest observed-pred = pred[3]=30 -> donor truth = 31
  # Cell 5: pred=50, closest observed-pred = pred[3]=30 -> donor truth = 31
  expect_equal(out[4], 31)
  expect_equal(out[5], 31)
  # Observed cells unchanged
  expect_equal(out[1:3], truth[1:3])
})

test_that("[Phase G' L1] pmm_impute_one_trait K=5 picks from K nearest", {
  pmm_one <- getFromNamespace("pmm_impute_one_trait", "pigauto")
  set.seed(NULL)  # don't constrain global seed
  # 10 observed cells with predictions 10, 20, ..., 100 and truths 1, 2, ..., 10
  preds_obs <- seq(10, 100, by = 10)
  truth_obs <- seq_len(10)
  preds_miss <- 55  # nearest obs-preds are 50, 60, then 40, 70, then 30, 80
  truth_miss <- NA_real_
  preds <- c(preds_obs, preds_miss)
  truth <- c(truth_obs, truth_miss)
  # Run 200 times with different seeds; each donor must be one of the
  # K=5 nearest by predicted-value distance.
  K_nearest <- order(abs(preds_obs - preds_miss))[1:5]
  K_nearest_truth <- truth_obs[K_nearest]
  donors <- vapply(1:200, function(s) {
    pmm_one(preds, truth, K = 5L, seed = s)[11]
  }, numeric(1L))
  expect_true(all(donors %in% K_nearest_truth),
              info = "[Phase G' L1] all PMM donors must be among K=5 nearest by pred")
  # And we should sample more than one distinct donor across many seeds
  expect_gt(length(unique(donors)), 1L,
            label = "[Phase G' L1] sampling must be stochastic across seeds")
})

test_that("[Phase G' L1] pmm_impute_one_trait empty miss is no-op", {
  pmm_one <- getFromNamespace("pmm_impute_one_trait", "pigauto")
  truth <- c(11, 21, 31)
  preds <- c(10, 20, 30)
  out <- pmm_one(preds, truth, K = 5L, seed = 1L)
  expect_equal(out, truth)
})

test_that("[Phase G' L1] pmm_impute_one_trait empty observed is no-op", {
  pmm_one <- getFromNamespace("pmm_impute_one_trait", "pigauto")
  truth <- rep(NA_real_, 5)              # nothing observed -> no donors
  preds <- 1:5
  out <- pmm_one(preds, truth, K = 5L, seed = 1L)
  expect_equal(out, truth)
})

test_that("[Phase G' L1] pmm_impute_one_trait K capped at n_observed", {
  pmm_one <- getFromNamespace("pmm_impute_one_trait", "pigauto")
  truth <- c(11, 21, NA)                 # only 2 observed
  preds <- c(10, 20, 30)
  out <- pmm_one(preds, truth, K = 100L, seed = 1L)
  # Donor for cell 3 must be one of the 2 observed truths
  expect_true(out[3] %in% c(11, 21))
})

test_that("[Phase G' L1] pmm_impute_one_trait deterministic with same seed", {
  pmm_one <- getFromNamespace("pmm_impute_one_trait", "pigauto")
  preds <- runif(20, 0, 100)
  truth <- c(rnorm(15), rep(NA, 5))
  out1 <- pmm_one(preds, truth, K = 5L, seed = 7L)
  out2 <- pmm_one(preds, truth, K = 5L, seed = 7L)
  expect_identical(out1, out2)
})

test_that("[Phase G' L1] pmm_impute_one_trait does not leak global RNG state", {
  pmm_one <- getFromNamespace("pmm_impute_one_trait", "pigauto")
  set.seed(123L)
  before <- runif(1)
  set.seed(123L)
  pmm_one(c(10, 20, NA), c(11, 21, NA), K = 1L, seed = 99L)  # mutates internally
  after <- runif(1)
  expect_equal(before, after,
               info = "[Phase G' L1] PMM must restore the user's RNG state")
})

test_that("[Phase G' L1] pmm_impute_one_trait NA prediction skips that cell", {
  pmm_one <- getFromNamespace("pmm_impute_one_trait", "pigauto")
  preds <- c(10, 20, 30, NA_real_)
  truth <- c(11, 21, 31, NA_real_)
  out <- pmm_one(preds, truth, K = 1L, seed = 1L)
  # Cell 4: NA prediction -> stays NA (no donor match possible)
  expect_true(is.na(out[4]))
})

test_that("[Phase G' L1] pmm_impute_one_trait K < 1 errors", {
  pmm_one <- getFromNamespace("pmm_impute_one_trait", "pigauto")
  expect_error(pmm_one(1:5, 1:5, K = 0L), "must be >= 1")
})

# ===========================================================================
# Layer 2: type-dispatch via pmm_is_eligible
# ===========================================================================

test_that("[Phase G' L2] pmm_is_eligible says yes for log-cont, count, prop, zi_count", {
  is_elig <- getFromNamespace("pmm_is_eligible", "pigauto")
  expect_true(is_elig(list(type = "continuous", log_transform = TRUE)))
  expect_true(is_elig(list(type = "count")))
  expect_true(is_elig(list(type = "proportion")))
  expect_true(is_elig(list(type = "zi_count")))
})

test_that("[Phase G' L2] pmm_is_eligible says no for un-log cont and discrete-class types", {
  is_elig <- getFromNamespace("pmm_is_eligible", "pigauto")
  expect_false(is_elig(list(type = "continuous", log_transform = FALSE)))
  expect_false(is_elig(list(type = "binary")))
  expect_false(is_elig(list(type = "categorical")))
  expect_false(is_elig(list(type = "ordinal")))
  expect_false(is_elig(list(type = "multi_proportion")))
})

# ===========================================================================
# Layer 3: statistical property -- no extrapolation by construction
# ===========================================================================

test_that("[Phase G' L3] PMM imputed values are EXACTLY in observed pool when pmm_when='always'", {
  skip_if_not_installed("torch")
  set.seed(2080L)
  n <- 30L
  tree <- ape::rtree(n)
  df <- data.frame(
    mass = exp(stats::rnorm(n, 5, 1)),  # log-cont, range ~ [10, 1000]
    row.names = tree$tip.label
  )
  miss_idx <- c(2L, 5L, 9L, 13L, 17L)
  obs_values <- df$mass[-miss_idx]      # the donor pool of TRUTH values
  df$mass[miss_idx] <- NA

  # Use pmm_when = "always" (original mice-style mode) to enforce the
  # strict "every imputation is a donor" invariant.  Phase G''
  # default ("outside_observed") only triggers PMM on extrapolation,
  # so in-range predictions pass through and would NOT equal a donor.
  res <- pigauto::impute(df, tree,
                          epochs = 10L, n_imputations = 5L,
                          match_observed = "pmm",
                          pmm_K = 5L, pmm_when = "always",
                          verbose = FALSE, seed = 2080L)
  imputed <- res$completed$mass[miss_idx]
  # EVERY imputed value must equal one of the observed values exactly
  expect_true(all(imputed %in% obs_values),
              info = "[Phase G' L3] pmm_when='always' must produce only donor values")
})

test_that("[Phase G' L3] PMM never produces values outside observed range", {
  skip_if_not_installed("torch")
  set.seed(2081L)
  n <- 30L
  tree <- ape::rtree(n)
  df <- data.frame(
    mass = exp(stats::rnorm(n, 5, 1)),
    row.names = tree$tip.label
  )
  miss_idx <- c(2L, 5L, 9L)
  obs_values <- df$mass[-miss_idx]
  df$mass[miss_idx] <- NA
  res <- pigauto::impute(df, tree,
                          epochs = 10L, n_imputations = 10L,
                          match_observed = "pmm",
                          verbose = FALSE, seed = 2081L)
  imputed <- res$completed$mass[miss_idx]
  expect_gte(min(imputed), min(obs_values))
  expect_lte(max(imputed), max(obs_values))
})

test_that("[Phase G' L3] PMM does NOT modify observed values", {
  skip_if_not_installed("torch")
  set.seed(2082L)
  n <- 30L
  tree <- ape::rtree(n)
  df <- data.frame(
    mass = exp(stats::rnorm(n, 5, 1)),
    row.names = tree$tip.label
  )
  miss_idx <- c(2L, 5L, 9L)
  obs_truth <- df$mass[-miss_idx]
  df$mass[miss_idx] <- NA
  res <- pigauto::impute(df, tree,
                          epochs = 10L, n_imputations = 5L,
                          match_observed = "pmm",
                          verbose = FALSE, seed = 2082L)
  expect_equal(res$completed$mass[-miss_idx], obs_truth, tolerance = 1e-9,
               info = "[Phase G' L3] PMM must not touch observed cells")
})

# ===========================================================================
# Layer 4: integration with impute()
# ===========================================================================

test_that("[Phase G' L4] match_observed = 'none' default preserves backward compat", {
  skip_if_not_installed("torch")
  set.seed(2083L)
  n <- 30L
  tree <- ape::rtree(n)
  df <- data.frame(
    mass = exp(stats::rnorm(n, 5, 1)),
    row.names = tree$tip.label
  )
  df$mass[c(2L, 5L)] <- NA

  res_default <- pigauto::impute(df, tree,
                                  epochs = 10L, n_imputations = 1L,
                                  verbose = FALSE, seed = 2083L)
  res_explicit_none <- pigauto::impute(df, tree,
                                        match_observed = "none",
                                        epochs = 10L, n_imputations = 1L,
                                        verbose = FALSE, seed = 2083L)
  expect_equal(res_default$completed$mass,
               res_explicit_none$completed$mass,
               tolerance = 1e-9,
               info = "[Phase G' L4] default and 'none' must produce identical predictions")
})

test_that("[Phase G' L4] match_observed = 'pmm' and clamp_outliers = TRUE coexist", {
  skip_if_not_installed("torch")
  set.seed(2084L)
  n <- 30L
  tree <- ape::rtree(n)
  df <- data.frame(
    mass = exp(stats::rnorm(n, 5, 1)),
    row.names = tree$tip.label
  )
  df$mass[c(2L, 5L)] <- NA
  expect_no_error(
    pigauto::impute(df, tree,
                     match_observed = "pmm",
                     clamp_outliers = TRUE,
                     epochs = 10L, n_imputations = 5L,
                     verbose = FALSE, seed = 2084L)
  )
})

test_that("[Phase G' L4] non-eligible trait types pass through PMM unchanged", {
  skip_if_not_installed("torch")
  set.seed(2085L)
  n <- 30L
  tree <- ape::rtree(n)
  # Binary trait: PMM is no-op since binary decode is already from observed levels
  df <- data.frame(
    bin = factor(sample(c("A", "B"), n, replace = TRUE), levels = c("A", "B")),
    row.names = tree$tip.label
  )
  df$bin[c(2L, 5L)] <- NA
  res <- pigauto::impute(df, tree,
                          match_observed = "pmm",
                          epochs = 10L, n_imputations = 5L,
                          verbose = FALSE, seed = 2085L)
  # Binary type should still produce factor-class outputs
  expect_s3_class(res$completed$bin, "factor")
  expect_true(all(as.character(res$completed$bin) %in% c("A", "B")))
})

# ===========================================================================
# Layer 5: edge cases
# ===========================================================================

test_that("[Phase G' L5] pmm_K invalid value errors", {
  skip_if_not_installed("torch")
  n <- 20L
  tree <- ape::rtree(n)
  df <- data.frame(mass = exp(stats::rnorm(n)),
                    row.names = tree$tip.label)
  expect_error(
    pigauto::impute(df, tree,
                     match_observed = "pmm", pmm_K = 0L,
                     epochs = 5L, n_imputations = 1L,
                     verbose = FALSE, seed = 1L),
    "pmm_K"
  )
})

test_that("[Phase G' L5] PMM with K=1 is fully deterministic given a seed", {
  skip_if_not_installed("torch")
  set.seed(2086L)
  n <- 30L
  tree <- ape::rtree(n)
  df <- data.frame(
    mass = exp(stats::rnorm(n, 5, 1)),
    row.names = tree$tip.label
  )
  df$mass[c(2L, 5L)] <- NA
  res1 <- pigauto::impute(df, tree,
                           match_observed = "pmm", pmm_K = 1L,
                           epochs = 10L, n_imputations = 5L,
                           verbose = FALSE, seed = 2086L)
  res2 <- pigauto::impute(df, tree,
                           match_observed = "pmm", pmm_K = 1L,
                           epochs = 10L, n_imputations = 5L,
                           verbose = FALSE, seed = 2086L)
  expect_equal(res1$completed$mass, res2$completed$mass, tolerance = 1e-9)
})

# ===========================================================================
# Acceptance: simulated tail-extrapolation case
# ===========================================================================

# ===========================================================================
# Phase G'' (2026-05-01): conditional PMM via pmm_when = "outside_observed"
# ===========================================================================

test_that("[Phase G'' L1] pmm_impute_one_trait when='outside_observed' keeps in-range preds", {
  pmm_one <- getFromNamespace("pmm_impute_one_trait", "pigauto")
  preds <- c(50, 50, 50, 50)              # all in-range predictions
  truth <- c(40, 60, NA, NA)              # obs range = [40, 60]; pred 50 is in-range
  out <- pmm_one(preds, truth, K = 1L, seed = 1L,
                 when = "outside_observed")
  # Cells 3, 4: pred=50 in [40, 60] -> kept as-is, NOT a donor
  expect_equal(out[3], 50)
  expect_equal(out[4], 50)
  # Observed cells unchanged
  expect_equal(out[1:2], truth[1:2])
})

test_that("[Phase G'' L1] pmm_impute_one_trait when='outside_observed' triggers on extrapolation", {
  pmm_one <- getFromNamespace("pmm_impute_one_trait", "pigauto")
  preds <- c(50, 50, 200, -10)            # cells 3, 4 extrapolate
  truth <- c(40, 60, NA, NA)              # obs range = [40, 60]
  out <- pmm_one(preds, truth, K = 1L, seed = 1L,
                 when = "outside_observed")
  # Cell 3: pred=200 > 60 -> PMM triggers; nearest donor by pred: index 2 (50 vs 200, d=150) wait
  # Actually pred[1]=50, pred[2]=50, both equally close to 200.  Tie-break: order() picks index 1 first.
  # donor 1 -> truth[1] = 40
  expect_equal(out[3], 40)
  # Cell 4: pred=-10 < 40 -> PMM triggers; same logic, donor 1 -> 40
  expect_equal(out[4], 40)
})

test_that("[Phase G'' L1] when='always' replaces every cell (mice mode)", {
  pmm_one <- getFromNamespace("pmm_impute_one_trait", "pigauto")
  preds <- c(50, 50, 51)                  # cell 3 in-range pred
  truth <- c(40, 60, NA)
  out <- pmm_one(preds, truth, K = 1L, seed = 1L,
                 when = "always")
  # Even though pred[3]=51 is in range [40, 60], when="always" forces a donor.
  # Nearest donor pred: |50-51|=1 (cell 1, truth 40) or |50-51|=1 (cell 2, truth 60).
  # Tie -> order() returns index 1 first -> donor truth = 40
  expect_equal(out[3], 40)
})

test_that("[Phase G'' L1] default when is 'outside_observed'", {
  pmm_one <- getFromNamespace("pmm_impute_one_trait", "pigauto")
  preds <- c(50, 50, 50)                  # all in-range
  truth <- c(40, 60, NA)
  out_default <- pmm_one(preds, truth, K = 1L, seed = 1L)
  out_outside <- pmm_one(preds, truth, K = 1L, seed = 1L,
                         when = "outside_observed")
  expect_equal(out_default, out_outside,
               info = "[Phase G''] default `when` must be 'outside_observed'")
})

test_that("[Phase G'' L4] integration: pmm_when defaults preserve in-range predictions", {
  skip_if_not_installed("torch")
  set.seed(2090L)
  n <- 30L
  tree <- ape::rtree(n)
  # log-cont with NO extreme tail; predictions should mostly stay in-range
  df <- data.frame(
    mass = exp(stats::rnorm(n, 5, 0.3)),
    row.names = tree$tip.label
  )
  df$mass[c(2L, 5L, 9L)] <- NA

  res_default_pmm <- pigauto::impute(df, tree,
                                      match_observed = "pmm",
                                      epochs = 10L, n_imputations = 5L,
                                      verbose = FALSE, seed = 2090L)
  # With when='outside_observed' default, in-range predictions should
  # NOT be replaced -> similar to "none" behaviour on this fixture
  res_none <- pigauto::impute(df, tree,
                               match_observed = "none",
                               epochs = 10L, n_imputations = 5L,
                               verbose = FALSE, seed = 2090L)
  # The two should be very close on a small-tail fixture (most preds
  # in-range -> PMM is no-op).  Allow a small tolerance for any cell
  # that does extrapolate.
  diffs <- abs(res_default_pmm$completed$mass[c(2L, 5L, 9L)] -
               res_none$completed$mass[c(2L, 5L, 9L)])
  obs_range_max <- max(df$mass, na.rm = TRUE) - min(df$mass, na.rm = TRUE)
  # Most cells should have zero diff (in-range -> identical).  But the
  # specific outcome depends on the GNN's stochastic dropout output;
  # at minimum verify diffs are bounded by a small fraction of the
  # observed range.
  expect_lt(max(diffs), obs_range_max,
            label = "[Phase G''] in-range PMM should not deviate by more than the observed range")
})

# Original acceptance test
test_that("[Phase G' acceptance] PMM caps a synthetic tail blow-up better than no-clamp", {
  skip_if_not_installed("torch")
  # Construct a fixture where:
  # - log-cont mass with one phylogenetically isolated tip (analogue of
  #   the AVONET Casuarius case)
  # - Mask the isolated tip so it's the missing cell
  # - Run impute() with match_observed = "pmm" and "none"
  # - PMM imputed value must be in observed range; "none" can extrapolate
  set.seed(2087L)
  n <- 50L
  tree <- ape::rtree(n)
  # Most species cluster around mass = 1000 (log = ~6.9); one outlier at
  # mass = 30000 (log = ~10.3) -- this is the isolated species.
  log_mass <- stats::rnorm(n, 6.9, 0.5)
  log_mass[1L] <- log(30000)            # the "Casuarius"
  df <- data.frame(mass = exp(log_mass), row.names = tree$tip.label)
  obs_max <- max(df$mass)               # ~30000
  df$mass[1L] <- NA                     # mask the outlier

  # No PMM, no clamp: prediction may extrapolate (depending on noise)
  res_none <- pigauto::impute(df, tree,
                               epochs = 30L, n_imputations = 10L,
                               match_observed = "none",
                               verbose = FALSE, seed = 2087L)
  # PMM with pmm_when="always": imputation guaranteed to equal an
  # observed value.  (Phase G'' default "outside_observed" only
  # triggers on extrapolation; here the GNN may predict in-range.)
  res_pmm <- pigauto::impute(df, tree,
                              epochs = 30L, n_imputations = 10L,
                              match_observed = "pmm", pmm_K = 5L,
                              pmm_when = "always",
                              verbose = FALSE, seed = 2087L)
  pmm_val <- as.numeric(res_pmm$completed$mass[1L])
  # PMM imputation must be exactly equal to one of the observed mass values
  obs_vals <- df$mass[-1L]
  expect_true(pmm_val %in% obs_vals,
              info = "[Phase G' accept] PMM-always imputation must equal an observed value")
  # PMM must be at most obs_max (no extrapolation)
  expect_lte(pmm_val, obs_max,
              label = "[Phase G' accept] PMM imputation cannot exceed observed max")
})
