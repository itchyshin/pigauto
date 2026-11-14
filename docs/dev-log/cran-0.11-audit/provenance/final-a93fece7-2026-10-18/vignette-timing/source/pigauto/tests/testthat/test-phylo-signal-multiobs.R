# P1-12: the phylo-signal gate was a silent no-op in multi-obs mode.
#
# X_scaled has one row per OBSERVATION; species_names has one entry per
# SPECIES. The old code indexed species_names with an observation-length
# logical, produced NA tip names, ape::keep.tip() errored, the caller's
# tryCatch swallowed it to NA_real_, and the gate never fired.

skip_if_not_installed("phytools")

make_multiobs <- function(n_sp = 40L, reps = 3L, seed = 11L) {
  set.seed(seed)
  tr <- ape::rtree(n_sp)
  sp <- tr$tip.label
  # strong-signal continuous trait, BM on the tree
  z <- as.numeric(ape::rTraitCont(tr, model = "BM", sigma = 1))
  names(z) <- sp
  df <- data.frame(
    species = rep(sp, each = reps),
    x = as.numeric(rep(z, each = reps)) + rnorm(n_sp * reps, 0, 0.05),
    stringsAsFactors = FALSE)
  list(tree = tr, df = df, n_sp = n_sp)
}

test_that("[P1-12] per-trait phylo signal is finite in multi-obs mode", {
  fx <- make_multiobs()
  d <- preprocess_traits(fx$df, fx$tree, species_col = "species")
  expect_true(!is.null(d$obs_to_species))
  expect_lt(length(d$species_names), nrow(d$X_scaled))   # genuinely multi-obs

  sig <- pigauto:::compute_phylo_signal_per_trait(d, fx$tree, method = "lambda",
                                                  min_tips = 20L)
  expect_length(sig, length(d$trait_map))
  # The defect: this was NA for every trait. A BM-simulated trait must
  # return a finite lambda.
  expect_false(is.na(sig[["x"]]))
  expect_gte(sig[["x"]], 0)
})

test_that("[P1-12] multi-obs signal recovers strong BM signal (lambda near 1)", {
  fx <- make_multiobs()
  d <- preprocess_traits(fx$df, fx$tree, species_col = "species")
  sig <- pigauto:::compute_phylo_signal_per_trait(d, fx$tree, method = "lambda",
                                                  min_tips = 20L)
  # BM-simulated with small observation noise -> lambda should be high.
  expect_gt(sig[["x"]], 0.5)
})

test_that("[P1-12] multi-obs and single-obs agree when reps = 1", {
  # With one observation per species the aggregation is the identity, so
  # the multi-obs path must reproduce the single-obs answer exactly.
  fx <- make_multiobs(reps = 1L)
  d_multi <- preprocess_traits(fx$df, fx$tree, species_col = "species")

  single <- data.frame(row.names = fx$df$species, x = fx$df$x)
  d_single <- preprocess_traits(single, fx$tree)

  s_multi <- pigauto:::compute_phylo_signal_per_trait(d_multi, fx$tree,
                                                      method = "lambda", min_tips = 20L)
  s_single <- pigauto:::compute_phylo_signal_per_trait(d_single, fx$tree,
                                                        method = "lambda", min_tips = 20L)
  expect_equal(unname(s_multi[["x"]]), unname(s_single[["x"]]), tolerance = 1e-6)
})

test_that("[P1-12] single-obs path is unchanged", {
  set.seed(5L)
  tr <- ape::rtree(40L)
  z <- as.numeric(ape::rTraitCont(tr, model = "BM", sigma = 1))
  df <- data.frame(row.names = tr$tip.label, x = z)
  d <- preprocess_traits(df, tr)
  expect_null(d$obs_to_species)
  sig <- pigauto:::compute_phylo_signal_per_trait(d, tr, method = "lambda",
                                                  min_tips = 20L)
  expect_false(is.na(sig[["x"]]))
})
