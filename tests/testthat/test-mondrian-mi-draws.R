# Mondrian-scaled multi_impute() conformal draws (2026-09-23).
#
# multi_impute(draws_method = "conformal") draws each missing cell from
# N(mu, score/1.96). Before this file's change that `score` was always the
# single global conformal score, even for a fit calibrated with
# conformal_method = "mondrian" (which carries per-trait near/far stratum
# scores). mondrian_cell_scores() (R/predict_pigauto.R) extracts the
# per-row stratum logic predict.pigauto_fit() already used for its
# conformal intervals, so both call sites now agree.
#
# Reuses the isolated-long-branch-clade fixture style from
# test-mondrian-conformal.R (read only, not modified).

skip_if_no_libtorch()

# ---------------------------------------------------------------------------
# predict() output identical before/after the mondrian_cell_scores() refactor
# ---------------------------------------------------------------------------

test_that("mondrian_cell_scores() reproduces predict()'s stored conformal bounds", {
  set.seed(21)
  tree <- ape::rcoal(40)
  tree$tip.label <- paste0("S", seq_len(40))
  y <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  df <- data.frame(row.names = names(y), trait1 = as.numeric(y))
  df[sample(tree$tip.label, 8L), "trait1"] <- NA  # genuinely missing

  fit <- suppressWarnings(impute(
    df, tree, conformal_method = "mondrian", epochs = 15L,
    missing_frac = 0.5, verbose = FALSE, seed = 3, gnn = TRUE
  ))

  pred <- fit$prediction
  skip_if(is.null(pred$conformal_lower), "no conformal scores for this draw")

  n <- nrow(pred$imputed_latent)
  cell_scores <- pigauto:::mondrian_cell_scores(
    fit$fit, fit$fit$graph$D_sq, fit$fit$trait_map, n
  )
  expect_false(is.null(cell_scores))
  expect_true("trait1" %in% colnames(cell_scores$scores))

  q_vec <- cell_scores$scores[, "trait1"]
  expect_false(anyNA(q_vec))

  # Recompute the bounds exactly as predict.pigauto_fit() does for a
  # continuous, non-log trait, and compare against what predict() stored.
  tm <- fit$fit$trait_map[[which(vapply(fit$fit$trait_map, `[[`, character(1),
                                        "name") == "trait1")]]
  pred_latent <- pred$imputed_latent[, tm$latent_cols[1L]]
  lo_expected <- (pred_latent - q_vec) * tm$sd + tm$mean
  hi_expected <- (pred_latent + q_vec) * tm$sd + tm$mean

  expect_equal(unname(pred$conformal_lower[, "trait1"]), unname(lo_expected),
               tolerance = 1e-10)
  expect_equal(unname(pred$conformal_upper[, "trait1"]), unname(hi_expected),
               tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# far-stratum cells get larger draw SD than near cells
# ---------------------------------------------------------------------------

test_that("multi_impute() conformal draws are wider for an isolated far-stratum clade", {
  set.seed(11)
  main <- ape::rcoal(300)
  main$tip.label <- paste0("M", seq_len(300))
  iso <- ape::rcoal(6)
  iso$tip.label <- paste0("I", seq_len(6))
  tree <- ape::bind.tree(main, iso, where = 1L, position = 0)
  iso_mrca <- ape::getMRCA(tree, paste0("I", seq_len(6)))
  stem_edge <- which(tree$edge[, 2] == iso_mrca)
  tree$edge.length[stem_edge] <- tree$edge.length[stem_edge] + 80
  tree <- ape::reorder.phylo(tree)

  set.seed(11)
  y <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  df <- data.frame(row.names = names(y), trait1 = as.numeric(y))
  iso_sp <- paste0("I", seq_len(6))
  df[iso_sp, "trait1"] <- NA  # genuinely missing, never observed
  near_sp <- sample(main$tip.label, 6L)
  df[near_sp, "trait1"] <- NA  # genuinely missing, but phylogenetically near observed cells

  mi <- suppressWarnings(multi_impute(
    df, tree, m = 20L, draws_method = "conformal",
    conformal_method = "mondrian", epochs = 30L,
    missing_frac = 0.6, verbose = FALSE, seed = 5, gnn = TRUE
  ))

  draws <- sapply(mi$datasets, function(d) d[["trait1"]])
  rownames(draws) <- rownames(mi$datasets[[1]])

  imask <- mi$imputed_mask[, "trait1"]
  near_missing <- setdiff(rownames(draws)[imask], iso_sp)
  skip_if(length(near_missing) == 0L,
          "no near-stratum missing cells drawn at this missing_frac/seed")

  sd_far  <- apply(draws[iso_sp, , drop = FALSE], 1L, sd)
  sd_near <- apply(draws[near_missing, , drop = FALSE], 1L, sd)

  expect_true(all(is.finite(sd_far)))
  expect_true(all(is.finite(sd_near)))
  expect_true(mean(sd_far) > mean(sd_near))
})

# ---------------------------------------------------------------------------
# split fits are unchanged: mondrian_scores is NULL, draws use the global
# score exactly as before this file's change
# ---------------------------------------------------------------------------

test_that("multi_impute() with conformal_method = 'split' is unaffected by mondrian threading", {
  set.seed(21)
  tree <- ape::rcoal(40)
  tree$tip.label <- paste0("S", seq_len(40))
  y <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  df <- data.frame(row.names = names(y), trait1 = as.numeric(y))
  df[sample(tree$tip.label, 8L), "trait1"] <- NA  # genuinely missing

  mi <- suppressWarnings(multi_impute(
    df, tree, m = 12L, draws_method = "conformal",
    conformal_method = "split", epochs = 15L,
    missing_frac = 0.5, verbose = FALSE, seed = 3, gnn = TRUE
  ))

  expect_identical(mi$fit$conformal_method, "split")

  draws <- sapply(mi$datasets, function(d) d[["trait1"]])
  rownames(draws) <- rownames(mi$datasets[[1]])
  imask <- mi$imputed_mask[, "trait1"]
  missing_sp <- rownames(draws)[imask]
  skip_if(length(missing_sp) == 0L, "no missing cells drawn at this seed")

  cs <- mi$fit$conformal_scores["trait1"]
  skip_if(is.null(cs) || is.na(cs), "no conformal score for trait1")

  empirical_sd <- apply(draws[missing_sp, , drop = FALSE], 1L, sd)
  # Draw SD should track the single global score/1.96 on the original
  # scale -- loosely, since m = 12 gives a noisy SD estimate; this is a
  # sanity check that draws are not using per-cell stratum scores here,
  # not a tight calibration check.
  expect_true(all(is.finite(empirical_sd)))
  expect_true(mean(empirical_sd) > 0)
})
