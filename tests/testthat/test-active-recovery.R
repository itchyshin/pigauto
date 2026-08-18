# Integrity checks for the Stage-A evidence harness.  The numerical active
# formulas themselves have independent tests in test-active-impute.R; these
# tests exercise the experimental treatment, split protection, and receipt.

ar_env <- new.env(parent = globalenv())
sys.source(testthat::test_path("..", "..", "script", "active_recovery",
                               "00_prepare_active_recovery.R"), envir = ar_env)

test_that("[active recovery] fixed test cells cannot be selected", {
  sp <- paste0("sp", seq_len(100L))
  split <- ar_env$active_recovery_split(sp, seed = 1L)
  expect_length(intersect(split$initial, split$candidate), 0L)
  expect_length(intersect(split$initial, split$test), 0L)
  expect_length(intersect(split$candidate, split$test), 0L)
  expect_equal(length(c(split$initial, split$candidate, split$test)), length(sp))
})

test_that("[active recovery] exactly one candidate value is restored", {
  draw <- ar_env$active_recovery_draw(30L, 1, "continuous", seed = 2L)
  split <- ar_env$active_recovery_split(draw$tree$tip.label, seed = 3L)
  initial <- ar_env$active_recovery_mask(draw$truth, split$initial)
  changed <- ar_env$active_recovery_restore_one(initial, draw$truth, split$candidate[[1L]])
  expect_equal(sum(!is.na(changed$trait)), sum(!is.na(initial$trait)) + 1L)
  expect_false(identical(ar_env$active_recovery_data_hash(initial),
                         ar_env$active_recovery_data_hash(changed)))
  expect_true(all(is.na(changed[split$test, "trait"])))
})

test_that("[active recovery] uncertainty policy uses only decision-time scores", {
  scores <- c(a = 0.2, b = 0.9, c = 0.4)
  fake <- list(data = list(species_names = names(scores)),
               prediction = list(se = matrix(scores, ncol = 1L)))
  class(fake) <- "fake"
  expect_equal(ar_env$active_recovery_select(fake, "continuous", names(scores),
                                              "uncertainty"), "b")
})

test_that("[active recovery] active selection is maximal within candidate set", {
  skip_if_no_libtorch()
  draw <- ar_env$active_recovery_draw(30L, 1, "continuous", seed = 4L)
  split <- ar_env$active_recovery_split(draw$tree$tip.label, seed = 5L)
  initial <- ar_env$active_recovery_mask(draw$truth, split$initial)
  fit <- suppressWarnings(ar_env$active_recovery_fit(initial, draw$tree, seed = 6L, epochs = 5L))
  candidate <- split$candidate
  chosen <- ar_env$active_recovery_select(fit, "continuous", candidate, "active")
  all_scores <- pigauto::suggest_next_observation(
    fit, top_n = length(draw$tree$tip.label), types = "continuous")
  all_scores <- all_scores[all_scores$species %in% candidate, , drop = FALSE]
  expect_equal(chosen, all_scores$species[which.max(all_scores$delta)])
  refit_data <- ar_env$active_recovery_restore_one(initial, draw$truth, chosen)
  refit <- suppressWarnings(ar_env$active_recovery_fit(refit_data, draw$tree, seed = 6L, epochs = 5L))
  pred <- refit$prediction$imputed$trait[match(split$test, refit$data$species_names)]
  expect_true(all(is.finite(pred)))
})

test_that("[active recovery] binary active score is maximal within candidate set", {
  skip_if_no_libtorch()
  draw <- ar_env$active_recovery_draw(30L, 1, "binary", seed = 7L)
  split <- ar_env$active_recovery_split(draw$tree$tip.label, seed = 8L)
  initial <- ar_env$active_recovery_mask(draw$truth, split$initial)
  fit <- suppressWarnings(ar_env$active_recovery_fit(initial, draw$tree, seed = 9L, epochs = 5L))
  chosen <- ar_env$active_recovery_select(fit, "binary", split$candidate, "active")
  ranked <- pigauto::suggest_next_observation(
    fit, top_n = length(draw$tree$tip.label), types = "binary")
  ranked <- ranked[ranked$species %in% split$candidate, , drop = FALSE]
  expect_equal(chosen, ranked$species[which.max(ranked$delta)])
  changed <- ar_env$active_recovery_restore_one(initial, draw$truth, chosen)
  refit <- suppressWarnings(ar_env$active_recovery_fit(changed, draw$tree, seed = 9L, epochs = 5L))
  p <- refit$prediction$probabilities$trait[match(split$test, refit$data$species_names)]
  expect_true(all(is.finite(p) & p >= 0 & p <= 1))
})
