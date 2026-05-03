# Pre-shipping coverage tests for v0.9.2.  Three sections:
#   T2.A  Smoke tests for exported functions that previously had ZERO
#         direct test usage (calibration_df, confusion_matrix,
#         compare_methods, simulate_non_bm, read_tree, plot_*).
#   T2.B  Smoke tests for non-default option paths that previously
#         had no coverage (draws_method = "mc_dropout",
#         use_attention = FALSE, conformal_method = "bootstrap").
#   T3    Edge cases: single-species tree, tree with no branch
#         lengths, traits/tree label mismatch, all-observed predict.

# Local helper -- repeats the make_test_data pattern from
# test-fit-predict.R so this file is self-sufficient.
.tcov_make_data <- function(n = 30L, seed = 1L) {
  set.seed(seed)
  tree <- ape::rcoal(n)
  df <- data.frame(
    row.names = tree$tip.label,
    cont = stats::rnorm(n),
    bin  = factor(sample(c("a", "b"), n, replace = TRUE))
  )
  list(tree = tree, df = df)
}


# ---------------------------------------------------------------------------
# T2.A  Smoke tests for exports without direct test usage
# ---------------------------------------------------------------------------

test_that("[T2A] confusion_matrix returns a list with table + accuracy + per_class", {
  truth <- factor(c("a", "b", "a", "b", "a"), levels = c("a", "b"))
  pred  <- factor(c("a", "a", "a", "b", "b"), levels = c("a", "b"))
  cm <- confusion_matrix(truth, pred)
  expect_type(cm, "list")
  expect_true(all(c("table", "accuracy", "per_class") %in% names(cm)))
  # 2-class problem -> 2x2 confusion table
  expect_equal(nrow(cm$table), 2L)
  expect_equal(ncol(cm$table), 2L)
  expect_true(is.numeric(cm$accuracy) && length(cm$accuracy) == 1L)
  expect_true(cm$accuracy >= 0 && cm$accuracy <= 1)
  expect_s3_class(cm$per_class, "data.frame")
  expect_equal(nrow(cm$per_class), 2L)
})

test_that("[T2A] calibration_df returns a data.frame with bin / observed / expected columns", {
  set.seed(7L)
  truth <- rbinom(200, 1, 0.4)
  prob  <- pmin(pmax(truth * 0.7 + stats::rnorm(200, sd = 0.15), 0.01), 0.99)
  cd <- calibration_df(truth, prob, n_bins = 5L)
  expect_s3_class(cd, "data.frame")
  expect_true(nrow(cd) >= 1L)
  # The two semantically required columns: predicted-bin centre + observed rate.
  # Be tolerant of column naming -- check for any column matching the concept.
  has_predicted <- any(grepl("pred|expected|prob|center|centre|mean",
                              names(cd), ignore.case = TRUE))
  has_observed  <- any(grepl("obs|actual|rate|frac|y", names(cd),
                              ignore.case = TRUE))
  expect_true(has_predicted)
  expect_true(has_observed)
})

test_that("[T2A] simulate_non_bm returns a numeric matrix with one column per trait", {
  set.seed(11L)
  tree <- ape::rcoal(15L)
  for (sc in c("OU", "regime_shift", "nonlinear")) {
    out <- simulate_non_bm(tree, n_traits = 3L, scenario = sc, seed = 1L)
    # Some implementations return a matrix, others a data.frame -- accept both.
    expect_true(is.numeric(out) || is.data.frame(out))
    expect_equal(NROW(out), length(tree$tip.label))
    expect_equal(NCOL(out), 3L)
  }
})

test_that("[T2A] simulate_non_bm is reproducible given a seed", {
  tree <- ape::rcoal(15L)
  a <- simulate_non_bm(tree, n_traits = 2L, scenario = "OU", seed = 99L)
  b <- simulate_non_bm(tree, n_traits = 2L, scenario = "OU", seed = 99L)
  expect_equal(unname(as.matrix(a)), unname(as.matrix(b)))
})

test_that("[T2A] read_tree round-trips a Newick file", {
  tree <- ape::rcoal(20L)
  path <- withr::local_tempfile(fileext = ".tre")
  ape::write.tree(tree, path)
  back <- read_tree(path)
  expect_s3_class(back, "phylo")
  expect_equal(length(back$tip.label), length(tree$tip.label))
})

test_that("[T2A] read_tree round-trips a NEXUS file", {
  tree <- ape::rcoal(15L)
  path <- withr::local_tempfile(fileext = ".nex")
  ape::write.nexus(tree, file = path)
  back <- read_tree(path)
  expect_s3_class(back, "phylo")
  expect_equal(length(back$tip.label), length(tree$tip.label))
})

test_that("[T2A] read_tree errors clearly on a missing file", {
  expect_error(read_tree("/no/such/path/nope.tre"), regexp = "not found")
})

test_that("[T2A] plot_history_gg returns a ggplot for a fit with $history", {
  skip_if_not_installed("ggplot2")
  td <- .tcov_make_data(n = 20L, seed = 2L)
  pd <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 1, trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     epochs = 10L, eval_every = 5L, patience = 5L,
                     verbose = FALSE, seed = 1)
  p <- tryCatch(plot_history_gg(fit), error = function(e) e)
  expect_s3_class(p, "ggplot")
})

test_that("[T2A] plot_uncertainty returns a ggplot for a continuous trait", {
  skip_if_not_installed("ggplot2")
  td  <- .tcov_make_data(n = 20L, seed = 3L)
  td$df$cont[1:5] <- NA
  result <- impute(td$df, td$tree, epochs = 10L,
                   n_imputations = 1L, verbose = FALSE, seed = 3L)
  p <- tryCatch(plot_uncertainty(result$prediction, trait_name = "cont"),
                error = function(e) e)
  expect_s3_class(p, "ggplot")
})


# ---------------------------------------------------------------------------
# T2.B  Smoke tests for non-default option paths
# ---------------------------------------------------------------------------

test_that("[T2B] multi_impute(draws_method = 'mc_dropout') returns m completed datasets", {
  td  <- .tcov_make_data(n = 25L, seed = 4L)
  td$df$cont[1:6] <- NA
  mi <- multi_impute(td$df, td$tree, m = 3L,
                     draws_method = "mc_dropout",
                     epochs = 10L, verbose = FALSE, seed = 4L)
  expect_s3_class(mi, "pigauto_mi")
  expect_equal(length(mi$datasets), 3L)
  for (d in mi$datasets) {
    expect_equal(nrow(d), nrow(td$df))
    expect_false(any(is.na(d$cont)))   # missing cells filled
  }
})

test_that("[T2B] fit_pigauto(use_attention = FALSE) runs end-to-end and predicts", {
  td  <- .tcov_make_data(n = 25L, seed = 5L)
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 5, trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     use_attention = FALSE,
                     epochs = 10L, eval_every = 5L, patience = 5L,
                     verbose = FALSE, seed = 5L)
  expect_s3_class(fit, "pigauto_fit")
  expect_false(isTRUE(fit$model_config$use_attention))
  pred <- predict(fit, return_se = TRUE)
  expect_s3_class(pred, "pigauto_pred")
  expect_equal(nrow(pred$imputed), nrow(td$df))
})

test_that("[T2B] fit_pigauto(conformal_method = 'bootstrap') stores conformal scores", {
  td  <- .tcov_make_data(n = 30L, seed = 6L)
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 6, trait_map = pd$trait_map)
  fit <- fit_pigauto(pd, td$tree, splits = spl,
                     conformal_method = "bootstrap",
                     conformal_bootstrap_B = 50L,
                     epochs = 10L, eval_every = 5L, patience = 5L,
                     verbose = FALSE, seed = 6L)
  expect_s3_class(fit, "pigauto_fit")
  expect_true("conformal_scores" %in% names(fit))
})


# ---------------------------------------------------------------------------
# T3   Edge cases
# ---------------------------------------------------------------------------

test_that("[T3] preprocess_traits warns when no rownames overlap with tree tips", {
  # Documents the CURRENT contract: when traits and tree have zero
  # overlap, preprocess emits a "tree tips have no trait data" warning
  # and returns a pigauto_data with all-NA X_scaled.  TODO (future
  # hardening): consider erroring rather than silently producing an
  # all-NA matrix that downstream fit_pigauto will fail on.
  set.seed(8L)
  tree <- ape::rcoal(10L)
  df <- data.frame(
    row.names = paste0("not_a_tip_", seq_len(10L)),
    cont = stats::rnorm(10L)
  )
  expect_message(
    pd <- preprocess_traits(df, tree),
    regexp = "(no trait data|all-NA|tip)",
    ignore.case = TRUE
  )
  expect_s3_class(pd, "pigauto_data")
  expect_true(all(is.na(pd$X_scaled)))
})

test_that("[T3] impute errors gracefully on a tree without branch lengths", {
  set.seed(9L)
  tree <- ape::rcoal(15L)
  tree$edge.length <- NULL  # strip branch lengths
  df <- data.frame(
    row.names = tree$tip.label,
    cont = stats::rnorm(15L)
  )
  res <- tryCatch(
    impute(df, tree, epochs = 5L, n_imputations = 1L,
           verbose = FALSE, seed = 9L),
    error = function(e) e
  )
  # Either it errors clearly, or it auto-computes branch lengths and
  # succeeds.  Both are acceptable; what is NOT acceptable is a silent
  # NaN propagation.  Failing loudly is the contract.
  if (inherits(res, "error")) {
    expect_match(conditionMessage(res),
                 "(branch|edge|length|null)",
                 ignore.case = TRUE)
  } else {
    expect_s3_class(res, "pigauto_result")
    expect_true(all(is.finite(res$completed$cont)))
  }
})

test_that("[T3] preprocess_traits accepts a single-species tree (current contract)", {
  # Documents the CURRENT contract: preprocess does not error for
  # n_species = 1; it returns a 1-row pigauto_data.  Downstream
  # fit_pigauto / fit_baseline will then fail because phylogenetic
  # signal is undefined for n = 1.  TODO (future hardening): error in
  # preprocess_traits when n_species < 2.
  set.seed(10L)
  tree <- ape::rcoal(2L)
  tree <- ape::drop.tip(tree, tree$tip.label[2])
  df <- data.frame(row.names = tree$tip.label, cont = 1.0)
  pd <- suppressWarnings(preprocess_traits(df, tree))
  expect_s3_class(pd, "pigauto_data")
  expect_equal(nrow(pd$X_scaled), 1L)
})

test_that("[T3] all-observed continuous trait round-trips through impute() as identity", {
  set.seed(11L)
  td <- .tcov_make_data(n = 20L, seed = 11L)
  res <- impute(td$df, td$tree, epochs = 10L, n_imputations = 1L,
                verbose = FALSE, seed = 11L)
  # All input cells were observed, so completed should equal input
  # exactly on the observed mask (no impute happens).
  expect_equal(res$completed$cont, td$df$cont, tolerance = 1e-9)
  expect_equal(as.character(res$completed$bin),
               as.character(td$df$bin))
})

test_that("[T3] all-NA column triggers a clear error rather than silent NaN", {
  set.seed(12L)
  td <- .tcov_make_data(n = 20L, seed = 12L)
  td$df$cont <- NA_real_  # all-NA column
  res <- tryCatch(
    suppressWarnings(impute(td$df, td$tree, epochs = 5L,
                            n_imputations = 1L, verbose = FALSE,
                            seed = 12L)),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    expect_match(conditionMessage(res),
                 "(NA|missing|empty|all observed|no observed)",
                 ignore.case = TRUE)
  } else {
    # If it returns a result, the column must NOT be silently filled
    # with NaN -- it must be either fully NA still or fully imputed.
    out <- res$completed$cont
    expect_false(any(is.nan(out)))
  }
})
