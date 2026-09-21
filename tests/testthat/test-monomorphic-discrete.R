# Regression tests for discrete traits that are monomorphic among the
# observed (unmasked) cells.
#
# Two shapes, both reported from the four-arm imputation campaign at
# lambda = 1, where maximum phylogenetic signal plus a threshold cut can
# leave every observed species on one side:
#
#   (A) the factor carries a SINGLE declared level (nlevels == 1).  Type
#       detection routes this to "categorical" with K = 1, so the one-hot
#       block is a single column.  torch's `[` drops the trailing axis on a
#       length-1 column index, so `$argmax(dim = 2L)` in the categorical
#       loss received a 1-D tensor and raised
#       "Dimension out of range (expected to be in range of [-1, 1], but got 2)".
#
#   (B) the factor declares K >= 2 levels but every observed cell falls in
#       one class.  This path never errored; it is pinned here so the
#       "predict the observed class" behaviour cannot regress.
#
# A monomorphic trait has an obvious answer -- predict that class -- so
# neither shape may error.

make_monomorphic <- function(n, levs, seed = 7L, ordered = FALSE) {
  set.seed(seed)
  tree <- ape::rtree(n)
  df <- data.frame(
    disc = factor(rep(levs[1], n), levels = levs, ordered = ordered),
    cont = as.numeric(ape::rTraitCont(tree)[tree$tip.label]),
    row.names = tree$tip.label
  )
  miss <- sample(n, size = round(0.3 * n))
  df$disc[miss] <- NA
  list(tree = tree, df = df, miss = miss)
}

test_that("compute_mixed_loss keeps the K axis for a single-level categorical", {
  # Tightest regression on the root cause: a K = 1 one-hot block must not
  # lose its trailing dimension before argmax(dim = 2L).
  n <- 6L
  trait_map <- list(disc = list(
    type = "categorical", n_latent = 1L, latent_cols = 1L, levels = "a"
  ))
  pred  <- torch::torch_randn(n, 1L)
  truth <- torch::torch_ones(n, 1L)
  corrupt_mask <- torch::torch_tensor(
    matrix(c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE), ncol = 1L)
  )

  loss <- compute_mixed_loss(pred, truth, corrupt_mask, trait_map)

  expect_true(inherits(loss, "torch_tensor"))
  expect_true(is.finite(as.numeric(loss$item())))
})

test_that("composite_val_loss keeps the K axis for a single-level categorical", {
  n <- 6L
  trait_map <- list(disc = list(
    type = "categorical", n_latent = 1L, latent_cols = 1L, levels = "a"
  ))
  pred  <- torch::torch_randn(n, 1L)
  truth <- torch::torch_ones(n, 1L)
  val_mask <- torch::torch_tensor(
    matrix(c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE), ncol = 1L)
  )

  expect_true(is.finite(composite_val_loss(pred, truth, val_mask, trait_map)))
})

test_that("impute() predicts the single class for a one-level factor", {
  fx <- make_monomorphic(60L, levs = "a")

  res <- expect_no_error(
    impute(fx$df, fx$tree, epochs = 20L, verbose = FALSE, seed = 1L)
  )

  filled <- res$completed$disc
  expect_false(anyNA(filled))
  expect_true(all(as.character(filled) == "a"))
  # K = 1: the only class carries all the mass.
  probs <- res$prediction$probabilities$disc
  expect_equal(dim(probs), c(60L, 1L))
  expect_equal(unname(probs[fx$miss, 1]), rep(1, length(fx$miss)),
               tolerance = 1e-6)
})

test_that("impute() predicts the single class for a one-level ordered factor", {
  fx <- make_monomorphic(60L, levs = "a", ordered = TRUE)

  res <- expect_no_error(
    impute(fx$df, fx$tree, epochs = 20L, verbose = FALSE, seed = 1L)
  )

  expect_false(anyNA(res$completed$disc))
  expect_true(all(as.character(res$completed$disc) == "a"))
})

test_that("impute() handles a K-level factor monomorphic among observed cells", {
  # Shape (B): both binary (K = 2) and categorical (K = 3) declare unseen
  # levels, so the observed class should win with high -- not certain --
  # probability, and the unseen levels must survive on the output factor.
  for (levs in list(c("a", "b"), c("a", "b", "c"))) {
    fx <- make_monomorphic(60L, levs = levs)

    res <- expect_no_error(
      impute(fx$df, fx$tree, epochs = 20L, verbose = FALSE, seed = 1L)
    )

    filled <- res$completed$disc
    expect_false(anyNA(filled))
    expect_true(all(as.character(filled) == "a"))
    expect_equal(levels(filled), levs)

    probs <- res$prediction$probabilities$disc
    # binary stores P(level 2) as a vector; categorical an n x K matrix.
    p_observed <- if (is.null(dim(probs))) 1 - probs else probs[, 1]
    expect_true(all(p_observed[fx$miss] > 0.5))
    expect_true(all(p_observed[fx$miss] <= 1))
  }
})

test_that("the whole pipeline survives a single-level categorical trait", {
  # fit_baseline / fit_pigauto / predict / evaluate each index the K-block
  # separately, so walk them explicitly rather than only through impute().
  fx <- make_monomorphic(60L, levs = "a")

  pd <- preprocess_traits(fx$df, fx$tree, log_transform = FALSE)
  expect_equal(pd$trait_map$disc$type, "categorical")
  expect_equal(pd$trait_map$disc$n_latent, 1L)

  graph  <- build_phylo_graph(fx$tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.25,
                                seed = 1L, trait_map = pd$trait_map)
  bl <- fit_baseline(pd, fx$tree, splits = splits, graph = graph)
  ft <- expect_no_error(
    fit_pigauto(pd, fx$tree, splits = splits, graph = graph, baseline = bl,
                epochs = 20L, verbose = FALSE, seed = 1L)
  )
  pr <- expect_no_error(predict(ft))
  ev <- expect_no_error(
    evaluate_imputation(pr, pd$X_scaled, splits, trait_map = pd$trait_map)
  )
  expect_s3_class(ev, "data.frame")
})
