test_that("predict.pigauto_fit accepts baseline_override and uses it in place of object$baseline", {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "libtorch not installed")

  set.seed(42L)
  tree <- ape::rcoal(30L)
  df   <- pigauto:::simulate_bm_traits(tree, n_traits = 2L, seed = 1L)

  pd    <- pigauto::preprocess_traits(df, tree, log_transform = FALSE)
  spl   <- pigauto::make_missing_splits(pd$X_scaled, missing_frac = 0.20,
                                         seed = 1L, trait_map = pd$trait_map)
  graph <- pigauto::build_phylo_graph(tree, k_eigen = 4L)
  bl    <- pigauto::fit_baseline(pd, tree, splits = spl, graph = graph)
  fit   <- pigauto::fit_pigauto(pd, tree, splits = spl, baseline = bl, graph = graph,
                                 epochs = 30L, verbose = FALSE, seed = 1L)

  # default: uses object$baseline
  pred_default <- stats::predict(fit, return_se = FALSE, n_imputations = 1L)
  expect_s3_class(pred_default, "pigauto_pred")

  # override with a shifted baseline -- predictions must differ
  bl2   <- bl
  bl2$mu <- bl$mu + 5.0
  pred_shifted <- stats::predict(fit, return_se = FALSE, n_imputations = 1L,
                                  baseline_override = bl2)
  expect_s3_class(pred_shifted, "pigauto_pred")

  # Imputed latent values should reflect the shift (gate is open at some level)
  expect_false(
    identical(pred_default$imputed_latent, pred_shifted$imputed_latent),
    info = "baseline_override should change predictions when gate is not fully closed"
  )
})

test_that("resolve_reference_tree returns user-supplied tree when given", {
  trees <- list(ape::rcoal(20L, tip.label = paste0("t", 1:20)),
                ape::rcoal(20L, tip.label = paste0("t", 1:20)),
                ape::rcoal(20L, tip.label = paste0("t", 1:20)))
  user_tree <- ape::rcoal(20L, tip.label = paste0("t", 1:20))

  out <- pigauto:::resolve_reference_tree(trees, reference_tree = user_tree)
  expect_identical(out, user_tree)
})

test_that("resolve_reference_tree uses phangorn MCC when phangorn available", {
  skip_if_not_installed("phangorn")

  set.seed(1L)
  trees <- lapply(1:5, function(i) ape::rcoal(20L, tip.label = paste0("t", 1:20)))
  class(trees) <- "multiPhylo"

  out <- pigauto:::resolve_reference_tree(trees, reference_tree = NULL)
  expect_s3_class(out, "phylo")
  # MCC must be identical (in structure) to one of the input trees
  any_match <- any(vapply(trees, function(tr) {
    isTRUE(all.equal(tr$edge, out$edge)) &&
    isTRUE(all.equal(tr$tip.label, out$tip.label))
  }, logical(1)))
  expect_true(any_match, info = "MCC tree must be one of the input trees")
})

test_that("resolve_reference_tree warns and falls back when phangorn missing", {
  trees <- list(ape::rcoal(15L), ape::rcoal(15L))
  local_mocked_bindings(
    requireNamespace = function(pkg, ...) !(identical(pkg, "phangorn")),
    .package = "base"
  )

  expect_warning(
    out <- pigauto:::resolve_reference_tree(trees, reference_tree = NULL),
    "phangorn"
  )
  expect_identical(out, trees[[1]])
})
