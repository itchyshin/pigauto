#!/usr/bin/env Rscript
# Closer reproducer of the multi_proportion crash, mirroring bench_multi_proportion
# more carefully.  Uses the bench's exact pipeline pieces.
suppressPackageStartupMessages({
  pkg_path <- "/Users/z3437171/Dropbox/Github Local/pigauto"
  ok <- tryCatch({ devtools::load_all(pkg_path, quiet = TRUE); TRUE },
                  error = function(e) {
                    cat("LOAD_ALL FAILED:", conditionMessage(e), "\n"); FALSE })
  if (!ok) quit(status = 1L, save = "no")
})

set.seed(122)  # 1 (rep) * 100 + 22 (signal_0.2 index in bench)
res <- tryCatch({
  tree <- ape::rtree(300)  # bench's n_species
  K <- 5L
  signal <- 0.2

  df <- simulate_multi_proportion_traits(tree, K = K, signal = signal, seed = 122)
  group_cols <- names(df)
  pd  <- preprocess_traits(df, tree,
                            multi_proportion_groups = setNames(list(group_cols), "comp"))
  spl <- make_missing_splits(pd$X_scaled, missing_frac = 0.25,
                              seed = 122, trait_map = pd$trait_map)
  graph <- build_phylo_graph(tree, k_eigen = "auto")
  bl <- fit_baseline(pd, tree, splits = spl, graph = graph)
  fit <- fit_pigauto(pd, tree, splits = spl, baseline = bl,
                      graph = graph, epochs = 30L, verbose = FALSE, seed = 122)
  pred <- stats::predict(fit, return_se = FALSE, n_imputations = 1L)
  TRUE
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n"); FALSE
})

if (isTRUE(res)) {
  cat("PASS\n")
  quit(status = 0L, save = "no")
} else {
  cat("FAIL\n")
  quit(status = 1L, save = "no")
}
