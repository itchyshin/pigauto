test_that("baseline-only package use does not load torch in a fresh session", {
  skip_if_not_installed("callr")
  if ("pkgload" %in% loadedNamespaces() && pkgload::is_dev_package("pigauto")) {
    skip("requires the installed package; pkgload eagerly loads Imports")
  }
  lib <- dirname(system.file(package = "pigauto"))
  state <- callr::r(function(lib) {
    library(pigauto, lib.loc = lib)
    after_load <- "torch" %in% loadedNamespaces()
    set.seed(19)
    tree <- ape::rtree(20)
    traits <- data.frame(mass = seq_len(20) + 0.1)
    rownames(traits) <- tree$tip.label
    traits$mass[c(2, 5)] <- NA_real_
    result <- impute(traits, tree, gnn = FALSE, verbose = FALSE, seed = 19)
    list(after_load = after_load, after_fit = "torch" %in% loadedNamespaces(),
         completed = all(is.finite(result$completed$mass)))
  }, args = list(lib = lib))
  expect_false(state$after_load)
  expect_false(state$after_fit)
  expect_true(state$completed)
})
