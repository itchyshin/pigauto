# ---------------------------------------------------------------------------
# fit_baseline_bace(final_imp =) -- opt-in proper-MI path
#
# The default path (final_imp = FALSE) summarises the bace_imp() chain
# datasets and must stay exactly as it was.  The opt-in path appends
# BACE::bace_final_imp() and summarises its n_final independent draws
# instead.  Everything that actually calls BACE is guarded twice: once on
# the package being installed, once on the call succeeding, because BACE
# is Suggests-only and its MCMCglmm fits can fail on tiny synthetic data.
# ---------------------------------------------------------------------------

.bfi_make_data <- function(n = 20L, seed = 7L) {
  set.seed(seed)
  tree <- ape::rcoal(n)
  df <- data.frame(
    row.names = tree$tip.label,
    cont = stats::rnorm(n),
    bin  = factor(sample(c("a", "b"), n, replace = TRUE))
  )
  list(tree = tree, df = df)
}

.bfi_try_baseline <- function(...) {
  tryCatch(fit_baseline_bace(...), error = function(e) e)
}


test_that("final_imp defaults to FALSE and n_final to 15L", {
  fmls <- formals(fit_baseline_bace)
  expect_false(eval(fmls$final_imp))
  expect_identical(eval(fmls$n_final), 15L)
})


test_that("final_imp rejects non-scalar / NA flags", {
  skip_if_not_installed("BACE")
  td <- .bfi_make_data()
  pd <- preprocess_traits(td$df, td$tree)
  expect_error(
    fit_baseline_bace(pd, td$tree, final_imp = NA),
    "single TRUE or FALSE"
  )
  expect_error(
    fit_baseline_bace(pd, td$tree, final_imp = c(TRUE, FALSE)),
    "single TRUE or FALSE"
  )
})


test_that("n_final is validated only when final_imp = TRUE", {
  skip_if_not_installed("BACE")
  td <- .bfi_make_data()
  pd <- preprocess_traits(td$df, td$tree)
  expect_error(
    fit_baseline_bace(pd, td$tree, final_imp = TRUE, n_final = 1L),
    "n_final"
  )
  # Ignored on the default path: an absurd n_final must not trip a check
  # before BACE is ever called.
  res <- .bfi_try_baseline(pd, td$tree, runs = 1L, nitt = 200L,
                           burnin = 50L, thin = 5L, n_final = -3L,
                           verbose = FALSE)
  if (inherits(res, "error")) {
    expect_false(grepl("n_final", conditionMessage(res)))
  } else {
    expect_true(all(c("mu", "se") %in% names(res)))
  }
})


test_that("default path is unchanged and reproducible", {
  skip_if_not_installed("BACE")
  td  <- .bfi_make_data()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 7L, trait_map = pd$trait_map)

  set.seed(101L)
  a <- .bfi_try_baseline(pd, td$tree, splits = spl, runs = 2L, nitt = 200L,
                         burnin = 50L, thin = 5L, verbose = FALSE)
  if (inherits(a, "error")) {
    skip(paste("BACE chain not runnable here:", conditionMessage(a)))
  }
  set.seed(101L)
  b <- .bfi_try_baseline(pd, td$tree, splits = spl, runs = 2L, nitt = 200L,
                         burnin = 50L, thin = 5L, verbose = FALSE)
  skip_if(inherits(b, "error"), "BACE chain not runnable here")

  # Same seed, same defaults -> identical output. This is the guard that
  # the final_imp branch did not perturb the default path.
  expect_equal(a$mu, b$mu)
  expect_equal(a$se, b$se)
  expect_equal(nrow(a$mu), nrow(pd$X_scaled))
  expect_equal(ncol(a$mu), ncol(pd$X_scaled))
})


test_that("final_imp = TRUE returns same-shaped finite mu / se", {
  skip_if_not_installed("BACE")
  skip_if_not(
    "bace_final_imp" %in% getNamespaceExports("BACE"),
    "installed BACE has no bace_final_imp()"
  )
  td  <- .bfi_make_data()
  pd  <- preprocess_traits(td$df, td$tree)
  spl <- make_missing_splits(pd$X_scaled, seed = 7L, trait_map = pd$trait_map)

  res <- .bfi_try_baseline(pd, td$tree, splits = spl, runs = 2L, nitt = 200L,
                           burnin = 50L, thin = 5L, final_imp = TRUE,
                           n_final = 3L, verbose = FALSE)
  if (inherits(res, "error")) {
    skip(paste("bace_final_imp not runnable here:", conditionMessage(res)))
  }

  expect_type(res, "list")
  expect_true(all(c("mu", "se") %in% names(res)))
  expect_equal(dim(res$mu), dim(pd$X_scaled))
  expect_equal(dim(res$se), dim(pd$X_scaled))
  expect_true(all(is.finite(res$mu)))
  expect_true(all(is.finite(res$se)))
  expect_true(all(res$se >= 0))
  expect_identical(rownames(res$mu), pd$species_names)
})


test_that("final_imp = TRUE errors clearly when BACE lacks bace_final_imp", {
  skip_if_not_installed("BACE")
  skip_if(
    "bace_final_imp" %in% getNamespaceExports("BACE"),
    "installed BACE does export bace_final_imp(); nothing to assert"
  )
  td <- .bfi_make_data()
  pd <- preprocess_traits(td$df, td$tree)
  expect_error(
    fit_baseline_bace(pd, td$tree, final_imp = TRUE),
    "bace_final_imp"
  )
})
