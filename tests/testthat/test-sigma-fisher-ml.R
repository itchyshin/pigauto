# Tests for `sigma_method = c("single_pass", "fisher_ml")` on
# fit_mvn_bm_inhouse() / fit_joint_solver() (R/joint_mvn_solver.R).
#
# Ported (adapted, not cherry-picked) from commit e7ca41c (2026-05-16,
# never merged). See docs/dev-log/2026-08-16-continuous-gap-diagnosis.md
# for why it was revived, and the header comment above
# `.mvn_par_to_chol()` in R/joint_mvn_solver.R for how the ported
# machinery coexists with the three post-e7ca41c "Bug fixes: in-house
# Sigma solver" (NEWS.md, v0.9.2).

make_correlated_liability <- function(seed, n = 30) {
  set.seed(seed)
  tree <- ape::rcoal(n)
  C <- ape::vcv(tree)
  Lc <- chol(C)
  x1 <- as.numeric(t(Lc) %*% stats::rnorm(n))
  x2 <- 0.8 * x1 + as.numeric(t(Lc) %*% stats::rnorm(n)) * 0.3
  L <- cbind(t1 = x1, t2 = x2)
  rownames(L) <- tree$tip.label
  L[sample(n, 4), "t1"] <- NA
  L[sample(n, 5), "t2"] <- NA
  list(tree = tree, L = L)
}

test_that("sigma_method defaults to 'single_pass' and is byte-identical whether passed explicitly or not", {
  d <- make_correlated_liability(seed = 2)

  fit_default  <- fit_mvn_bm_inhouse(d$L, tree = d$tree)
  fit_explicit <- fit_mvn_bm_inhouse(d$L, tree = d$tree, sigma_method = "single_pass")

  expect_identical(fit_default, fit_explicit)
})

test_that("sigma_method = 'single_pass' numeric regression on a fixed-seed L/tree", {
  # Deterministic linear algebra (no torch RNG involved) -- cross-platform
  # stable at generous tolerance. Values computed at write-time on the
  # authoring machine via fit_mvn_bm_inhouse(L, tree = tree) with
  # seed = 42, n = 25 (see make_correlated_liability()).
  d <- make_correlated_liability(seed = 42, n = 25)

  fit <- fit_mvn_bm_inhouse(d$L, tree = d$tree)

  expect_equal(
    fit$pars$phylocov,
    matrix(c(1.710944020304429, 0.936012345856949,
             0.936012345856949, 1.008967782985983),
           nrow = 2, ncol = 2),
    tolerance = 1e-6
  )
  expect_equal(sum(fit$anc_recon), -36.3778205620947, tolerance = 1e-6)
  expect_equal(sum(fit$anc_var), 0.124196977768639, tolerance = 1e-6)
  expect_identical(fit$n_iter, 0L)
  expect_true(fit$converged)
})

test_that("sigma_method = 'fisher_ml' returns a valid SPD phylocov and finite outputs, and differs from single_pass on correlated data", {
  d <- make_correlated_liability(seed = 3, n = 40)

  fit_fm <- fit_mvn_bm_inhouse(d$L, tree = d$tree, sigma_method = "fisher_ml")
  fit_sp <- fit_mvn_bm_inhouse(d$L, tree = d$tree, sigma_method = "single_pass")

  Sig <- fit_fm$pars$phylocov
  expect_true(isSymmetric(Sig, tol = 1e-6))
  evals <- eigen(Sig, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(evals > 0))

  expect_true(all(is.finite(fit_fm$anc_recon)))
  expect_true(all(is.finite(fit_fm$anc_var)))

  # Different Sigma-estimation strategies (species-iid Fisher-ML vs
  # R-credited Kronecker MLE) -- should differ measurably on data with
  # real phylogenetic signal, not coincide to machine precision.
  expect_gt(max(abs(Sig - fit_sp$pars$phylocov)), 1e-4)
})

test_that("K=1 short-circuit is identical for single_pass and fisher_ml", {
  set.seed(5)
  tree <- ape::rcoal(20)
  C <- ape::vcv(tree)
  Lc <- chol(C)
  x <- as.numeric(t(Lc) %*% stats::rnorm(20))
  L <- matrix(x, ncol = 1, dimnames = list(tree$tip.label, "t1"))
  L[sample(20, 3), 1] <- NA_real_

  fit_sp <- fit_mvn_bm_inhouse(L, tree = tree, sigma_method = "single_pass")
  fit_fm <- fit_mvn_bm_inhouse(L, tree = tree, sigma_method = "fisher_ml")

  expect_identical(fit_sp, fit_fm)
})

test_that("sigma_method = 'fisher_ml' falls back to the single_pass Sigma on optim() non-convergence, with a warning", {
  d <- make_correlated_liability(seed = 4, n = 25)

  testthat::local_mocked_bindings(
    optim = function(...) list(par = numeric(0), convergence = 1L),
    .package = "stats"
  )

  expect_warning(
    fit_fallback <- fit_mvn_bm_inhouse(d$L, tree = d$tree, sigma_method = "fisher_ml"),
    "falling back to the single_pass Sigma estimate"
  )

  fit_single_pass <- fit_mvn_bm_inhouse(d$L, tree = d$tree, sigma_method = "single_pass")
  expect_equal(fit_fallback$pars$phylocov, fit_single_pass$pars$phylocov)
})

test_that("sigma_method = 'fisher_ml' falls back on a hard optim() error too, with a warning", {
  d <- make_correlated_liability(seed = 7, n = 25)

  testthat::local_mocked_bindings(
    optim = function(...) stop("boom"),
    .package = "stats"
  )

  expect_warning(
    fit_fallback <- fit_mvn_bm_inhouse(d$L, tree = d$tree, sigma_method = "fisher_ml"),
    "falling back to the single_pass Sigma estimate"
  )

  fit_single_pass <- fit_mvn_bm_inhouse(d$L, tree = d$tree, sigma_method = "single_pass")
  expect_equal(fit_fallback$pars$phylocov, fit_single_pass$pars$phylocov)
})

test_that("fit_joint_solver() threads sigma_method through to fit_mvn_bm_inhouse()", {
  d <- make_correlated_liability(seed = 6)

  out_sp <- fit_joint_solver(d$L, d$tree, joint_solver = "inhouse", sigma_method = "single_pass")
  out_fm <- fit_joint_solver(d$L, d$tree, joint_solver = "inhouse", sigma_method = "fisher_ml")

  expect_identical(out_sp, fit_mvn_bm_inhouse(d$L, tree = d$tree, sigma_method = "single_pass"))
  expect_gt(max(abs(out_sp$pars$phylocov - out_fm$pars$phylocov)), 1e-4)
})

test_that("fit_joint_solver() default sigma_method matches the fit_mvn_bm_inhouse() default", {
  d <- make_correlated_liability(seed = 8)

  out_default <- fit_joint_solver(d$L, d$tree, joint_solver = "inhouse")
  out_explicit <- fit_mvn_bm_inhouse(d$L, tree = d$tree, sigma_method = "single_pass")

  expect_identical(out_default, out_explicit)
})
