test_that("pool_mi aligns named coefficients and covariance by term", {
  make_fit <- function(beta, vars) {
    V <- diag(vars)
    dimnames(V) <- list(names(beta), names(beta))
    list(beta = beta, V = V)
  }
  fits <- list(
    make_fit(c("mu:(Intercept)" = 1, "mu:x" = 2), c(0.04, 0.09)),
    make_fit(c("mu:x" = 4, "mu:(Intercept)" = 3), c(0.16, 0.25))
  )
  out <- suppressWarnings(pool_mi(
    fits, coef_fun = function(f) f$beta, vcov_fun = function(f) f$V
  ))
  expect_identical(out$term, c("mu:(Intercept)", "mu:x"))
  expect_equal(out$estimate, c(2, 3))
  expect_equal(out$std.error, sqrt(c(3.145, 3.125)))
})

test_that("pool_mi reports absent and aliased terms explicitly", {
  absent <- list(list(beta = c(a = 1, b = 2)), list(beta = c(a = 2)))
  expect_error(pool_mi(absent, coef_fun = function(f) f$beta,
                       vcov_fun = function(f) diag(length(f$beta))),
               "Coefficient names differ.*fit 2.*Missing terms: b")

  aliased <- list(list(beta = c(a = 1, b = 2)),
                  list(beta = c(a = 2, b = NA_real_)))
  expect_error(pool_mi(aliased, coef_fun = function(f) f$beta,
                       vcov_fun = function(f) diag(length(f$beta))),
               "finite, uniquely named")

  missing_name <- list(list(beta = c(a = 1, b = 2)),
                       list(beta = stats::setNames(c(2, 3), c("a", NA))))
  expect_error(pool_mi(missing_name, coef_fun = function(f) f$beta,
                       vcov_fun = function(f) diag(length(f$beta))),
               "finite, uniquely named")
})

test_that("pool_mi rejects bad tidy standard errors", {
  fits <- list(list(i = 1), list(i = 2))
  bad_se <- function(f) data.frame(
    term = "x", estimate = f$i, std.error = if (f$i == 1) Inf else 0.2
  )
  expect_error(pool_mi(fits, tidy_fun = bad_se), "finite numeric estimates")
})

test_that("optional package adapters preserve drmTMB component names", {
  blocks <- list(mu = c("(Intercept)" = 1, x = 2),
                 sigma = c("(Intercept)" = 3, x = 4))
  testthat::local_mocked_bindings(
    .pool_mi_drm_coef_blocks = function(fit) blocks,
    .package = "pigauto"
  )
  fit <- structure(list(opt = list(convergence = 0L), sdr = list(pdHess = TRUE)),
                   class = "drmTMB")
  expect_identical(unname(.pool_mi_auto_coef(fit)), c(1, 2, 3, 4))
  expect_identical(names(.pool_mi_auto_coef(fit)), c(
    "mu:(Intercept)", "mu:x", "sigma:(Intercept)", "sigma:x"
  ))
})

test_that("drmTMB adapter rejects duplicate component names", {
  testthat::local_mocked_bindings(
    .pool_mi_drm_coef_blocks = function(fit) {
      list(mu = c(x = 1), mu = c(x = 2))
    },
    .package = "pigauto"
  )
  fit <- structure(list(opt = list(convergence = 0L), sdr = list(pdHess = TRUE)),
                   class = "drmTMB")
  expect_error(.pool_mi_auto_coef(fit), "uniquely named")
})

test_that("gllvmTMB adapter recognizes fits with a leading wrapper class", {
  testthat::local_mocked_bindings(
    .pool_mi_gllvm_tidy = function(fit) {
      data.frame(term = "beta", estimate = 1, std.error = 0.2)
    },
    .package = "pigauto"
  )
  fit <- structure(list(), class = c("saved_fit_wrapper", "gllvmTMB_multi"))
  expect_identical(.pool_mi_auto_coef(fit), c(beta = 1))
})

test_that("automatic drmTMB and gllvmTMB adapters reject failed convergence or Hessians", {
  drm_bad_convergence <- structure(list(opt = list(convergence = 1L)),
                                   class = "drmTMB")
  drm_bad_hessian <- structure(list(opt = list(convergence = 0L),
                                   sdr = list(pdHess = FALSE)),
                               class = "drmTMB")
  gllvm_bad_convergence <- structure(list(opt = list(convergence = 1L)),
                                     class = "gllvmTMB_multi")
  gllvm_bad_hessian <- structure(list(opt = list(convergence = 0L),
                                     sd_report = list(pdHess = FALSE)),
                                 class = "gllvmTMB_multi")

  expect_error(.pool_mi_auto_coef(drm_bad_convergence), "convergence code")
  expect_error(.pool_mi_auto_coef(drm_bad_hessian), "positive-definite Hessian")
  expect_error(.pool_mi_auto_coef(gllvm_bad_convergence), "convergence code")
  expect_error(.pool_mi_auto_coef(gllvm_bad_hessian), "positive-definite Hessian")
})

test_that("optional gllvm adapter reports its runtime namespace requirement", {
  expect_error(.pool_mi_require_namespace("pigauto_no_such_adapter_pkg"),
               "requires the package.*installed and loadable at runtime")
})

test_that("gllvmTMB fixture pools only fixed terms by name", {
  testthat::local_mocked_bindings(
    .pool_mi_gllvm_tidy = function(fit) {
      i <- fit$i
      out <- data.frame(term = c("traitB:x", "traitA:x"),
                        estimate = c(i + 2, i),
                        std.error = c(0.3, 0.2))
      if (i == 2L) out <- out[2:1, , drop = FALSE]
      out
    },
    .package = "pigauto"
  )
  make_fit <- function(i) structure(
    list(i = i, opt = list(convergence = 0L),
         sd_report = list(pdHess = TRUE)),
    class = "gllvmTMB_multi"
  )
  out <- suppressWarnings(pool_mi(list(make_fit(1L), make_fit(2L))))
  expect_identical(out$term, c("traitB:x", "traitA:x"))
  expect_equal(out$estimate, c(3.5, 1.5))
  expect_equal(out$std.error, sqrt(c(0.09, 0.04) + 0.75))
})

test_that("drmTMB fixture rejects invalid component terms before pooling", {
  testthat::local_mocked_bindings(
    .pool_mi_drm_coef_blocks = function(fit) {
      list(mu = c(x = 1), sigma = stats::setNames(2, NA_character_))
    },
    .package = "pigauto"
  )
  fit <- structure(list(opt = list(convergence = 0L),
                        sdr = list(pdHess = TRUE)), class = "drmTMB")
  expect_error(.pool_mi_auto_coef(fit), "uniquely named")
})
