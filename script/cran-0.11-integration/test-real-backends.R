# End-to-end pool_mi() backend tests: multi_impute(draws_method = "posterior")
# -> with_imputations() -> pool_mi(), so the provenance path is exercised
# (the backend tests in test-multi-impute.R pass bare lists of fits).

pmb_mi <- function(n = 60L, m = 3L) {
  set.seed(21L)
  x <- stats::rnorm(n)
  z <- stats::rnorm(n)
  alpha <- c(0.5, -0.5, 1, 0)
  beta <- c(0.7, -0.4, 0.3, 0.5)
  loading <- c(0.6, 0.4, -0.5, 0.7)
  Y <- matrix(alpha, n, 4L, byrow = TRUE) + outer(x, beta) +
    outer(z, loading) + 0.8 * matrix(stats::rnorm(n * 4L), n)
  df <- data.frame(x = x, Y, row.names = paste0("site", seq_len(n)))
  names(df) <- c("x", paste0("y", seq_len(4L)))
  df$x[c(2L, 5L, 9L, 14L)] <- NA_real_
  multi_impute_analysis(
    data = df, formula = y1 ~ x, missing = "x", model = "lm",
    m = m, seed = 1L
  )
}

pmb_check <- function(pooled, terms) {
  expect_s3_class(pooled, "data.frame")
  expect_true(all(terms %in% pooled$term))
  expect_true(all(is.finite(pooled$estimate)))
  expect_true(all(is.finite(pooled$std.error)))
  expect_true(all(pooled$std.error > 0))
  # Between-imputation variance must enter: pooling only one fit, or fits on
  # identical datasets, would give riv = 0.
  expect_true(all(pooled$riv[pooled$term %in% terms] > 0))
}

pmb_reload_check <- function(fits, backend) {
  fit_path <- tempfile(fileext = ".rds")
  saveRDS(fits, fit_path)
  check_path <- normalizePath("reload-fit.R")
  audit_lib <- Sys.getenv("PIGAUTO_AUDIT_INSTALLED_LIB", unset = "")
  installed <- nzchar(audit_lib)
  location <- if (installed) audit_lib else normalizePath(file.path("..", ".."))
  reload_output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(check_path), shQuote(fit_path),
      shQuote(backend), shQuote(location),
      shQuote(if (installed) "installed" else "source")),
    stdout = TRUE, stderr = TRUE
  ))
  reload_status <- attr(reload_output, "status")
  if (is.null(reload_status)) reload_status <- 0L
  expect_identical(reload_status, 0L,
                   info = paste(reload_output, collapse = "\n"))
  expect_true(any(grepl(paste("RELOAD_POOL_OK", backend), reload_output,
                        fixed = TRUE)),
              info = paste(reload_output, collapse = "\n"))
}

test_that("pool_mi() pools drmTMB fits from with_imputations()", {
  skip_on_cran()
  skip_if_not_installed("drmTMB")
  drm_fit <- getExportedValue("drmTMB", "drmTMB")
  bf <- getExportedValue("drmTMB", "bf")
  mi <- pmb_mi()
  fits <- with_imputations(mi, function(d) drm_fit(bf(y1 ~ x, sigma ~ x), data = d))
  pmb_check(pool_mi(fits),
            c("mu:(Intercept)", "mu:x", "sigma:(Intercept)", "sigma:x"))
  pmb_reload_check(fits, "drmTMB")
})

test_that("pool_mi() pools gllvmTMB fixed effects from with_imputations()", {
  skip_on_cran()
  skip_if_not_installed("gllvmTMB")
  mi <- pmb_mi()
  # Four genuine Gaussian traits are represented by four rows per unit.
  # The latent factor is a shared unit effect, matching the declared model
  # structure rather than manufacturing repeated jittered observations.
  fits <- with_imputations(mi, function(d) {
    tr <- paste0("y", seq_len(4L))
    long <- data.frame(
      unit = factor(rep(rownames(d), times = 4L)),
      trait = factor(rep(tr, each = nrow(d)), levels = tr),
      x = rep(d$x, times = 4L),
      value = unlist(d[tr], use.names = FALSE)
    )
    stopifnot(nrow(long) == 4L * nrow(d), all(table(long$unit) == 4L))
    suppressMessages(suppressWarnings(gllvmTMB::gllvmTMB(
      value ~ 0 + trait + trait:x + latent(0 + trait | unit, d = 1,
                                          unique = FALSE),
      data = long, trait = "trait", unit = "unit", family = gaussian(),
      REML = FALSE, engine = "tmb", silent = TRUE
    )))
  }, .on_error = "stop")
  pmb_reload_check(fits, "gllvmTMB")
  pooled <- pool_mi(fits)
  pmb_check(pooled, c(paste0("traity", 1:4), paste0("traity", 1:4, ":x")))
})
