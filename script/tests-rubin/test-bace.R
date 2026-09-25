# script/tests-rubin/test-bace.R
# Gate G-S3 for the two BACE arms in script/rubin_bace.R: mi_bace_shipped() (BACE's own n_final
# datasets, untouched) and mi_bace_resid() (BACE's shipped datasets plus an independent per-cell
# residual draw on gaussian-modelled missing cells, from that final fit's own VCV[,"units"]
# posterior). See script/rubin_bace.R's header for why "BACE + residual draw" is a materially
# different claim against the INSTALLED BACE package than against the stale on-disk BACE/R
# source -- this test exercises the installed package (the one BACE::bace() actually is here),
# not the source comments.
#
# ONE tiny shared BACE fit for all three expectations (n = 40, seed = 1, lambda = 0.7, rho = 0.5,
# nitt = 3000, burnin = 1000, thin = 10, runs = 2, M = 5) -- BACE itself dominates runtime, so
# everything downstream reuses the single fit rather than re-fitting per expectation.

source(normalizePath(file.path(testthat::test_path(), "..", "campaign_gnn_off_lib.R")))
source(normalizePath(file.path(testthat::test_path(), "..", "rubin_bace.R")))

testthat::skip_if_not_installed("BACE")

suppressPackageStartupMessages({ library(ape) })

cell <- make_cell("types_mixed", 40, seed = 1, lambda = 0.7, rho = 0.5)
df_miss <- cell$df_miss

t_fit0 <- Sys.time()
fit <- suppressWarnings(fit_bace_mi(cell, M = 5L, nitt = 3000, burnin = 1000, thin = 10, runs = 2))
cat(sprintf("[test-bace] BACE tiny fit (n=40, nitt=3000, burnin=1000, thin=10, runs=2, M=5): %.1fs "
           , as.numeric(difftime(Sys.time(), t_fit0, units = "secs"))))
cat(sprintf("(fit$wall_s = %.1fs, converged = %s, ess_med = %.0f)\n",
           fit$wall_s, fit$converged, fit$ess_med %||% NA_real_))

shipped <- mi_bace_shipped(fit$outb, df_miss)
M <- length(fit$outb$imputed_datasets)

testthat::test_that("fit_bace_mi returns the bace object, wall time, and BACE's own diagnostics", {
  testthat::expect_s3_class(fit$outb, "bace_complete")
  testthat::expect_length(fit$outb$imputed_datasets, 5L)
  testthat::expect_true(is.numeric(fit$wall_s) && fit$wall_s > 0)
  testthat::expect_identical(fit$converged, isTRUE(fit$outb$converged))
  bd <- bace_diagnostics(fit$outb)
  testthat::expect_identical(fit$diag[names(bd)], bd)   # plus the input-cleaning record
  testthat::expect_true(is.character(fit$diag$input_fix))
})

testthat::test_that("mi_bace_shipped is outb$imputed_datasets, aligned to df_miss row order/columns only", {
  testthat::expect_type(shipped, "list")
  testthat::expect_length(shipped, M)
  for (i in seq_len(M)) {
    raw <- fit$outb$imputed_datasets[[i]]
    raw <- raw[rownames(df_miss), names(df_miss), drop = FALSE]
    s <- shipped[[i]]
    testthat::expect_identical(rownames(s), rownames(df_miss))
    testthat::expect_identical(names(s), names(df_miss))
    for (v in names(df_miss)) {
      # compare on values only -- alignment is allowed to change storage class (e.g. BACE's
      # rounded-double count predictions vs df_miss's integer class), never the value
      testthat::expect_equal(as.character(s[[v]]), as.character(raw[[v]]),
                             info = sprintf("dataset %d, trait %s", i, v))
    }
  }
})

testthat::test_that("mi_bace_resid leaves observed cells untouched and changes missing gaussian cells", {
  resid <- mi_bace_resid(fit$outb, df_miss, seed = 20260924L)
  testthat::expect_length(resid, M)

  types <- fit$outb$final_results$types
  gaussian_traits <- intersect(names(df_miss),
                               names(types)[vapply(types, identical, logical(1), y = "gaussian")])
  testthat::expect_true(all(c("c1", "c2", "prp") %in% gaussian_traits))
  testthat::expect_false("cnt" %in% gaussian_traits)   # count trait, poisson-classified by BACE

  for (i in seq_len(M)) {
    for (v in names(df_miss)) {
      obs_idx <- which(!is.na(df_miss[[v]]))
      # observed cells: bit-identical to df_miss, for every trait, in every dataset
      testthat::expect_equal(as.character(resid[[i]][[v]][obs_idx]), as.character(df_miss[[v]][obs_idx]),
                             info = sprintf("dataset %d, trait %s (observed cells)", i, v))
      miss_idx <- which(is.na(df_miss[[v]]))
      if (!length(miss_idx)) next
      if (v %in% gaussian_traits) {
        # missing gaussian cells: resid != shipped (an independent continuous draw was added)
        testthat::expect_false(isTRUE(all.equal(resid[[i]][[v]][miss_idx], shipped[[i]][[v]][miss_idx])),
                               info = sprintf("dataset %d, trait %s (missing, gaussian)", i, v))
      } else {
        # missing non-gaussian cells: untouched, identical to the shipped BACE prediction
        testthat::expect_equal(as.character(resid[[i]][[v]][miss_idx]), as.character(shipped[[i]][[v]][miss_idx]),
                               info = sprintf("dataset %d, trait %s (missing, non-gaussian)", i, v))
      }
    }
  }
})

testthat::test_that("residual draw variance matches mean(VCV[,'units']) * sd_val^2 within 5%", {
  # One dataset/trait (i = 1, v = "c1"). Replicate the draw at ONE missing cell 4000 times, with the
  # sampled VCV iteration drawn fresh each replicate (a fresh seed per call to mi_bace_resid()), and
  # compare the empirical variance of the added noise to the theoretical target computed from the
  # SAME sd_val reconstruction mi_bace_resid() itself uses (.bace_gaussian_sd(), exported by
  # rubin_bace.R) -- this checks the sampling mechanics (right VCV column, sd not variance, no
  # off-by-scale error), not BACE's own internals.
  i <- 1L; v <- "c1"
  miss_idx <- which(is.na(df_miss[[v]]))
  testthat::expect_true(length(miss_idx) > 0)
  row1 <- miss_idx[1]

  n_rep <- 4000L
  noise <- numeric(n_rep)
  base <- shipped[[i]][[v]][row1]
  for (r in seq_len(n_rep)) {
    draw <- mi_bace_resid(fit$outb, df_miss, seed = r)
    noise[r] <- draw[[i]][[v]][row1] - base
  }

  sd_val <- .bace_gaussian_sd(df_miss, v)
  units <- as.matrix(fit$outb$final_results$all_models[[i]][[v]]$VCV)[, "units"]
  theoretical_var <- mean(units) * sd_val^2
  empirical_var <- stats::var(noise)

  testthat::expect_equal(empirical_var, theoretical_var, tolerance = 0.05)
})
