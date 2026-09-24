# script/tests-rubin/test-freq-condmean.R
# Gate G-S2a for cond_draw()'s conditional-mean machinery (fit_block() + joint_cov()). On a
# make_cell("types_mixed", n = 60, seed = 1, lambda = 0.7, rho = 0.5) cell, the conditional mean of
# the missing block cells implied by joint_cov()'s kron(Sigma_p, C_lambda) covariance should
# reproduce Rphylopars' own fit$anc_recon for those same cells (on the fitted/transformed scale) --
# proof that the covariance builder (trait-major vec ordering, C_lambda = lambda*C +
# (1-lambda)*diag(diag(C)) on a height-1 ultrametric tree, no separate phenocov term for this
# single-obs-per-species data) matches what Rphylopars itself uses internally.

source(normalizePath(file.path(testthat::test_path(), "..", "rubin_freq.R")))
source(normalizePath(file.path(testthat::test_path(), "..", "campaign_gnn_off_lib.R")))

testthat::test_that("cond_draw's conditional mean reproduces Rphylopars anc_recon", {
  cell <- make_cell("types_mixed", n = 60, seed = 1, lambda = 0.7, rho = 0.5)
  block_traits <- default_block_traits(cell)
  testthat::expect_setequal(block_traits, c("c1", "c2", "prp"))

  pars <- fit_block(cell$df_miss, cell$tree, block_traits, cell$trait_types, phylo_model = "lambda")
  Yt <- transform_block(cell$df_miss, block_traits, pars$is_prp)
  cd <- cond_draw(Yt, pars, cell$tree, M = 1L)

  testthat::expect_gt(length(cd$mis_idx), 0L)

  n <- length(pars$species); p <- length(block_traits)
  mis_trait_idx <- ((cd$mis_idx - 1L) %/% n) + 1L
  mis_sp_idx <- ((cd$mis_idx - 1L) %% n) + 1L
  sp <- pars$species

  rec <- pars$fit$anc_recon[sp, block_traits, drop = FALSE]
  anc_at_missing <- vapply(seq_along(cd$mis_idx), function(i) {
    rec[mis_sp_idx[i], mis_trait_idx[i]]
  }, numeric(1))

  diffs <- cd$cond_mean - anc_at_missing
  max_diff <- max(abs(diffs))
  cat(sprintf("G-S2a max abs diff (cond_mean vs anc_recon): %.3e over %d missing cells\n",
              max_diff, length(cd$mis_idx)))
  testthat::expect_lt(max_diff, 1e-6)
})
