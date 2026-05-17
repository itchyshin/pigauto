#!/usr/bin/env Rscript
# script/exp6_synth_cat.R
#
# Fast synthetic K=4 + K=5 categorical recovery test for the Sigma M-step
# variants. Bypasses the full multi_impute / OVR pipeline and tests the
# joint-MVN solver directly: matrix-normal data is generated under known
# BM cov, cells are masked, fit_mvn_bm_inhouse() imputes, argmax of the
# imputed K liability columns decodes the predicted class. Runs in <60 s
# per (n, K) on n=100 vs the ~10+ minute full pipeline cost.
#
# What it does:
#   1. Build a 100-tip coalescent tree.
#   2. Simulate K liability columns under matrix-normal BM, decode the
#      observed K-class category via winner-take-all on the K liabilities.
#   3. Apply 30% MCAR to the K liability cells (whole-row mask: when a row
#      is masked, ALL K liability columns become NA so the OVR semantics
#      hold).
#   4. Fit pigauto's in-house joint MVN solver on the K-column liability
#      matrix (uses Henderson sparse R^{-1} + Kronecker M-step).
#   5. Decode the imputed L_hat via argmax for masked rows.
#   6. Report decoded classification accuracy.
#
# Reports both accuracy numbers (K=4 and K=5) plus delta vs random
# baseline (acc - 1/K). Variant winners are those that improve acc on K=5
# without degrading K=4.

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(ape)
})

set.seed(2026L)

simulate_one <- function(n_tips, K) {
  tree <- ape::rcoal(n_tips)
  R <- cov2cor(ape::vcv(tree))
  R <- R[tree$tip.label, tree$tip.label]
  # K liability columns from a single matrix-normal BM draw.
  Sigma_true <- matrix(rnorm(K^2), K)
  Sigma_true <- crossprod(Sigma_true) / K
  Sigma_true <- (Sigma_true + t(Sigma_true)) / 2 + diag(0.1, K)
  A <- chol(R); B <- chol(Sigma_true)
  Z <- matrix(rnorm(n_tips * K), n_tips, K)
  L <- t(A) %*% Z %*% B
  rownames(L) <- tree$tip.label
  class_obs <- apply(L, 1L, which.max)
  list(L = L, R = R, tree = tree, true_class = class_obs)
}

run_one_variant <- function(K, n_tips = 100L, miss_rate = 0.30, seed = 2026L) {
  set.seed(seed)
  sim <- simulate_one(n_tips = n_tips, K = K)
  L <- sim$L
  miss_idx <- sample(n_tips, floor(miss_rate * n_tips))
  L_masked <- L
  L_masked[miss_idx, ] <- NA  # whole-row mask matches OVR semantics
  t0 <- Sys.time()
  fit <- tryCatch(
    fit_mvn_bm_inhouse(L_masked, tree = sim$tree),
    error = function(e) NULL
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(fit)) return(list(acc = NA_real_, K = K, err = "fit failed",
                                 time = elapsed))
  pred_class <- apply(fit$anc_recon[miss_idx, , drop = FALSE], 1L, which.max)
  acc <- mean(pred_class == sim$true_class[miss_idx])
  list(acc = acc, K = K, time = elapsed, n_iter = fit$n_iter,
       converged = fit$converged)
}

cat("=== Synthetic K=4 + K=5 categorical recovery test (solver-direct) ===\n")
cat("Tree: n_tips=100, miss_rate=0.30, seed=2026\n\n")

r4 <- run_one_variant(K = 4L)
cat(sprintf("K=4: accuracy = %.4f (random baseline = %.4f)  [%.1f s, %d iter, conv=%s]\n",
            r4$acc, 1/4, r4$time, r4$n_iter %||% NA, r4$converged %||% NA))

r5 <- run_one_variant(K = 5L)
cat(sprintf("K=5: accuracy = %.4f (random baseline = %.4f)  [%.1f s, %d iter, conv=%s]\n",
            r5$acc, 1/5, r5$time, r5$n_iter %||% NA, r5$converged %||% NA))

cat(sprintf("\nDelta-vs-random: K=4: %+.4f  |  K=5: %+.4f\n",
            r4$acc - 1/4, r5$acc - 1/5))

dir.create("dev", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(K4 = r4, K5 = r5, ts = Sys.time()),
        file.path("dev", "exp6_synth_results.rds"))
