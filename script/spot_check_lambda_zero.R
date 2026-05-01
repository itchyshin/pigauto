#!/usr/bin/env Rscript
# script/spot_check_lambda_zero.R
#
# Quick spot-check: monkey-patch phylo_cor_matrix to return identity (= Pagel's
# lambda = 0 extreme).  If pigauto's multi-obs RMSE on a BACE-sim cell drops to
# ~column-mean level, the lambda hypothesis is confirmed.  If not, there's
# another bug to find.
#
# Single cell, n_species=100, multi_obs=4, beta=0, ps=0.4, ~3 min wall.

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  pkg_path <- "/Users/z3437171/Dropbox/Github Local/pigauto"
  devtools::load_all(pkg_path, quiet = TRUE)
  source(file.path(pkg_path, "script", "sim_bace_dgp.R"))
})

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - t0),
                              ..., "\n", sep = "")

# Generate one cell --------------------------------------------------------
seed <- 9001L
log_line("Generating BACE-sim multi_obs cell (n=100, mo=4, ps=0.4, beta=0)...")
d <- sim_bace_dgp(n_species = 100L, multi_obs_ratio = 4L,
                   phylo_signal = 0.4, beta_resp_strength = 0,
                   response_type = "gaussian", n_predictors = 3L,
                   miss_frac = 0.30, seed = seed)
log_line(sprintf("  n_cases=%d, mask=%d", nrow(d$df_observed), sum(d$mask)))

truth <- d$df_complete$y[d$mask]

# Method 1: column_mean baseline -------------------------------------------
cm <- mean(d$df_observed$y, na.rm = TRUE)
rmse_cm <- sqrt(mean((truth - cm)^2))
log_line(sprintf("  column_mean RMSE = %.4f", rmse_cm))

# Method 2: pigauto NO-COV with default lambda=1 (current behaviour) -------
log_line("Fitting pigauto no_cov, default (lambda=1)...")
res1 <- impute(traits = d$df_observed[, c("species", "y"), drop = FALSE],
                tree = d$tree, species_col = "species",
                epochs = 200L, verbose = FALSE, seed = seed,
                missing_frac = 0.0)
pred1 <- res1$completed$y[d$mask]
rmse_pig_lam1 <- sqrt(mean((truth - pred1)^2))
log_line(sprintf("  pigauto (lambda=1) RMSE = %.4f", rmse_pig_lam1))

# Method 3: pigauto NO-COV with phylo_cor_matrix monkey-patched to identity (lambda=0)
log_line("Monkey-patching phylo_cor_matrix to return identity (lambda=0 extreme)...")
orig_phylo_cor_matrix <- pigauto:::phylo_cor_matrix
identity_phylo_cor <- function(tree) {
  n <- ape::Ntip(tree)
  R <- diag(1, n)
  rownames(R) <- colnames(R) <- tree$tip.label
  R
}
# Patch in the package namespace
unlockBinding("phylo_cor_matrix", asNamespace("pigauto"))
assign("phylo_cor_matrix", identity_phylo_cor, envir = asNamespace("pigauto"))
lockBinding("phylo_cor_matrix", asNamespace("pigauto"))

log_line("Fitting pigauto no_cov with lambda=0 (no phylo borrowing)...")
res2 <- impute(traits = d$df_observed[, c("species", "y"), drop = FALSE],
                tree = d$tree, species_col = "species",
                epochs = 200L, verbose = FALSE, seed = seed,
                missing_frac = 0.0)
pred2 <- res2$completed$y[d$mask]
rmse_pig_lam0 <- sqrt(mean((truth - pred2)^2))
log_line(sprintf("  pigauto (lambda=0) RMSE = %.4f", rmse_pig_lam0))

# Restore
unlockBinding("phylo_cor_matrix", asNamespace("pigauto"))
assign("phylo_cor_matrix", orig_phylo_cor_matrix, envir = asNamespace("pigauto"))
lockBinding("phylo_cor_matrix", asNamespace("pigauto"))

# Method 4: per-species mean baseline (the optimal naive predictor) --------
sp_means <- tapply(d$df_observed$y, d$df_observed$species,
                    function(v) { v <- v[!is.na(v)]; if (length(v)) mean(v) else NA })
sp_pred <- sp_means[as.character(d$df_observed$species[d$mask])]
sp_pred[is.na(sp_pred)] <- cm
rmse_sm <- sqrt(mean((truth - sp_pred)^2))
log_line(sprintf("  species_mean RMSE = %.4f", rmse_sm))

cat("\n=== Verdict ===\n")
cat(sprintf("column_mean       : %.4f\n", rmse_cm))
cat(sprintf("species_mean      : %.4f\n", rmse_sm))
cat(sprintf("pigauto (lambda=1): %.4f  (ratio to column_mean: %.3f)\n",
            rmse_pig_lam1, rmse_pig_lam1 / rmse_cm))
cat(sprintf("pigauto (lambda=0): %.4f  (ratio to column_mean: %.3f)\n",
            rmse_pig_lam0, rmse_pig_lam0 / rmse_cm))
cat("\nIf lambda=0 pigauto >> lambda=1 pigauto AND drops to species_mean level,\n")
cat("the lambda hypothesis is confirmed.  Otherwise there's another bug.\n")
log_line("done")
