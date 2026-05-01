#!/usr/bin/env Rscript
# script/spot_check_predictions.R
#
# Print pigauto predictions vs truth for a BACE-sim multi_obs cell so we can
# see what the bias looks like.

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  pkg_path <- "/Users/z3437171/Dropbox/Github Local/pigauto"
  devtools::load_all(pkg_path, quiet = TRUE)
  source(file.path(pkg_path, "script", "sim_bace_dgp.R"))
})

t0 <- proc.time()[["elapsed"]]

seed <- 9001L
d <- sim_bace_dgp(n_species = 100L, multi_obs_ratio = 4L,
                   phylo_signal = 0.4, beta_resp_strength = 0,
                   response_type = "gaussian", n_predictors = 3L,
                   miss_frac = 0.30, seed = seed)

cat("=== Data summary ===\n")
cat(sprintf("n_cases    : %d\n", nrow(d$df_observed)))
cat(sprintf("mask       : %d cells held out\n", sum(d$mask)))
cat(sprintf("y range    : [%.2f, %.2f]\n",
            min(d$df_complete$y), max(d$df_complete$y)))
cat(sprintf("y mean     : %.4f\n", mean(d$df_complete$y)))
cat(sprintf("y sd       : %.4f\n", sd(d$df_complete$y)))
cat(sprintf("y_obs mean : %.4f (after mask)\n",
            mean(d$df_observed$y, na.rm = TRUE)))
cat(sprintf("y_obs sd   : %.4f\n", sd(d$df_observed$y, na.rm = TRUE)))
cat("\n")

# Per-species summary
sp_obs_count <- table(d$df_observed$species)
cat(sprintf("Species: n=%d, obs/sp range [%d, %d], mean %.1f\n",
            length(sp_obs_count), min(sp_obs_count), max(sp_obs_count),
            mean(sp_obs_count)))

# How many held-out cells have species fully missing? Should be ~zero in BACE.
sp_obs_after_mask <- tapply(d$df_observed$y, d$df_observed$species,
                             function(v) sum(!is.na(v)))
sp_fully_missing <- names(sp_obs_after_mask)[sp_obs_after_mask == 0]
cat(sprintf("Species with all observations masked: %d\n",
            length(sp_fully_missing)))
cat("\n")

# Run pigauto
res <- impute(traits = d$df_observed[, c("species", "y"), drop = FALSE],
               tree = d$tree, species_col = "species",
               epochs = 200L, verbose = FALSE, seed = seed,
               missing_frac = 0.0)

cat(sprintf("Time to fit: %.1fs\n\n", proc.time()[["elapsed"]] - t0))

held_out_idx <- which(d$mask)
truth <- d$df_complete$y[held_out_idx]
pred  <- res$completed$y[held_out_idx]

cat("=== Prediction summary ===\n")
cat(sprintf("Predictions range: [%.2f, %.2f]\n", min(pred), max(pred)))
cat(sprintf("Predictions mean : %.4f (truth mean: %.4f, diff: %.4f)\n",
            mean(pred), mean(truth), mean(pred) - mean(truth)))
cat(sprintf("Predictions sd   : %.4f (truth sd:  %.4f)\n", sd(pred), sd(truth)))
cat(sprintf("Pearson r        : %.4f\n", cor(pred, truth)))
cat(sprintf("RMSE             : %.4f\n", sqrt(mean((pred - truth)^2))))
cat("\n")

# First 20 held-out cells
cat("=== First 20 held-out cells (truth, pred, residual) ===\n")
for (i in seq_len(min(20L, length(held_out_idx)))) {
  sp <- as.character(d$df_observed$species[held_out_idx[i]])
  obs_for_sp <- d$df_observed$y[d$df_observed$species == sp]
  obs_count <- sum(!is.na(obs_for_sp))
  obs_mean <- mean(obs_for_sp, na.rm = TRUE)
  cat(sprintf("  sp=%s  n_obs=%d  sp_mean=%.3f  truth=%.3f  pred=%.3f  resid=%.3f\n",
              sp, obs_count, obs_mean, truth[i], pred[i], pred[i] - truth[i]))
}
cat("\n")

# Compare: what would species_mean give?
sp_means_obs <- tapply(d$df_observed$y, d$df_observed$species,
                        function(v) mean(v, na.rm = TRUE))
sp_pred_baseline <- sp_means_obs[as.character(d$df_observed$species[held_out_idx])]
sp_pred_baseline[is.na(sp_pred_baseline)] <- mean(d$df_observed$y, na.rm = TRUE)

cat("=== species_mean RMSE comparison ===\n")
cat(sprintf("species_mean RMSE: %.4f\n",
            sqrt(mean((sp_pred_baseline - truth)^2))))

# Diagnostic: is pigauto's pred close to species_mean? Or way off?
cat("\n")
cat(sprintf("Correlation between pigauto pred and species_mean: %.4f\n",
            cor(pred, as.numeric(sp_pred_baseline))))
cat(sprintf("Mean abs diff (pigauto - species_mean): %.4f\n",
            mean(abs(pred - sp_pred_baseline))))
cat(sprintf("RMSE pigauto - species_mean: %.4f\n",
            sqrt(mean((pred - sp_pred_baseline)^2))))
