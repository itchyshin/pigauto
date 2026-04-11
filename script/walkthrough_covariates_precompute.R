#!/usr/bin/env Rscript
#
# script/walkthrough_covariates_precompute.R
#
# Precompute imputation results for the Delhey covariates walkthrough.
# Runs impute() with and without environmental covariates on the full
# 5,809-species Delhey plumage dataset. Downstream regression is NOT
# included -- at 5,809 species a phylogenetic mixed model takes hours
# per fit. The walkthrough (Option B) shows the MI workflow as code
# but does not execute it.
#
# Output: script/walkthrough_covariates.rds
# Runtime: ~15 min on a laptop.

suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto",
    quiet = TRUE
  )
})

log_line <- function(...) {
  cat(sprintf("[%6.1fs] %s\n", proc.time()["elapsed"], paste0(...)))
  flush.console()
}

log_line("Loading data...")
data(delhey5809, tree_delhey, package = "pigauto")
df <- delhey5809
n  <- nrow(df)
log_line(sprintf("Species: %d, traits: 2, covariates: 4", n))

# ---- Introduce 20% MCAR missingness ----------------------------------------
set.seed(2026)
traits_full <- df[, c("lightness_male", "lightness_female")]
rownames(traits_full) <- df$Species_Key   # Species_Key matches tree tip labels
miss_idx  <- sample(seq_len(n * 2L), round(0.20 * n * 2L))
traits_mat <- as.matrix(traits_full)
traits_mat[miss_idx] <- NA
traits_miss <- as.data.frame(traits_mat)
log_line(sprintf("Missing: %d / %d cells (%.1f%%)",
                 sum(is.na(traits_miss)), n * 2L,
                 100 * mean(is.na(traits_miss))))

covs <- df[, c("annual_mean_temperature", "annual_precipitation",
               "percent_tree_cover", "midLatitude")]
rownames(covs) <- df$Species_Key

# ---- Impute without covariates ---------------------------------------------
log_line("Imputing without covariates...")
t0 <- proc.time()["elapsed"]
res_no_cov <- tryCatch(
  impute(traits_miss, tree_delhey,
         covariates = NULL,
         epochs = 500L, verbose = FALSE, seed = 1L,
         eval_every = 50L, patience = 10L),
  error = function(e) { log_line("ERROR: ", conditionMessage(e)); NULL }
)
time_no_cov <- round(proc.time()["elapsed"] - t0, 1)
log_line(sprintf("  Done in %.1f min", time_no_cov / 60))

# ---- Impute with covariates ------------------------------------------------
log_line("Imputing with covariates...")
t0 <- proc.time()["elapsed"]
res_cov <- tryCatch(
  impute(traits_miss, tree_delhey,
         covariates = covs,
         epochs = 500L, verbose = FALSE, seed = 1L,
         eval_every = 50L, patience = 10L),
  error = function(e) { log_line("ERROR: ", conditionMessage(e)); NULL }
)
time_cov <- round(proc.time()["elapsed"] - t0, 1)
log_line(sprintf("  Done in %.1f min", time_cov / 60))

# ---- RMSE comparison -------------------------------------------------------
mask_mat <- matrix(FALSE, n, 2L,
                   dimnames = list(rownames(traits_full), colnames(traits_full)))
mask_mat[miss_idx] <- TRUE

rmse_row <- function(res, label, j) {
  if (is.null(res)) return(data.frame(trait = NA, condition = label,
                                      rmse = NA, pearson_r = NA))
  preds <- res$prediction$imputed[, colnames(traits_full)[j]]
  truth <- traits_full[, j]
  ok    <- mask_mat[, j] & is.finite(preds) & is.finite(truth)
  data.frame(
    trait     = colnames(traits_full)[j],
    condition = label,
    rmse      = round(sqrt(mean((truth[ok] - preds[ok])^2)), 3),
    pearson_r = round(cor(truth[ok], preds[ok]), 3),
    stringsAsFactors = FALSE
  )
}

rmse_table <- rbind(
  rmse_row(res_no_cov, "no_covariates",   1L),
  rmse_row(res_cov,    "with_covariates", 1L),
  rmse_row(res_no_cov, "no_covariates",   2L),
  rmse_row(res_cov,    "with_covariates", 2L)
)
log_line("RMSE table:")
print(rmse_table)

# ---- Save ------------------------------------------------------------------
out <- list(
  rmse_table    = rmse_table,
  time_no_cov   = time_no_cov,
  time_cov      = time_cov,
  n_species     = n,
  n_missing     = sum(is.na(traits_miss)),
  n_total_cells = n * 2L,
  timestamp     = format(Sys.time(), "%Y-%m-%d %H:%M")
)
saveRDS(out, "script/walkthrough_covariates.rds")
log_line("Saved script/walkthrough_covariates.rds")
log_line("Done.")
