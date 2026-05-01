#!/usr/bin/env Rscript
#
# script/diag_avonet_mass.R
#
# Diagnose the seed-2030 Mass instability uncovered by the
# 2026-05-01 multi-seed N_IMP=20 AVONET re-verify.  Seed 2030
# produces RMSE ~19960 on Mass while seeds 2031/2032 produce
# 293 / 2715 (both beat the column-mean baseline of ~1820).
#
# This script:
#   1. Reproduces the seed-2030 AVONET n=1500 setup
#   2. Runs impute() with N_IMP=20
#   3. Per imputation m in 1..20, computes residual on masked
#      Mass cells: truth - imputed
#   4. Identifies the top-K worst residuals across all draws
#   5. Tests two hypotheses:
#         H1: a SINGLE rogue imputation draw drives the RMSE
#             (median pooling would fix it)
#         H2: ALL draws have a few extreme tail-back-transform
#             cells (would need a value clamp instead)
#
# Output:
#   script/diag_avonet_mass.rds  -- per-imputation residuals
#   useful/MEMO_2026-05-01_avonet_mass_diag.md  -- findings memo
#
# Wall: ~9 min on Apple MPS.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH", unset = "")
  if (nzchar(pkg_path) && dir.exists(pkg_path) &&
      file.exists(file.path(pkg_path, "DESCRIPTION"))) {
    devtools::load_all(pkg_path, quiet = TRUE)
  } else {
    library(pigauto)
  }
})

here       <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds    <- file.path(here, "script", "diag_avonet_mass.rds")
out_memo   <- file.path(here, "useful",
                        "MEMO_2026-05-01_avonet_mass_diag.md")

SEED      <- 2030L
MISS_FRAC <- 0.30
N_SUB     <- 1500L
N_IMP     <- 20L

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - t0),
      ..., "\n", sep = "")
  flush.console()
}

log_line("Loading avonet_full + tree_full ...")
e <- new.env(parent = emptyenv())
utils::data("avonet_full", package = "pigauto", envir = e)
utils::data("tree_full",   package = "pigauto", envir = e)
df   <- e$avonet_full
tree <- e$tree_full
rownames(df) <- df$Species_Key
df$Species_Key <- NULL

set.seed(SEED)
keep <- sample(rownames(df), N_SUB)
df   <- df[keep, , drop = FALSE]
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
df   <- df[tree$tip.label, , drop = FALSE]
stopifnot(all(rownames(df) == tree$tip.label))
log_line(sprintf("Subset: %d species x %d traits", nrow(df), ncol(df)))

# Mirror bench mask exactly
set.seed(SEED)
df_truth  <- df
mask_test <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df),
                    dimnames = list(rownames(df), names(df)))
for (v in names(df)) {
  obs_idx <- which(!is.na(df[[v]]))
  to_hide <- sample(obs_idx, ceiling(length(obs_idx) * MISS_FRAC))
  mask_test[to_hide, v] <- TRUE
}
df_miss <- df
for (v in names(df)) df_miss[[v]][mask_test[, v]] <- NA

log_line(sprintf("Mass test cells: %d", sum(mask_test[, "Mass"])))

log_line(sprintf("impute(seed=%d, N_IMP=%d) ...", SEED, N_IMP))
res <- pigauto::impute(df_miss, tree,
                         epochs = 200L,
                         n_imputations = N_IMP,
                         verbose = FALSE,
                         seed = SEED)

# ---- Per-imputation Mass residuals ----------------------------------------

mass_idx <- which(mask_test[, "Mass"])
mass_truth <- df_truth$Mass[mass_idx]
species_idx <- rownames(df_truth)[mass_idx]
log_line(sprintf("Mass: held-out cells = %d, truth range = [%.2f, %.2f], "
                 ,
                  length(mass_idx),
                  min(mass_truth, na.rm = TRUE),
                  max(mass_truth, na.rm = TRUE)))

# imputed_datasets is a list of M data.frames (each n_species x n_traits)
mi <- res$prediction$imputed_datasets
if (is.null(mi) || length(mi) == 0L) {
  stop("No imputed_datasets in result -- check N_IMP > 1?")
}
log_line(sprintf("Got %d imputation draws", length(mi)))

# Per-imputation: extract Mass at masked rows
per_imp <- vapply(mi, function(d) as.numeric(d$Mass[mass_idx]),
                  numeric(length(mass_idx)))
# per_imp is now (n_held x M); each column is one imputation draw

# Pooled point estimate (median per cell, what predict.pigauto_fit returns
# under pool_method = "median")
pooled_median <- apply(per_imp, 1, stats::median, na.rm = TRUE)
pooled_mean   <- rowMeans(per_imp, na.rm = TRUE)

rmse <- function(x) sqrt(mean((mass_truth - x)^2, na.rm = TRUE))
log_line(sprintf("Median-pooled Mass RMSE = %.2f", rmse(pooled_median)))
log_line(sprintf("Mean-pooled   Mass RMSE = %.2f", rmse(pooled_mean)))

per_imp_rmse <- apply(per_imp, 2, rmse)
log_line(sprintf("Per-imputation Mass RMSE: min=%.2f, q25=%.2f, median=%.2f, q75=%.2f, max=%.2f",
                  min(per_imp_rmse), stats::quantile(per_imp_rmse, 0.25),
                  stats::median(per_imp_rmse),
                  stats::quantile(per_imp_rmse, 0.75),
                  max(per_imp_rmse)))

# Top-5 single-cell residuals at the median pool
res_pool <- mass_truth - pooled_median
top5_pool_idx <- order(abs(res_pool), decreasing = TRUE)[1:5]

log_line("Top-5 worst residuals at median-pooled Mass:")
for (k in top5_pool_idx) {
  log_line(sprintf("  species=%s truth=%.2f pooled=%.2f resid=%.2f draws_range=[%.2f, %.2f]",
                    species_idx[k], mass_truth[k], pooled_median[k],
                    res_pool[k], min(per_imp[k, ]), max(per_imp[k, ])))
}

# Per-imputation top-1 cell
worst_per_imp <- vapply(seq_len(ncol(per_imp)), function(m) {
  r <- mass_truth - per_imp[, m]
  k <- which.max(abs(r))
  c(species = species_idx[k], truth = mass_truth[k],
    imp = per_imp[k, m], resid = r[k])
}, character(4L))

# H1 vs H2 verdict
n_extreme_imp <- sum(per_imp_rmse > 5 * stats::median(per_imp_rmse))

log_line("---- Verdict ----")
log_line(sprintf("Per-imputation RMSE max/median ratio = %.2f",
                  max(per_imp_rmse) / stats::median(per_imp_rmse)))
log_line(sprintf("Imputations with RMSE > 5x median: %d / %d",
                  n_extreme_imp, length(per_imp_rmse)))
if (n_extreme_imp <= 2 && rmse(pooled_median) < rmse(pooled_mean) * 0.5) {
  verdict <- "H1: single rogue draw -- median pool fixes it"
} else if (rmse(pooled_median) > 2 * stats::median(per_imp_rmse)) {
  verdict <- "H2: all draws have tail outliers -- needs value clamp"
} else {
  verdict <- "MIXED: see per-imputation table"
}
log_line(verdict)

# ---- Persist ---------------------------------------------------------------

saveRDS(list(
  per_imp        = per_imp,
  pooled_median  = pooled_median,
  pooled_mean    = pooled_mean,
  per_imp_rmse   = per_imp_rmse,
  mass_truth     = mass_truth,
  species_idx    = species_idx,
  top5_pool_idx  = top5_pool_idx,
  worst_per_imp  = worst_per_imp,
  verdict        = verdict
), out_rds)
log_line("DONE -- ", out_rds)
