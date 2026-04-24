#!/usr/bin/env Rscript
# script/bench_delhey_covariates.R
#
# Covariate-lift test on the bundled delhey5809 dataset (Delhey et al.
# 2019 Ecology Letters): 5,809 bird species × plumage lightness (male +
# female) × 5 pre-computed environmental covariates.
#
# Plumage lightness follows Gloger's rule: darker in warm/humid
# climates.  Classic textbook climate-driven trait.  Covariates already
# present -- no GBIF + WorldClim pipeline needed.  Clean, curated,
# species-level.
#
# This is the CLEAN environmental-covariate test that BIEN couldn't be.
# Three fits: baseline / cov + safety_floor = TRUE / cov + safety_floor
# = FALSE.
#
# Invocation:
#   PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_delhey_covariates.R

suppressPackageStartupMessages({
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH", unset = "")
  if (nzchar(pkg_path)) devtools::load_all(pkg_path, quiet = TRUE)
  else library(pigauto)
  library(ape)
})
options(warn = -1L)

SEED      <- 2026L
MISS_FRAC <- 0.30
N_IMP     <- 20L
EPOCHS    <- as.integer(Sys.getenv("PIGAUTO_BENCH_EPOCHS", "150"))
N_TARGET  <- as.integer(Sys.getenv("PIGAUTO_BENCH_N", "5809"))

cat_line <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

data(delhey5809)
data(tree_delhey)
cat_line(sprintf("delhey5809: %d rows x %d cols", nrow(delhey5809), ncol(delhey5809)))

# Traits: lightness_male + lightness_female
trait_cols <- c("lightness_male", "lightness_female")
cov_cols   <- c("annual_mean_temperature", "annual_precipitation",
                 "percent_tree_cover",
                 "mean_temperature_of_warmest_quarter",
                 "precipitation_of_warmest_quarter",
                 "midLatitude")

# Align rownames to species + tree tips
df_full <- delhey5809
rownames(df_full) <- df_full$Species_Key
matched <- intersect(df_full$Species_Key, tree_delhey$tip.label)
cat_line(sprintf("matched to tree: %d species", length(matched)))

set.seed(SEED)
sp <- sample(matched, min(N_TARGET, length(matched)))
df_aligned <- df_full[sp, , drop = FALSE]
tr <- ape::keep.tip(tree_delhey, sp)
cat_line(sprintf("using n = %d species", length(sp)))

# Build trait-only df
df_traits <- df_aligned[, trait_cols, drop = FALSE]
# Covariates (fully observed, no NA handling needed for these pre-computed
# environmental features)
cov_df <- df_aligned[, cov_cols, drop = FALSE]
# Drop rows where ANY covariate is NA (small handful if any)
ok_cov <- stats::complete.cases(cov_df)
if (any(!ok_cov)) {
  cat_line(sprintf("dropping %d species with any NA covariate", sum(!ok_cov)))
  df_traits <- df_traits[ok_cov, , drop = FALSE]
  cov_df    <- cov_df[ok_cov, , drop = FALSE]
  tr        <- ape::keep.tip(tr, rownames(df_traits))
}
cat_line(sprintf("final n = %d species, %d traits, %d covariates",
                  nrow(df_traits), ncol(df_traits), ncol(cov_df)))

# 30% MCAR mask on each trait
mask <- matrix(FALSE, nrow = nrow(df_traits), ncol = length(trait_cols),
               dimnames = list(NULL, trait_cols))
df_obs <- df_traits
for (v in trait_cols) {
  ok <- which(!is.na(df_traits[[v]]))
  if (length(ok) < 20L) next
  idx <- sample(ok, round(MISS_FRAC * length(ok)))
  mask[idx, v] <- TRUE
  df_obs[[v]][idx] <- NA_real_
}
for (v in trait_cols) {
  cat_line(sprintf("  %-20s: %d held-out cells", v, sum(mask[, v])))
}

# 3 fits
cat_line("=============== fit baseline (no covariates) ===============")
t0 <- proc.time()[["elapsed"]]
res_none <- pigauto::impute(df_obs, tr,
                               epochs = EPOCHS, n_imputations = N_IMP,
                               verbose = FALSE, seed = SEED,
                               safety_floor = TRUE)
w_none <- proc.time()[["elapsed"]] - t0
cat_line(sprintf("baseline done in %.0fs", w_none))

cat_line("=============== fit cov + safety_floor=TRUE ===============")
t0 <- proc.time()[["elapsed"]]
res_bio <- pigauto::impute(df_obs, tr, covariates = cov_df,
                              epochs = EPOCHS, n_imputations = N_IMP,
                              verbose = FALSE, seed = SEED,
                              safety_floor = TRUE)
w_bio <- proc.time()[["elapsed"]] - t0
cat_line(sprintf("cov sf=TRUE done in %.0fs", w_bio))

cat_line("=============== fit cov + safety_floor=FALSE ===============")
t0 <- proc.time()[["elapsed"]]
res_off <- pigauto::impute(df_obs, tr, covariates = cov_df,
                              epochs = EPOCHS, n_imputations = N_IMP,
                              verbose = FALSE, seed = SEED,
                              safety_floor = FALSE)
w_off <- proc.time()[["elapsed"]] - t0
cat_line(sprintf("cov sf=FALSE done in %.0fs", w_off))

# Report
cat("\n\n============ RESULTS (n =", nrow(df_traits), ") ============\n\n")
cat(sprintf("%-20s %5s %10s %10s %10s %10s %8s %8s %8s %8s %8s\n",
             "trait", "n", "mean", "none", "cov_on", "cov_off",
             "none_r", "on_r", "off_r", "rat_on", "rat_off"))

rows <- list()
for (v in trait_cols) {
  if (!any(mask[, v])) next
  truth <- df_traits[[v]][mask[, v]]
  ok <- is.finite(truth)
  if (sum(ok) < 5L) next
  mean_pred <- mean(df_obs[[v]], na.rm = TRUE)
  rmse_mean <- sqrt(mean((mean_pred - truth[ok])^2, na.rm = TRUE))
  rmse_none <- sqrt(mean((res_none$completed[[v]][mask[, v]][ok] - truth[ok])^2,
                          na.rm = TRUE))
  rmse_bio  <- sqrt(mean((res_bio$completed[[v]][mask[, v]][ok]  - truth[ok])^2,
                          na.rm = TRUE))
  rmse_off  <- sqrt(mean((res_off$completed[[v]][mask[, v]][ok]  - truth[ok])^2,
                          na.rm = TRUE))
  r_none <- tryCatch({
    p <- res_none$completed[[v]][mask[, v]][ok]
    if (!is.finite(sd(p, na.rm = TRUE)) || sd(p, na.rm = TRUE) < 1e-10) NA_real_
    else suppressWarnings(stats::cor(p, truth[ok], use = "complete.obs"))
  }, error = function(e) NA_real_)
  r_bio <- tryCatch({
    p <- res_bio$completed[[v]][mask[, v]][ok]
    if (!is.finite(sd(p, na.rm = TRUE)) || sd(p, na.rm = TRUE) < 1e-10) NA_real_
    else suppressWarnings(stats::cor(p, truth[ok], use = "complete.obs"))
  }, error = function(e) NA_real_)
  r_off <- tryCatch({
    p <- res_off$completed[[v]][mask[, v]][ok]
    if (!is.finite(sd(p, na.rm = TRUE)) || sd(p, na.rm = TRUE) < 1e-10) NA_real_
    else suppressWarnings(stats::cor(p, truth[ok], use = "complete.obs"))
  }, error = function(e) NA_real_)
  cat(sprintf("%-20s %5d %10.3g %10.3g %10.3g %10.3g %+8.3f %+8.3f %+8.3f %8.3f %8.3f\n",
               v, sum(ok), rmse_mean, rmse_none, rmse_bio, rmse_off,
               ifelse(is.na(r_none), NA, r_none),
               ifelse(is.na(r_bio), NA, r_bio),
               ifelse(is.na(r_off), NA, r_off),
               rmse_bio / rmse_none, rmse_off / rmse_none))
  rows[[length(rows) + 1L]] <- data.frame(
    trait = v, n_held = sum(ok),
    mean_RMSE = rmse_mean, none_RMSE = rmse_none,
    cov_on_RMSE = rmse_bio, cov_off_RMSE = rmse_off,
    none_r = r_none, cov_on_r = r_bio, cov_off_r = r_off,
    ratio_on = rmse_bio / rmse_none, ratio_off = rmse_off / rmse_none)
}
all_res <- do.call(rbind, rows)

saveRDS(list(results = all_res, n = nrow(df_traits),
              config = list(seed = SEED, miss_frac = MISS_FRAC,
                             epochs = EPOCHS, n_imp = N_IMP)),
         "script/bench_delhey_covariates.rds", compress = "xz")
cat_line("wrote script/bench_delhey_covariates.rds")
