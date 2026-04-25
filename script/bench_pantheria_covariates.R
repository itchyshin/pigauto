#!/usr/bin/env Rscript
# script/bench_pantheria_covariates.R
#
# Covariate-lift test on PanTHERIA mammal traits with the
# bundled climate columns serving as covariates.  Real species-
# level mammal phylogeny + real climate values from species ranges.
# Mirrors bench_delhey_covariates.R but for mammals.
#
# Targets (continuous, BM-eligible):
#   - 5-1_AdultBodyMass_g          (Bergmann; expect strong phylo signal)
#   - 9-1_GestationLen_d           (life-history)
#   - 17-1_MaxLongevity_m          (life-history; covaries with mass)
#   - 15-1_LitterSize              (life-history; covaries inverse with mass)
#   - 21-1_PopulationDensity_n/km2 (geographic + life-history)
#
# Covariates (climate from species ranges):
#   - 28-1_Precip_Mean_mm          (mean annual precipitation in range)
#   - 28-2_Temp_Mean_01degC        (mean annual temperature in range)
#   - 26-2_GR_MaxLat_dd            (max latitude in range)
#   - 26-3_GR_MinLat_dd            (min latitude in range)
#   - 30-2_PET_Mean_mm             (potential evapotranspiration mean)
#
# Three fits: baseline / cov + sf=TRUE / cov + sf=FALSE.
#
# Invocation:
#   PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_pantheria_covariates.R

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
N_TARGET  <- as.integer(Sys.getenv("PIGAUTO_BENCH_N", "1000"))

cat_line <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "),
                                ..., "\n", sep = "")

here     <- "/Users/z3437171/Dropbox/Github Local/pigauto"
pan_path <- file.path(here, "script", "data-cache", "pantheria.txt")
tree_path <- file.path(here, "script", "data-cache", "mammal_tree.tre")
stopifnot(file.exists(pan_path), file.exists(tree_path))

cat_line("loading PanTHERIA + mammal tree ...")
pan <- utils::read.delim(pan_path, sep = "\t", stringsAsFactors = FALSE,
                          na.strings = c("-999", "-999.00", "NA"))
cat_line(sprintf("PanTHERIA: %d rows x %d cols", nrow(pan), ncol(pan)))

mt <- ape::read.tree(tree_path)
cat_line(sprintf("mammal tree: %d tips, binary=%s",
                  length(mt$tip.label), ape::is.binary(mt)))

# Resolve polytomies + ensure positive branch lengths
set.seed(SEED)
mt <- ape::multi2di(mt, random = TRUE)
if (is.null(mt$edge.length) || any(mt$edge.length <= 0)) {
  if (is.null(mt$edge.length)) {
    mt <- ape::compute.brlen(mt, method = "Grafen")
  } else {
    mt$edge.length[mt$edge.length <= 0] <- 1e-6
  }
}
cat_line(sprintf("mammal tree resolved: %d tips, binary=%s",
                  length(mt$tip.label), ape::is.binary(mt)))

# Align species names: PanTHERIA uses "Genus Species" with space; tree uses underscore
pan$Species_Key <- gsub(" ", "_", trimws(pan$MSW93_Binomial))
matched <- intersect(pan$Species_Key, mt$tip.label)
cat_line(sprintf("PanTHERIA x tree intersect: %d species", length(matched)))

# Subsample
set.seed(SEED)
sp <- sample(matched, min(N_TARGET, length(matched)))
pan_aligned <- pan[match(sp, pan$Species_Key), , drop = FALSE]
rownames(pan_aligned) <- pan_aligned$Species_Key
tr <- ape::keep.tip(mt, sp)
cat_line(sprintf("subsample: n = %d species", length(sp)))

# Trait & covariate columns
trait_cols <- c("X5.1_AdultBodyMass_g",
                 "X9.1_GestationLen_d",
                 "X17.1_MaxLongevity_m",
                 "X15.1_LitterSize",
                 "X21.1_PopulationDensity_n.km2")
cov_cols <- c("X28.1_Precip_Mean_mm",
               "X28.2_Temp_Mean_01degC",
               "X26.2_GR_MaxLat_dd",
               "X26.3_GR_MinLat_dd",
               "X30.2_PET_Mean_mm")

# Validate columns are present (PanTHERIA naming uses dashes that R replaces)
trait_cols <- intersect(trait_cols, colnames(pan_aligned))
cov_cols   <- intersect(cov_cols,   colnames(pan_aligned))
stopifnot(length(trait_cols) >= 3L, length(cov_cols) >= 3L)
cat_line(sprintf("traits in: %s", paste(trait_cols, collapse = ", ")))
cat_line(sprintf("covs   in: %s", paste(cov_cols, collapse = ", ")))

df_traits <- pan_aligned[, trait_cols, drop = FALSE]
cov_df    <- pan_aligned[, cov_cols,   drop = FALSE]

# Coerce to numeric + log-transform very-skew traits (mass, longevity, density)
for (v in trait_cols) df_traits[[v]] <- suppressWarnings(as.numeric(df_traits[[v]]))
for (v in cov_cols)   cov_df[[v]]    <- suppressWarnings(as.numeric(cov_df[[v]]))

skew_cols <- c("X5.1_AdultBodyMass_g", "X9.1_GestationLen_d",
                "X17.1_MaxLongevity_m", "X21.1_PopulationDensity_n.km2")
for (v in intersect(skew_cols, trait_cols)) {
  ok <- !is.na(df_traits[[v]]) & df_traits[[v]] > 0
  df_traits[[v]] <- ifelse(ok, log(df_traits[[v]]), NA_real_)
}

# Drop species with all-NA covariates
cov_complete <- stats::complete.cases(cov_df)
cat_line(sprintf("species with all 5 covariates non-NA: %d / %d",
                  sum(cov_complete), nrow(cov_df)))
if (sum(cov_complete) < 100L)
  stop("Too few species with full covariates; widen pool or relax cov_complete.")

df_traits <- df_traits[cov_complete, , drop = FALSE]
cov_df    <- cov_df[cov_complete, , drop = FALSE]
tr        <- ape::keep.tip(tr, rownames(df_traits))
cat_line(sprintf("final n = %d species, %d traits, %d covariates",
                  nrow(df_traits), ncol(df_traits), ncol(cov_df)))

# 30% MCAR mask
mask <- matrix(FALSE, nrow = nrow(df_traits), ncol = length(trait_cols),
                dimnames = list(NULL, trait_cols))
df_obs <- df_traits
set.seed(SEED + 1L)
for (v in trait_cols) {
  ok <- which(!is.na(df_traits[[v]]))
  if (length(ok) < 20L) next
  idx <- sample(ok, round(MISS_FRAC * length(ok)))
  mask[idx, v] <- TRUE
  df_obs[[v]][idx] <- NA_real_
}
for (v in trait_cols) {
  cat_line(sprintf("  %-30s: %d held-out cells", v, sum(mask[, v])))
}

# Three fits
cat_line("=============== fit baseline (no covariates) ===============")
t0 <- proc.time()[["elapsed"]]
res_none <- tryCatch(
  pigauto::impute(df_obs, tr,
                    epochs = EPOCHS, n_imputations = N_IMP,
                    verbose = FALSE, seed = SEED,
                    safety_floor = TRUE),
  error = function(e) { cat_line("baseline ERROR: ", conditionMessage(e)); NULL })
w_none <- proc.time()[["elapsed"]] - t0
cat_line(sprintf("baseline done in %.0fs", w_none))

cat_line("=============== fit cov + safety_floor=TRUE ===============")
t0 <- proc.time()[["elapsed"]]
res_bio <- tryCatch(
  pigauto::impute(df_obs, tr, covariates = cov_df,
                    epochs = EPOCHS, n_imputations = N_IMP,
                    verbose = FALSE, seed = SEED,
                    safety_floor = TRUE),
  error = function(e) { cat_line("cov sf=TRUE ERROR: ", conditionMessage(e)); NULL })
w_bio <- proc.time()[["elapsed"]] - t0
cat_line(sprintf("cov sf=TRUE done in %.0fs", w_bio))

cat_line("=============== fit cov + safety_floor=FALSE ===============")
t0 <- proc.time()[["elapsed"]]
res_off <- tryCatch(
  pigauto::impute(df_obs, tr, covariates = cov_df,
                    epochs = EPOCHS, n_imputations = N_IMP,
                    verbose = FALSE, seed = SEED,
                    safety_floor = FALSE),
  error = function(e) { cat_line("cov sf=FALSE ERROR: ", conditionMessage(e)); NULL })
w_off <- proc.time()[["elapsed"]] - t0
cat_line(sprintf("cov sf=FALSE done in %.0fs", w_off))

# Report
cat("\n\n============ RESULTS (n =", nrow(df_traits), ") ============\n\n")
cat(sprintf("%-32s %5s %10s %10s %10s %10s %8s %8s %8s %8s %8s\n",
             "trait", "n", "mean", "none", "cov_on", "cov_off",
             "none_r", "on_r", "off_r", "rat_on", "rat_off"))

safe_get <- function(res, v, idx) {
  if (is.null(res) || is.null(res$completed)) return(rep(NA_real_, length(idx)))
  res$completed[[v]][idx]
}
safe_cor <- function(p, t) {
  ok <- is.finite(p) & is.finite(t)
  if (sum(ok) < 5L) return(NA_real_)
  if (sd(p[ok], na.rm = TRUE) < 1e-10) return(NA_real_)
  suppressWarnings(stats::cor(p[ok], t[ok]))
}

rows <- list()
for (v in trait_cols) {
  if (!any(mask[, v])) next
  truth <- df_traits[[v]][mask[, v]]
  ok <- is.finite(truth)
  if (sum(ok) < 5L) next
  m_pred <- mean(df_obs[[v]], na.rm = TRUE)
  rmse_mean <- sqrt(mean((m_pred - truth[ok])^2, na.rm = TRUE))
  rmse_none <- sqrt(mean((safe_get(res_none, v, which(mask[, v]))[ok] - truth[ok])^2,
                          na.rm = TRUE))
  rmse_bio  <- sqrt(mean((safe_get(res_bio,  v, which(mask[, v]))[ok] - truth[ok])^2,
                          na.rm = TRUE))
  rmse_off  <- sqrt(mean((safe_get(res_off,  v, which(mask[, v]))[ok] - truth[ok])^2,
                          na.rm = TRUE))
  r_none <- safe_cor(safe_get(res_none, v, which(mask[, v]))[ok], truth[ok])
  r_bio  <- safe_cor(safe_get(res_bio,  v, which(mask[, v]))[ok], truth[ok])
  r_off  <- safe_cor(safe_get(res_off,  v, which(mask[, v]))[ok], truth[ok])
  cat(sprintf("%-32s %5d %10.3g %10.3g %10.3g %10.3g %+8.3f %+8.3f %+8.3f %8.3f %8.3f\n",
               v, sum(ok), rmse_mean, rmse_none, rmse_bio, rmse_off,
               ifelse(is.na(r_none), NA, r_none),
               ifelse(is.na(r_bio),  NA, r_bio),
               ifelse(is.na(r_off),  NA, r_off),
               rmse_bio / rmse_none, rmse_off / rmse_none))
  rows[[length(rows) + 1L]] <- data.frame(
    trait = v, n_held = sum(ok),
    mean_RMSE = rmse_mean, none_RMSE = rmse_none,
    cov_on_RMSE = rmse_bio, cov_off_RMSE = rmse_off,
    none_r = r_none, cov_on_r = r_bio, cov_off_r = r_off,
    ratio_on  = rmse_bio / rmse_none,
    ratio_off = rmse_off / rmse_none)
}
all_res <- do.call(rbind, rows)

saveRDS(list(results = all_res, n = nrow(df_traits),
              config = list(seed = SEED, miss_frac = MISS_FRAC,
                              epochs = EPOCHS, n_imp = N_IMP,
                              n_target = N_TARGET)),
         "script/bench_pantheria_covariates.rds", compress = "xz")
cat_line("wrote script/bench_pantheria_covariates.rds")
