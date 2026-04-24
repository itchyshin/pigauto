#!/usr/bin/env Rscript
# script/bench_plants_bioclim_n240.R
#
# Mid-size plants bench using the 300-species per-occurrence bioclim
# fixture (240 of which have valid bioclim). Compares:
#   1. baseline         : no covariates, safety_floor = TRUE
#   2. per-occurrence   : bioclim covariates, safety_floor = TRUE
#
# Uses ONLY already-cached species -- zero new GBIF fetches.
#
# If per-occurrence lifts sla or leaf_area by >= 10% RMSE here, the
# n=4,745 full bench is likely to deliver the paper claim.
# If it's flat here, the paper claim is in trouble and we need
# to debug before spending 2+ hr on full-scale.
#
# Invocation:
#   PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_plants_bioclim_n240.R

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
EPOCHS    <- 150L

# -------------------------------------------------------------------------
# Load fixture + BIEN cache
# -------------------------------------------------------------------------

cat_line <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

wc_fx_path <- "tests/testthat/fixtures/worldclim_plants_300.rds"
gbif_fx_path <- "tests/testthat/fixtures/gbif_plants_300.rds"
stopifnot(file.exists(wc_fx_path), file.exists(gbif_fx_path))
wc_all <- readRDS(wc_fx_path)
cat_line(sprintf("loaded worldclim fixture: %d species x %d cols",
                  nrow(wc_all), ncol(wc_all)))

bien_cache_t <- "script/data-cache/bien_trait_means.rds"
bien_cache_tr <- "script/data-cache/bien_tree.rds"
stopifnot(file.exists(bien_cache_t), file.exists(bien_cache_tr))
trait_means <- readRDS(bien_cache_t)
tree_raw    <- readRDS(bien_cache_tr)
tree_all <- if (is.list(tree_raw) && !inherits(tree_raw, "phylo")) {
  t <- tree_raw$scenario.3; t$tip.label <- gsub("_", " ", t$tip.label); t
} else tree_raw

all_species <- Reduce(union, lapply(trait_means,
  function(d) if (!is.null(d)) d$species else character(0)))
wide <- data.frame(species = all_species, stringsAsFactors = FALSE)
for (nm in names(trait_means)) {
  d <- trait_means[[nm]]
  if (is.null(d)) { wide[[nm]] <- NA_real_; next }
  m <- match(wide$species, d$species)
  wide[[nm]] <- suppressWarnings(as.numeric(d$mean_value[m]))
}

# -------------------------------------------------------------------------
# Intersect species: fixture with valid bioclim + BIEN traits + tree tips
# -------------------------------------------------------------------------

iqr_cols <- grep("_iqr$", colnames(wc_all), value = TRUE)
sp_with_iqr <- rownames(wc_all)[rowSums(wc_all[, iqr_cols] > 0, na.rm = TRUE) > 0]
matched <- Reduce(intersect, list(wide$species, tree_all$tip.label, sp_with_iqr))
cat_line(sprintf("matched species (bioclim + BIEN + tree): %d", length(matched)))
stopifnot(length(matched) >= 100L)

set.seed(SEED)
sp_s <- matched  # use all of them
cat_line(sprintf("using n = %d species", length(sp_s)))

wide_s <- wide[wide$species %in% sp_s, , drop = FALSE]
rownames(wide_s) <- wide_s$species; wide_s$species <- NULL
cov_bio <- wc_all[sp_s, grepl("^bio", colnames(wc_all)), drop = FALSE]
tree_s <- ape::keep.tip(tree_all, sp_s)

# -------------------------------------------------------------------------
# 30% MCAR mask
# -------------------------------------------------------------------------

df <- wide_s
cont_cols <- colnames(df)
mask <- matrix(FALSE, nrow = nrow(df), ncol = length(cont_cols),
               dimnames = list(NULL, cont_cols))
for (v in cont_cols) {
  ok <- which(!is.na(wide_s[[v]]))
  if (length(ok) < 20L) next
  idx <- sample(ok, round(MISS_FRAC * length(ok)))
  mask[idx, v] <- TRUE
  df[[v]][idx] <- NA_real_
}
for (v in cont_cols) {
  n_held <- sum(mask[, v])
  cat_line(sprintf("trait %-15s: %d held-out cells", v, n_held))
}

# -------------------------------------------------------------------------
# Fit baseline (no covariates)
# -------------------------------------------------------------------------

cat_line("=============== baseline (no covariates) ===============")
t0 <- proc.time()[["elapsed"]]
res_base <- pigauto::impute(df, tree_s,
                               epochs = EPOCHS, n_imputations = N_IMP,
                               verbose = FALSE, seed = SEED,
                               safety_floor = TRUE)
w_base <- proc.time()[["elapsed"]] - t0
cat_line(sprintf("baseline done in %.1fs", w_base))

cat_line("=============== per-occurrence bioclim covariates ===============")
t0 <- proc.time()[["elapsed"]]
res_bio <- pigauto::impute(df, tree_s, covariates = cov_bio,
                              epochs = EPOCHS, n_imputations = N_IMP,
                              verbose = FALSE, seed = SEED,
                              safety_floor = TRUE)
w_bio <- proc.time()[["elapsed"]] - t0
cat_line(sprintf("per-occur done in %.1fs", w_bio))

# -------------------------------------------------------------------------
# Report
# -------------------------------------------------------------------------

cat("\n\n=============== RESULTS (n =", length(sp_s), ") ===============\n\n")
rows <- list()
for (v in cont_cols) {
  if (!any(mask[, v])) next
  truth <- wide_s[[v]][mask[, v]]
  ok <- is.finite(truth)
  if (sum(ok) < 10L) next
  mean_pred <- mean(df[[v]], na.rm = TRUE)
  rmse_mean <- sqrt(mean((mean_pred - truth[ok])^2))
  rmse_base <- sqrt(mean((res_base$completed[[v]][mask[, v]][ok] - truth[ok])^2))
  rmse_bio  <- sqrt(mean((res_bio$completed[[v]][mask[, v]][ok]  - truth[ok])^2))
  r_base <- if (sd(res_base$completed[[v]][mask[, v]][ok]) < 1e-10) NA_real_
            else cor(res_base$completed[[v]][mask[, v]][ok], truth[ok])
  r_bio  <- if (sd(res_bio$completed[[v]][mask[, v]][ok]) < 1e-10) NA_real_
            else cor(res_bio$completed[[v]][mask[, v]][ok], truth[ok])
  row <- data.frame(
    trait      = v,
    n_held     = sum(ok),
    mean_RMSE  = rmse_mean,
    base_RMSE  = rmse_base,
    bio_RMSE   = rmse_bio,
    base_vs_mean = rmse_base / rmse_mean - 1,
    bio_vs_mean  = rmse_bio  / rmse_mean - 1,
    bio_vs_base  = rmse_bio  / rmse_base - 1,
    base_r = r_base,
    bio_r  = r_bio)
  rows[[length(rows) + 1L]] <- row
  cat(sprintf("%-15s  n=%3d  mean_RMSE=%10.3g  base_RMSE=%10.3g (r=%+.3f)  bio_RMSE=%10.3g (r=%+.3f)\n",
               v, sum(ok), rmse_mean, rmse_base, r_base, rmse_bio, r_bio))
  cat(sprintf("                  base vs mean: %+0.3f   bio vs mean: %+0.3f   bio vs base: %+0.3f\n\n",
               row$base_vs_mean, row$bio_vs_mean, row$bio_vs_base))
}
all_res <- do.call(rbind, rows)

# Persist
saveRDS(list(results = all_res,
              config  = list(n = length(sp_s), seed = SEED,
                              miss_frac = MISS_FRAC, epochs = EPOCHS,
                              n_imp = N_IMP),
              wall    = list(baseline_s = w_base, per_occur_s = w_bio)),
         "script/bench_plants_bioclim_n240.rds", compress = "xz")
cat_line("wrote script/bench_plants_bioclim_n240.rds")
