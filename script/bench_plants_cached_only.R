#!/usr/bin/env Rscript
# script/bench_plants_cached_only.R
#
# Plants bench using ONLY species we have in the local GBIF cache
# with valid per-occurrence points.  No new GBIF API calls -- avoids
# the rate-limit lockout that broke earlier n=3,400 and n=4,745 runs.
#
# Intersects cached species × BIEN traits × tree tips → runs at
# whatever n that gives us (expect ~2,500-3,000 on current cache).
#
# Three fits: baseline (no cov) / bio + safety_floor=TRUE / bio +
# safety_floor=FALSE.
#
# Invocation:
#   PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_plants_cached_only.R

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

cat_line <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

gbif_dir <- "script/data-cache/gbif"
wc_dir   <- "script/data-cache/worldclim"

# -------------------------------------------------------------------------
# Step 1: enumerate cached species with valid per-occurrence points
# -------------------------------------------------------------------------

cat_line("enumerating cached GBIF species with valid points ...")
cache_files <- list.files(gbif_dir, full.names = TRUE, pattern = "\\.rds$")
cat_line(sprintf("  %d total cache files", length(cache_files)))

cached_sp <- character(0)
for (f in cache_files) {
  x <- tryCatch(readRDS(f), error = function(e) NULL)
  if (!is.null(x) && !is.null(x$points) && nrow(x$points) > 0L) {
    cached_sp <- c(cached_sp, x$species)
  }
}
cat_line(sprintf("  %d species with valid per-occurrence points", length(cached_sp)))

# -------------------------------------------------------------------------
# Step 2: load BIEN + tree, intersect with cached species
# -------------------------------------------------------------------------

bien_cache_t <- "script/data-cache/bien_trait_means.rds"
bien_cache_tr <- "script/data-cache/bien_tree.rds"
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

# Intersect cached × BIEN × tree
sp <- Reduce(intersect, list(cached_sp, wide$species, tree_all$tip.label))
cat_line(sprintf("intersection (cached × BIEN × tree): %d species", length(sp)))
stopifnot(length(sp) >= 100L)

set.seed(SEED)
wide_s <- wide[wide$species %in% sp, , drop = FALSE]
rownames(wide_s) <- wide_s$species; wide_s$species <- NULL
tr <- ape::keep.tip(tree_all, sp)

# -------------------------------------------------------------------------
# Step 3: bioclim extract (fully cached, no live fetch)
# -------------------------------------------------------------------------

cat_line("extracting bioclim for", length(sp), "species (all from cache) ...")
t_bio <- proc.time()[["elapsed"]]
wc_df <- pull_worldclim_per_species(sp,
                                      gbif_cache_dir = gbif_dir,
                                      worldclim_cache_dir = wc_dir,
                                      verbose = FALSE)
cat_line(sprintf("bioclim extract done in %.0fs", proc.time()[["elapsed"]] - t_bio))
cov_bio <- wc_df[sp, grepl("^bio", colnames(wc_df)), drop = FALSE]

# Drop species with all-NA bioclim (should be 0 given cached-points filter)
na_rows <- rowSums(is.na(cov_bio)) == ncol(cov_bio)
if (any(na_rows)) {
  cat_line(sprintf("dropping %d species with all-NA bioclim (unexpected)", sum(na_rows)))
  keep <- !na_rows
  cov_bio <- cov_bio[keep, , drop = FALSE]
  wide_s  <- wide_s[keep, , drop = FALSE]
  sp      <- sp[keep]
  tr      <- ape::keep.tip(tr, sp)
}
cat_line("final n =", length(sp), "species")

# -------------------------------------------------------------------------
# Step 4: 30% MCAR mask
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
  cat_line(sprintf("  %-15s: %d held-out cells", v, sum(mask[, v])))
}

# -------------------------------------------------------------------------
# Step 5: three pigauto fits
# -------------------------------------------------------------------------

cat_line("=============== fit baseline (no covariates) ===============")
t0 <- proc.time()[["elapsed"]]
res_none <- pigauto::impute(df, tr,
                               epochs = EPOCHS, n_imputations = N_IMP,
                               verbose = FALSE, seed = SEED,
                               safety_floor = TRUE)
w_none <- proc.time()[["elapsed"]] - t0
cat_line(sprintf("baseline done in %.0fs", w_none))

cat_line("=============== fit bio + safety_floor=TRUE ===============")
t0 <- proc.time()[["elapsed"]]
res_bio <- pigauto::impute(df, tr, covariates = cov_bio,
                              epochs = EPOCHS, n_imputations = N_IMP,
                              verbose = FALSE, seed = SEED,
                              safety_floor = TRUE)
w_bio <- proc.time()[["elapsed"]] - t0
cat_line(sprintf("bio sf=TRUE done in %.0fs", w_bio))

cat_line("=============== fit bio + safety_floor=FALSE ===============")
t0 <- proc.time()[["elapsed"]]
res_off <- pigauto::impute(df, tr, covariates = cov_bio,
                              epochs = EPOCHS, n_imputations = N_IMP,
                              verbose = FALSE, seed = SEED,
                              safety_floor = FALSE)
w_off <- proc.time()[["elapsed"]] - t0
cat_line(sprintf("bio sf=FALSE done in %.0fs", w_off))

# -------------------------------------------------------------------------
# Step 6: report per-trait
# -------------------------------------------------------------------------

cat("\n\n============ RESULTS (n =", length(sp), ") ============\n\n")
cat(sprintf("%-15s %5s %10s %10s %10s %10s %7s %7s %7s %7s\n",
             "trait", "n", "mean", "none", "bio_on", "bio_off", "r_on", "r_off", "rat_on", "rat_off"))
rows <- list()
for (v in cont_cols) {
  if (!any(mask[, v])) next
  truth <- wide_s[[v]][mask[, v]]
  ok <- is.finite(truth)
  if (sum(ok) < 5L) next
  mean_pred <- mean(df[[v]], na.rm = TRUE)
  rmse_mean <- sqrt(mean((mean_pred - truth[ok])^2, na.rm = TRUE))
  # Use na.rm = TRUE throughout to guard against NA in predictions
  rmse_none <- sqrt(mean((res_none$completed[[v]][mask[, v]][ok] - truth[ok])^2,
                          na.rm = TRUE))
  rmse_bio  <- sqrt(mean((res_bio$completed[[v]][mask[, v]][ok]  - truth[ok])^2,
                          na.rm = TRUE))
  rmse_off  <- sqrt(mean((res_off$completed[[v]][mask[, v]][ok]  - truth[ok])^2,
                          na.rm = TRUE))
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
  r_none <- tryCatch({
    p <- res_none$completed[[v]][mask[, v]][ok]
    if (!is.finite(sd(p, na.rm = TRUE)) || sd(p, na.rm = TRUE) < 1e-10) NA_real_
    else suppressWarnings(stats::cor(p, truth[ok], use = "complete.obs"))
  }, error = function(e) NA_real_)
  cat(sprintf("%-15s %5d %10.3g %10.3g %10.3g %10.3g %+7.3f %+7.3f %7.3f %7.3f\n",
               v, sum(ok), rmse_mean, rmse_none, rmse_bio, rmse_off,
               ifelse(is.na(r_bio), NA, r_bio),
               ifelse(is.na(r_off), NA, r_off),
               rmse_bio / rmse_none, rmse_off / rmse_none))
  rows[[length(rows) + 1L]] <- data.frame(
    trait = v, n_held = sum(ok),
    mean_RMSE = rmse_mean, none_RMSE = rmse_none,
    bio_on_RMSE = rmse_bio, bio_off_RMSE = rmse_off,
    none_r = r_none, bio_on_r = r_bio, bio_off_r = r_off,
    ratio_on = rmse_bio / rmse_none, ratio_off = rmse_off / rmse_none)
}
all_res <- do.call(rbind, rows)

saveRDS(list(results = all_res, n = length(sp),
              config = list(seed = SEED, miss_frac = MISS_FRAC,
                             epochs = EPOCHS, n_imp = N_IMP)),
         "script/bench_plants_cached_only.rds", compress = "xz")
cat_line("wrote script/bench_plants_cached_only.rds")
