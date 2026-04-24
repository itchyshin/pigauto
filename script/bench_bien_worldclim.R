#!/usr/bin/env Rscript
# script/bench_bien_worldclim.R
# Full plants bench rerun with WorldClim bioclim covariates.
# Run manually after PR merge to produce paper-ready numbers at
# n = 4,745 species.
#
# Invocation:
#   PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_bien_worldclim.R

suppressPackageStartupMessages({
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH", unset = "")
  if (nzchar(pkg_path)) devtools::load_all(pkg_path, quiet = TRUE)
  else library(pigauto)
  library(ape)
})
options(warn = -1L)

SEED <- 2026L; MISS_FRAC <- 0.30; N_IMP <- 20L; EPOCHS <- 300L

cache_trait <- "script/data-cache/bien_trait_means.rds"
cache_tree  <- "script/data-cache/bien_tree.rds"
stopifnot(file.exists(cache_trait), file.exists(cache_tree))

gbif_dir <- "script/data-cache/gbif"
wc_dir   <- "script/data-cache/worldclim"
stopifnot(dir.exists(gbif_dir), dir.exists(wc_dir))

trait_means <- readRDS(cache_trait)
tree_raw    <- readRDS(cache_tree)
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
matched <- intersect(wide$species, tree_all$tip.label)
set.seed(SEED)
# Parameterised by PIGAUTO_BENCH_N env var (default 4745)
N_TARGET <- as.integer(Sys.getenv("PIGAUTO_BENCH_N", "4745"))
sp <- sample(matched, min(N_TARGET, length(matched)))
wide_s <- wide[wide$species %in% sp, , drop = FALSE]
rownames(wide_s) <- wide_s$species; wide_s$species <- NULL
tr <- ape::keep.tip(tree_all, sp)

# Pre-fetch GBIF occurrences (with points) for any species not yet cached
# so that pull_worldclim_per_species() can do per-occurrence extraction.
# Cached species are instant; un-cached ones hit the API.
cat("[bench] pre-fetching GBIF occurrences for", length(sp), "species ...\n")
cat("[bench] (cached species instant; ~2-4 sec/species for uncached)\n")
t_gbif <- proc.time()[["elapsed"]]
gbif_df <- pull_gbif_centroids(sp,
                                 cache_dir = gbif_dir,
                                 occurrence_limit = 300L,
                                 sleep_ms = 50L,
                                 verbose = TRUE,
                                 store_points = TRUE,
                                 refresh_cache = FALSE)
cat(sprintf("[bench] GBIF pre-fetch done in %.0fs (%d species with valid centroids)\n",
             proc.time()[["elapsed"]] - t_gbif,
             sum(!is.na(gbif_df$centroid_lat))))

# Fetch bioclim for these species
cat("[bench] extracting bioclim for", length(sp), "species ...\n")
wc_df <- pull_worldclim_per_species(sp,
                                      gbif_cache_dir = gbif_dir,
                                      worldclim_cache_dir = wc_dir,
                                      verbose = TRUE)
cov_bio <- wc_df[sp, grepl("^bio", colnames(wc_df)), drop = FALSE]
# Drop species with all-NA bioclim (no GBIF hit) so pigauto sees clean covariates
na_rows <- rowSums(is.na(cov_bio)) == ncol(cov_bio)
if (any(na_rows)) {
  cat(sprintf("[bench] dropping %d species without any bioclim coverage\n", sum(na_rows)))
  keep <- !na_rows
  cov_bio <- cov_bio[keep, , drop = FALSE]
  wide_s  <- wide_s[keep, , drop = FALSE]
  sp      <- sp[keep]
  tr      <- ape::keep.tip(tr, sp)
}
cat("[bench] final n =", length(sp), "species after covariate filter\n")

# 30% MCAR mask
df <- wide_s
cont_cols <- colnames(df)
mask <- matrix(FALSE, nrow = nrow(df), ncol = length(cont_cols),
               dimnames = list(NULL, cont_cols))
for (v in cont_cols) {
  ok <- which(!is.na(wide_s[[v]]))
  idx <- sample(ok, round(MISS_FRAC * length(ok)))
  mask[idx, v] <- TRUE
  df[[v]][idx] <- NA_real_
}

cat("[bench] fit without covariates ...\n")
res_none <- pigauto::impute(df, tr, phylo_signal_gate = FALSE,
                               epochs = EPOCHS, n_imputations = N_IMP,
                               verbose = FALSE, seed = SEED)
cat("[bench] fit with bioclim ...\n")
res_bio  <- pigauto::impute(df, tr, covariates = cov_bio,
                               phylo_signal_gate = FALSE,
                               epochs = EPOCHS, n_imputations = N_IMP,
                               verbose = FALSE, seed = SEED)

for (v in cont_cols) {
  if (!any(mask[, v])) next
  truth <- wide_s[[v]][mask[, v]]
  ok <- is.finite(truth)
  r_none <- sqrt(mean((res_none$completed[[v]][mask[, v]][ok] - truth[ok])^2))
  r_bio  <- sqrt(mean((res_bio$completed[[v]][mask[, v]][ok]  - truth[ok])^2))
  corr_bio <- stats::cor(res_bio$completed[[v]][mask[, v]][ok], truth[ok])
  cat(sprintf("  %-15s: no-cov RMSE = %.3g, with-bio = %.3g (r = %.3f, ratio = %.3f)\n",
               v, r_none, r_bio, corr_bio, r_bio / r_none))
}
