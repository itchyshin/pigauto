#!/usr/bin/env Rscript
# script/bench_globtherm_covariates.R
#
# Covariate-lift test on the GlobTherm thermal-tolerance database
# (Bennett et al. 2018, Sci. Data 5: 180022).
#
# Why this dataset is the strongest candidate yet:
#   - Tmax (CTmax) and tmin (CTmin) have a textbook latitudinal
#     gradient.  Tropical species tolerate hotter, polar species
#     tolerate colder.  ~0.5 deg C per degree latitude.
#   - Sister species often live at different latitudes (range
#     shifts during diversification), so latitude can carry
#     information that phylogeny alone does not encode.
#   - This is the dataset where species-level climate covariates
#     have the best chance of breaking the phylo-redundancy that
#     killed lift on BIEN (n=3,450) and Delhey (n=5,809).
#
# Subset: ectotherms (Lepidosauria + Insecta + Amphibia + Actinopteri)
#   ~793 species with Tmax + lat_max non-NA.  Ectotherms are where
#   thermal tolerance is most directly tied to ambient conditions.
#
# Phylogeny: cross-class taxonomic tree from Class/Order/Family/
#   Genus/Species via ape::as.phylo.formula() + Grafen brlen.
#   This mirrors the AmphiBIO + LepTraits patterns since a
#   comprehensive cross-phylum molecular tree is not available.
#
# Three fits: baseline / cov + sf=TRUE / cov + sf=FALSE.
#
# Invocation:
#   PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_globtherm_covariates.R

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

cat_line <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "),
                                ..., "\n", sep = "")

here <- "/Users/z3437171/Dropbox/Github Local/pigauto"
gt_path <- file.path(here, "script", "data-cache", "external", "globTherm.rda")
stopifnot(file.exists(gt_path))

cat_line("loading globTherm ...")
load(gt_path)
# globTherm ships as a tibble; convert to data.frame so rownames stick
globTherm <- as.data.frame(globTherm, stringsAsFactors = FALSE)
cat_line(sprintf("globTherm: %d rows x %d cols", nrow(globTherm),
                  ncol(globTherm)))

# Subset: ectotherms with Tmax + lat
ECTO_CLASSES <- c("Lepidosauria", "Insecta", "Amphibia", "Actinopteri",
                   "Arachnida", "Malacostraca", "Bivalvia", "Gastropoda")
gt <- globTherm
gt <- gt[gt$Class %in% ECTO_CLASSES, , drop = FALSE]
gt <- gt[!is.na(gt$Tmax) & !is.na(gt$lat_max), , drop = FALSE]
cat_line(sprintf("ectotherm subset (Tmax + lat_max known): %d rows",
                  nrow(gt)))

# Drop rows missing taxonomy
ok_tax <- !is.na(gt$Class) & !is.na(gt$Order) & !is.na(gt$Family) &
           nzchar(gt$Class) & nzchar(gt$Order) & nzchar(gt$Family) &
           !is.na(gt$Genus) & nzchar(gt$Genus) &
           !is.na(gt$scientificNameStd) & nzchar(gt$scientificNameStd)
gt <- gt[ok_tax, , drop = FALSE]
cat_line(sprintf("after taxonomy filter: %d rows", nrow(gt)))

# Build species names (use scientificNameStd)
gt$Species_Key <- gsub(" ", "_", gt$scientificNameStd)
gt <- gt[!duplicated(gt$Species_Key), , drop = FALSE]
cat_line(sprintf("after dedup: %d species", nrow(gt)))

# Targets
trait_cols <- c("Tmax", "tmin")
# Use both since some species have only one
trait_cols <- intersect(trait_cols, colnames(gt))
df_traits  <- gt[, trait_cols, drop = FALSE]
for (v in trait_cols) df_traits[[v]] <- suppressWarnings(as.numeric(df_traits[[v]]))

# Covariates: latitude, longitude, elevation, |lat| (abs lat = thermal niche),
# and lat-difference between max and min collection (ecological breadth)
gt$abs_lat_max <- abs(gt$lat_max)
gt$lat_diff    <- gt$lat_max - gt$lat_min  # NA-safe; both used at the same site
cov_cols <- c("lat_max", "long_max", "elevation_max", "abs_lat_max")
cov_cols <- intersect(cov_cols, colnames(gt))
cov_df   <- gt[, cov_cols, drop = FALSE]
for (v in cov_cols) cov_df[[v]] <- suppressWarnings(as.numeric(cov_df[[v]]))

# Drop species missing any covariate
ok_cov <- stats::complete.cases(cov_df)
cat_line(sprintf("species with all covariates non-NA: %d / %d",
                  sum(ok_cov), nrow(cov_df)))
df_traits <- df_traits[ok_cov, , drop = FALSE]
cov_df    <- cov_df[ok_cov, , drop = FALSE]
gt        <- gt[ok_cov, , drop = FALSE]
rownames(df_traits) <- gt$Species_Key
rownames(cov_df)    <- gt$Species_Key
cat_line(sprintf("final n = %d species, %d traits, %d covariates",
                  nrow(df_traits), ncol(df_traits), ncol(cov_df)))

# Build taxonomic tree
cat_line("building cross-class taxonomic tree (Class/Order/Family/Genus/Species) ...")
tax_df <- gt[, c("Class", "Order", "Family", "Genus", "Species_Key"),
              drop = FALSE]
tax_df[] <- lapply(tax_df, factor)
tree <- as.phylo(~Class/Order/Family/Genus/Species_Key,
                  data = tax_df, collapse = FALSE)
set.seed(SEED)
tree <- ape::collapse.singles(tree)
if (!ape::is.rooted(tree))
  tree <- ape::root.phylo(tree, outgroup = 1L, resolve.root = TRUE)
tree <- ape::multi2di(tree, random = TRUE)
tree <- compute.brlen(tree, method = "Grafen")
tree$edge.length[tree$edge.length <= 0] <- 1e-8
cat_line(sprintf("tree: %d tips, binary=%s, rooted=%s",
                  length(tree$tip.label), ape::is.binary(tree),
                  ape::is.rooted(tree)))

# Align: keep only species present in tree
matched <- intersect(rownames(df_traits), tree$tip.label)
df_traits <- df_traits[matched, , drop = FALSE]
cov_df    <- cov_df[matched, , drop = FALSE]
tree      <- ape::keep.tip(tree, matched)
cat_line(sprintf("after tree alignment: n = %d species", length(matched)))

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
  cat_line(sprintf("  %-15s: %d held-out cells", v, sum(mask[, v])))
}

# Three fits
cat_line("=============== fit baseline (no covariates) ===============")
t0 <- proc.time()[["elapsed"]]
res_none <- tryCatch(
  pigauto::impute(df_obs, tree,
                    epochs = EPOCHS, n_imputations = N_IMP,
                    verbose = FALSE, seed = SEED,
                    safety_floor = TRUE),
  error = function(e) { cat_line("baseline ERROR: ", conditionMessage(e)); NULL })
w_none <- proc.time()[["elapsed"]] - t0
cat_line(sprintf("baseline done in %.0fs", w_none))

cat_line("=============== fit cov + safety_floor=TRUE ===============")
t0 <- proc.time()[["elapsed"]]
res_bio <- tryCatch(
  pigauto::impute(df_obs, tree, covariates = cov_df,
                    epochs = EPOCHS, n_imputations = N_IMP,
                    verbose = FALSE, seed = SEED,
                    safety_floor = TRUE),
  error = function(e) { cat_line("cov sf=TRUE ERROR: ", conditionMessage(e)); NULL })
w_bio <- proc.time()[["elapsed"]] - t0
cat_line(sprintf("cov sf=TRUE done in %.0fs", w_bio))

cat_line("=============== fit cov + safety_floor=FALSE ===============")
t0 <- proc.time()[["elapsed"]]
res_off <- tryCatch(
  pigauto::impute(df_obs, tree, covariates = cov_df,
                    epochs = EPOCHS, n_imputations = N_IMP,
                    verbose = FALSE, seed = SEED,
                    safety_floor = FALSE),
  error = function(e) { cat_line("cov sf=FALSE ERROR: ", conditionMessage(e)); NULL })
w_off <- proc.time()[["elapsed"]] - t0
cat_line(sprintf("cov sf=FALSE done in %.0fs", w_off))

# Report
cat("\n\n============ RESULTS (n =", nrow(df_traits), ") ============\n\n")
cat(sprintf("%-15s %5s %10s %10s %10s %10s %8s %8s %8s %8s %8s\n",
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
  cat(sprintf("%-15s %5d %10.3g %10.3g %10.3g %10.3g %+8.3f %+8.3f %+8.3f %8.3f %8.3f\n",
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
                              ecto_classes = ECTO_CLASSES)),
         "script/bench_globtherm_covariates.rds", compress = "xz")
cat_line("wrote script/bench_globtherm_covariates.rds")
