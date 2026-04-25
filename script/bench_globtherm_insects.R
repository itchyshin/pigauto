#!/usr/bin/env Rscript
# script/bench_globtherm_insects.R
#
# Insecta-only subset of the GlobTherm bench.  Tests whether the
# Tmax lift survives within a single class (insects), independent
# of cross-class taxonomic-tree artefacts.
#
# Insects are the cleanest CTmax-vs-latitude test: large within-
# class latitudinal range, well-known thermal biology, and 213
# species in GlobTherm with both Tmax and lat known.
#
# Three fits: baseline / cov + sf=TRUE / cov + sf=FALSE.

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
globTherm <- as.data.frame(globTherm, stringsAsFactors = FALSE)
gt <- globTherm[globTherm$Class == "Insecta", , drop = FALSE]
gt <- gt[!is.na(gt$Tmax) & !is.na(gt$lat_max), , drop = FALSE]
ok_tax <- !is.na(gt$Order) & !is.na(gt$Family) &
           nzchar(gt$Order) & nzchar(gt$Family) &
           !is.na(gt$Genus) & nzchar(gt$Genus) &
           !is.na(gt$scientificNameStd) & nzchar(gt$scientificNameStd)
gt <- gt[ok_tax, , drop = FALSE]
gt$Species_Key <- gsub(" ", "_", gt$scientificNameStd)
gt <- gt[!duplicated(gt$Species_Key), , drop = FALSE]
gt$abs_lat_max <- abs(gt$lat_max)
cat_line(sprintf("Insecta subset: %d species", nrow(gt)))

trait_cols <- "Tmax"
df_traits <- gt[, trait_cols, drop = FALSE]
df_traits$Tmax <- as.numeric(df_traits$Tmax)

cov_cols <- c("lat_max", "long_max", "elevation_max", "abs_lat_max")
cov_df <- gt[, cov_cols, drop = FALSE]
for (v in cov_cols) cov_df[[v]] <- as.numeric(cov_df[[v]])

ok_cov <- stats::complete.cases(cov_df)
df_traits <- df_traits[ok_cov, , drop = FALSE]
cov_df    <- cov_df[ok_cov, , drop = FALSE]
gt        <- gt[ok_cov, , drop = FALSE]
rownames(df_traits) <- gt$Species_Key
rownames(cov_df) <- gt$Species_Key
cat_line(sprintf("with covariates: %d species", nrow(df_traits)))

# Build Insecta taxonomic tree (Order/Family/Genus/Species)
tax_df <- gt[, c("Order", "Family", "Genus", "Species_Key"), drop = FALSE]
tax_df[] <- lapply(tax_df, factor)
tree <- as.phylo(~Order/Family/Genus/Species_Key, data = tax_df, collapse = FALSE)
set.seed(SEED)
tree <- ape::collapse.singles(tree)
if (!ape::is.rooted(tree))
  tree <- ape::root.phylo(tree, outgroup = 1L, resolve.root = TRUE)
tree <- ape::multi2di(tree, random = TRUE)
tree <- compute.brlen(tree, method = "Grafen")
tree$edge.length[tree$edge.length <= 0] <- 1e-8
cat_line(sprintf("tree: %d tips", length(tree$tip.label)))

matched <- intersect(rownames(df_traits), tree$tip.label)
df_traits <- df_traits[matched, , drop = FALSE]
cov_df    <- cov_df[matched, , drop = FALSE]
tree      <- ape::keep.tip(tree, matched)
cat_line(sprintf("after tree alignment: n = %d species", length(matched)))

# 30% MCAR mask
mask <- matrix(FALSE, nrow = nrow(df_traits), ncol = 1L,
                dimnames = list(NULL, "Tmax"))
df_obs <- df_traits
set.seed(SEED + 1L)
ok <- which(!is.na(df_traits$Tmax))
idx <- sample(ok, round(MISS_FRAC * length(ok)))
mask[idx, "Tmax"] <- TRUE
df_obs$Tmax[idx] <- NA_real_
cat_line(sprintf("Tmax: %d held-out cells", sum(mask[, "Tmax"])))

cat_line("=============== fit baseline (no covariates) ===============")
t0 <- proc.time()[["elapsed"]]
res_none <- tryCatch(
  pigauto::impute(df_obs, tree,
                    epochs = EPOCHS, n_imputations = N_IMP,
                    verbose = FALSE, seed = SEED,
                    safety_floor = TRUE),
  error = function(e) { cat_line("baseline ERROR: ", conditionMessage(e)); NULL })
cat_line(sprintf("done in %.0fs", proc.time()[["elapsed"]] - t0))

cat_line("=============== fit cov + safety_floor=TRUE ===============")
t0 <- proc.time()[["elapsed"]]
res_bio <- tryCatch(
  pigauto::impute(df_obs, tree, covariates = cov_df,
                    epochs = EPOCHS, n_imputations = N_IMP,
                    verbose = FALSE, seed = SEED,
                    safety_floor = TRUE),
  error = function(e) { cat_line("cov sf=TRUE ERROR: ", conditionMessage(e)); NULL })
cat_line(sprintf("done in %.0fs", proc.time()[["elapsed"]] - t0))

cat_line("=============== fit cov + safety_floor=FALSE ===============")
t0 <- proc.time()[["elapsed"]]
res_off <- tryCatch(
  pigauto::impute(df_obs, tree, covariates = cov_df,
                    epochs = EPOCHS, n_imputations = N_IMP,
                    verbose = FALSE, seed = SEED,
                    safety_floor = FALSE),
  error = function(e) { cat_line("cov sf=FALSE ERROR: ", conditionMessage(e)); NULL })
cat_line(sprintf("done in %.0fs", proc.time()[["elapsed"]] - t0))

# Report
cat("\n\n============ RESULTS (Insecta-only, n =", nrow(df_traits), ") ============\n\n")

safe_get <- function(res, idx) {
  if (is.null(res)) return(rep(NA_real_, length(idx)))
  res$completed$Tmax[idx]
}
safe_cor <- function(p, t) {
  ok <- is.finite(p) & is.finite(t)
  if (sum(ok) < 5L) return(NA_real_)
  if (sd(p[ok], na.rm = TRUE) < 1e-10) return(NA_real_)
  suppressWarnings(stats::cor(p[ok], t[ok]))
}
truth <- df_traits$Tmax[mask[, "Tmax"]]
ok <- is.finite(truth)
held_idx <- which(mask[, "Tmax"])

mean_pred <- mean(df_obs$Tmax, na.rm = TRUE)
rmse_mean <- sqrt(mean((mean_pred - truth[ok])^2, na.rm = TRUE))
rmse_none <- sqrt(mean((safe_get(res_none, held_idx)[ok] - truth[ok])^2, na.rm = TRUE))
rmse_bio  <- sqrt(mean((safe_get(res_bio,  held_idx)[ok] - truth[ok])^2, na.rm = TRUE))
rmse_off  <- sqrt(mean((safe_get(res_off,  held_idx)[ok] - truth[ok])^2, na.rm = TRUE))
r_none <- safe_cor(safe_get(res_none, held_idx)[ok], truth[ok])
r_bio  <- safe_cor(safe_get(res_bio,  held_idx)[ok], truth[ok])
r_off  <- safe_cor(safe_get(res_off,  held_idx)[ok], truth[ok])

cat(sprintf("trait %15s n_held %d\n", "Tmax", sum(ok)))
cat(sprintf("  mean   baseline  RMSE = %8.3f\n", rmse_mean))
cat(sprintf("  pigauto no-cov   RMSE = %8.3f, r = %+0.3f\n", rmse_none, r_none))
cat(sprintf("  pigauto cov sf=T RMSE = %8.3f, r = %s, ratio = %0.3f\n",
             rmse_bio, ifelse(is.na(r_bio), "  NA  ", sprintf("%+0.3f", r_bio)),
             rmse_bio / rmse_none))
cat(sprintf("  pigauto cov sf=F RMSE = %8.3f, r = %s, ratio = %0.3f\n",
             rmse_off, ifelse(is.na(r_off), "  NA  ", sprintf("%+0.3f", r_off)),
             rmse_off / rmse_none))

results <- data.frame(
  trait = "Tmax", n_held = sum(ok),
  mean_RMSE = rmse_mean, none_RMSE = rmse_none,
  cov_on_RMSE = rmse_bio, cov_off_RMSE = rmse_off,
  none_r = r_none, cov_on_r = r_bio, cov_off_r = r_off,
  ratio_on  = rmse_bio / rmse_none,
  ratio_off = rmse_off / rmse_none)

saveRDS(list(results = results, n = nrow(df_traits),
              config = list(seed = SEED, miss_frac = MISS_FRAC,
                              epochs = EPOCHS, n_imp = N_IMP,
                              class = "Insecta")),
         "script/bench_globtherm_insects.rds", compress = "xz")
cat_line("wrote script/bench_globtherm_insects.rds")
