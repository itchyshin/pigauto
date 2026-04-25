#!/usr/bin/env Rscript
# script/bench_amphibio_covariates.R
#
# Covariate-lift test on AmphiBIO with the bundled climate-zone
# binary occupancy indicators serving as covariates.  Real species-
# level amphibian phylogeny (taxonomic, Grafen branch lengths) +
# climate zones from species ranges.
#
# AmphiBIO climate columns are 0/1 occupancy indicators:
#   Wet_warm, Wet_cold, Dry_warm, Dry_cold
# These differ from PanTHERIA's continuous Precip_Mean / Temp_Mean
# (which are scaled climate values).  Test whether binary climate-
# zone indicators carry predictive signal beyond phylogeny.
#
# Targets: Body_size_mm, Body_mass_g, Longevity_max_y,
#           Age_at_maturity_min_y, Litter_size_min_n
#
# Three fits: baseline / cov + sf=TRUE / cov + sf=FALSE.
#
# Invocation:
#   PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_amphibio_covariates.R

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
N_TARGET  <- as.integer(Sys.getenv("PIGAUTO_BENCH_N", "1500"))

cat_line <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "),
                                ..., "\n", sep = "")

here     <- "/Users/z3437171/Dropbox/Github Local/pigauto"
csv_path <- file.path(here, "script", "data-cache", "AmphiBIO_v1.csv")
stopifnot(file.exists(csv_path))

cat_line("loading AmphiBIO ...")
amphi <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
cat_line(sprintf("AmphiBIO: %d rows x %d cols", nrow(amphi), ncol(amphi)))

# Resolve species + tax columns
amphi <- amphi[nzchar(amphi$Species) & !is.na(amphi$Species), , drop = FALSE]
amphi <- amphi[!duplicated(amphi$Species), , drop = FALSE]
tax_ok <- nzchar(amphi$Order) & nzchar(amphi$Family) & nzchar(amphi$Genus) &
           nzchar(amphi$Species)
amphi <- amphi[tax_ok, , drop = FALSE]
cat_line(sprintf("after dedup + tax filter: %d species", nrow(amphi)))

# Targets (continuous, log-skew where appropriate)
trait_cols <- c("Body_size_mm", "Body_mass_g", "Longevity_max_y",
                 "Age_at_maturity_min_y", "Litter_size_min_n")
cov_cols   <- c("Wet_warm", "Wet_cold", "Dry_warm", "Dry_cold")
trait_cols <- intersect(trait_cols, colnames(amphi))
cov_cols   <- intersect(cov_cols,   colnames(amphi))
stopifnot(length(trait_cols) >= 3L, length(cov_cols) >= 2L)

for (v in trait_cols) amphi[[v]] <- suppressWarnings(as.numeric(amphi[[v]]))
for (v in cov_cols)   amphi[[v]] <- suppressWarnings(as.numeric(amphi[[v]]))

# Log-transform skewed traits
for (v in c("Body_size_mm", "Body_mass_g", "Longevity_max_y",
             "Age_at_maturity_min_y", "Litter_size_min_n")) {
  if (v %in% trait_cols) {
    ok <- !is.na(amphi[[v]]) & amphi[[v]] > 0
    amphi[[v]] <- ifelse(ok, log(amphi[[v]]), NA_real_)
  }
}

# Drop species with all-NA covariates (a row of all zeros is OK -- means
# species occurs in none of the four zones, which is unusual but possible)
cov_complete <- stats::complete.cases(amphi[, cov_cols, drop = FALSE])
cat_line(sprintf("species with all 4 climate zones non-NA: %d / %d",
                  sum(cov_complete), nrow(amphi)))
amphi <- amphi[cov_complete, , drop = FALSE]

# Subsample
set.seed(SEED)
if (nrow(amphi) > N_TARGET) {
  amphi <- amphi[sample.int(nrow(amphi), N_TARGET), , drop = FALSE]
  cat_line(sprintf("sub-sampled to %d species", nrow(amphi)))
}

# Build taxonomic tree from Order/Family/Genus/Species
cat_line("building taxonomic tree ...")
tax_df <- amphi[, c("Order", "Family", "Genus", "Species"), drop = FALSE]
tax_df[] <- lapply(tax_df, factor)
tree <- as.phylo(~Order/Family/Genus/Species, data = tax_df, collapse = FALSE)
set.seed(SEED)
tree <- ape::collapse.singles(tree)
if (!ape::is.rooted(tree)) {
  tree <- ape::root.phylo(tree, outgroup = 1L, resolve.root = TRUE)
}
tree <- ape::multi2di(tree, random = TRUE)
tree <- compute.brlen(tree, method = "Grafen")
tree$edge.length[tree$edge.length <= 0] <- 1e-8
cat_line(sprintf("tree: %d tips, binary=%s, rooted=%s",
                  length(tree$tip.label), ape::is.binary(tree),
                  ape::is.rooted(tree)))

# Align species to tree
rownames(amphi) <- amphi$Species
matched <- intersect(amphi$Species, tree$tip.label)
amphi <- amphi[matched, , drop = FALSE]
tree  <- ape::keep.tip(tree, matched)
cat_line(sprintf("matched: n = %d species", length(matched)))

df_traits <- amphi[, trait_cols, drop = FALSE]
cov_df    <- amphi[, cov_cols,   drop = FALSE]

# Drop traits with too few observations
keep_trait <- vapply(df_traits, function(x) sum(!is.na(x)) >= 50L, logical(1))
trait_cols <- trait_cols[keep_trait]
df_traits  <- df_traits[, trait_cols, drop = FALSE]
cat_line(sprintf("kept traits with >= 50 non-NA: %d (%s)",
                  length(trait_cols), paste(trait_cols, collapse = ", ")))

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
cat(sprintf("%-30s %5s %10s %10s %10s %10s %8s %8s %8s %8s %8s\n",
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
  cat(sprintf("%-30s %5d %10.3g %10.3g %10.3g %10.3g %+8.3f %+8.3f %+8.3f %8.3f %8.3f\n",
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
         "script/bench_amphibio_covariates.rds", compress = "xz")
cat_line("wrote script/bench_amphibio_covariates.rds")
