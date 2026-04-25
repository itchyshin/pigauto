#!/usr/bin/env Rscript
# script/bench_leptraits_covariates.R
#
# Covariate-lift test on LepTraits 1.0 (Shirey et al. 2022, Sci.
# Data 9: 382), the global butterfly trait dataset.
#
# Targets (continuous):
#   - WS_L  (wingspan lower bound, mm)
#   - FW_L  (forewing length lower bound, mm)
#   - FlightDuration (number of flight months per year)
#   - NumberOfHostplantFamilies (count, log-transformed)
#
# Covariates: 12 monthly flight indicators (Jan, Feb, ..., Dec),
#   each binary 0/1 ("does this species fly in this month?").
#   Treats per-month phenology as a 12-D climate-correlated
#   feature -- tropical species fly year-round, temperate species
#   seasonally, polar species briefly in summer.  Sister species
#   often diverge in flight phenology, so this carries information
#   beyond phylogeny.
#
# Three fits: baseline / cov + sf=TRUE / cov + sf=FALSE.
#
# Invocation:
#   PIGAUTO_PKG_PATH="$(pwd)" Rscript script/bench_leptraits_covariates.R

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
csv_path <- file.path(here, "script", "data-cache", "external", "LepTraits",
                      "consensus", "consensus.csv")
stopifnot(file.exists(csv_path))

cat_line("loading LepTraits consensus ...")
lt <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
cat_line(sprintf("LepTraits: %d rows x %d cols", nrow(lt), ncol(lt)))

# Targets + covariates
trait_cols <- c("WS_L", "FW_L", "FlightDuration", "NumberOfHostplantFamilies")
cov_cols   <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
trait_cols <- intersect(trait_cols, colnames(lt))
cov_cols   <- intersect(cov_cols,   colnames(lt))
stopifnot(length(trait_cols) >= 2L, length(cov_cols) == 12L)

for (v in trait_cols) lt[[v]] <- suppressWarnings(as.numeric(lt[[v]]))
for (v in cov_cols)   lt[[v]] <- suppressWarnings(as.numeric(lt[[v]]))

# Log-transform skewed counts/sizes
for (v in c("WS_L", "FW_L", "NumberOfHostplantFamilies")) {
  if (v %in% trait_cols) {
    ok <- !is.na(lt[[v]]) & lt[[v]] > 0
    lt[[v]] <- ifelse(ok, log(lt[[v]]), NA_real_)
  }
}

# Drop species with all-NA covariates (those have no recorded flight months)
cov_complete <- stats::complete.cases(lt[, cov_cols, drop = FALSE])
cat_line(sprintf("species with all 12 monthly indicators non-NA: %d / %d",
                  sum(cov_complete), nrow(lt)))

# Also require at least 1 month active (zero across all 12 means no info)
has_any_month <- rowSums(lt[, cov_cols, drop = FALSE] > 0, na.rm = TRUE) >= 1L
keep <- cov_complete & has_any_month
cat_line(sprintf("species with cov + >= 1 active month: %d / %d",
                  sum(keep), nrow(lt)))
lt <- lt[keep, , drop = FALSE]

# Drop species with incomplete tax
tax_ok <- !is.na(lt$Family) & !is.na(lt$Genus) & !is.na(lt$Species) &
           nzchar(lt$Family) & nzchar(lt$Genus) & nzchar(lt$Species)
lt <- lt[tax_ok, , drop = FALSE]
cat_line(sprintf("after tax filter: %d species", nrow(lt)))

# Build species key + dedup
lt$Species_Key <- gsub(" ", "_", trimws(lt$Species))
lt <- lt[!duplicated(lt$Species_Key), , drop = FALSE]
cat_line(sprintf("after dedup: %d species", nrow(lt)))

# Subsample
set.seed(SEED)
if (nrow(lt) > N_TARGET) {
  lt <- lt[sample.int(nrow(lt), N_TARGET), , drop = FALSE]
  cat_line(sprintf("sub-sampled to %d species", nrow(lt)))
}

# Build taxonomic tree (Family / Genus / Species)
cat_line("building taxonomic tree ...")
tax_df <- lt[, c("Family", "Genus", "Species_Key"), drop = FALSE]
tax_df[] <- lapply(tax_df, factor)
tree <- as.phylo(~Family/Genus/Species_Key, data = tax_df, collapse = FALSE)
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

# Convert lt to data.frame (drop tibble class if present) so rownames stick
lt <- as.data.frame(lt, stringsAsFactors = FALSE)
rownames(lt) <- lt$Species_Key

matched <- intersect(rownames(lt), tree$tip.label)
lt   <- lt[matched, , drop = FALSE]
tree <- ape::keep.tip(tree, matched)
cat_line(sprintf("after tree alignment: n = %d species", length(matched)))

df_traits <- lt[, trait_cols, drop = FALSE]
cov_df    <- lt[, cov_cols,   drop = FALSE]

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
         "script/bench_leptraits_covariates.rds", compress = "xz")
cat_line("wrote script/bench_leptraits_covariates.rds")
