#!/usr/bin/env Rscript
#
# script/walkthrough_covariates_precompute.R
#
# Pre-compute all results for the "Comparative study with covariates"
# walkthrough. Saves an RDS consumed by make_walkthrough_covariates_html.R.
#
# Output: script/walkthrough_covariates.rds
#
# Runtime: ~15–25 min on a laptop (5,809 species)

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  library(glmmTMB)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto",
    quiet = TRUE
  )
})

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "walkthrough_covariates.rds")

script_start <- proc.time()[["elapsed"]]

log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start), ..., "\n",
      sep = "")
  flush.console()
}

# -------------------------------------------------------------------------
# 1. Data setup
# -------------------------------------------------------------------------

log_line("Loading data...")
data(delhey5809, package = "pigauto")
data(tree_delhey, package = "pigauto")

df <- delhey5809
rownames(df) <- df$Species_Key
n <- nrow(df)

# Response traits
traits <- df[, c("lightness_male", "lightness_female")]

# Covariates (4, dropping redundant warmest-quarter variants)
covs <- df[, c("annual_mean_temperature", "annual_precipitation",
               "percent_tree_cover", "midLatitude")]

log_line(sprintf("Species: %d, traits: %d, covariates: %d", n,
                 ncol(traits), ncol(covs)))

# -------------------------------------------------------------------------
# 2. Introduce 20% MCAR missingness
# -------------------------------------------------------------------------

set.seed(42)
miss_frac <- 0.20
traits_miss <- traits
mask <- matrix(FALSE, nrow = n, ncol = ncol(traits))
colnames(mask) <- names(traits)

for (j in seq_len(ncol(traits_miss))) {
  n_miss <- round(n * miss_frac)
  idx <- sample.int(n, n_miss)
  traits_miss[[j]][idx] <- NA
  mask[idx, j] <- TRUE
}

n_na <- sum(is.na(as.matrix(traits_miss)))
log_line(sprintf("Missing: %d / %d cells (%.1f%%)",
                 n_na, length(as.matrix(traits_miss)),
                 100 * n_na / length(as.matrix(traits_miss))))

# -------------------------------------------------------------------------
# 3. Impute WITHOUT covariates
# -------------------------------------------------------------------------

log_line("Imputing without covariates...")
t0 <- proc.time()[["elapsed"]]
res_no_cov <- impute(traits_miss, tree_delhey,
                     epochs = 500L, verbose = FALSE,
                     eval_every = 25L, patience = 15L,
                     seed = 42L)
t1 <- proc.time()[["elapsed"]]
log_line(sprintf("  Done in %.1f min", (t1 - t0) / 60))

# -------------------------------------------------------------------------
# 4. Impute WITH covariates
# -------------------------------------------------------------------------

log_line("Imputing with covariates...")
t0 <- proc.time()[["elapsed"]]
res_cov <- impute(traits_miss, tree_delhey, covariates = covs,
                  epochs = 500L, verbose = FALSE,
                  eval_every = 25L, patience = 15L,
                  seed = 42L)
t1 <- proc.time()[["elapsed"]]
log_line(sprintf("  Done in %.1f min", (t1 - t0) / 60))

# -------------------------------------------------------------------------
# 5. RMSE comparison on held-out cells
# -------------------------------------------------------------------------

truth_mat <- as.matrix(traits)
completed_no_cov <- as.matrix(res_no_cov$completed[, names(traits)])
completed_cov    <- as.matrix(res_cov$completed[, names(traits)])

rmse_compare <- data.frame(
  trait = character(), method = character(),
  rmse = numeric(), pearson_r = numeric(),
  stringsAsFactors = FALSE
)

for (j in names(traits)) {
  m_j <- mask[, j]
  t_j <- truth_mat[m_j, j]

  for (mname in c("no_covariates", "with_covariates")) {
    p_j <- if (mname == "no_covariates") completed_no_cov[m_j, j]
           else completed_cov[m_j, j]
    ok <- is.finite(t_j) & is.finite(p_j)
    rmse_compare <- rbind(rmse_compare, data.frame(
      trait     = j,
      method    = mname,
      rmse      = sqrt(mean((t_j[ok] - p_j[ok])^2)),
      pearson_r = cor(t_j[ok], p_j[ok]),
      stringsAsFactors = FALSE
    ))
  }
}

log_line("RMSE comparison:")
print(rmse_compare)

# -------------------------------------------------------------------------
# 6. Multiple imputation with covariates (m = 50)
# -------------------------------------------------------------------------

log_line("Running multi_impute (m = 50, with covariates)...")
t0 <- proc.time()[["elapsed"]]
mi <- multi_impute(traits_miss, tree_delhey, m = 50L,
                   covariates = covs,
                   epochs = 500L, verbose = FALSE,
                   eval_every = 25L, patience = 15L,
                   seed = 42L)
t1 <- proc.time()[["elapsed"]]
log_line(sprintf("  Done in %.1f min", (t1 - t0) / 60))

# -------------------------------------------------------------------------
# 7. Downstream regression: lightness_male ~ environment + phylogeny
# -------------------------------------------------------------------------

log_line("Computing Vphy (correlation matrix)...")
Vphy <- cov2cor(ape::vcv(tree_delhey))

log_line("Fitting 50 glmmTMB models...")
t0 <- proc.time()[["elapsed"]]

fits <- with_imputations(mi, function(d) {
  d$temperature   <- covs$annual_mean_temperature
  d$precipitation <- covs$annual_precipitation
  d$tree_cover    <- covs$percent_tree_cover
  d$latitude      <- covs$midLatitude
  d$species       <- factor(rownames(d), levels = rownames(Vphy))
  d$dummy         <- factor(1)

  glmmTMB(
    lightness_male ~ temperature + precipitation +
      tree_cover + latitude +
      propto(0 + species | dummy, Vphy),
    data = d
  )
})

t1 <- proc.time()[["elapsed"]]
log_line(sprintf("  Done in %.1f min", (t1 - t0) / 60))

# Check how many fits succeeded
n_ok <- sum(vapply(fits$fits, function(f) !inherits(f, "pigauto_mi_error"),
                   logical(1)))
log_line(sprintf("  Successful fits: %d / %d", n_ok, length(fits$fits)))

# -------------------------------------------------------------------------
# 8. Pool with Rubin's rules
# -------------------------------------------------------------------------

log_line("Pooling with Rubin's rules...")
pooled <- pool_mi(
  fits,
  coef_fun = function(f) fixef(f)$cond,
  vcov_fun = function(f) vcov(f)$cond
)

log_line("Pooled results:")
print(pooled)

# -------------------------------------------------------------------------
# 9. Save everything
# -------------------------------------------------------------------------

meta <- list(
  n_species     = n,
  n_traits      = ncol(traits),
  n_covs        = ncol(covs),
  miss_frac     = miss_frac,
  epochs        = 500L,
  m_imputations = 50L,
  seed          = 42L,
  timestamp     = Sys.time(),
  wall_time     = proc.time()[["elapsed"]] - script_start
)

saveRDS(list(
  rmse_compare = rmse_compare,
  pooled       = pooled,
  meta         = meta,
  n_fits_ok    = n_ok,
  data_summary = list(
    traits_summary  = summary(traits),
    covs_summary    = summary(covs),
    n_missing       = colSums(is.na(traits_miss))
  )
), out_rds)

log_line(sprintf("Wrote %s (%.1f KB)", out_rds, file.size(out_rds) / 1024))
log_line(sprintf("Total wall time: %.1f min", meta$wall_time / 60))
log_line("done")
