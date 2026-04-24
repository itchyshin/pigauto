#!/usr/bin/env Rscript
# script/bench_plants_sf_off_diagnostic.R
# Diagnostic: does safety_floor = FALSE unlock covariate lift?
# Reuses the exact n=271 data setup from bench_plants_bioclim_n240.R.

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

cat_line <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

prev <- readRDS("script/bench_plants_bioclim_n240.rds")
cat_line("loaded existing bench rds")

wc_all <- readRDS("tests/testthat/fixtures/worldclim_plants_300.rds")
trait_means <- readRDS("script/data-cache/bien_trait_means.rds")
tree_raw    <- readRDS("script/data-cache/bien_tree.rds")
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

iqr_cols <- grep("_iqr$", colnames(wc_all), value = TRUE)
sp_with_iqr <- rownames(wc_all)[rowSums(wc_all[, iqr_cols] > 0, na.rm = TRUE) > 0]
matched <- Reduce(intersect, list(wide$species, tree_all$tip.label, sp_with_iqr))
set.seed(SEED)
sp_s <- matched
wide_s <- wide[wide$species %in% sp_s, , drop = FALSE]
rownames(wide_s) <- wide_s$species; wide_s$species <- NULL
cov_bio <- wc_all[sp_s, grepl("^bio", colnames(wc_all)), drop = FALSE]
tree_s <- ape::keep.tip(tree_all, sp_s)

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

cat_line("using n =", length(sp_s), "species")
cat_line("=== FIT: safety_floor = FALSE + bioclim (force gate open) ===")
t0 <- proc.time()[["elapsed"]]
res_off <- pigauto::impute(df, tree_s, covariates = cov_bio,
                              epochs = EPOCHS, n_imputations = N_IMP,
                              verbose = FALSE, seed = SEED,
                              safety_floor = FALSE)
wall <- proc.time()[["elapsed"]] - t0
cat_line(sprintf("done in %.1fs", wall))

cat("\n\n=============== DIAGNOSTIC: safety_floor = FALSE + bioclim ===============\n\n")
cat(sprintf("%-15s %6s %10s %10s %14s %10s\n",
             "trait", "n_held", "mean_RMSE", "base_RMSE", "bio_off_RMSE", "bio_vs_base"))
for (v in cont_cols) {
  if (!any(mask[, v])) next
  truth <- wide_s[[v]][mask[, v]]
  ok <- is.finite(truth)
  if (sum(ok) < 10L) next
  mean_pred <- mean(df[[v]], na.rm = TRUE)
  rmse_mean <- sqrt(mean((mean_pred - truth[ok])^2))
  rmse_off  <- sqrt(mean((res_off$completed[[v]][mask[, v]][ok] - truth[ok])^2))
  prev_row <- prev$results[prev$results$trait == v, ]
  base_RMSE <- if (nrow(prev_row) > 0) prev_row$base_RMSE else NA
  ratio <- rmse_off / base_RMSE - 1
  r_off <- if (sd(res_off$completed[[v]][mask[, v]][ok]) < 1e-10) NA_real_
           else cor(res_off$completed[[v]][mask[, v]][ok], truth[ok])
  cat(sprintf("%-15s %6d %10.3g %10.3g %14.3g %+10.3f (r=%+0.3f)\n",
               v, sum(ok), rmse_mean, base_RMSE, rmse_off, ratio, r_off))
}
saveRDS(list(res = res_off, wall = wall),
         "script/bench_plants_bioclim_sf_off.rds", compress = "xz")
cat_line("wrote script/bench_plants_bioclim_sf_off.rds")
