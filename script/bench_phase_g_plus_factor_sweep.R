#!/usr/bin/env Rscript
#
# script/bench_phase_g_plus_factor_sweep.R
#
# Phase G+: empirically tune `clamp_factor` for the AVONET Mass
# tail-extrapolation mode.
#
# Phase G v1 shipped with `clamp_factor = 5` (Tukey-style heuristic).
# That achieved -74 % RMSE reduction on seed 2030 but the Casuarius
# prediction (167,847 g) still landed just under the 175,000-g cap.
# This sweep checks whether a tighter factor extracts more juice
# without harming the seeds where Phase G's effect is small (2031,
# 2032).
#
# Pre-registered question:
#   For factor in {2, 3, 4, 5, 7, 10, Inf}, what does Mass RMSE look
#   like across all 3 seeds at N_IMP = 20?
#
# Default-flip rule (proposed):
#   Flip the default clamp_factor in v0.9.2 to the value F* that
#   minimises mean(RMSE) across 3 seeds AND has at most a 5 % RMSE
#   degradation on any seed compared to factor = Inf.
#   (Inf = no clamp, equivalent to clamp_outliers = FALSE.)
#
# Bench design:
#   3 seeds  x  1 impute() run each (clamp OFF, capture imputed_datasets)
#   x  post-hoc clamp at 7 factors  =  3 fits + 21 re-pools.
#   Wall: 3 x ~9 min = ~27 min.
#
# Output: script/bench_phase_g_plus_factor_sweep.{rds,md}

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH", unset = "")
  if (nzchar(pkg_path) && dir.exists(pkg_path) &&
      file.exists(file.path(pkg_path, "DESCRIPTION"))) {
    devtools::load_all(pkg_path, quiet = TRUE)
  } else {
    library(pigauto)
  }
})

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "bench_phase_g_plus_factor_sweep.rds")
out_md  <- file.path(here, "script", "bench_phase_g_plus_factor_sweep.md")

SEEDS     <- c(2030L, 2031L, 2032L)
MISS_FRAC <- 0.30
N_SUB     <- 1500L
N_IMP     <- 20L
FACTORS   <- c(2, 3, 4, 5, 7, 10, Inf)   # Inf = no clamp

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - t0),
      ..., "\n", sep = "")
  flush.console()
}

per_seed <- function(seed) {
  log_line(sprintf("Seed %d: data + mask ...", seed))
  e <- new.env(parent = emptyenv())
  utils::data("avonet_full", package = "pigauto", envir = e)
  utils::data("tree_full",   package = "pigauto", envir = e)
  df   <- e$avonet_full
  tree <- e$tree_full
  rownames(df) <- df$Species_Key
  df$Species_Key <- NULL

  set.seed(seed)
  keep <- sample(rownames(df), N_SUB)
  df   <- df[keep, , drop = FALSE]
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
  df   <- df[tree$tip.label, , drop = FALSE]

  set.seed(seed)
  df_truth  <- df
  mask_test <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df),
                      dimnames = list(rownames(df), names(df)))
  for (v in names(df)) {
    obs_idx <- which(!is.na(df[[v]]))
    to_hide <- sample(obs_idx, ceiling(length(obs_idx) * MISS_FRAC))
    mask_test[to_hide, v] <- TRUE
  }
  df_miss <- df
  for (v in names(df)) df_miss[[v]][mask_test[, v]] <- NA

  log_line(sprintf("Seed %d: impute(N_IMP=%d, clamp_outliers=FALSE) ...",
                   seed, N_IMP))
  res <- pigauto::impute(df_miss, tree,
                          epochs = 200L,
                          n_imputations = N_IMP,
                          clamp_outliers = FALSE,        # capture uncapped draws
                          verbose = FALSE,
                          seed = seed)

  # imputed_datasets: list of M data.frames (each n_species x n_traits)
  mi <- res$prediction$imputed_datasets
  if (is.null(mi) || length(mi) == 0L) stop("No imputed_datasets")

  mass_idx   <- which(mask_test[, "Mass"])
  mass_truth <- df_truth$Mass[mass_idx]

  # Per-imputation Mass at the held-out cells (n_held x M)
  per_imp <- vapply(mi, function(d) as.numeric(d$Mass[mass_idx]),
                    numeric(length(mass_idx)))

  # Recover obs_max from the trait_map of the fit (Phase G recorded it)
  tm <- NULL
  for (t in res$fit$trait_map) {
    if (identical(t$name, "Mass")) { tm <- t; break }
  }
  obs_max <- tm$obs_max
  log_line(sprintf("Seed %d: training obs_max(Mass) = %g g", seed, obs_max))

  # ---- Post-hoc clamp sweep ---------------------------------------------
  rmse_per_factor <- numeric(length(FACTORS))
  names(rmse_per_factor) <- as.character(FACTORS)
  top1_per_factor <- numeric(length(FACTORS))
  names(top1_per_factor) <- as.character(FACTORS)
  for (i in seq_along(FACTORS)) {
    f   <- FACTORS[i]
    cap <- if (is.finite(f)) obs_max * f else Inf
    clamped <- pmin(per_imp, cap)
    pooled  <- apply(clamped, 1L, stats::median, na.rm = TRUE)
    rmse_per_factor[i] <- sqrt(mean((mass_truth - pooled)^2,
                                     na.rm = TRUE))
    top1_per_factor[i] <- max(pooled, na.rm = TRUE)
    log_line(sprintf("  Seed %d factor=%5g cap=%9g  RMSE=%.0f  top1_pool=%.0f",
                     seed, f, cap, rmse_per_factor[i], top1_per_factor[i]))
  }

  data.frame(
    seed       = seed,
    factor     = FACTORS,
    rmse       = rmse_per_factor,
    top1_pool  = top1_per_factor,
    obs_max    = obs_max,
    n_test     = length(mass_idx),
    stringsAsFactors = FALSE
  )
}

# ---- Sweep -----------------------------------------------------------------

results <- do.call(rbind, lapply(SEEDS, per_seed))
saveRDS(results, out_rds)

# ---- Summary -------------------------------------------------------------

# Per-factor mean and worst RMSE across seeds; identify the F that
# minimises mean RMSE while never being more than 5 % worse than no-clamp
# on any seed.
inf_rmse <- subset(results, factor == Inf)$rmse
names(inf_rmse) <- as.character(subset(results, factor == Inf)$seed)

per_factor_summary <- do.call(rbind, lapply(FACTORS, function(f) {
  rows <- subset(results, factor == f)
  inf_per_seed <- inf_rmse[as.character(rows$seed)]
  worst_pct_change <- max(100 * (rows$rmse - inf_per_seed) / inf_per_seed,
                          na.rm = TRUE)
  data.frame(
    factor          = f,
    mean_rmse       = mean(rows$rmse, na.rm = TRUE),
    worst_pct_above_no_clamp = worst_pct_change,
    stringsAsFactors = FALSE
  )
}))

eligible <- subset(per_factor_summary, worst_pct_above_no_clamp <= 5)
best <- if (nrow(eligible) > 0L) {
  eligible[which.min(eligible$mean_rmse), "factor"]
} else {
  per_factor_summary[which.min(per_factor_summary$mean_rmse), "factor"]
}

# ---- Markdown report -----------------------------------------------------

md <- c(
  "# Phase G+ smoke bench: clamp_factor sweep on AVONET Mass",
  "",
  sprintf("Run: %s", format(Sys.time())),
  sprintf("AVONET n=%d, miss_frac=%.2f, seeds=%s, N_IMP=%d",
          N_SUB, MISS_FRAC, paste(SEEDS, collapse = ","), N_IMP),
  sprintf("Factors: %s", paste(FACTORS, collapse = ", ")),
  "",
  "## Per-seed Mass RMSE by factor",
  "",
  "```",
  capture.output(print(results, row.names = FALSE)),
  "```",
  "",
  "## Aggregated across seeds",
  "",
  "```",
  capture.output(print(per_factor_summary, row.names = FALSE)),
  "```",
  "",
  sprintf("Best factor (minimises mean RMSE while staying <= 5 %% above no-clamp on every seed): **%s**",
          format(best)),
  ""
)
writeLines(md, out_md)
log_line("DONE -- ", out_rds)
log_line("DONE -- ", out_md)
