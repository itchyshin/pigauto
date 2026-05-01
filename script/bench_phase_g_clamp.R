#!/usr/bin/env Rscript
#
# script/bench_phase_g_clamp.R
#
# Phase G acceptance bench: does opt-in clamp_outliers close the
# AVONET Mass tail-extrapolation mode?
#
# Pre-registered acceptance criterion (from PR #53 memo):
#   - Seed-2030 Mass RMSE drops below 5,000 (currently 23,723 at
#     N_IMP=20 without clamp; the diagnostic showed Casuarius
#     bennetti predicted at 538,166 g vs truth 35,000 g).
#   - No regression on seeds 2031, 2032 (Mass RMSE change < 5%).
#
# Bench: 3 seeds (2030 / 2031 / 2032)  x  clamp in {OFF, ON} at
# clamp_factor = 5.  N_IMP = 20 (the regime where the blow-up is
# largest).  Wall: ~9 min per impute() call x 6 = ~54 min.
#
# Output: script/bench_phase_g_clamp.{rds,md}

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
out_rds <- file.path(here, "script", "bench_phase_g_clamp.rds")
out_md  <- file.path(here, "script", "bench_phase_g_clamp.md")

SEEDS        <- c(2030L, 2031L, 2032L)
MISS_FRAC    <- 0.30
N_SUB        <- 1500L
N_IMP        <- 20L
CLAMP_FACTOR <- 5

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - t0),
      ..., "\n", sep = "")
  flush.console()
}

run_one <- function(seed, clamp_on) {
  log_line(sprintf("Seed %d clamp=%s: loading + masking ...",
                   seed, if (clamp_on) "ON" else "OFF"))
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

  log_line(sprintf("Seed %d clamp=%s: impute() N_IMP=%d ...",
                   seed, if (clamp_on) "ON" else "OFF", N_IMP))
  res <- pigauto::impute(df_miss, tree,
                          epochs = 200L,
                          n_imputations = N_IMP,
                          clamp_outliers = clamp_on,
                          clamp_factor   = CLAMP_FACTOR,
                          verbose = FALSE,
                          seed = seed)
  mass_idx   <- which(mask_test[, "Mass"])
  mass_truth <- df_truth$Mass[mass_idx]
  mass_pred  <- as.numeric(res$completed$Mass[mass_idx])
  rmse <- sqrt(mean((mass_truth - mass_pred)^2, na.rm = TRUE))

  baseline_rmse <- sqrt(mean(
    (mass_truth - mean(df_miss$Mass, na.rm = TRUE))^2, na.rm = TRUE))

  # Also report top-5 worst residuals to confirm the Casuarius-class
  # outliers are catched
  resid_abs <- abs(mass_truth - mass_pred)
  ord <- order(resid_abs, decreasing = TRUE)
  top5_residual_max <- if (length(ord) >= 1L) mass_pred[ord[1L]] else NA_real_

  log_line(sprintf("Seed %d clamp=%s: Mass RMSE = %.0f (baseline %.0f); top-1 pred = %.0f",
                   seed, if (clamp_on) "ON" else "OFF",
                   rmse, baseline_rmse, top5_residual_max))

  data.frame(
    seed         = seed,
    clamp        = if (clamp_on) "ON" else "OFF",
    rmse         = rmse,
    baseline     = baseline_rmse,
    top1_pred    = top5_residual_max,
    n_test       = length(mass_idx),
    stringsAsFactors = FALSE
  )
}

# ---- Sweep ---------------------------------------------------------------

results <- list()
for (s in SEEDS) {
  for (cl in c(FALSE, TRUE)) {
    results[[length(results) + 1L]] <- run_one(s, cl)
  }
}
all_rows <- do.call(rbind, results)
saveRDS(all_rows, out_rds)

# ---- Acceptance check ----------------------------------------------------

off  <- subset(all_rows, clamp == "OFF")
on   <- subset(all_rows, clamp == "ON")
off  <- off[order(off$seed), ]
on   <- on[order(on$seed), ]

# Criterion 1: seed-2030 with clamp ON, RMSE < 5000
crit1 <- on$rmse[on$seed == 2030L] < 5000
# Criterion 2: seeds 2031/2032, RMSE change with clamp ON < 5% of OFF RMSE
delta_2031 <- abs(on$rmse[on$seed == 2031L] - off$rmse[off$seed == 2031L]) /
              off$rmse[off$seed == 2031L]
delta_2032 <- abs(on$rmse[on$seed == 2032L] - off$rmse[off$seed == 2032L]) /
              off$rmse[off$seed == 2032L]
crit2 <- delta_2031 < 0.05 && delta_2032 < 0.05

verdict <- if (crit1 && crit2) {
  "**PHASE G PASSES**"
} else if (crit1 && !crit2) {
  sprintf("**PHASE G PARTIAL**: seed-2030 closure OK, but Δ on 2031/2032 = (%.1f%%, %.1f%%) > 5%%",
          100 * delta_2031, 100 * delta_2032)
} else if (!crit1) {
  sprintf("**PHASE G FAILS criterion 1**: seed-2030 RMSE = %.0f >= 5000",
          on$rmse[on$seed == 2030L])
} else {
  "**PHASE G FAILS** (both criteria)"
}

# ---- Markdown report -----------------------------------------------------

md <- c(
  "# Phase G acceptance bench: clamp_outliers on AVONET Mass",
  "",
  sprintf("Run: %s", format(Sys.time())),
  sprintf("AVONET n=%d, miss_frac=%.2f, seeds=%s, N_IMP=%d, clamp_factor=%g",
          N_SUB, MISS_FRAC,
          paste(SEEDS, collapse = ","), N_IMP, CLAMP_FACTOR),
  "",
  "## Pre-registered acceptance criteria",
  "",
  "1. Seed-2030 Mass RMSE with clamp ON drops below 5,000 (currently 23,723).",
  "2. Seeds 2031 and 2032 Mass RMSE change with clamp ON < 5 % of OFF RMSE.",
  "",
  "## Per-seed results",
  "",
  "```",
  capture.output(print(all_rows[order(all_rows$seed, all_rows$clamp),
                                 c("seed", "clamp", "rmse",
                                   "baseline", "top1_pred")],
                       row.names = FALSE)),
  "```",
  "",
  sprintf("Seed-2030 RMSE drop: %.0f -> %.0f (%.1f%% reduction)",
          off$rmse[off$seed == 2030L], on$rmse[on$seed == 2030L],
          100 * (off$rmse[off$seed == 2030L] - on$rmse[on$seed == 2030L]) /
            off$rmse[off$seed == 2030L]),
  sprintf("Seed-2031 RMSE Δ: %.0f -> %.0f (%.2f%%)",
          off$rmse[off$seed == 2031L], on$rmse[on$seed == 2031L],
          100 * delta_2031),
  sprintf("Seed-2032 RMSE Δ: %.0f -> %.0f (%.2f%%)",
          off$rmse[off$seed == 2032L], on$rmse[on$seed == 2032L],
          100 * delta_2032),
  "",
  "## Verdict",
  "",
  verdict,
  ""
)
writeLines(md, out_md)
log_line("DONE -- ", out_rds)
log_line("DONE -- ", out_md)
