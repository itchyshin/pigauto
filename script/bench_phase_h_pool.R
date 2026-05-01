#!/usr/bin/env Rscript
#
# script/bench_phase_h_pool.R
#
# Phase H: discrete-trait pooling investigation.
#
# Pre-registered hypothesis (from useful/MEMO_2026-05-01_phase_f_smoke_results.md):
#   AVONET Migration accuracy regresses from -3.3 +/- 1.0 pp (N_IMP=1)
#   to -10.9 +/- 3.1 pp (N_IMP=20).  The baseline does not change with
#   N_IMP, so this cannot be baseline misspecification.  The current
#   ordinal pooling (mean of integer class indices, then round) likely
#   biases toward middle classes when dropout spreads predictions across
#   adjacent classes.
#
# Pre-registered acceptance criterion:
#   With pool_method = "mode" (per-cell majority vote), Migration acc
#   at N_IMP=20 >= Migration acc at N_IMP=1 across at least 2 of 3 seeds.
#
# Bench design:
#   3 seeds (2030, 2031, 2032)  x  N_IMP in {1, 5, 20}  x  pool_method in
#   {"median" (= current default for ordinal: mean+round),
#    "mode"   (= Phase H: per-cell majority vote)}
#   = 18 cells, ~3-9 min wall per cell, ~2 hr total.
#
# We use a SINGLE FIT per (seed, N_IMP) and re-pool the same M draws
# under both pool methods to keep the comparison apples-to-apples
# (no model retraining noise between pool methods).
#
# Output:
#   script/bench_phase_h_pool.{rds,md}
#   useful/MEMO_2026-05-01_phase_h_results.md  (verdict)

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

here     <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds  <- file.path(here, "script", "bench_phase_h_pool.rds")
out_md   <- file.path(here, "script", "bench_phase_h_pool.md")

SEEDS     <- c(2030L, 2031L, 2032L)
MISS_FRAC <- 0.30
N_SUB     <- 1500L
N_IMP_LIST <- c(1L, 20L)
POOLS      <- c("mode")     # mode-only; compare vs existing default-pool data:
#   N_IMP=20 default ("median"): Migration acc 0.713 +/- 0.032
#     (useful/MEMO_2026-05-01_multiseed_n20_and_default_flip.md)
#   N_IMP=1 default ("median"):  Migration acc 0.767 +/- 0.011
#     (useful/MEMO_2026-05-01_phase_f_smoke_results.md)

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - t0),
      ..., "\n", sep = "")
  flush.console()
}

# ---- Per-seed, per-N_IMP run -------------------------------------------------

run_one <- function(seed, n_imp) {
  log_line(sprintf("Seed %d N_IMP=%d: loading + masking ...", seed, n_imp))
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

  out <- list()
  for (pm in POOLS) {
    log_line(sprintf("Seed %d N_IMP=%d pool=%s: impute() ...",
                      seed, n_imp, pm))
    res <- pigauto::impute(df_miss, tree,
                            epochs = 200L,
                            n_imputations = n_imp,
                            pool_method = pm,
                            verbose = FALSE,
                            seed = seed)
    mig_idx <- which(mask_test[, "Migration"])
    truth   <- df_truth$Migration[mig_idx]
    pred    <- res$completed$Migration[mig_idx]
    acc     <- mean(as.character(pred) == as.character(truth), na.rm = TRUE)

    mode_class   <- names(sort(table(df_miss$Migration), decreasing = TRUE))[1]
    acc_baseline <- mean(as.character(truth) == mode_class, na.rm = TRUE)
    out[[pm]] <- data.frame(
      seed         = seed,
      n_imp        = n_imp,
      pool_method  = pm,
      pigauto_acc  = acc,
      baseline_acc = acc_baseline,
      lift_pp      = 100 * (acc - acc_baseline),
      n_test       = length(mig_idx),
      stringsAsFactors = FALSE
    )
    log_line(sprintf("  -> Migration pigauto = %.3f baseline = %.3f lift = %+.1f pp",
                      acc, acc_baseline, 100 * (acc - acc_baseline)))
  }
  do.call(rbind, out)
}

# ---- Sweep ----------------------------------------------------------------

results <- list()
i <- 0L
total <- length(SEEDS) * length(N_IMP_LIST)
for (s in SEEDS) {
  for (m in N_IMP_LIST) {
    i <- i + 1L
    log_line(sprintf("==== cell %d/%d: seed=%d N_IMP=%d ====",
                      i, total, s, m))
    results[[length(results) + 1L]] <- run_one(s, m)
  }
}
all_rows <- do.call(rbind, results)
saveRDS(all_rows, out_rds)

# ---- Summary -------------------------------------------------------------

agg <- stats::aggregate(
  cbind(pigauto_acc, lift_pp) ~ pool_method + n_imp,
  data = all_rows,
  FUN = function(x) c(mean = mean(x), sd = stats::sd(x))
)

# Acceptance: with pool_method = "mode", N_IMP=20 acc >= N_IMP=1 acc on
# at least 2 of 3 seeds.
mode_n1   <- subset(all_rows, pool_method == "mode" & n_imp == 1L)
mode_n20  <- subset(all_rows, pool_method == "mode" & n_imp == 20L)
mode_n1   <- mode_n1[order(mode_n1$seed), ]
mode_n20  <- mode_n20[order(mode_n20$seed), ]
n_pass    <- sum(mode_n20$pigauto_acc >= mode_n1$pigauto_acc, na.rm = TRUE)
verdict   <- if (n_pass >= 2L) {
  sprintf("PHASE H PASSES: with pool_method = 'mode', N_IMP=20 acc >= N_IMP=1 acc on %d / 3 seeds",
          n_pass)
} else {
  sprintf("PHASE H FAILS: with pool_method = 'mode', N_IMP=20 acc >= N_IMP=1 acc on only %d / 3 seeds",
          n_pass)
}

# ---- Markdown report -----------------------------------------------------

md <- c(
  "# Phase H smoke bench: ordinal pooling (mean+round vs mode)",
  "",
  sprintf("Run: %s", format(Sys.time())),
  sprintf("AVONET n=%d random subset, miss_frac=%.2f, seeds=%s, N_IMP=%s, pools=%s",
          N_SUB, MISS_FRAC,
          paste(SEEDS, collapse = ","),
          paste(N_IMP_LIST, collapse = ","),
          paste(POOLS, collapse = ",")),
  "",
  "## Per-seed Migration accuracy",
  "",
  "```",
  capture.output(print(all_rows[order(all_rows$pool_method, all_rows$seed,
                                       all_rows$n_imp),
                                 c("pool_method", "seed", "n_imp",
                                   "pigauto_acc", "baseline_acc", "lift_pp")],
                       row.names = FALSE)),
  "```",
  "",
  "## Aggregated (mean +/- SD over 3 seeds)",
  "",
  "```",
  capture.output(print(agg, row.names = FALSE)),
  "```",
  "",
  "## Verdict",
  "",
  verdict,
  ""
)
writeLines(md, out_md)
log_line("DONE -- ", out_rds)
log_line("DONE -- ", out_md)
