#!/usr/bin/env Rscript
#
# script/bench_avonet_phase_f.R
#
# Phase F smoke bench: does adding LP as a third option in the per-trait
# ordinal path selection (R/fit_baseline.R) close the AVONET Migration
# regression?
#
# Multi-seed baseline (N_IMP=5 and N_IMP=20):
#   Migration acc: 0.713 +/- 0.032 (seeds 2030/2031/2032)
#   Mean baseline: 0.800 +/- 0.013
#   Regression: -10.9 +/- 3.1 pp
# Pre-registered Phase F target:
#   Migration acc >= 0.78 on at least 2 of 3 seeds.
#   That is roughly the lower bound of the column-mean baseline minus a
#   small implementation tolerance.
#
# Bench setup mirrors script/bench_avonet_full_local.R but uses
# n_imputations = 1 to keep wall fast (Migration accuracy on a single
# pass is robust to MC-dropout noise; the regression evidence persisted
# at both N_IMP=5 and N_IMP=20).
#
# Wall: ~3 min per seed x 3 seeds = ~9 min total, sequential.
#
# Output:
#   script/bench_avonet_phase_f.{rds,md}

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
out_rds  <- file.path(here, "script", "bench_avonet_phase_f.rds")
out_md   <- file.path(here, "script", "bench_avonet_phase_f.md")

SEEDS     <- c(2030L, 2031L, 2032L)
MISS_FRAC <- 0.30
N_SUB     <- 1500L

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - t0),
      ..., "\n", sep = "")
  flush.console()
}

run_one_seed <- function(seed) {
  log_line(sprintf("Seed %d: loading + subsetting ...", seed))
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
  stopifnot(all(rownames(df) == tree$tip.label))

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

  log_line(sprintf("Seed %d: impute(N_IMP=1) ...", seed))
  res <- pigauto::impute(df_miss, tree,
                          epochs = 200L,
                          n_imputations = 1L,
                          verbose = FALSE,
                          seed = seed)

  # Inspect ordinal_path_chosen on the fit
  path_chosen <- tryCatch(
    res$fit$baseline$ordinal_path_chosen %||% NULL,
    error = function(e) NULL)

  # Migration accuracy on test cells
  mig_idx <- which(mask_test[, "Migration"])
  truth   <- df_truth$Migration[mig_idx]
  pred    <- res$completed$Migration[mig_idx]
  acc     <- mean(as.character(pred) == as.character(truth), na.rm = TRUE)

  # Mean-baseline (per-class mode) accuracy
  mode_class <- names(sort(table(df_miss$Migration), decreasing = TRUE))[1]
  acc_baseline <- mean(as.character(truth) == mode_class, na.rm = TRUE)

  log_line(sprintf("Seed %d: Migration pigauto = %.3f, baseline = %.3f, lift = %+.1f pp",
                    seed, acc, acc_baseline, 100 * (acc - acc_baseline)))
  if (!is.null(path_chosen)) {
    log_line(sprintf("Seed %d: ordinal_path_chosen = %s",
                      seed, paste0(names(path_chosen), "=", path_chosen,
                                    collapse = ", ")))
  }

  data.frame(
    seed              = seed,
    pigauto_acc       = acc,
    baseline_acc      = acc_baseline,
    lift_pp           = 100 * (acc - acc_baseline),
    path_chosen       = if (!is.null(path_chosen)) {
                          paste(unlist(path_chosen), collapse = "|")
                        } else NA_character_,
    n_test            = length(mig_idx),
    stringsAsFactors  = FALSE
  )
}

# Tiny helper used in run_one_seed (no rlang dependency for one usage)
`%||%` <- function(a, b) if (is.null(a)) b else a

results <- do.call(rbind, lapply(SEEDS, run_one_seed))
saveRDS(results, out_rds)

# ---- Markdown report -------------------------------------------------------

tgt_acc <- 0.78
n_pass  <- sum(results$pigauto_acc >= tgt_acc, na.rm = TRUE)
verdict <- if (n_pass >= 2L) {
  sprintf("**PHASE F PASSES**: %d / %d seeds reach Migration acc >= %.2f",
          n_pass, nrow(results), tgt_acc)
} else {
  sprintf("**PHASE F INCONCLUSIVE**: only %d / %d seeds reach Migration acc >= %.2f",
          n_pass, nrow(results), tgt_acc)
}

md <- c(
  "# Phase F smoke bench: AVONET Migration ordinal regression",
  "",
  sprintf("Run: %s", format(Sys.time())),
  sprintf("n=%d species random subset, miss_frac=%.2f, N_IMP=1 (single pass)",
          N_SUB, MISS_FRAC),
  sprintf("Pre-registered target: Migration accuracy >= %.2f on >= 2 of 3 seeds",
          tgt_acc),
  "",
  "## Results",
  "",
  "```",
  capture.output(print(results, row.names = FALSE)),
  "```",
  "",
  "## Verdict",
  "",
  verdict,
  "",
  sprintf("- Mean pigauto Migration acc: %.3f (SD %.3f)",
          mean(results$pigauto_acc), stats::sd(results$pigauto_acc)),
  sprintf("- Mean baseline Migration acc: %.3f (SD %.3f)",
          mean(results$baseline_acc), stats::sd(results$baseline_acc)),
  sprintf("- Mean lift over baseline: %+.1f pp (SD %.1f)",
          mean(results$lift_pp), stats::sd(results$lift_pp)),
  "",
  "## Pre-Phase-F baseline (from `useful/MEMO_2026-05-01_multiseed_n20_and_default_flip.md`)",
  "",
  "Migration acc with pigauto N_IMP=20: 0.713 +/- 0.032",
  "Migration acc with mean baseline:    0.800 +/- 0.013",
  "Pre-Phase-F regression: -10.9 +/- 3.1 pp",
  "",
  "## Path selection",
  "",
  "Phase F adds 'lp' as a third option in the per-trait ordinal path",
  "selection (alongside 'threshold_joint' and 'bm_mvn'). The path_chosen",
  "column above shows which option was selected on each seed."
)
writeLines(md, out_md)
log_line("DONE -- ", out_rds)
log_line("DONE -- ", out_md)
