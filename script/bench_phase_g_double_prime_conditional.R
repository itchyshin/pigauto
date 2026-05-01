#!/usr/bin/env Rscript
#
# script/bench_phase_g_double_prime_conditional.R
#
# Phase G'' acceptance bench: pmm_when = "outside_observed" (default)
# vs the Phase G' "always" mode + reference numbers from PR #61.
#
# Pre-registered hypothesis (from useful/MEMO_2026-05-01_phase_g_prime_results.md):
#   The Phase G' (PMM-always) seed-2032 Mass regression of +805 %
#   was driven by donor-mismatch noise on accurate predictions.
#   Phase G'' default (`pmm_when = "outside_observed"`) only triggers
#   PMM on extrapolating predictions; in-range predictions are
#   trusted as-is.  Expected:
#     * Seed 2030 Casuarius win preserved (predictions extrapolate;
#       PMM still triggers).
#     * Seed 2032 Mass regression eliminated (predictions in-range;
#       PMM does NOT trigger; result = "none" baseline).
#     * Seeds 2031 + non-Mass traits: no harm.
#
# Pre-registered acceptance:
#   1. PMM-outside seed-2030 Mass RMSE <= PMM-always seed-2030 RMSE
#      (don't lose the Phase G' win).
#   2. PMM-outside doesn't regress > 5 % vs `none` on any cell
#      (the criterion that PMM-always failed).
#   3. PMM-outside imputed values stay in observed range when in-range
#      predictions are kept (verified at the test layer).
#
# Bench design:
#   3 AVONET seeds (2030/2031/2032) x 1 new config
#   ("pmm_when = outside_observed").  Reference numbers for "none",
#   "clamp", and "pmm_when = always" come from
#   script/bench_phase_g_prime_pmm.rds (PR #61).
#   Plus 3 synthetic heavy-tail reps for completeness.
#   Wall: ~30 min AVONET + ~3 min synthetic.
#
# Output: script/bench_phase_g_double_prime_conditional.{rds,md}

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
out_rds <- file.path(here, "script", "bench_phase_g_double_prime_conditional.rds")
out_md  <- file.path(here, "script", "bench_phase_g_double_prime_conditional.md")
ref_rds <- file.path(here, "script", "bench_phase_g_prime_pmm.rds")

SEEDS     <- c(2030L, 2031L, 2032L)
MISS_FRAC <- 0.30
N_SUB     <- 1500L
N_IMP     <- 20L
CONT_TRAITS <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - t0),
      ..., "\n", sep = "")
  flush.console()
}

run_avonet_one <- function(seed) {
  log_line(sprintf("AVONET seed %d (pmm + outside_observed) ...", seed))
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

  res <- pigauto::impute(df_miss, tree,
                          epochs = 200L, n_imputations = N_IMP,
                          match_observed = "pmm",
                          pmm_when = "outside_observed",
                          pmm_K = 5L,
                          verbose = FALSE, seed = seed)

  out <- vector("list", length(CONT_TRAITS))
  for (i in seq_along(CONT_TRAITS)) {
    tr <- CONT_TRAITS[i]
    idx <- which(mask_test[, tr])
    truth <- as.numeric(df_truth[[tr]][idx])
    pred  <- as.numeric(res$completed[[tr]][idx])
    rmse  <- sqrt(mean((truth - pred)^2, na.rm = TRUE))
    obs_max <- max(as.numeric(df_truth[[tr]]), na.rm = TRUE)
    obs_min <- min(as.numeric(df_truth[[tr]]), na.rm = TRUE)
    in_range <- mean(pred >= obs_min & pred <= obs_max, na.rm = TRUE)
    out[[i]] <- data.frame(
      dataset = "AVONET",
      seed    = seed,
      config  = "pmm_outside",
      trait   = tr,
      rmse    = rmse,
      pct_in_observed_range = 100 * in_range,
      max_pred = max(pred, na.rm = TRUE),
      obs_max  = obs_max,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

run_synthetic_one <- function(rep_id) {
  log_line(sprintf("Synthetic rep %d (pmm + outside_observed) ...", rep_id))
  set.seed(rep_id * 100L + 1L)
  n <- 200L
  tree <- ape::rtree(n)
  log_mass <- stats::rnorm(n, 7, 0.5)
  outlier_idx <- sample.int(n, ceiling(0.05 * n))
  log_mass[outlier_idx] <- log_mass[outlier_idx] + stats::rnorm(length(outlier_idx),
                                                                 mean = 2,
                                                                 sd = 0.5)
  df <- data.frame(mass = exp(log_mass), row.names = tree$tip.label)
  df_truth <- df
  obs_max <- max(df$mass)
  obs_min <- min(df$mass)
  set.seed(rep_id * 100L + 2L)
  mask_idx <- sample.int(n, ceiling(MISS_FRAC * n))
  df$mass[mask_idx] <- NA

  res <- pigauto::impute(df, tree,
                          epochs = 100L, n_imputations = 10L,
                          match_observed = "pmm",
                          pmm_when = "outside_observed",
                          pmm_K = 5L,
                          verbose = FALSE,
                          seed = rep_id * 100L + 3L)
  truth <- df_truth$mass[mask_idx]
  pred  <- as.numeric(res$completed$mass[mask_idx])
  rmse  <- sqrt(mean((truth - pred)^2, na.rm = TRUE))
  in_range <- mean(pred >= obs_min & pred <= obs_max, na.rm = TRUE)
  data.frame(
    dataset = "synthetic_heavy_tail",
    seed    = rep_id,
    config  = "pmm_outside",
    trait   = "mass",
    rmse    = rmse,
    pct_in_observed_range = 100 * in_range,
    max_pred = max(pred, na.rm = TRUE),
    obs_max  = obs_max,
    stringsAsFactors = FALSE
  )
}

# ---- New runs --------------------------------------------------------------

new_results <- list()
for (s in SEEDS) {
  new_results[[length(new_results) + 1L]] <- run_avonet_one(s)
}
for (rp in 1:3) {
  new_results[[length(new_results) + 1L]] <- run_synthetic_one(rp)
}
new_rows <- do.call(rbind, new_results)

# ---- Combine with reference (PR #61) ---------------------------------------

if (file.exists(ref_rds)) {
  ref_rows <- readRDS(ref_rds)
  # Rename ref's "pmm" config -> "pmm_always" for clarity vs new "pmm_outside"
  ref_rows$config[ref_rows$config == "pmm"] <- "pmm_always"
  all_rows <- rbind(ref_rows, new_rows)
} else {
  all_rows <- new_rows
}
saveRDS(all_rows, out_rds)

# ---- Acceptance check ------------------------------------------------------

mass_2030 <- subset(all_rows, dataset == "AVONET" & trait == "Mass" &
                              seed == 2030L)
crit_1 <- mass_2030$rmse[mass_2030$config == "pmm_outside"] <=
          mass_2030$rmse[mass_2030$config == "pmm_always"]

# Criterion 2: pmm_outside doesn't regress > 5 % vs none on any cell
none_rows <- subset(all_rows, config == "none")
out_rows  <- subset(all_rows, config == "pmm_outside")
key  <- paste(none_rows$dataset, none_rows$seed, none_rows$trait, sep = "/")
okey <- paste(out_rows$dataset,  out_rows$seed,  out_rows$trait,  sep = "/")
ord  <- match(key, okey)
pct_change <- 100 * (out_rows$rmse[ord] - none_rows$rmse) / none_rows$rmse
crit_2_violators <- which(pct_change > 5)
crit_2 <- length(crit_2_violators) == 0L

verdict <- if (crit_1 && crit_2) {
  "**PHASE G'' PASSES** -- conditional PMM preserves seed-2030 Casuarius win AND eliminates the Phase G' regressions on accurate-prediction cells."
} else {
  parts <- c(
    if (!crit_1) sprintf("FAIL 1: pmm_outside seed-2030 Mass %.0f > pmm_always %.0f",
                          mass_2030$rmse[mass_2030$config == "pmm_outside"],
                          mass_2030$rmse[mass_2030$config == "pmm_always"]),
    if (!crit_2) sprintf("FAIL 2: pmm_outside regresses > 5 %% on %d cells (worst delta = %.0f%%)",
                          length(crit_2_violators), max(pct_change, na.rm = TRUE))
  )
  paste("**PHASE G'' PARTIAL** --", paste(parts, collapse = "; "))
}

# ---- Markdown report -------------------------------------------------------

md <- c(
  "# Phase G'' acceptance bench: conditional PMM (pmm_when = 'outside_observed')",
  "",
  sprintf("Run: %s", format(Sys.time())),
  sprintf("AVONET n=%d, miss=%.2f, seeds=%s, N_IMP=%d.  Synthetic heavy-tail: 3 reps.",
          N_SUB, MISS_FRAC, paste(SEEDS, collapse = ","), N_IMP),
  "Reference numbers for `none`, `clamp`, `pmm_always` come from PR #61's bench (script/bench_phase_g_prime_pmm.rds).",
  "",
  "## Pre-registered acceptance criteria",
  "",
  "1. pmm_outside seed-2030 Mass RMSE <= pmm_always seed-2030 RMSE",
  "   (don't lose the Phase G' win).",
  "2. pmm_outside doesn't regress > 5 % vs `none` on any cell",
  "   (the criterion that pmm_always failed).",
  "",
  "## Per-cell results across all configs",
  "",
  "```",
  capture.output(print(all_rows[order(all_rows$dataset, all_rows$trait,
                                        all_rows$seed, all_rows$config),
                                  c("dataset", "trait", "seed", "config",
                                    "rmse", "pct_in_observed_range",
                                    "max_pred", "obs_max")],
                        row.names = FALSE)),
  "```",
  "",
  sprintf("Crit 1 (pmm_outside seed-2030 Mass <= pmm_always): %s",
          if (crit_1) "PASS" else "FAIL"),
  sprintf("Crit 2 (pmm_outside no-regress vs none, threshold 5%%): %s",
          if (crit_2) "PASS"
          else sprintf("FAIL on %d cells (worst %.1f%%)",
                       length(crit_2_violators),
                       max(pct_change, na.rm = TRUE))),
  "",
  "## Verdict",
  "",
  verdict,
  ""
)
writeLines(md, out_md)
log_line("DONE -- ", out_rds)
log_line("DONE -- ", out_md)
