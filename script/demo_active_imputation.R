#!/usr/bin/env Rscript
#
# script/demo_active_imputation.R
#
# Demonstration / validation of suggest_next_observation():
# does following its top-N recommendation actually reduce
# imputation error more than random selection?
#
# Setup:
#   * AVONET 300 (4 continuous traits: Mass, Beak.Length_Culmen,
#     Tarsus.Length, Wing.Length)
#   * 30 % MCAR mask -> baseline imputation + variance estimate
#   * Compare three strategies for "next 10 cells to observe":
#      (a) ACTIVE: top-10 by suggest_next_observation() (by = "cell")
#      (b) RANDOM: 10 cells drawn uniformly from missing cells
#      (c) HIGH_SE: 10 cells with the highest individual posterior SE
#         (the naive heuristic)
#   * For each strategy, "observe" those 10 cells (use the truth
#     from the masked-out values), refit, evaluate test-set RMSE on
#     the OTHER missing cells.
#
# The hypothesis: the closed-form expected variance reduction
# correlates with actual realised improvement, so ACTIVE strategy
# should reduce test RMSE more than RANDOM or HIGH_SE.
#
# Run: Rscript script/demo_active_imputation.R
# Output: script/demo_active_imputation.md (table)

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  here_path <- "/Users/z3437171/Dropbox/Github Local/pigauto"
  devtools::load_all(here_path, quiet = TRUE)
})

out_md  <- file.path(here_path, "script", "demo_active_imputation.md")

script_start <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start),
      ..., "\n", sep = "")
  flush.console()
}

# Load AVONET 300
data("avonet300", package = "pigauto")
data("tree300",   package = "pigauto")
df <- avonet300
rownames(df) <- df$Species_Key
df$Species_Key <- NULL
cont_cols <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")

# Mask 30% of cells across the 4 continuous traits, deterministically
set.seed(2026)
mask <- matrix(FALSE, nrow = nrow(df), ncol = length(cont_cols),
               dimnames = list(NULL, cont_cols))
df_truth <- df
for (v in cont_cols) {
  ok <- which(!is.na(df_truth[[v]]))
  n_mask <- round(0.30 * length(ok))
  idx <- sample(ok, n_mask)
  mask[idx, v] <- TRUE
  df[[v]][idx] <- NA_real_
}

n_miss_total <- sum(mask)
log_line("Total masked cells: ", n_miss_total)

# ---- Step 1: baseline impute on the masked dataset ----------------------
log_line("Baseline impute (drops 30% of cells, no oracle)...")
res0 <- pigauto::impute(df, tree300, epochs = 60L, n_imputations = 1L,
                          verbose = FALSE, seed = 2026L)

baseline_rmse <- numeric(length(cont_cols)); names(baseline_rmse) <- cont_cols
for (v in cont_cols) {
  truth_v <- df_truth[[v]][mask[, v]]
  pred_v  <- res0$completed[[v]][mask[, v]]
  baseline_rmse[v] <- sqrt(mean((truth_v - pred_v)^2, na.rm = TRUE))
}
log_line("Baseline RMSE per trait:")
log_line(paste(sprintf("%s = %.3g", cont_cols, baseline_rmse), collapse = " | "))

# ---- Step 2: collect ACTIVE / RANDOM / HIGH_SE candidate cells ---------
log_line("Collecting top-10 ACTIVE candidates...")
active_top <- pigauto::suggest_next_observation(res0, top_n = 10, by = "cell")
log_line("Active top-3 cells:")
print(utils::head(active_top, 3))

# RANDOM: sample 10 cells uniformly from masked cells
miss_lin <- which(mask)
set.seed(2027)
random_idx <- sample(miss_lin, 10)
# Convert linear -> (row, col)
n <- nrow(mask); p <- ncol(mask)
random_rows <- ((random_idx - 1L) %% n) + 1L
random_cols <- ((random_idx - 1L) %/% n) + 1L
random_cells <- data.frame(species = rownames(df_truth)[random_rows],
                            trait = colnames(mask)[random_cols],
                            stringsAsFactors = FALSE)

# HIGH_SE: cells with highest predicted individual SE
pred0 <- res0$prediction
se_mat <- pred0$se   # may be on user-scale; same monotone for ranking
# Get SE at masked cells, rank descending
high_se_records <- list()
for (v in cont_cols) {
  rows <- which(mask[, v])
  for (r in rows) {
    se_v <- if (is.matrix(se_mat)) se_mat[r, v] else NA_real_
    if (!is.finite(se_v)) next
    high_se_records[[length(high_se_records) + 1L]] <- data.frame(
      species = rownames(df_truth)[r],
      trait = v,
      se = as.numeric(se_v),
      stringsAsFactors = FALSE)
  }
}
high_se_df <- do.call(rbind, high_se_records)
high_se_top <- utils::head(high_se_df[order(-high_se_df$se), ], 10)

# ---- Step 3: simulate observing each strategy's top-10 ----------------
simulate_observe <- function(strategy_cells) {
  df_new <- df  # current masked
  remaining_mask <- mask
  rownames(remaining_mask) <- rownames(df_truth)
  for (i in seq_len(nrow(strategy_cells))) {
    sp <- strategy_cells$species[i]; tr <- strategy_cells$trait[i]
    df_new[sp, tr] <- df_truth[sp, tr]
    remaining_mask[sp, tr] <- FALSE
  }
  res_new <- pigauto::impute(df_new, tree300, epochs = 60L,
                                n_imputations = 1L, verbose = FALSE,
                                seed = 2026L)
  per_trait <- numeric(length(cont_cols)); names(per_trait) <- cont_cols
  for (v in cont_cols) {
    if (!any(remaining_mask[, v])) {
      per_trait[v] <- NA_real_
    } else {
      truth_v <- df_truth[[v]][remaining_mask[, v]]
      pred_v  <- res_new$completed[[v]][remaining_mask[, v]]
      per_trait[v] <- sqrt(mean((truth_v - pred_v)^2, na.rm = TRUE))
    }
  }
  per_trait
}

log_line("Simulating ACTIVE observation...")
rmse_active <- simulate_observe(active_top)
log_line("Simulating RANDOM observation...")
rmse_random <- simulate_observe(random_cells)
log_line("Simulating HIGH_SE observation...")
rmse_high_se <- simulate_observe(high_se_top)

# ---- Step 4: report -----------------------------------------------------
md <- c(
  "# Active imputation demo: AVONET 300, 4 continuous traits",
  "",
  sprintf("Run: %s", format(Sys.time())),
  sprintf("n=%d species, mask=30%%, %d masked cells; observing 10 of them",
          nrow(df_truth), n_miss_total),
  "",
  "## Test-set RMSE on REMAINING masked cells after observing 10 cells",
  "",
  "Lower RMSE = better strategy.  Baseline = RMSE on ALL masked cells",
  "(no observation) for reference.",
  "",
  sprintf("| Trait | Baseline (no obs) | RANDOM | HIGH_SE | ACTIVE |"),
  sprintf("|---|---|---|---|---|"),
  paste(
    "| ",
    cont_cols,
    " | ",
    sprintf("%.3g", baseline_rmse),
    " | ",
    sprintf("%.3g", rmse_random),
    " | ",
    sprintf("%.3g", rmse_high_se),
    " | ",
    sprintf("%.3g", rmse_active),
    " |",
    sep = "", collapse = "\n"
  ),
  "",
  "## Interpretation",
  "",
  sprintf("- **Mean RMSE across traits**: baseline=%.3g, RANDOM=%.3g, HIGH_SE=%.3g, ACTIVE=%.3g",
          mean(baseline_rmse), mean(rmse_random, na.rm = TRUE),
          mean(rmse_high_se, na.rm = TRUE), mean(rmse_active, na.rm = TRUE)),
  "",
  "If suggest_next_observation()'s closed-form variance-reduction",
  "formula correlates with actual realised improvement, the ACTIVE",
  "strategy should produce lower test RMSE than RANDOM and",
  "(usually) than HIGH_SE.  The HIGH_SE heuristic is a sensible-",
  "looking but suboptimal strategy: cells with high individual SE",
  "may not have high TOTAL variance reduction (their SE is high",
  "because they are far from observed species, but observing them",
  "doesn't necessarily inform many other cells).",
  "",
  "Top-3 ACTIVE recommendations (highest variance reduction):",
  "",
  "```",
  capture.output(print(utils::head(active_top, 3), row.names = FALSE)),
  "```"
)
writeLines(md, out_md)
log_line("DONE -- ", out_md)
