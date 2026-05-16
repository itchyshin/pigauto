# script/gha/_ci_config.R
#
# Shared CI config for the BACE head-to-head workflow.
# Sourced by every script/gha/run_bench_*.R wrapper.
#
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md.

PIGAUTO_CI_CONFIG <- list(
  subset_n          = 2000L,     # cap; smaller datasets use all available species
  n_imputations     = 10L,
  missing_frac      = 0.30,
  seed              = 2026L,
  pool_method       = "median",
  clamp_outliers    = FALSE,
  phylo_signal_gate = TRUE,
  mc_cores          = 4L
)

# Canonical schema for cross-method comparison.
# Pigauto's evaluate_imputation() output gets passed through here.
# BACE snapshots also conform to this schema (see script/gha/snapshot_bace.R).
.normalize_eval <- function(df, dataset, method) {
  stopifnot(is.data.frame(df), is.character(dataset), is.character(method))
  expected_cols <- c("trait", "type", "imputation_idx", "rmse", "mae",
                     "pearson_r", "accuracy", "brier", "time_sec")
  for (col in expected_cols) {
    if (!col %in% colnames(df)) df[[col]] <- NA
  }
  data.frame(
    dataset        = rep(dataset, nrow(df)),
    trait          = as.character(df$trait),
    type           = as.character(df$type),
    method         = rep(method, nrow(df)),
    imputation_idx = as.integer(df$imputation_idx),
    rmse           = as.numeric(df$rmse),
    mae            = as.numeric(df$mae),
    pearson_r      = as.numeric(df$pearson_r),
    accuracy       = as.numeric(df$accuracy),
    brier          = as.numeric(df$brier),
    time_sec       = as.numeric(df$time_sec),
    stringsAsFactors = FALSE
  )
}
