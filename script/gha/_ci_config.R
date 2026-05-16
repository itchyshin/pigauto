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

# Per-imputation evaluation in user-data scale (not latent).
# Operates on the list `mi$datasets` returned by `multi_impute()`.
# Returns a long-format data.frame with one row per (trait, imputation_idx).
.trait_type_of <- function(trait_map, user_col) {
  if (is.null(trait_map)) return(NA_character_)
  for (tm in trait_map) {
    user_cols <- if (!is.null(tm$input_cols)) tm$input_cols else tm$name
    if (!is.null(user_cols) && user_col %in% user_cols) return(tm$type)
  }
  NA_character_
}

.eval_per_imputation <- function(completed_list, truth_df, mask, trait_map,
                                  t_fit_sec) {
  stopifnot(is.list(completed_list), is.data.frame(truth_df), is.matrix(mask))
  rows <- list()
  for (i in seq_along(completed_list)) {
    pred_df <- completed_list[[i]]
    for (j in colnames(mask)) {
      masked <- which(mask[, j])
      if (length(masked) == 0L) next
      t_v <- truth_df[[j]][masked]
      p_v <- pred_df[[j]][masked]
      type_j <- .trait_type_of(trait_map, j)
      if (is.factor(t_v) || is.character(t_v) || is.logical(t_v)) {
        ok  <- !is.na(t_v) & !is.na(p_v)
        acc <- if (any(ok)) mean(as.character(p_v[ok]) == as.character(t_v[ok])) else NA_real_
        rows[[length(rows) + 1L]] <- data.frame(
          trait = j, type = type_j, imputation_idx = i,
          rmse = NA_real_, mae = NA_real_, pearson_r = NA_real_,
          accuracy = acc, brier = NA_real_, time_sec = t_fit_sec,
          stringsAsFactors = FALSE
        )
      } else {
        t_num <- as.numeric(t_v); p_num <- as.numeric(p_v)
        ok <- !is.na(t_num) & !is.na(p_num)
        rmse <- if (any(ok)) sqrt(mean((t_num[ok] - p_num[ok])^2)) else NA_real_
        mae  <- if (any(ok)) mean(abs(t_num[ok] - p_num[ok])) else NA_real_
        pear <- if (sum(ok) > 2L) {
          suppressWarnings(stats::cor(t_num[ok], p_num[ok]))
        } else NA_real_
        rows[[length(rows) + 1L]] <- data.frame(
          trait = j, type = type_j, imputation_idx = i,
          rmse = rmse, mae = mae, pearson_r = pear,
          accuracy = NA_real_, brier = NA_real_, time_sec = t_fit_sec,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(
      trait = character(0), type = character(0), imputation_idx = integer(0),
      rmse = numeric(0), mae = numeric(0), pearson_r = numeric(0),
      accuracy = numeric(0), brier = numeric(0), time_sec = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}
