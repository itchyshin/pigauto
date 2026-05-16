# script/gha/_bace_compat.R
#
# Port of daniel1noble/BACE's benchmark_engine.R subset + mask
# procedure, so pigauto's CI wrappers see the EXACT SAME masked
# data BACE saw. Without this we'd be running on a different DGP
# than BACE's snapshot results, and any head-to-head comparison
# is meaningless.
#
# Source of truth: daniel1noble/BACE/dev/benchmark_engine.R
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md.

# Per-dataset trait subsets + log-transformed traits, copied verbatim
# from BACE's dev/0[0-7]_benchmark_*.R scripts.
BACE_DATASET_CONFIG <- list(
  avonet = list(
    trait_subset = c("mass_g", "wing_length_mm", "beak_length_culmen_mm",
                      "tarsus_length_mm", "trophic_level",
                      "primary_lifestyle", "migration"),
    log_traits   = c("mass_g", "wing_length_mm",
                      "beak_length_culmen_mm", "tarsus_length_mm")
  ),
  pantheria = list(
    trait_subset = NULL,   # default: all non-tax cols
    log_traits   = c("body_mass_g", "head_body_length_mm",
                      "gestation_d", "max_longevity_m")
  ),
  amphibio = list(
    trait_subset = c("body_size_mm", "body_mass_g"),
    log_traits   = c("body_size_mm", "body_mass_g")
  ),
  bien = list(
    trait_subset = NULL,
    log_traits   = c("height_m", "leaf_area", "sla",
                      "seed_mass", "wood_density")
  ),
  globtherm = list(
    trait_subset = NULL,
    log_traits   = character(0)
  ),
  leptraits = list(
    trait_subset = c("wingspan_lower", "forewing_length_lower",
                      "flight_duration", "n_hostplant_families"),
    log_traits   = c("wingspan_lower", "forewing_length_lower",
                      "flight_duration", "n_hostplant_families")
  )
)

# Mirror BACE's tax-col drop default (BACE skips taxonomy cols when
# trait_subset is NULL).
.bace_default_subset <- function(traits_df) {
  tax_cols <- intersect(c("order_name", "family_name", "genus_name",
                          "class_name"), colnames(traits_df))
  setdiff(colnames(traits_df), tax_cols)
}

# Identical to BACE's set.seed + sample(rownames(traits_df), subset_n)
# subset step. Returns indices, not data, so the caller can re-use the
# same indices for the truth lookup.
.bace_subset_species <- function(traits_df, subset_n, seed) {
  set.seed(seed)
  n_full <- nrow(traits_df)
  if (is.na(subset_n) || subset_n >= n_full) {
    return(rownames(traits_df))
  }
  sample(rownames(traits_df), subset_n)
}

# Port of BACE's .apply_mask(): per-trait MCAR at cont_miss_rate or
# cat_miss_rate. Categoricals/ordinals use cat_miss_rate; continuous /
# count use cont_miss_rate. Returns (masked_df, truth) where truth is
# a long-format data.frame with one row per masked cell.
#
# Important: BACE applies the masks in column order with intervening
# sample() calls. To match BACE bit-for-bit, the caller must seed
# AGAIN before this (BACE re-seeds inside benchmark_dataset right
# after .apply_mask's sample calls move the RNG).
.bace_apply_mask <- function(traits_df, trait_cols,
                              cont_miss_rate = 0.30,
                              cat_miss_rate  = 0.30) {
  masked_df <- traits_df
  truth_rows <- list()
  for (v in trait_cols) {
    vals <- masked_df[[v]]
    is_cat <- is.factor(vals) || is.ordered(vals) || is.character(vals)
    rate <- if (is_cat) cat_miss_rate else cont_miss_rate
    obs <- which(!is.na(vals))
    if (length(obs) == 0L) next
    n_mask <- floor(length(obs) * rate)
    if (n_mask == 0L) next
    mask_idx <- sample(obs, n_mask)
    truth_rows[[length(truth_rows) + 1L]] <- data.frame(
      species_tip = rownames(masked_df)[mask_idx],
      trait       = v,
      true_value  = if (is_cat) as.character(vals[mask_idx])
                    else        as.numeric(vals[mask_idx]),
      stringsAsFactors = FALSE
    )
    masked_df[[v]][mask_idx] <- NA
  }
  truth <- if (length(truth_rows) > 0L)
    do.call(rbind, truth_rows)
  else
    data.frame(species_tip = character(0), trait = character(0),
               true_value  = character(0), stringsAsFactors = FALSE)
  list(masked = masked_df, truth = truth)
}

# Apply log on the per-dataset LOG_TRAITS. Replaces 0 and negative
# values with NA (mirrors BACE's behaviour: log(0) and log(-1) are
# disallowed; the affected cells become missing).
.bace_log_transform <- function(df, log_traits) {
  for (v in intersect(log_traits, colnames(df))) {
    x <- df[[v]]
    if (!is.numeric(x)) next
    bad <- !is.finite(x) | x <= 0
    if (any(bad, na.rm = TRUE)) x[bad] <- NA
    df[[v]] <- log(x)
  }
  df
}

# Reverse the log: predictions are on log scale, truth is on raw
# scale; we evaluate on raw scale via exp().
.bace_unlog_pred <- function(pred_df, log_traits) {
  for (v in intersect(log_traits, colnames(pred_df))) {
    if (is.numeric(pred_df[[v]])) {
      pred_df[[v]] <- exp(pred_df[[v]])
    }
  }
  pred_df
}

# Evaluate per imputation at the masked cells using BACE's truth
# long-format. Returns canonical-schema rows compatible with
# .normalize_eval().
.bace_eval_per_imputation <- function(completed_list, truth_long,
                                       trait_types_lookup, t_fit_sec) {
  # truth_long: data.frame(species_tip, trait, true_value)
  if (nrow(truth_long) == 0L) {
    return(data.frame(
      trait = character(0), type = character(0), imputation_idx = integer(0),
      rmse = numeric(0), mae = numeric(0), pearson_r = numeric(0),
      accuracy = numeric(0), brier = numeric(0), time_sec = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  rows <- list()
  for (i in seq_along(completed_list)) {
    pred_df <- completed_list[[i]]
    for (trait_name in unique(truth_long$trait)) {
      trait_truth <- truth_long[truth_long$trait == trait_name, , drop = FALSE]
      sp <- trait_truth$species_tip
      sp_keep <- intersect(sp, rownames(pred_df))
      if (length(sp_keep) == 0L) next
      t_v <- trait_truth$true_value[match(sp_keep, trait_truth$species_tip)]
      p_v <- pred_df[[trait_name]][match(sp_keep, rownames(pred_df))]
      type_t <- trait_types_lookup[[trait_name]]
      if (is.null(type_t)) type_t <- NA_character_
      # Branch on declared trait type (not on observed value class). When
      # .bace_apply_mask rbinds truth rows across mixed-type traits, the
      # true_value column gets coerced to character for ALL rows — so the
      # storage class of t_v / p_v can't be trusted; trust type_t.
      is_discrete <- isTRUE(type_t %in% c("categorical", "binary", "ordinal"))
      if (is_discrete) {
        ok <- !is.na(p_v) & !is.na(t_v)
        acc <- if (any(ok)) mean(as.character(p_v[ok]) ==
                                  as.character(t_v[ok])) else NA_real_
        rows[[length(rows) + 1L]] <- data.frame(
          trait = trait_name, type = type_t, imputation_idx = i,
          rmse = NA_real_, mae = NA_real_, pearson_r = NA_real_,
          accuracy = acc, brier = NA_real_, time_sec = t_fit_sec,
          stringsAsFactors = FALSE
        )
      } else {
        # Continuous / count / proportion / zi_count: coerce truth to
        # numeric even if storage is character (after rbind coercion).
        t_num <- suppressWarnings(as.numeric(as.character(t_v)))
        p_num <- suppressWarnings(as.numeric(as.character(p_v)))
        ok <- is.finite(t_num) & is.finite(p_num)
        rmse <- if (any(ok)) sqrt(mean((t_num[ok] - p_num[ok])^2)) else NA_real_
        mae  <- if (any(ok)) mean(abs(t_num[ok] - p_num[ok])) else NA_real_
        pear <- if (sum(ok) > 2L) {
          suppressWarnings(stats::cor(t_num[ok], p_num[ok]))
        } else NA_real_
        rows[[length(rows) + 1L]] <- data.frame(
          trait = trait_name, type = type_t, imputation_idx = i,
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

# Top-level BACE-compat bench. Mirrors BACE's benchmark_dataset()
# flow exactly: load, set.seed(seed), subset to subset_n, mask,
# log-transform, fit, evaluate. Output schema matches the rest of
# the gha harness.
.run_bench_bace_compat <- function(traits_df, tree, dataset, out_dir, cfg,
                                    log_transform = TRUE) {
  stopifnot(is.data.frame(traits_df), inherits(tree, "phylo"),
            is.character(dataset), is.list(cfg))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ds_cfg <- BACE_DATASET_CONFIG[[dataset]]
  if (is.null(ds_cfg)) stop(sprintf("No BACE config for dataset '%s'", dataset))

  # 1) Subset species via BACE's exact set.seed(2026); sample(rownames, 2000)
  keep <- .bace_subset_species(traits_df, subset_n = cfg$subset_n,
                                seed = cfg$seed)
  traits_df <- traits_df[keep, , drop = FALSE]
  tree <- ape::keep.tip(tree, keep)
  traits_df <- traits_df[tree$tip.label, , drop = FALSE]

  # 2) Restrict columns to trait_subset (default: all non-tax)
  trait_cols <- if (!is.null(ds_cfg$trait_subset))
    intersect(ds_cfg$trait_subset, colnames(traits_df))
  else
    .bace_default_subset(traits_df)
  traits_df <- traits_df[, trait_cols, drop = FALSE]
  cat(sprintf("[%s] subset to n=%d species x %d traits (%s)\n",
              dataset, nrow(traits_df), ncol(traits_df),
              paste(trait_cols, collapse = ",")))

  # 3) Apply BACE's exact .apply_mask() (continues from the same RNG
  # state set by .bace_subset_species; do NOT re-seed)
  masked_out <- .bace_apply_mask(traits_df, trait_cols,
                                  cont_miss_rate = cfg$missing_frac,
                                  cat_miss_rate  = cfg$missing_frac)
  masked_df <- masked_out$masked
  truth_long <- masked_out$truth
  cat(sprintf("[%s] %d masked cells across %d traits\n",
              dataset, nrow(truth_long),
              length(unique(truth_long$trait))))

  # 4) Log-transform the specified traits
  masked_df_logged <- .bace_log_transform(masked_df,
                                           log_traits = ds_cfg$log_traits)
  # Truth on raw scale (BACE convention)
  truth_long$true_value_chr <- truth_long$true_value
  # For numeric truth, coerce; for categoricals leave as character
  # (.bace_apply_mask already stored as character when factor).
  is_log <- truth_long$trait %in% ds_cfg$log_traits

  # 5) Fit pigauto on the masked + log-transformed data
  t0_fit <- Sys.time()
  mi <- multi_impute(
    masked_df_logged, tree,
    m              = cfg$n_imputations,
    pool_method    = cfg$pool_method,
    clamp_outliers = cfg$clamp_outliers,
    log_transform  = FALSE,   # we already log-transformed selected cols
    seed           = cfg$seed
  )
  t_fit <- as.numeric(difftime(Sys.time(), t0_fit, units = "secs"))

  # 6) Back-transform log preds to raw scale and evaluate
  t0_eval <- Sys.time()
  completed_raw <- lapply(mi$datasets, function(d) {
    .bace_unlog_pred(d, log_traits = ds_cfg$log_traits)
  })

  # Trait-type lookup from the original (pre-log) traits_df
  trait_types_lookup <- setNames(
    vapply(trait_cols, function(v) {
      x <- traits_df[[v]]
      if (is.ordered(x))       "ordinal"
      else if (is.factor(x))   "categorical"
      else if (is.character(x)) "categorical"
      else if (is.integer(x) && all(x[!is.na(x)] >= 0)) "count"
      else                     "continuous"
    }, character(1)),
    trait_cols
  )

  ev <- .bace_eval_per_imputation(completed_raw, truth_long,
                                   trait_types_lookup, t_fit)
  t_eval <- as.numeric(difftime(Sys.time(), t0_eval, units = "secs"))

  # 7) Normalise + persist
  out_tbl <- .normalize_eval(ev, dataset = dataset, method = "pigauto_ci")
  saveRDS(out_tbl, file.path(out_dir, "results.rds"))
  jsonlite::write_json(
    list(fit_sec = t_fit, eval_sec = t_eval,
         total_sec = t_fit + t_eval,
         n_species = nrow(traits_df),
         n_traits  = length(trait_cols),
         n_masked_cells = nrow(truth_long),
         n_imputations = cfg$n_imputations),
    file.path(out_dir, "timings.json"),
    auto_unbox = TRUE, pretty = TRUE
  )

  md <- c(
    sprintf("# %s — pigauto CI bench (BACE-compat)", dataset),
    sprintf("Seed=%d, subset_n=%d, missing_frac=%.2f, n_imputations=%d",
            cfg$seed, cfg$subset_n, cfg$missing_frac, cfg$n_imputations),
    sprintf("N species used: %d", nrow(traits_df)),
    sprintf("Traits: %s", paste(trait_cols, collapse = ", ")),
    sprintf("Log-transformed: %s",
            if (length(ds_cfg$log_traits) > 0L)
              paste(ds_cfg$log_traits, collapse = ", ") else "(none)"),
    sprintf("Masked cells: %d", nrow(truth_long)),
    sprintf("Wall time: %.1f s (fit) + %.1f s (eval)", t_fit, t_eval),
    ""
  )
  if (any(!is.na(out_tbl$rmse))) {
    agg <- stats::aggregate(rmse ~ trait + type,
                             data = out_tbl[!is.na(out_tbl$rmse), ],
                             FUN = function(x) stats::median(x, na.rm = TRUE))
    md <- c(md, "## Continuous-family (RMSE, raw scale)", "",
            knitr::kable(agg, format = "markdown", digits = 4), "")
  }
  if (any(!is.na(out_tbl$accuracy))) {
    agg <- stats::aggregate(accuracy ~ trait + type,
                             data = out_tbl[!is.na(out_tbl$accuracy), ],
                             FUN = function(x) stats::median(x, na.rm = TRUE))
    md <- c(md, "## Discrete-family (accuracy)", "",
            knitr::kable(agg, format = "markdown", digits = 4), "")
  }
  writeLines(md, file.path(out_dir, "results.md"))

  cat(sprintf("[%s] done in %.1f s total\n",
              dataset, t_fit + t_eval))
  invisible(out_tbl)
}
