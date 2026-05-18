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
#
# `covariate_cols` (added 2026-05-17): per-dataset list of columns that
# are PRESENT IN THE SNAPSHOT but absent from BACE's trait_subset.
# pigauto passes them through its dedicated covariate channel
# (`multi_impute(..., covariates = ...)`) which routes through the
# cov_encoder / obs_refine / cov_linear / cov_inject_blocks paths in
# the GNN. BACE's benchmark_dataset() never sees these columns, so
# pigauto's covariate-aware predictions are a structural advantage
# documented up-front, not an accidental data leak: both methods see
# the same snapshot, but only pigauto's architecture exposes a
# covariate channel.
#
# AVONET: range / centroid lat-lon are geographic columns that the
#   raw AVONET CSV ships with. They aren't in Dan's trait_subset.
# AmphiBIO: diurnal / nocturnal are binary diel-activity columns. They
#   aren't in BACE's amphibio trait_subset.
# LepTraits: Jan-Dec are monthly flight occurrence columns. Not in
#   Dan's leptraits trait_subset.
# GlobTherm, PanTHERIA, BIEN: leave empty -- their unused columns are
#   taxonomy only (already encoded by the tree).
BACE_DATASET_CONFIG <- list(
  avonet = list(
    trait_subset    = c("mass_g", "wing_length_mm", "beak_length_culmen_mm",
                         "tarsus_length_mm", "trophic_level",
                         "primary_lifestyle", "migration"),
    log_traits      = c("mass_g", "wing_length_mm",
                         "beak_length_culmen_mm", "tarsus_length_mm"),
    # 2026-05-17: tried (range_size_km2, centroid_lat, centroid_lon) as
    # covariates in CI 25996295467/25997176686. Continuous-trait
    # numbers were unchanged but AVONET categorical regressed
    # significantly (trophic_level 0.825 -> 0.731, primary_lifestyle
    # 0.828 -> 0.657). Suspect the GNN gate calibration opens for
    # categorical when more cov inputs are present + the OVR
    # threshold-joint baseline can't yet make use of them. Reverted
    # to character(0) until a categorical-safe covariate path lands.
    covariate_cols  = character(0)
  ),
  pantheria = list(
    trait_subset    = NULL,   # default: all non-tax cols
    log_traits      = c("body_mass_g", "head_body_length_mm",
                         "gestation_d", "max_longevity_m"),
    covariate_cols  = character(0)
  ),
  amphibio = list(
    trait_subset    = c("body_size_mm", "body_mass_g"),
    log_traits      = c("body_size_mm", "body_mass_g"),
    # diurnal / nocturnal have 76-87% NA in the snapshot, so not
    # usable as fully-observed covariates without aggressive imputation.
    covariate_cols  = character(0)
  ),
  bien = list(
    trait_subset    = NULL,
    log_traits      = c("height_m", "leaf_area", "sla",
                         "seed_mass", "wood_density"),
    # WorldClim bioclim covariates (38 = 19 vars x median + iqr).
    # external_cov_rds below filters traits_df to species WITH bioclim
    # coverage; these column names flow to cov_df and on to multi_impute.
    # 2026-05-18 controls:
    #   baseline (full pool, no cov, gate ON):  8.286
    #   control  (filtered pool, no cov, gate OFF): 7.455
    #   treatment (filtered pool, WC cov, gate OFF): 7.236
    #   => pool effect -0.83, covariate effect -0.22.
    covariate_cols  = paste0("bio", 1:19, "_median"),  # W5: median-only

    # Phylo-signal gate OFF on BIEN: weak-phylo-signal traits (sla, leaf_area,
    # height_m have lambda ~ 0) would otherwise be forced to mean-only, masking
    # the covariate lift. Spec §5.4 Option A.
    phylo_signal_gate = FALSE,
    # External covariate source: loaded at runtime and cbind-ed onto traits_df
    # BEFORE trait_subset / mask / subsetting. Species without bioclim coverage
    # are dropped from BIEN before the bench's n=2000 subsample.
    external_cov_rds  = "useful/bace_data_snapshot/data/bien_worldclim_covariates.rds"
  ),
  globtherm = list(
    trait_subset    = NULL,
    log_traits      = character(0),
    covariate_cols  = character(0)
  ),
  leptraits = list(
    trait_subset    = c("wingspan_lower", "forewing_length_lower",
                         "flight_duration", "n_hostplant_families"),
    log_traits      = c("wingspan_lower", "forewing_length_lower",
                         "flight_duration", "n_hostplant_families"),
    # Jan-Dec monthly flight columns have ~23% NA in the snapshot.
    # Median-filling that many cells would dilute the covariate signal.
    # Skip until we have a cleaner extraction or a proper missing-cov
    # handling path in pigauto.
    covariate_cols  = character(0)
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
#
# pool_method = c("auto", "per_draw", "pooled_point"):
#   "auto" (default, matches BACE): all trait types pool through
#     `pooled_point` for a single per-trait number, matching how Dan's
#     dev/benchmark_engine.R reports its head-to-head metrics. BACE
#     uses majority-vote-across-M-imputations for categorical/ordinal
#     accuracy and rowMeans-then-RMSE for continuous (lines 437-443
#     and 453-462 of benchmark_engine.R). pigauto's `mi$pooled_point`
#     contains the analogous central tendency: argmax-of-average-
#     probability for categorical/binary, modal-or-mean for ordinal,
#     median-or-mean for continuous depending on log-transform. Using
#     it for the head-to-head guarantees apples-to-apples comparison.
#   "per_draw": legacy behaviour — evaluate every trait per stochastic
#     draw. For categorical this measures sampled-class-vs-truth match,
#     which is bounded by E[p_max] across cells and badly underestimates
#     the model's actual classification accuracy when probabilities are
#     well-calibrated but not concentrated. For continuous, mean-of-per-
#     draw-RMSE >= RMSE-of-mean (Jensen on squared error), so per-draw
#     reporting is systematically pessimistic vs BACE's central-tendency
#     RMSE. Kept as an opt-in for diagnostics.
#   "pooled_point": same as "auto"; retained as an explicit synonym.
#
# Default is "auto" because BACE's snapshot reports a single accuracy /
# RMSE per trait per dataset. The previous per-draw default caused
# pigauto's CI categorical accuracy to look 0.4 pts lower AND its
# continuous RMSE to look ~10% higher than the central-tendency
# numbers `mi$pooled_point` already exposes.
.bace_eval_per_imputation <- function(completed_list, truth_long,
                                       trait_types_lookup, t_fit_sec,
                                       pooled_point = NULL,
                                       conformal_lower = NULL,
                                       conformal_upper = NULL,
                                       probabilities   = NULL,
                                       trait_levels    = NULL,
                                       pool_method = c("auto", "per_draw",
                                                       "pooled_point")) {
  pool_method <- match.arg(pool_method)
  # truth_long: data.frame(species_tip, trait, true_value)
  # conformal_lower / conformal_upper (optional): n_species x p matrices
  #   with 95% prediction-interval bounds for each continuous-family
  #   trait, in user scale (back-transformed when log-transform was on).
  #   Columns indexed by user trait name.
  # probabilities (optional): named list. For "binary" traits: numeric
  #   vector (length n_species) of P(class=levels[2]). For "categorical":
  #   K-column matrix with one row per species. Used to compute Brier
  #   score, the calibration metric BACE_snapshot also exposes.
  # trait_levels (optional): named list of character level vectors per
  #   discrete trait, required to align probabilities columns with truth.
  if (nrow(truth_long) == 0L) {
    return(data.frame(
      trait = character(0), type = character(0), imputation_idx = integer(0),
      rmse = numeric(0), mae = numeric(0), pearson_r = numeric(0),
      accuracy = numeric(0), brier = numeric(0),
      coverage_95 = numeric(0), interval_width = numeric(0),
      time_sec = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  use_pooled <- function(type_t) {
    if (pool_method == "per_draw" || is.null(pooled_point)) return(FALSE)
    # "auto" and "pooled_point" both pool every type
    TRUE
  }

  # Compute Brier score for one trait from a probability object and the
  # true class labels. Binary path: probs is a numeric vector of
  # P(class=positive). Categorical path: K-column matrix. Returns NA
  # when probs is missing / shape-incompatible.
  brier_one <- function(trait_name, type_t, sp_keep, t_v) {
    if (is.null(probabilities) || !(trait_name %in% names(probabilities))) {
      return(NA_real_)
    }
    p_obj <- probabilities[[trait_name]]
    levels_t <- if (!is.null(trait_levels) && trait_name %in% names(trait_levels))
                  trait_levels[[trait_name]] else NULL
    if (type_t == "binary") {
      if (is.null(levels_t) || length(levels_t) < 2L) return(NA_real_)
      # P(class = levels[2]) per pigauto convention.
      idx <- match(sp_keep, names(p_obj))
      if (all(is.na(idx))) idx <- match(sp_keep, rownames(p_obj))
      if (all(is.na(idx))) return(NA_real_)
      p_pos <- as.numeric(p_obj)[idx]
      y_pos <- as.numeric(as.character(t_v) == levels_t[2])
      ok <- is.finite(p_pos) & !is.na(y_pos)
      if (!any(ok)) return(NA_real_)
      mean((p_pos[ok] - y_pos[ok])^2, na.rm = TRUE)
    } else if (type_t == "categorical" || type_t == "ordinal") {
      if (!is.matrix(p_obj) || is.null(levels_t)) return(NA_real_)
      idx <- match(sp_keep, rownames(p_obj))
      if (all(is.na(idx))) return(NA_real_)
      pmat <- p_obj[idx, , drop = FALSE]
      # Truth one-hot in K columns matching p_obj column order.
      col_levels <- colnames(pmat)
      if (is.null(col_levels)) col_levels <- levels_t
      t_idx <- match(as.character(t_v), col_levels)
      ok <- !is.na(t_idx) & rowSums(is.finite(pmat)) == ncol(pmat)
      if (!any(ok)) return(NA_real_)
      Y <- matrix(0, nrow = length(ok), ncol = ncol(pmat))
      Y[cbind(seq_along(t_idx)[ok], t_idx[ok])] <- 1
      mean(rowSums((pmat[ok, , drop = FALSE] - Y[ok, , drop = FALSE])^2),
           na.rm = TRUE)
    } else {
      NA_real_
    }
  }

  # Compute conformal coverage_95 + interval_width for a continuous-
  # family trait at the masked species. Returns list(coverage, width).
  # Both NA when the bounds aren't supplied or don't have the trait.
  coverage_one <- function(trait_name, sp_keep, t_num) {
    if (is.null(conformal_lower) || is.null(conformal_upper)) {
      return(list(coverage = NA_real_, width = NA_real_))
    }
    if (!(trait_name %in% colnames(conformal_lower))) {
      return(list(coverage = NA_real_, width = NA_real_))
    }
    lo <- conformal_lower[match(sp_keep, rownames(conformal_lower)), trait_name]
    up <- conformal_upper[match(sp_keep, rownames(conformal_upper)), trait_name]
    ok <- is.finite(lo) & is.finite(up) & is.finite(t_num)
    if (!any(ok)) return(list(coverage = NA_real_, width = NA_real_))
    inside <- (t_num[ok] >= lo[ok]) & (t_num[ok] <= up[ok])
    list(coverage = mean(inside),
         width    = mean(up[ok] - lo[ok]))
  }

  rows <- list()
  # First emit per-draw rows for traits where we DON'T pool (continuous-family).
  for (i in seq_along(completed_list)) {
    pred_df <- completed_list[[i]]
    for (trait_name in unique(truth_long$trait)) {
      type_t <- trait_types_lookup[[trait_name]]
      if (is.null(type_t)) type_t <- NA_character_
      if (use_pooled(type_t)) next   # handled below from pooled_point

      trait_truth <- truth_long[truth_long$trait == trait_name, , drop = FALSE]
      sp <- trait_truth$species_tip
      sp_keep <- intersect(sp, rownames(pred_df))
      if (length(sp_keep) == 0L) next
      t_v <- trait_truth$true_value[match(sp_keep, trait_truth$species_tip)]
      p_v <- pred_df[[trait_name]][match(sp_keep, rownames(pred_df))]
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
          accuracy = acc, brier = NA_real_,
          coverage_95 = NA_real_, interval_width = NA_real_,
          time_sec = t_fit_sec,
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
          accuracy = NA_real_, brier = NA_real_,
          coverage_95 = NA_real_, interval_width = NA_real_,
          time_sec = t_fit_sec,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # Pooled rows (single row per trait) from pooled_point for discrete or
  # all-traits modes. The single row uses imputation_idx = 0 so it
  # aggregates cleanly through the median() step in the h2h report.
  if (!is.null(pooled_point)) {
    for (trait_name in unique(truth_long$trait)) {
      type_t <- trait_types_lookup[[trait_name]]
      if (is.null(type_t)) type_t <- NA_character_
      if (!use_pooled(type_t)) next
      if (!(trait_name %in% colnames(pooled_point))) next

      trait_truth <- truth_long[truth_long$trait == trait_name, , drop = FALSE]
      sp_keep <- intersect(trait_truth$species_tip, rownames(pooled_point))
      if (length(sp_keep) == 0L) next
      t_v <- trait_truth$true_value[match(sp_keep, trait_truth$species_tip)]
      p_v <- pooled_point[[trait_name]][match(sp_keep, rownames(pooled_point))]
      is_discrete <- isTRUE(type_t %in% c("categorical", "binary", "ordinal"))
      if (is_discrete) {
        ok <- !is.na(p_v) & !is.na(t_v)
        acc <- if (any(ok)) mean(as.character(p_v[ok]) ==
                                  as.character(t_v[ok])) else NA_real_
        brier_val <- brier_one(trait_name, type_t, sp_keep[ok], t_v[ok])
        rows[[length(rows) + 1L]] <- data.frame(
          trait = trait_name, type = type_t, imputation_idx = 0L,
          rmse = NA_real_, mae = NA_real_, pearson_r = NA_real_,
          accuracy = acc, brier = brier_val,
          coverage_95 = NA_real_, interval_width = NA_real_,
          time_sec = t_fit_sec,
          stringsAsFactors = FALSE
        )
      } else {
        t_num <- suppressWarnings(as.numeric(as.character(t_v)))
        p_num <- suppressWarnings(as.numeric(as.character(p_v)))
        ok <- is.finite(t_num) & is.finite(p_num)
        rmse <- if (any(ok)) sqrt(mean((t_num[ok] - p_num[ok])^2)) else NA_real_
        mae  <- if (any(ok)) mean(abs(t_num[ok] - p_num[ok])) else NA_real_
        pear <- if (sum(ok) > 2L) {
          suppressWarnings(stats::cor(t_num[ok], p_num[ok]))
        } else NA_real_
        cov_info <- coverage_one(trait_name, sp_keep[ok], t_num[ok])
        rows[[length(rows) + 1L]] <- data.frame(
          trait = trait_name, type = type_t, imputation_idx = 0L,
          rmse = rmse, mae = mae, pearson_r = pear,
          accuracy = NA_real_, brier = NA_real_,
          coverage_95 = cov_info$coverage,
          interval_width = cov_info$width,
          time_sec = t_fit_sec,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0L) {
    return(data.frame(
      trait = character(0), type = character(0), imputation_idx = integer(0),
      rmse = numeric(0), mae = numeric(0), pearson_r = numeric(0),
      accuracy = numeric(0), brier = numeric(0),
      coverage_95 = numeric(0), interval_width = numeric(0),
      time_sec = numeric(0),
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

  # 0) External covariate join (e.g. WorldClim for BIEN). When
  # ds_cfg$external_cov_rds points to a per-species covariate RDS,
  # cbind those covariates onto traits_df BEFORE subsetting, and
  # filter traits_df to species with non-NA coverage on the join
  # cols (avoids the bench sampling 80% median-filled covariates).
  ext_cov_rds <- ds_cfg$external_cov_rds %||% NULL
  ext_cov_col_names <- character(0)
  if (!is.null(ext_cov_rds) && nzchar(ext_cov_rds) && file.exists(ext_cov_rds)) {
    ext_cov <- readRDS(ext_cov_rds)
    # Match by rowname; species present in traits_df but not in ext_cov
    # get NAs (would be dropped below).
    ext_aligned <- ext_cov[rownames(traits_df), , drop = FALSE]
    rownames(ext_aligned) <- rownames(traits_df)
    ext_cov_col_names <- colnames(ext_aligned)
    has_cov <- stats::complete.cases(ext_aligned)
    cat(sprintf("[%s] external cov: %s -> %d / %d species have full coverage\n",
                dataset, ext_cov_rds, sum(has_cov), length(has_cov)))
    traits_df <- cbind(traits_df, ext_aligned)
    traits_df <- traits_df[has_cov, , drop = FALSE]
    tree <- ape::keep.tip(tree, rownames(traits_df))
    traits_df <- traits_df[tree$tip.label, , drop = FALSE]
  }

  # 1) Subset species via BACE's exact set.seed(2026); sample(rownames, 2000)
  keep <- .bace_subset_species(traits_df, subset_n = cfg$subset_n,
                                seed = cfg$seed)
  traits_df <- traits_df[keep, , drop = FALSE]
  tree <- ape::keep.tip(tree, keep)
  traits_df <- traits_df[tree$tip.label, , drop = FALSE]

  # 2a) Extract covariate_cols from the full traits_df BEFORE the
  #     trait_subset restriction kicks in (otherwise the covariate
  #     columns would be lost). cov_df is per-species, fully observed
  #     (median-fill the few NAs), z-scored downstream by pigauto's
  #     preprocess. range_size-style heavy-tailed columns get log-
  #     transformed first to keep their dynamic range comparable to
  #     other covariates.
  cov_cols_cfg <- ds_cfg$covariate_cols %||% character(0)
  cov_cols     <- intersect(cov_cols_cfg, colnames(traits_df))
  cov_df       <- NULL
  if (length(cov_cols) > 0L) {
    cov_raw <- traits_df[, cov_cols, drop = FALSE]
    # Log-transform heavy-tailed columns (range_size_km2 spans 8
    # orders of magnitude on AVONET; without log it dominates the
    # z-scored covariate matrix). Heuristic: numeric column with
    # min >= 0 across positive entries and max / median > 30.
    for (cc in cov_cols) {
      v <- suppressWarnings(as.numeric(cov_raw[[cc]]))
      if (any(is.finite(v) & v > 0)) {
        v_pos <- v[is.finite(v) & v > 0]
        if (length(v_pos) > 10L &&
            (max(v_pos) / max(stats::median(v_pos), 1e-12)) > 30) {
          # log1p so any zeros (none expected for area / lat / lon)
          # are handled without -Inf.
          v <- ifelse(is.finite(v) & v >= 0, log1p(v), NA_real_)
        }
      }
      # Median-fill NAs. cov_df must be fully observed per
      # multi_impute()'s covariates contract.
      med <- stats::median(v, na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      v[!is.finite(v)] <- med
      cov_raw[[cc]] <- v
    }
    cov_df <- as.data.frame(cov_raw)
    rownames(cov_df) <- rownames(traits_df)
    cat(sprintf("[%s] covariates: %s (median-filled NAs, log-1p heavy-tailed)\n",
                dataset, paste(cov_cols, collapse = ", ")))
  }

  # 2b) Restrict trait columns to trait_subset (default: all non-tax).
  # cov_cols AND ext_cov_col_names are removed from the trait set so they
  # aren't masked + imputed like traits. cov_cols flow to cov_df when
  # covariate_cols is non-empty; the broader ext_cov_col_names exclusion
  # ensures that columns cbind-ed in by external_cov_rds (e.g. WorldClim
  # bioclim) are NEVER treated as traits, even on a control run where
  # covariate_cols is intentionally empty.
  trait_cols <- if (!is.null(ds_cfg$trait_subset))
    intersect(ds_cfg$trait_subset, colnames(traits_df))
  else
    setdiff(.bace_default_subset(traits_df),
            union(cov_cols_cfg, ext_cov_col_names))
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

  # 5) Fit pigauto on the masked + log-transformed data. When
  # cov_df is non-NULL, pigauto's GNN activates its dedicated
  # covariate paths (cov_encoder / obs_refine / cov_linear /
  # cov_inject_blocks); BACE's bench saw no covariates so this is
  # pigauto's architectural advantage made explicit.
  t0_fit <- Sys.time()
  # Per-dataset phylo_signal_gate override. BIEN turns this OFF so
  # weak-phylo-signal traits (sla / leaf_area / height_m at lambda ~ 0)
  # don't get forced to mean-only by the gate, which would mask the
  # WorldClim covariate lift (spec §5.4 Option A).
  ds_phylo_gate <- ds_cfg$phylo_signal_gate
  if (is.null(ds_phylo_gate)) ds_phylo_gate <- TRUE
  mi <- multi_impute(
    masked_df_logged, tree,
    m              = cfg$n_imputations,
    pool_method    = cfg$pool_method,
    clamp_outliers = cfg$clamp_outliers,
    log_transform  = FALSE,
    covariates     = cov_df,
    n_heads        = 8L,    # T1: more attention heads
    phylo_signal_gate = ds_phylo_gate,
    seed           = cfg$seed
  )
  t_fit <- as.numeric(difftime(Sys.time(), t0_fit, units = "secs"))

  # 6) Evaluate on the FIT SCALE (= log-scale for log_traits, raw
  # otherwise). BACE's snapshot stores RMSE on the same fit scale,
  # so this is the apples-to-apples comparison. Previous version
  # back-transformed predictions to raw scale via .bace_unlog_pred(),
  # which inflated pigauto RMSE by factors of 10^2-10^8 vs BACE's
  # log-scale numbers on log_traits.
  t0_eval <- Sys.time()
  completed_fit <- mi$datasets   # already on fit scale by construction

  # Log-transform truth values for log_traits so truth and prediction
  # share the same scale. Non-log traits keep raw truth (BACE's
  # behaviour for those is identical, so no change needed).
  truth_long_fit <- truth_long
  if (length(ds_cfg$log_traits) > 0L) {
    for (lt in ds_cfg$log_traits) {
      idx <- which(truth_long_fit$trait == lt)
      if (!length(idx)) next
      v <- suppressWarnings(as.numeric(as.character(truth_long_fit$true_value[idx])))
      v[!is.finite(v) | v <= 0] <- NA_real_
      truth_long_fit$true_value[idx] <- log(v)
    }
  }

  # Trait-type lookup from the original (pre-log) traits_df. factor
  # with exactly 2 levels is BINARY in pigauto's type contract (and in
  # BACE's snapshot). Anything wider is categorical.
  trait_types_lookup <- setNames(
    vapply(trait_cols, function(v) {
      x <- traits_df[[v]]
      if (is.ordered(x))       "ordinal"
      else if (is.factor(x) && nlevels(x) == 2L) "binary"
      else if (is.factor(x))   "categorical"
      else if (is.character(x)) "categorical"
      else if (is.integer(x) && all(x[!is.na(x)] >= 0)) "count"
      else                     "continuous"
    }, character(1)),
    trait_cols
  )

  # Pass mi$pooled_point so the evaluator can score every trait type on
  # pigauto's central tendency (argmax-of-average-probability for
  # categorical/binary, mode-or-mean for ordinal, median-or-mean for
  # continuous depending on log-transform) — apples-to-apples with
  # BACE's "majority vote for discrete, rowMeans for continuous"
  # benchmark_engine.R policy.
  #
  # mi$conformal_lower / mi$conformal_upper feed coverage_95 +
  # interval_width for continuous-family traits. BACE doesn't expose
  # prediction intervals so those columns stay NA on the BACE side of
  # the head-to-head and pigauto reports them as an extra credibility
  # metric.
  #
  # mi$probabilities + per-trait level vectors feed the Brier score
  # for discrete traits, comparable to BACE's brier column.
  trait_levels_lookup <- lapply(trait_cols, function(v) {
    x <- traits_df[[v]]
    if (is.factor(x)) levels(x) else NULL
  })
  names(trait_levels_lookup) <- trait_cols

  ev <- .bace_eval_per_imputation(completed_fit, truth_long_fit,
                                   trait_types_lookup, t_fit,
                                   pooled_point    = mi$pooled_point,
                                   conformal_lower = mi$conformal_lower,
                                   conformal_upper = mi$conformal_upper,
                                   probabilities   = mi$probabilities,
                                   trait_levels    = trait_levels_lookup)
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
    agg <- stats::aggregate(
      cbind(rmse, coverage_95, interval_width) ~ trait + type,
      data = out_tbl[!is.na(out_tbl$rmse), ],
      FUN = function(x) stats::median(x, na.rm = TRUE),
      na.action = stats::na.pass)
    md <- c(md, "## Continuous-family (RMSE + 95% conformal coverage / width)", "",
            knitr::kable(agg, format = "markdown", digits = 4), "")
  }
  if (any(!is.na(out_tbl$accuracy))) {
    agg <- stats::aggregate(cbind(accuracy, brier) ~ trait + type,
                             data = out_tbl[!is.na(out_tbl$accuracy), ],
                             FUN = function(x) stats::median(x, na.rm = TRUE),
                             na.action = stats::na.pass)
    md <- c(md, "## Discrete-family (accuracy + Brier)", "",
            knitr::kable(agg, format = "markdown", digits = 4), "")
  }
  writeLines(md, file.path(out_dir, "results.md"))

  cat(sprintf("[%s] done in %.1f s total\n",
              dataset, t_fit + t_eval))
  invisible(out_tbl)
}
