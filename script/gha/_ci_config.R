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

# Build a taxonomic phylo from columns like Order / Family / Genus / Species.
# Mirrors the pattern used in script/bench_amphibio.R and
# script/bench_globtherm_covariates.R (no molecular tree dependency).
.build_tax_tree <- function(tax_df, seed) {
  stopifnot(is.data.frame(tax_df),
            all(c("Order", "Family", "Genus", "Species") %in% colnames(tax_df)))
  tax_df[] <- lapply(tax_df, factor)
  tree <- ape::as.phylo(~Order/Family/Genus/Species, data = tax_df,
                         collapse = FALSE)
  set.seed(seed)
  tree <- ape::collapse.singles(tree)
  if (!ape::is.rooted(tree)) {
    tree <- ape::root.phylo(tree, outgroup = 1L, resolve.root = TRUE)
  }
  tree <- ape::multi2di(tree, random = TRUE)
  tree <- ape::compute.brlen(tree, method = "Grafen")
  tree$edge.length[tree$edge.length <= 0] <- 1e-8
  tree
}

# Read a file from a URL into a local cache; return the cached path.
# CI uses /tmp; local dev uses script/data-cache/.
.cache_path <- function(name) {
  root <- Sys.getenv("PIGAUTO_CI_CACHE", unset = file.path("script", "data-cache"))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  file.path(root, name)
}

.download_to_cache <- function(url, name, mode = "wb") {
  dst <- .cache_path(name)
  if (!file.exists(dst)) {
    utils::download.file(url, dst, mode = mode, quiet = TRUE)
  }
  dst
}

# End-to-end CI bench loop. Called by each script/gha/run_bench_<dataset>.R
# after that wrapper has populated (df, tree). Writes results.{rds,md} +
# timings.json under script/gha/results/<dataset>/.
.run_bench <- function(df, tree, dataset, out_dir, cfg,
                       trait_types = NULL,
                       multi_proportion_groups = NULL,
                       log_transform = TRUE) {
  stopifnot(is.data.frame(df), inherits(tree, "phylo"),
            is.character(dataset), is.list(cfg))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Subset to cap
  n_target <- min(cfg$subset_n, nrow(df))
  if (n_target < nrow(df)) {
    set.seed(cfg$seed)
    keep <- sort(sample.int(nrow(df), n_target))
    df   <- df[keep, , drop = FALSE]
    tree <- ape::keep.tip(tree, rownames(df))
  }
  cat(sprintf("[%s] n=%d species x %d traits\n",
              dataset, nrow(df), ncol(df)))

  # Preprocess + split
  t0_split <- Sys.time()
  pd0 <- preprocess_traits(df, tree, trait_types = trait_types,
                            multi_proportion_groups = multi_proportion_groups,
                            log_transform = log_transform)
  splits <- make_missing_splits(pd0$X_scaled,
                                 missing_frac = cfg$missing_frac,
                                 val_frac     = 0.5,
                                 seed         = cfg$seed,
                                 trait_map    = pd0$trait_map)
  mask_latent <- matrix(FALSE, nrow = nrow(pd0$X_scaled),
                         ncol = ncol(pd0$X_scaled))
  mask_latent[splits$test_idx] <- TRUE
  user_mask_test <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df),
                            dimnames = list(rownames(df), colnames(df)))
  for (k in seq_along(pd0$trait_map)) {
    tm <- pd0$trait_map[[k]]
    user_cols   <- if (!is.null(tm$input_cols)) tm$input_cols else tm$name
    latent_cols <- tm$latent_cols
    hit_rows <- apply(mask_latent[, latent_cols, drop = FALSE], 1L, any)
    for (uc in user_cols) {
      user_mask_test[, uc] <- user_mask_test[, uc] | hit_rows
    }
  }
  df_masked <- df
  df_masked[user_mask_test] <- NA
  t_split <- as.numeric(difftime(Sys.time(), t0_split, units = "secs"))

  # Multi-impute
  t0_fit <- Sys.time()
  mi <- multi_impute(
    df_masked, tree,
    m              = cfg$n_imputations,
    pool_method    = cfg$pool_method,
    clamp_outliers = cfg$clamp_outliers,
    trait_types    = trait_types,
    multi_proportion_groups = multi_proportion_groups,
    log_transform  = log_transform,
    seed           = cfg$seed
  )
  t_fit <- as.numeric(difftime(Sys.time(), t0_fit, units = "secs"))

  # Evaluate
  t0_eval <- Sys.time()
  ev <- .eval_per_imputation(mi$datasets, df, user_mask_test, pd0$trait_map,
                              t_fit)
  t_eval <- as.numeric(difftime(Sys.time(), t0_eval, units = "secs"))

  # Persist
  out_tbl <- .normalize_eval(ev, dataset = dataset, method = "pigauto_ci")
  saveRDS(out_tbl, file.path(out_dir, "results.rds"))
  jsonlite::write_json(
    list(split_sec = t_split, fit_sec = t_fit, eval_sec = t_eval,
         total_sec = t_split + t_fit + t_eval,
         n_species = nrow(df), n_imputations = cfg$n_imputations),
    file.path(out_dir, "timings.json"),
    auto_unbox = TRUE, pretty = TRUE)

  md_lines <- c(
    sprintf("# %s — pigauto CI bench", dataset),
    sprintf("Run config: seed=%d, missing_frac=%.2f, n_imputations=%d",
            cfg$seed, cfg$missing_frac, cfg$n_imputations),
    sprintf("N species used: %d", nrow(df)),
    sprintf("Wall time: %.1f s (fit) + %.1f s (eval)", t_fit, t_eval),
    "",
    "## Per-trait medians across imputations",
    ""
  )
  if (any(!is.na(out_tbl$rmse))) {
    agg_cont <- stats::aggregate(rmse ~ trait + type,
                                  data = out_tbl[!is.na(out_tbl$rmse), ],
                                  FUN = function(x) stats::median(x, na.rm = TRUE))
    md_lines <- c(md_lines, "### Continuous-family (RMSE)", "",
                  knitr::kable(agg_cont, format = "markdown", digits = 4), "")
  }
  if (any(!is.na(out_tbl$accuracy))) {
    agg_disc <- stats::aggregate(accuracy ~ trait + type,
                                  data = out_tbl[!is.na(out_tbl$accuracy), ],
                                  FUN = function(x) stats::median(x, na.rm = TRUE))
    md_lines <- c(md_lines, "### Discrete-family (accuracy)", "",
                  knitr::kable(agg_disc, format = "markdown", digits = 4), "")
  }
  writeLines(md_lines, file.path(out_dir, "results.md"))

  cat(sprintf("[%s] done in %.1f s total\n",
              dataset, t_split + t_fit + t_eval))
  invisible(out_tbl)
}

# Graceful failure marker: write a results.md/results.rds saying the
# dataset couldn't be loaded, so the aggregator job doesn't fall over.
.write_skip_marker <- function(dataset, out_dir, reason) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(.normalize_eval(
    data.frame(trait = character(0), type = character(0),
               imputation_idx = integer(0), rmse = numeric(0),
               mae = numeric(0), pearson_r = numeric(0),
               accuracy = numeric(0), brier = numeric(0),
               time_sec = numeric(0)),
    dataset = dataset, method = "pigauto_ci"
  ), file.path(out_dir, "results.rds"))
  writeLines(c(sprintf("# %s — SKIPPED", dataset), "",
               sprintf("Data load failed: %s", reason)),
             file.path(out_dir, "results.md"))
  jsonlite::write_json(list(error = reason),
                       file.path(out_dir, "timings.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  cat(sprintf("[%s] SKIPPED: %s\n", dataset, reason))
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
