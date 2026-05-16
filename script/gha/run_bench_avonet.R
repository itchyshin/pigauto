#!/usr/bin/env Rscript
# script/gha/run_bench_avonet.R
#
# CI wrapper for the AVONET head-to-head bench.
# Loads pigauto via devtools::load_all() (CI doesn't install),
# runs the documented-default config, writes a normalized
# .rds + a short .md + timings.json under script/gha/results/avonet/.
#
# Spec: specs/2026-05-16-bace-headtohead-ci-design.md.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  here <- "."  # workflow runs from repo root
  source(file.path(here, "script", "gha", "_ci_config.R"))
  devtools::load_all(here, quiet = TRUE)
})

cfg     <- PIGAUTO_CI_CONFIG
out_dir <- file.path("script", "gha", "results", "avonet")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Data — bundled with the package
# -------------------------------------------------------------------------
e <- new.env()
utils::data("avonet300", package = "pigauto", envir = e)
utils::data("tree300",   package = "pigauto", envir = e)
df <- e$avonet300
rownames(df) <- df$Species_Key
df$Species_Key <- NULL
tree <- e$tree300

# Subset to cap (no-op for AVONET300 which is already 300 species).
n_target <- min(cfg$subset_n, nrow(df))
if (n_target < nrow(df)) {
  set.seed(cfg$seed)
  keep <- sort(sample.int(nrow(df), n_target))
  df   <- df[keep, , drop = FALSE]
  tree <- ape::keep.tip(tree, rownames(df))
}
cat(sprintf("[avonet] n=%d species x %d traits\n", nrow(df), ncol(df)))

# -------------------------------------------------------------------------
# Hold-out mask (30% MCAR, seeded)
# -------------------------------------------------------------------------
t0_split <- Sys.time()
pd0 <- preprocess_traits(df, tree, log_transform = TRUE)
splits <- make_missing_splits(pd0$X_scaled,
                              missing_frac = cfg$missing_frac,
                              val_frac     = 0.5,
                              seed         = cfg$seed,
                              trait_map    = pd0$trait_map)

# Mark which LATENT cells are in the test split.
mask_latent <- matrix(FALSE, nrow = nrow(pd0$X_scaled), ncol = ncol(pd0$X_scaled))
mask_latent[splits$test_idx] <- TRUE

# Map latent cols -> user cols via trait_map (matches the pattern
# used in script/bench_bace_avonet_head_to_head.R).
user_mask_test <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df),
                          dimnames = list(rownames(df), colnames(df)))
for (k in seq_along(pd0$trait_map)) {
  tm <- pd0$trait_map[[k]]
  # Per trait_map convention: single-col traits use tm$name; multi_proportion
  # uses tm$input_cols (a K-vector). Always fall back to tm$name.
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

# -------------------------------------------------------------------------
# Multi-impute with CI config
# -------------------------------------------------------------------------
t0_fit <- Sys.time()
mi <- multi_impute(
  df_masked, tree,
  m              = cfg$n_imputations,
  pool_method    = cfg$pool_method,
  clamp_outliers = cfg$clamp_outliers,
  seed           = cfg$seed
)
t_fit <- as.numeric(difftime(Sys.time(), t0_fit, units = "secs"))

# -------------------------------------------------------------------------
# Evaluate per imputation in user-data scale (long-format canonical rows)
# -------------------------------------------------------------------------
t0_eval <- Sys.time()
ev <- .eval_per_imputation(mi$datasets, df, user_mask_test, pd0$trait_map, t_fit)
t_eval <- as.numeric(difftime(Sys.time(), t0_eval, units = "secs"))

# -------------------------------------------------------------------------
# Normalize + write outputs
# -------------------------------------------------------------------------
out_tbl <- .normalize_eval(ev, dataset = "avonet", method = "pigauto_ci")
saveRDS(out_tbl, file.path(out_dir, "results.rds"))

timings <- list(split_sec = t_split, fit_sec = t_fit, eval_sec = t_eval,
                total_sec = t_split + t_fit + t_eval,
                n_species = nrow(df), n_imputations = cfg$n_imputations)
jsonlite::write_json(timings, file.path(out_dir, "timings.json"),
                     auto_unbox = TRUE, pretty = TRUE)

# Markdown summary
md_lines <- c(
  "# avonet — pigauto CI bench",
  sprintf("Run config: seed=%d, missing_frac=%.2f, n_imputations=%d",
          cfg$seed, cfg$missing_frac, cfg$n_imputations),
  sprintf("N species used: %d", nrow(df)),
  sprintf("Wall time: %.1f s (fit) + %.1f s (eval)", t_fit, t_eval),
  "",
  "## Per-trait medians across imputations",
  ""
)
agg_cont <- stats::aggregate(rmse ~ trait + type, data = out_tbl[!is.na(out_tbl$rmse), ],
                              FUN = function(x) stats::median(x, na.rm = TRUE))
agg_disc <- stats::aggregate(accuracy ~ trait + type,
                              data = out_tbl[!is.na(out_tbl$accuracy), ],
                              FUN = function(x) stats::median(x, na.rm = TRUE))
if (nrow(agg_cont) > 0L) {
  md_lines <- c(md_lines, "### Continuous-family (RMSE)", "",
                knitr::kable(agg_cont, format = "markdown", digits = 4), "")
}
if (nrow(agg_disc) > 0L) {
  md_lines <- c(md_lines, "### Discrete-family (accuracy)", "",
                knitr::kable(agg_disc, format = "markdown", digits = 4), "")
}
writeLines(md_lines, file.path(out_dir, "results.md"))

cat(sprintf("[avonet] done in %.1f s total\n", t_split + t_fit + t_eval))
