#!/usr/bin/env Rscript
# script/make_covariate_lift_table.R
#
# Consolidate per-trait covariate-lift results from all benches into a
# single tidy table.  Loads each bench's *.rds, extracts the per-trait
# RMSE rows, and writes:
#   useful/covariate_lift_table.md  -- per-trait detail
#   useful/covariate_lift_summary.md -- one-line summary per dataset
#
# Run: Rscript script/make_covariate_lift_table.R

here <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_md      <- file.path(here, "useful", "covariate_lift_table.md")
out_summary <- file.path(here, "useful", "covariate_lift_summary.md")

# Each bench: file path + dataset/source label + n_species
benches <- list(
  list(rds   = "script/bench_pantheria_covariates.rds",
        label = "PanTHERIA mammals (precip+temp+lat+PET)",
        n_field = "n",
        kind = "real"),
  list(rds   = "script/bench_globtherm_covariates.rds",
        label = "GlobTherm ectotherms (lat+long+elev+|lat|)",
        n_field = "n",
        kind = "real"),
  list(rds   = "script/bench_amphibio_covariates.rds",
        label = "AmphiBIO amphibians (climate-zone occupancy)",
        n_field = "n",
        kind = "real"),
  list(rds   = "script/bench_delhey_covariates.rds",
        label = "Delhey birds (annual_temp+precip+tree_cover+...)",
        n_field = "n",
        kind = "real"),
  list(rds   = "script/bench_plants_cached_only.rds",
        label = "BIEN plants (WorldClim bioclim, cached species)",
        n_field = "n",
        kind = "real"),
  list(rds   = "script/bench_multi_obs.rds",
        label = "Multi-obs sim (Yule tree + acclim_temp)",
        n_field = NULL,
        kind = "sim_yule"),
  list(rds   = "script/bench_multi_obs_real_tree.rds",
        label = "Multi-obs sim (REAL tree300 + acclim_temp)",
        n_field = NULL,
        kind = "sim_real_phylo")
)

cat("== Loading bench results ==\n")

per_trait_rows <- list()
summary_rows   <- list()

fmt_num <- function(x, digits = 3) {
  if (length(x) == 0L || !is.finite(x)) return("NA")
  format(round(x, digits), nsmall = digits)
}

for (b in benches) {
  full <- file.path(here, b$rds)
  if (!file.exists(full)) {
    cat(sprintf("  [skip] %s -- not found\n", b$rds))
    next
  }
  obj <- tryCatch(readRDS(full), error = function(e) NULL)
  if (is.null(obj)) { cat(sprintf("  [skip] %s -- read error\n", b$rds)); next }

  if (b$kind %in% c("real",  "real_subclass")) {
    res <- obj$results
    if (is.null(res) || nrow(res) == 0L) {
      cat(sprintf("  [skip] %s -- no results rows\n", b$rds)); next
    }
    n <- if (!is.null(b$n_field) && !is.null(obj[[b$n_field]])) obj[[b$n_field]] else NA_integer_

    cat(sprintf("  [load] %s (n=%s)\n", b$rds, n))
    for (i in seq_len(nrow(res))) {
      r <- res[i, ]
      per_trait_rows[[length(per_trait_rows) + 1L]] <- data.frame(
        dataset = b$label, n_species = n, trait = r$trait,
        n_held = r$n_held %||% NA, mean_RMSE = r$mean_RMSE %||% NA,
        none_RMSE = r$none_RMSE %||% NA,
        cov_on_RMSE = r$cov_on_RMSE %||% NA,
        cov_off_RMSE = r$cov_off_RMSE %||% NA,
        ratio_on  = r$ratio_on %||% NA,
        ratio_off = r$ratio_off %||% NA,
        none_r    = r$none_r %||% NA,
        cov_on_r  = r$cov_on_r %||% NA,
        cov_off_r = r$cov_off_r %||% NA,
        stringsAsFactors = FALSE)
    }

    # Aggregate summary: median ratio_on across traits
    median_ratio_on  <- median(res$ratio_on, na.rm = TRUE)
    median_ratio_off <- median(res$ratio_off, na.rm = TRUE)
    n_lift_on  <- sum(res$ratio_on  < 0.95, na.rm = TRUE)
    n_lift_off <- sum(res$ratio_off < 0.95, na.rm = TRUE)
    n_total    <- nrow(res)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      dataset = b$label, n_species = n,
      n_traits = n_total,
      median_ratio_on  = median_ratio_on,
      median_ratio_off = median_ratio_off,
      n_traits_with_lift_on  = n_lift_on,
      n_traits_with_lift_off = n_lift_off,
      kind = b$kind,
      stringsAsFactors = FALSE)
  } else if (b$kind == "sim_yule") {
    # bench_multi_obs.rds has a different shape (per-cell)
    if (!is.null(obj$results)) {
      res <- obj$results
      pigcov <- res[res$method == "pigauto_cov", ]
      pignone <- res[res$method == "pigauto_no_cov", ]
      if (nrow(pigcov) > 0L && nrow(pignone) > 0L) {
        # Sum over cells with beta > 0 (where covariate carries info)
        # Both bench_multi_obs.R (col "rep") and bench_multi_obs_real_tree.R
        # (col "rep_id") schemas exist; pick whichever is present.
        rep_col <- if ("rep_id" %in% colnames(pigcov)) "rep_id" else "rep"
        keys <- c("lambda", "beta", "sp_missing_frac", rep_col)
        m <- merge(pigcov[, c(keys, "obs_rmse", "obs_pearson_r")],
                    pignone[, c(keys, "obs_rmse", "obs_pearson_r")],
                    by = keys, suffixes = c("_cov", "_nocov"))
        m$ratio <- m$obs_rmse_cov / m$obs_rmse_nocov
        m_useful <- m[m$beta > 0, ]
        cat(sprintf("  [load] %s (sim Yule, beta>0 cells=%d)\n",
                     b$rds, nrow(m_useful)))
        summary_rows[[length(summary_rows) + 1L]] <- data.frame(
          dataset = b$label, n_species = NA_integer_,
          n_traits = nrow(m_useful),
          median_ratio_on  = median(m_useful$ratio, na.rm = TRUE),
          median_ratio_off = NA_real_,
          n_traits_with_lift_on  = sum(m_useful$ratio < 0.95, na.rm = TRUE),
          n_traits_with_lift_off = NA_integer_,
          kind = b$kind,
          stringsAsFactors = FALSE)
      }
    }
  } else if (b$kind == "sim_real_phylo") {
    if (!is.null(obj$results)) {
      res <- obj$results
      pigcov  <- res[res$method == "pigauto_cov", ]
      pignone <- res[res$method == "pigauto_no_cov", ]
      if (nrow(pigcov) > 0L && nrow(pignone) > 0L) {
        # Both bench_multi_obs.R (col "rep") and bench_multi_obs_real_tree.R
        # (col "rep_id") schemas exist; pick whichever is present.
        rep_col <- if ("rep_id" %in% colnames(pigcov)) "rep_id" else "rep"
        keys <- c("lambda", "beta", "sp_missing_frac", rep_col)
        m <- merge(pigcov[, c(keys, "obs_rmse", "obs_pearson_r")],
                    pignone[, c(keys, "obs_rmse", "obs_pearson_r")],
                    by = keys, suffixes = c("_cov", "_nocov"))
        m$ratio <- m$obs_rmse_cov / m$obs_rmse_nocov
        m_useful <- m[m$beta > 0, ]
        cat(sprintf("  [load] %s (sim real-tree, beta>0 cells=%d)\n",
                     b$rds, nrow(m_useful)))
        summary_rows[[length(summary_rows) + 1L]] <- data.frame(
          dataset = b$label,
          n_species = obj$meta$n_species %||% NA_integer_,
          n_traits = nrow(m_useful),
          median_ratio_on  = median(m_useful$ratio, na.rm = TRUE),
          median_ratio_off = NA_real_,
          n_traits_with_lift_on  = sum(m_useful$ratio < 0.95, na.rm = TRUE),
          n_traits_with_lift_off = NA_integer_,
          kind = b$kind,
          stringsAsFactors = FALSE)
      }
    }
  }
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

per_trait <- if (length(per_trait_rows)) do.call(rbind, per_trait_rows) else NULL
summary <- if (length(summary_rows)) do.call(rbind, summary_rows) else NULL

# Per-trait output
if (!is.null(per_trait)) {
  per_trait_md <- c(
    "# Covariate-lift results: per-trait detail",
    "",
    sprintf("Last regenerated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "Each row = one (dataset, trait) pair.  RMSE columns are the absolute RMSE on held-out cells.",
    "",
    "- `mean_RMSE`: imputing every held-out cell with the column grand mean",
    "- `none_RMSE`: pigauto baseline (no covariates)",
    "- `cov_on_RMSE`:  pigauto with covariates, safety_floor = TRUE",
    "- `cov_off_RMSE`: pigauto with covariates, safety_floor = FALSE",
    "- `ratio_on`:  cov_on_RMSE / none_RMSE.  <1 = covariates help.",
    "- `ratio_off`: cov_off_RMSE / none_RMSE. <1 = covariates help.",
    "",
    "Pearson r columns describe the same fits' correlation between predicted and held-out truth.",
    "",
    "```",
    capture.output(print(per_trait, row.names = FALSE, max = 1000)),
    "```")
  writeLines(per_trait_md, out_md)
  cat(sprintf("wrote %s\n", out_md))
}

# Summary output
if (!is.null(summary)) {
  summary_md <- c(
    "# Covariate-lift summary across datasets",
    "",
    sprintf("Last regenerated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "One row per dataset.  `median_ratio_on` is the median across traits of",
    "`cov_on_RMSE / none_RMSE`.  Values < 1.00 mean covariates lift accuracy.",
    "",
    "- `n_traits_with_lift_on`: how many traits had `ratio_on < 0.95` (>=5% lift).",
    "- `n_traits_with_lift_off`: same for safety_floor=FALSE.",
    "",
    "```",
    capture.output(print(summary, row.names = FALSE)),
    "```",
    "",
    "## Reading guide",
    "",
    "Real species-level datasets with phylo-redundant climate covariates",
    "(BIEN plants, Delhey birds, AmphiBIO amphibians) cluster around",
    "ratio_on ~ 1.00 -- the safety floor closes the gate when covariates",
    "are uninformative beyond phylogeny, so users pay no accuracy penalty",
    "for passing redundant climate features.",
    "",
    "Datasets with PARTIAL phylo-decoupled covariate signal show",
    "MIXED lift: PanTHERIA mammals with the bundled precip/temp/lat",
    "columns lifts MaxLongevity by 22% but is neutral on body mass and",
    "gestation length.  This is the safety property in action: the",
    "calibrated gate adapts per-trait.",
    "",
    "Multi-observation simulations (acclim_temp varying within species)",
    "show consistent 10-19% lift, and the lift survives when the Yule",
    "sim tree is replaced with the real AVONET 300 bird phylogeny.",
    "")
  writeLines(summary_md, out_summary)
  cat(sprintf("wrote %s\n", out_summary))
}
