#!/usr/bin/env Rscript
#
# script/bench_delhey.R
#
# Benchmark: environmental covariate support on Delhey 2019 plumage data
#
# Purpose
#   Demonstrate that environmental covariates (temperature, precipitation,
#   tree cover, latitude) improve imputation of plumage lightness traits
#   compared to phylogeny alone. This is the key comparison for closing
#   the feature gap with TrEvol.
#
# Design
#   - Dataset: delhey5809 (5,809 passerine species, 2 continuous traits)
#   - Tree:    tree_delhey (Hackett MCC, pruned)
#   - Covariates: 6 environmental variables from Delhey 2019
#   - Missingness: artificially remove 20%, 40%, 60% of lightness values
#   - Methods:
#       (a) mean:            column mean imputation
#       (b) BM_only:         phylogenetic baseline (Rphylopars BM)
#       (c) pigauto:         pigauto without covariates
#       (d) pigauto_covs:    pigauto with 6 environmental covariates
#   - Metrics: RMSE, Pearson r (on test split)
#   - Replicates: 3 per missingness level
#
# Output
#   script/bench_delhey.rds    tidy results + metadata
#   script/bench_delhey.md     human-readable summary
#
# Run with
#   cd pigauto && /usr/local/bin/Rscript script/bench_delhey.R

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto",
    quiet = TRUE
  )
})

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "bench_delhey.rds")
out_md  <- file.path(here, "script", "bench_delhey.md")

script_start <- proc.time()[["elapsed"]]

log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start), ..., "\n",
      sep = "")
  flush.console()
}

# -------------------------------------------------------------------------
# Data
# -------------------------------------------------------------------------

data(delhey5809, tree_delhey, package = "pigauto")

# Traits: only lightness columns (to be imputed)
trait_cols <- c("lightness_male", "lightness_female")

# Covariates: all environmental variables
cov_cols <- c("annual_mean_temperature", "annual_precipitation",
              "percent_tree_cover", "mean_temperature_of_warmest_quarter",
              "precipitation_of_warmest_quarter", "midLatitude")

df_full <- delhey5809
rownames(df_full) <- df_full$Species_Key
df_full$Species_Key <- NULL
df_full$family <- NULL  # not a trait, not a covariate for this benchmark

# Separate traits and covariates
covs <- df_full[, cov_cols, drop = FALSE]
traits_only <- df_full[, trait_cols, drop = FALSE]

log_line("Data: ", nrow(traits_only), " species, ",
         length(trait_cols), " traits, ", length(cov_cols), " covariates")

# -------------------------------------------------------------------------
# Parameters
# -------------------------------------------------------------------------

miss_fracs <- c(0.20, 0.40, 0.60)
n_reps     <- 3L
epochs     <- 1000L

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

punch_holes <- function(df, frac, seed) {
  set.seed(seed)
  df_out <- df
  n <- nrow(df_out)
  for (col in names(df_out)) {
    n_miss <- round(n * frac)
    idx <- sample.int(n, n_miss)
    df_out[[col]][idx] <- NA
  }
  df_out
}

mean_impute <- function(X_latent, splits) {
  X_train <- X_latent
  X_train[c(splits$val_idx, splits$test_idx)] <- NA
  col_means <- colMeans(X_train, na.rm = TRUE)
  col_means[!is.finite(col_means)] <- 0
  matrix(col_means,
         nrow = nrow(X_latent), ncol = ncol(X_latent),
         byrow = TRUE, dimnames = dimnames(X_latent))
}

tag_rows <- function(ev_df, method, frac, rep_id) {
  if (is.null(ev_df) || nrow(ev_df) == 0L) return(NULL)
  ev_df$method       <- method
  ev_df$missing_frac <- frac
  ev_df$rep          <- rep_id
  ev_df
}

# -------------------------------------------------------------------------
# Main loop
# -------------------------------------------------------------------------

all_results <- list()
result_idx  <- 0L

for (frac in miss_fracs) {
  for (rep_id in seq_len(n_reps)) {
    rep_seed <- rep_id * 1000L + round(frac * 100)

    log_line(sprintf("=== frac=%.2f, rep=%d (seed=%d) ===", frac, rep_id, rep_seed))

    df_miss <- punch_holes(traits_only, frac, rep_seed)
    n_na <- sum(is.na(df_miss))
    log_line(sprintf("  Missing cells: %d / %d (%.1f%%)",
                     n_na, nrow(df_miss) * ncol(df_miss),
                     100 * n_na / (nrow(df_miss) * ncol(df_miss))))

    # Common preprocessing (no covariates)
    pd <- tryCatch(
      preprocess_traits(df_miss, tree_delhey, log_transform = FALSE),
      error = function(e) { log_line("  ERROR preprocess: ", e$message); NULL }
    )
    if (is.null(pd)) next

    spl <- make_missing_splits(pd$X_scaled, missing_frac = 0.25,
                               seed = rep_seed, trait_map = pd$trait_map)
    graph <- build_phylo_graph(tree_delhey, k_eigen = "auto")

    # ---- (a) Mean imputation -----------------------------------------------
    log_line("  mean...")
    mean_pred <- mean_impute(pd$X_scaled, spl)
    ev_mean <- evaluate_imputation(mean_pred, pd$X_scaled, spl,
                                   trait_map = pd$trait_map)
    result_idx <- result_idx + 1L
    all_results[[result_idx]] <- tag_rows(ev_mean, "mean", frac, rep_id)

    # ---- (b) BM baseline ---------------------------------------------------
    log_line("  BM baseline...")
    bl <- tryCatch(
      fit_baseline(pd, tree_delhey, splits = spl, graph = graph),
      error = function(e) { log_line("    ERROR: ", e$message); NULL }
    )
    if (!is.null(bl)) {
      ev_bm <- evaluate_imputation(bl$mu, pd$X_scaled, spl,
                                   trait_map = pd$trait_map)
      result_idx <- result_idx + 1L
      all_results[[result_idx]] <- tag_rows(ev_bm, "BM_only", frac, rep_id)
    }

    # Free D before GNN training
    graph$D <- NULL
    invisible(gc(full = TRUE, verbose = FALSE))

    # ---- (c) pigauto without covariates ------------------------------------
    log_line("  pigauto (no covariates)...")
    t0 <- proc.time()[["elapsed"]]
    fit_nocov <- tryCatch(
      fit_pigauto(pd, tree_delhey, splits = spl, graph = graph,
                  baseline = bl, epochs = epochs, verbose = FALSE,
                  seed = rep_seed, eval_every = 100L, patience = 10L),
      error = function(e) { log_line("    ERROR: ", e$message); NULL }
    )
    if (!is.null(fit_nocov)) {
      pred_nocov <- predict(fit_nocov, return_se = TRUE, n_imputations = 1L)
      ev_nocov <- evaluate_imputation(pred_nocov, pd$X_scaled, spl)
      result_idx <- result_idx + 1L
      all_results[[result_idx]] <- tag_rows(ev_nocov, "pigauto", frac, rep_id)
      t1 <- proc.time()[["elapsed"]]
      log_line(sprintf("    Done in %.1fs", t1 - t0))
    }

    # ---- (d) pigauto WITH covariates ---------------------------------------
    log_line("  pigauto (with covariates)...")
    t0 <- proc.time()[["elapsed"]]

    # Re-preprocess with covariates
    pd_cov <- tryCatch(
      preprocess_traits(df_miss, tree_delhey, log_transform = FALSE,
                        covariates = covs),
      error = function(e) { log_line("    ERROR preprocess_cov: ", e$message); NULL }
    )
    if (!is.null(pd_cov)) {
      # Reuse splits and graph from before (same masked data, same tree)
      graph_cov <- build_phylo_graph(tree_delhey, k_eigen = "auto")
      bl_cov <- tryCatch(
        fit_baseline(pd_cov, tree_delhey, splits = spl, graph = graph_cov),
        error = function(e) { log_line("    ERROR bl_cov: ", e$message); NULL }
      )
      graph_cov$D <- NULL
      invisible(gc(full = TRUE, verbose = FALSE))

      if (!is.null(bl_cov)) {
        fit_cov <- tryCatch(
          fit_pigauto(pd_cov, tree_delhey, splits = spl, graph = graph_cov,
                      baseline = bl_cov, epochs = epochs, verbose = FALSE,
                      seed = rep_seed, eval_every = 100L, patience = 10L),
          error = function(e) { log_line("    ERROR fit_cov: ", e$message); NULL }
        )
        if (!is.null(fit_cov)) {
          pred_cov <- predict(fit_cov, return_se = TRUE, n_imputations = 1L)
          ev_cov <- evaluate_imputation(pred_cov, pd_cov$X_scaled, spl)
          result_idx <- result_idx + 1L
          all_results[[result_idx]] <- tag_rows(ev_cov, "pigauto_covs", frac, rep_id)
          t1 <- proc.time()[["elapsed"]]
          log_line(sprintf("    Done in %.1fs", t1 - t0))
        }
      }
    }

    log_line("  Cell done.\n")
  }
}

# -------------------------------------------------------------------------
# Assemble results
# -------------------------------------------------------------------------

results <- do.call(rbind, all_results)
rownames(results) <- NULL

log_line("Total rows: ", nrow(results))

# -------------------------------------------------------------------------
# Save RDS
# -------------------------------------------------------------------------

meta <- list(
  n_species   = nrow(delhey5809),
  n_traits    = length(trait_cols),
  n_covariates = length(cov_cols),
  trait_cols   = trait_cols,
  cov_cols     = cov_cols,
  epochs       = epochs,
  miss_fracs   = miss_fracs,
  n_reps       = n_reps,
  timestamp    = Sys.time(),
  wall_time    = proc.time()[["elapsed"]] - script_start
)

saveRDS(list(results = results, meta = meta), out_rds)
log_line("Wrote ", out_rds)

# -------------------------------------------------------------------------
# Markdown summary
# -------------------------------------------------------------------------

md <- character()
md <- c(md, "# Delhey plumage lightness benchmark (covariates)\n")
md <- c(md, sprintf("- Species: %d (Delhey et al. 2019)", meta$n_species))
md <- c(md, sprintf("- Traits: %s", paste(trait_cols, collapse = ", ")))
md <- c(md, sprintf("- Covariates: %s", paste(cov_cols, collapse = ", ")))
md <- c(md, sprintf("- Methods: mean, BM_only, pigauto, pigauto_covs"))
md <- c(md, sprintf("- Missingness: %s", paste(miss_fracs, collapse = ", ")))
md <- c(md, sprintf("- Replicates: %d", n_reps))
md <- c(md, sprintf("- Total wall time: %.1f min\n", meta$wall_time / 60))

if (nrow(results) > 0) {
  # RMSE comparison
  md <- c(md, "## Test-set RMSE\n")
  rmse_agg <- aggregate(rmse ~ method + missing_frac + trait,
                        data = results[results$split == "test", ],
                        FUN = mean, na.rm = TRUE)
  rmse_agg <- rmse_agg[order(rmse_agg$missing_frac, rmse_agg$trait, rmse_agg$method), ]

  md <- c(md, "| method | miss_frac | trait | mean_RMSE |")
  md <- c(md, "|--------|-----------|-------|-----------|")
  for (i in seq_len(nrow(rmse_agg))) {
    md <- c(md, sprintf("| %s | %.2f | %s | %.4f |",
                        rmse_agg$method[i], rmse_agg$missing_frac[i],
                        rmse_agg$trait[i], rmse_agg$rmse[i]))
  }

  # Pearson r comparison
  md <- c(md, "\n## Test-set Pearson r\n")
  r_agg <- aggregate(pearson_r ~ method + missing_frac + trait,
                     data = results[results$split == "test", ],
                     FUN = mean, na.rm = TRUE)
  r_agg <- r_agg[order(r_agg$missing_frac, r_agg$trait, r_agg$method), ]

  md <- c(md, "| method | miss_frac | trait | mean_r |")
  md <- c(md, "|--------|-----------|-------|--------|")
  for (i in seq_len(nrow(r_agg))) {
    md <- c(md, sprintf("| %s | %.2f | %s | %.4f |",
                        r_agg$method[i], r_agg$missing_frac[i],
                        r_agg$trait[i], r_agg$pearson_r[i]))
  }

  # RMSE improvement: pigauto_covs vs pigauto
  md <- c(md, "\n## Covariate lift (pigauto_covs / pigauto RMSE ratio)\n")
  rmse_pg <- rmse_agg[rmse_agg$method == "pigauto", ]
  rmse_pc <- rmse_agg[rmse_agg$method == "pigauto_covs", ]
  if (nrow(rmse_pg) > 0 && nrow(rmse_pc) > 0) {
    merged <- merge(rmse_pg, rmse_pc,
                    by = c("missing_frac", "trait"),
                    suffixes = c("_pg", "_pc"))
    merged$ratio <- merged$rmse_pc / merged$rmse_pg

    md <- c(md, "| miss_frac | trait | RMSE_pigauto | RMSE_covs | ratio |")
    md <- c(md, "|-----------|-------|--------------|-----------|-------|")
    for (i in seq_len(nrow(merged))) {
      md <- c(md, sprintf("| %.2f | %s | %.4f | %.4f | %.3f |",
                          merged$missing_frac[i], merged$trait[i],
                          merged$rmse_pg[i], merged$rmse_pc[i],
                          merged$ratio[i]))
    }
  }
}

md <- c(md, "\n---\nGenerated:", format(Sys.time(), "%Y-%m-%d %H:%M"))
writeLines(md, out_md)
log_line("Wrote ", out_md)
log_line("done")
