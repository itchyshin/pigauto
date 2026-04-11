#!/usr/bin/env Rscript
#
# script/bench_tree_uncertainty.R
#
# Benchmark: tree uncertainty propagation via Rubin's rules
#
# Purpose
#   Demonstrate that pooling across multiple posterior trees via
#   multi_impute_trees() produces wider (and more honest) standard errors
#   than single-tree MI, because the between-imputation variance now
#   includes phylogenetic uncertainty in addition to imputation noise.
#
# Design
#   - Dataset:    avonet300 (300 species, 7 mixed-type traits)
#   - Trees:      trees300 (50 BirdTree Hackett posterior trees)
#   - Missingness: Artificially remove 20%, 40%, 60% of continuous trait
#                  values (MCAR) to create substantial imputation uncertainty.
#   - Methods:
#       (a) single_tree:  multi_impute(tree300, m = 50)
#       (b) multi_tree:   multi_impute_trees(trees300[1:10], m_per_tree = 5)
#                         = 10 trees x 5 = 50 datasets
#   - Downstream model: lm(log(Mass) ~ log(Wing.Length) + log(Beak.Length_Culmen))
#   - Metrics: pooled estimate, SE, df, FMI, RIV from pool_mi()
#   - Replicates: 3 per missingness level (different random masks)
#
# Output
#   script/bench_tree_uncertainty.rds    tidy results + metadata
#   script/bench_tree_uncertainty.md     human-readable summary
#
# Run with
#   cd pigauto && /usr/local/bin/Rscript script/bench_tree_uncertainty.R

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
out_rds <- file.path(here, "script", "bench_tree_uncertainty.rds")
out_md  <- file.path(here, "script", "bench_tree_uncertainty.md")

script_start <- proc.time()[["elapsed"]]

log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start), ..., "\n",
      sep = "")
  flush.console()
}

# -------------------------------------------------------------------------
# Data
# -------------------------------------------------------------------------

data(avonet300, tree300, trees300, package = "pigauto")

# Prepare traits
df_full <- avonet300
rownames(df_full) <- df_full$Species_Key
df_full$Species_Key <- NULL

# Downstream model formula: use only continuous traits that we'll mask
model_formula <- log(Mass) ~ log(Wing.Length) + log(Beak.Length_Culmen)

# -------------------------------------------------------------------------
# Parameters
# -------------------------------------------------------------------------

miss_fracs   <- c(0.20, 0.40, 0.60)
n_reps       <- 3L
n_trees_use  <- 10L          # 10 posterior trees (feasible runtime)
m_per_tree   <- 5L           # 5 MC dropout per tree = 50 total
m_single     <- 50L          # 50 MC dropout for single-tree MI (same total budget)
epochs       <- 300L
cont_cols    <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")

# -------------------------------------------------------------------------
# Helper: introduce MCAR missingness in continuous traits only
# -------------------------------------------------------------------------

punch_holes <- function(df, frac, seed) {
  set.seed(seed)
  df_out <- df
  for (col in cont_cols) {
    n <- nrow(df_out)
    n_miss <- round(n * frac)
    idx <- sample.int(n, n_miss)
    df_out[[col]][idx] <- NA
  }
  df_out
}

# -------------------------------------------------------------------------
# Helper: fit downstream model on one dataset
# -------------------------------------------------------------------------

fit_downstream <- function(d) {
  stats::lm(model_formula, data = d)
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

    df_miss <- punch_holes(df_full, frac, rep_seed)
    n_na <- sum(is.na(df_miss[, cont_cols]))
    log_line(sprintf("  Missing cells in continuous traits: %d / %d (%.1f%%)",
                     n_na, nrow(df_miss) * length(cont_cols),
                     100 * n_na / (nrow(df_miss) * length(cont_cols))))

    # ---- (a) Single-tree MI ------------------------------------------------
    log_line("  Single-tree MI (m=", m_single, ")...")
    t0 <- proc.time()[["elapsed"]]

    mi_single <- tryCatch({
      multi_impute(
        traits       = df_miss,
        tree         = tree300,
        m            = m_single,
        epochs       = epochs,
        missing_frac = 0.25,
        verbose      = FALSE,
        seed         = rep_seed,
        eval_every   = 50L,
        patience     = 20L
      )
    }, error = function(e) {
      log_line("    ERROR: ", conditionMessage(e))
      NULL
    })

    if (!is.null(mi_single)) {
      t1 <- proc.time()[["elapsed"]]
      log_line(sprintf("    Done in %.1fs", t1 - t0))

      fits_single <- with_imputations(mi_single, fit_downstream,
                                      .progress = FALSE,
                                      .on_error = "continue")
      ok_single <- !vapply(fits_single, inherits, logical(1), "pigauto_mi_error")
      log_line(sprintf("    Fits: %d/%d succeeded", sum(ok_single), length(fits_single)))

      if (sum(ok_single) >= 2L) {
        pooled_single <- pool_mi(fits_single[ok_single])
        pooled_single$method       <- "single_tree"
        pooled_single$missing_frac <- frac
        pooled_single$rep          <- rep_id
        pooled_single$n_datasets   <- sum(ok_single)
        pooled_single$wall_time    <- t1 - t0

        result_idx <- result_idx + 1L
        all_results[[result_idx]] <- pooled_single
      }
    }

    # ---- (b) Multi-tree MI -------------------------------------------------
    log_line("  Multi-tree MI (", n_trees_use, " trees x ", m_per_tree, ")...")
    t0 <- proc.time()[["elapsed"]]

    mi_multi <- tryCatch({
      multi_impute_trees(
        traits       = df_miss,
        trees        = trees300[seq_len(n_trees_use)],
        m_per_tree   = m_per_tree,
        epochs       = epochs,
        missing_frac = 0.25,
        verbose      = FALSE,
        seed         = rep_seed,
        eval_every   = 50L,
        patience     = 20L
      )
    }, error = function(e) {
      log_line("    ERROR: ", conditionMessage(e))
      NULL
    })

    if (!is.null(mi_multi)) {
      t1 <- proc.time()[["elapsed"]]
      log_line(sprintf("    Done in %.1fs (%.1f min)", t1 - t0, (t1 - t0) / 60))

      fits_multi <- with_imputations(mi_multi, fit_downstream,
                                     .progress = FALSE,
                                     .on_error = "continue")
      ok_multi <- !vapply(fits_multi, inherits, logical(1), "pigauto_mi_error")
      log_line(sprintf("    Fits: %d/%d succeeded", sum(ok_multi), length(fits_multi)))

      if (sum(ok_multi) >= 2L) {
        pooled_multi <- pool_mi(fits_multi[ok_multi])
        pooled_multi$method       <- "multi_tree"
        pooled_multi$missing_frac <- frac
        pooled_multi$rep          <- rep_id
        pooled_multi$n_datasets   <- sum(ok_multi)
        pooled_multi$wall_time    <- t1 - t0

        result_idx <- result_idx + 1L
        all_results[[result_idx]] <- pooled_multi
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
  n_species    = nrow(avonet300),
  n_trees      = n_trees_use,
  m_per_tree   = m_per_tree,
  m_single     = m_single,
  epochs       = epochs,
  miss_fracs   = miss_fracs,
  n_reps       = n_reps,
  model_formula = deparse(model_formula),
  timestamp    = Sys.time(),
  wall_time    = proc.time()[["elapsed"]] - script_start
)

saveRDS(list(results = results, meta = meta), out_rds)
log_line("Wrote ", out_rds)

# -------------------------------------------------------------------------
# Markdown summary
# -------------------------------------------------------------------------

md <- character()
md <- c(md, "# Tree uncertainty benchmark\n")
md <- c(md, sprintf("- Species: %d (avonet300)", meta$n_species))
md <- c(md, sprintf("- Trees: %d posterior (BirdTree Hackett)", n_trees_use))
md <- c(md, sprintf("- Methods: single_tree (m=%d) vs multi_tree (%d trees x %d = %d)",
                     m_single, n_trees_use, m_per_tree, n_trees_use * m_per_tree))
md <- c(md, sprintf("- Downstream model: %s", meta$model_formula))
md <- c(md, sprintf("- Missingness fracs: %s", paste(miss_fracs, collapse = ", ")))
md <- c(md, sprintf("- Replicates: %d", n_reps))
md <- c(md, sprintf("- Total wall time: %.1f min\n", meta$wall_time / 60))

# Summary table: average SE and FMI by method and frac
if (nrow(results) > 0) {
  md <- c(md, "## Pooled SE comparison\n")
  md <- c(md, "| method | miss_frac | term | mean_estimate | mean_SE | mean_FMI | mean_df |")
  md <- c(md, "|--------|-----------|------|---------------|---------|----------|---------|")

  agg <- aggregate(
    cbind(estimate, std.error, fmi, df) ~ method + missing_frac + term,
    data = results, FUN = mean, na.rm = TRUE
  )
  agg <- agg[order(agg$missing_frac, agg$term, agg$method), ]

  for (i in seq_len(nrow(agg))) {
    md <- c(md, sprintf("| %s | %.2f | %s | %.4f | %.4f | %.3f | %.1f |",
                        agg$method[i], agg$missing_frac[i], agg$term[i],
                        agg$estimate[i], agg$std.error[i],
                        agg$fmi[i], agg$df[i]))
  }

  # SE ratio: multi / single
  md <- c(md, "\n## SE ratio (multi_tree / single_tree)\n")
  single <- results[results$method == "single_tree", ]
  multi  <- results[results$method == "multi_tree", ]

  if (nrow(single) > 0 && nrow(multi) > 0) {
    se_single <- aggregate(std.error ~ missing_frac + term, data = single,
                           FUN = mean, na.rm = TRUE)
    se_multi  <- aggregate(std.error ~ missing_frac + term, data = multi,
                           FUN = mean, na.rm = TRUE)
    merged <- merge(se_single, se_multi, by = c("missing_frac", "term"),
                    suffixes = c("_single", "_multi"))
    merged$ratio <- merged$std.error_multi / merged$std.error_single

    md <- c(md, "| miss_frac | term | SE_single | SE_multi | ratio |")
    md <- c(md, "|-----------|------|-----------|----------|-------|")
    for (i in seq_len(nrow(merged))) {
      md <- c(md, sprintf("| %.2f | %s | %.4f | %.4f | %.2f |",
                          merged$missing_frac[i], merged$term[i],
                          merged$std.error_single[i],
                          merged$std.error_multi[i],
                          merged$ratio[i]))
    }
  }
}

md <- c(md, "\n---\nGenerated:", format(Sys.time(), "%Y-%m-%d %H:%M"))
writeLines(md, out_md)
log_line("Wrote ", out_md)
log_line("done")
