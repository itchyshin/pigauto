#!/usr/bin/env Rscript
#
# script/bench_scaling_v031.R
#
# Post-fix rerun of the Part 3 scaling benchmark, after landing Fix A
# (RSpectra sparse Lanczos eigensolver) and Fix B (cophenetic caching
# across the pipeline). Same workload and measurement harness as
# script/bench_scaling.R; only the output filenames differ, so that
# both before and after artifacts can coexist on disk for the pkgdown
# scaling article.
#
# Design
#   * Tree: ape::rcoal(n), seeded by n for reproducibility (identical
#     to the pre-fix baseline).
#   * Traits: 3 continuous (BM) + 1 binary (thresholded BM) +
#             1 ordinal (3-level binned BM), 15% MCAR on observed cells.
#   * Stages (sequential, dependent): preprocess, graph, splits,
#             baseline, train (epochs = 100), predict (5 imputations).
#   * Per-stage wall budget: 30 minutes.
#   * Overall wall budget: 3 hours.
#   * Checkpointing: results are saved after every stage.
#
# Output
#   script/bench_scaling_v031.rds           per-stage data.frame
#   script/bench_scaling_v031.png           wall-clock plot (log y)
#   script/bench_scaling_v031_summary.md    machine info + verdict
#
# Run with
#   cd pigauto && /usr/local/bin/Rscript script/bench_scaling_v031.R

options(warn = 1, stringsAsFactors = FALSE)

# IMPORTANT: load the dev tree, not the installed package. This script
# lives on the fix/scaling-rspectra-cophenetic-cache branch and depends
# on Fix A + Fix B which are not yet in the installed v0.3.0 build.
suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto",
    quiet = TRUE
  )
})

# -------------------------------------------------------------------------
# Paths and state
# -------------------------------------------------------------------------

here       <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds    <- file.path(here, "script", "bench_scaling_v031.rds")
out_png    <- file.path(here, "script", "bench_scaling_v031.png")
out_md     <- file.path(here, "script", "bench_scaling_v031_summary.md")

n_grid             <- c(300L, 1000L, 2000L, 3000L, 5000L, 7500L, 10000L)
stage_budget_sec   <- 30 * 60      # 30 minutes per stage
overall_budget_sec <- 3 * 3600     # 3 hours total

# Sanity check: RSpectra should be available so Fix A fires.
has_rspectra <- requireNamespace("RSpectra", quietly = TRUE)

results <- data.frame(
  n          = integer(),
  stage      = character(),
  wall_sec   = numeric(),
  max_mb     = numeric(),
  status     = character(),
  error_msg  = character(),
  stringsAsFactors = FALSE
)

script_start <- proc.time()[["elapsed"]]

log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start), ..., "\n",
      sep = "")
  flush.console()
}

append_result <- function(n, stage, wall_sec, max_mb, status,
                          error_msg = NA_character_) {
  row <- data.frame(
    n         = n, stage = stage,
    wall_sec  = wall_sec, max_mb = max_mb,
    status    = status, error_msg = error_msg,
    stringsAsFactors = FALSE
  )
  results <<- rbind(results, row)
  saveRDS(results, out_rds)
}

# -------------------------------------------------------------------------
# Per-stage timing wrapper
# -------------------------------------------------------------------------

time_stage <- function(n, stage_label, expr_fun) {
  gc(reset = TRUE, full = TRUE)
  gc_before <- gc(reset = TRUE, full = TRUE)

  t0 <- proc.time()[["elapsed"]]
  status    <- "ok"
  error_msg <- NA_character_
  value     <- NULL

  tryCatch({
    setTimeLimit(elapsed = stage_budget_sec, transient = TRUE)
    value <- expr_fun()
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (grepl("reached elapsed time limit", msg, fixed = TRUE)) {
      status    <<- "timeout"
      error_msg <<- sprintf("exceeded %.0fs per-stage budget", stage_budget_sec)
    } else {
      status    <<- "error"
      error_msg <<- msg
    }
  })
  setTimeLimit(elapsed = Inf, transient = FALSE)

  wall <- proc.time()[["elapsed"]] - t0
  gc_after <- gc(full = TRUE)
  mb <- tryCatch(as.numeric(gc_after["Vcells", 7L]),
                 error = function(e) NA_real_)

  append_result(n, stage_label, wall, mb, status, error_msg)

  log_line(sprintf("[n=%5d] %-10s %-7s wall=%7.1fs  maxR=%6.0f MB%s",
                   n, stage_label, status, wall, mb,
                   if (!is.na(error_msg)) paste0("  -- ", substr(error_msg, 1, 120)) else ""))

  list(status = status, value = value, wall = wall)
}

# -------------------------------------------------------------------------
# Mixed-type dataset simulator (same as pre-fix baseline)
# -------------------------------------------------------------------------

simulate_dataset <- function(n) {
  set.seed(n)
  tree <- ape::rcoal(n)

  x1 <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  x2 <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  x3 <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  z1 <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  z2 <- ape::rTraitCont(tree, model = "BM", sigma = 1)

  binary  <- factor(ifelse(z1 > stats::median(z1), "high", "low"),
                    levels = c("low", "high"))
  ordinal <- cut(z2,
                 breaks = stats::quantile(z2, probs = c(0, 1/3, 2/3, 1)),
                 include.lowest = TRUE,
                 labels = c("small", "medium", "large"))
  ordinal <- factor(ordinal, levels = c("small", "medium", "large"),
                    ordered = TRUE)

  traits <- data.frame(
    cont1   = as.numeric(x1),
    cont2   = as.numeric(x2),
    cont3   = as.numeric(x3),
    bin     = binary,
    ord     = ordinal,
    row.names = tree$tip.label,
    stringsAsFactors = FALSE
  )

  set.seed(1000L + n)
  mask <- matrix(stats::runif(nrow(traits) * ncol(traits)) < 0.15,
                 nrow = nrow(traits), ncol = ncol(traits))
  for (j in seq_len(ncol(traits))) {
    traits[mask[, j], j] <- NA
  }

  list(traits = traits, tree = tree)
}

# -------------------------------------------------------------------------
# Per-size runner
# -------------------------------------------------------------------------

run_one <- function(n) {
  log_line(sprintf("===== n = %d =====", n))

  sim <- tryCatch(simulate_dataset(n), error = function(e) {
    append_result(n, "simulate", NA_real_, NA_real_, "error", conditionMessage(e))
    NULL
  })
  if (is.null(sim)) return(invisible(NULL))
  traits <- sim$traits
  tree   <- sim$tree

  st <- time_stage(n, "preprocess", function() {
    preprocess_traits(traits, tree, log_transform = FALSE)
  })
  if (st$status != "ok") return(invisible(NULL))
  pd <- st$value

  st <- time_stage(n, "graph", function() {
    build_phylo_graph(tree, k_eigen = 8L, sigma_mult = 0.5)
  })
  if (st$status != "ok") return(invisible(NULL))
  graph <- st$value

  st <- time_stage(n, "splits", function() {
    make_missing_splits(pd$X_scaled, missing_frac = 0.25,
                        val_frac = 0.25, seed = 42L,
                        trait_map = pd$trait_map)
  })
  if (st$status != "ok") return(invisible(NULL))
  splits <- st$value

  # IMPORTANT: pass graph to fit_baseline so Fix B fires (cophenetic cache).
  st <- time_stage(n, "baseline", function() {
    fit_baseline(pd, tree, splits = splits, model = "BM", graph = graph)
  })
  if (st$status != "ok") return(invisible(NULL))
  baseline <- st$value

  st <- time_stage(n, "train", function() {
    fit_pigauto(
      data            = pd,
      tree            = tree,
      splits          = splits,
      graph           = graph,
      baseline        = baseline,
      hidden_dim      = 64L,
      k_eigen         = 8L,
      dropout         = 0.10,
      lr              = 3e-3,
      epochs          = 100L,
      corruption_rate = 0.55,
      lambda_shrink   = 0.03,
      eval_every      = 25L,
      patience        = 3L,
      verbose         = FALSE,
      seed            = 1L
    )
  })
  if (st$status != "ok") return(invisible(NULL))
  fit <- st$value

  st <- time_stage(n, "predict", function() {
    stats::predict(fit, return_se = TRUE, n_imputations = 5L)
  })
  invisible(NULL)
}

# -------------------------------------------------------------------------
# Main loop
# -------------------------------------------------------------------------

log_line(sprintf("bench_scaling_v031 starting at %s", format(Sys.time())))
log_line(sprintf("RSpectra available: %s", has_rspectra))
log_line(sprintf("n_grid = %s", paste(n_grid, collapse = ", ")))
log_line(sprintf("stage budget = %ds, overall budget = %ds",
                 stage_budget_sec, overall_budget_sec))

for (n in n_grid) {
  elapsed <- proc.time()[["elapsed"]] - script_start
  if (elapsed > overall_budget_sec) {
    log_line(sprintf("overall budget exhausted at %.0fs, skipping n = %d and later",
                     elapsed, n))
    append_result(n, "budget", elapsed, NA_real_, "budget-exhausted")
    break
  }
  run_one(n)
}

log_line(sprintf("loop finished at %s (total elapsed %.0fs)",
                 format(Sys.time()),
                 proc.time()[["elapsed"]] - script_start))

# -------------------------------------------------------------------------
# Plot
# -------------------------------------------------------------------------

if (requireNamespace("ggplot2", quietly = TRUE) && nrow(results) > 0) {
  plot_df <- results[results$status == "ok" & !is.na(results$wall_sec), ]
  if (nrow(plot_df) > 0) {
    plot_df$stage <- factor(
      plot_df$stage,
      levels = c("preprocess", "graph", "splits", "baseline", "train", "predict")
    )
    p <- ggplot2::ggplot(plot_df,
                         ggplot2::aes(x = n, y = wall_sec, colour = stage)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_y_log10() +
      ggplot2::scale_x_log10(breaks = n_grid) +
      ggplot2::labs(
        title    = "pigauto v0.3.1 scaling (Fix A + Fix B)",
        subtitle = "Wall time per pipeline stage vs tree size (CPU-only)",
        x        = "Number of tips (log)",
        y        = "Wall-clock seconds (log)",
        colour   = "Stage"
      ) +
      ggplot2::theme_minimal(base_size = 12)
    ggplot2::ggsave(out_png, plot = p, width = 10, height = 6, dpi = 120)
    log_line(sprintf("plot written to %s", out_png))
  }
} else {
  log_line("ggplot2 not available or no ok rows, skipping plot")
}

# -------------------------------------------------------------------------
# Markdown summary
# -------------------------------------------------------------------------

machine <- tryCatch(
  sprintf("%s %s (%s), R %s",
          Sys.info()[["sysname"]],
          Sys.info()[["release"]],
          Sys.info()[["machine"]],
          paste(R.version$major, R.version$minor, sep = ".")),
  error = function(e) "machine info unavailable"
)

md_lines <- c(
  "# pigauto scaling diagnostic (Part 4 rerun, post Fix A + B)",
  "",
  sprintf("Run on: %s", format(Sys.time())),
  sprintf("Machine: %s", machine),
  sprintf("Commit:  %s", system("git rev-parse HEAD", intern = TRUE)),
  sprintf("RSpectra available: %s", has_rspectra),
  "",
  "## Per-stage results",
  "",
  "| N | stage | status | wall (s) | R heap max (MB) | error |",
  "|---|---|---|---|---|---|"
)

if (nrow(results) > 0) {
  for (i in seq_len(nrow(results))) {
    r <- results[i, ]
    md_lines <- c(md_lines,
      sprintf("| %d | %s | %s | %s | %s | %s |",
              r$n, r$stage, r$status,
              if (is.na(r$wall_sec)) "-" else sprintf("%.1f", r$wall_sec),
              if (is.na(r$max_mb))   "-" else sprintf("%.0f",  r$max_mb),
              if (is.na(r$error_msg)) "" else substr(r$error_msg, 1, 80))
    )
  }
}

writeLines(md_lines, out_md)
log_line(sprintf("summary written to %s", out_md))
log_line("done")
