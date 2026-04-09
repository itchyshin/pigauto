#!/usr/bin/env Rscript
#
# script/bench_scaling.R
#
# Part 3 of the CRAN warm-up + scaling plan
# (/Users/z3437171/.claude/plans/fizzy-baking-treasure.md).
#
# Purpose
#   Measure wall-clock time and R-heap memory for each pigauto pipeline
#   stage at N = 300, 1000, 2000, 3000, 5000, 7500, 10000 tips on a
#   simulated mixed-type trait dataset, on the current (v0.3.0) code.
#   This is a pre-fix baseline: the results determine which of the Part 4
#   fixes in the plan are mandatory versus optional.
#
# Design
#   * Tree: ape::rcoal(n), seeded by n for reproducibility.
#   * Traits: 3 continuous (BM) + 1 binary (thresholded BM) +
#             1 ordinal (3-level binned BM), 15% MCAR on observed cells.
#   * Stages (sequential, dependent): preprocess, graph, splits,
#             baseline, train (epochs = 100), predict (5 imputations).
#   * Per-stage wall budget: 30 minutes. If a stage errors or exceeds the
#             budget, subsequent stages for that n are skipped.
#   * Overall wall budget: 3 hours. If exceeded, remaining n are
#             recorded as "budget-exhausted" and the script exits.
#   * Checkpointing: results are saved to script/bench_scaling.rds after
#             every stage, so a crash at the last size still leaves
#             usable data for smaller sizes.
#
# Output
#   script/bench_scaling.rds          per-stage data.frame
#   script/bench_scaling.png          wall-clock plot (log y)
#   script/bench_scaling_summary.md   machine info + verdict
#   script/bench_scaling.log          this script's stdout/stderr
#
# Run with
#   cd pigauto && /usr/local/bin/Rscript script/bench_scaling.R

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(pigauto)
  library(ape)
})

# -------------------------------------------------------------------------
# Paths and state
# -------------------------------------------------------------------------

here       <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds    <- file.path(here, "script", "bench_scaling.rds")
out_png    <- file.path(here, "script", "bench_scaling.png")
out_md     <- file.path(here, "script", "bench_scaling_summary.md")

n_grid             <- c(300L, 1000L, 2000L, 3000L, 5000L, 7500L, 10000L)
stage_budget_sec   <- 30 * 60      # 30 minutes per stage
overall_budget_sec <- 3 * 3600     # 3 hours total

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
  # gc() returns a matrix whose columns are c("used", "(Mb)", "gc trigger",
  # "(Mb)", "limit (Mb)", "max used", "(Mb)"). Column 7 is max-used in Mb.
  # Duplicate "(Mb)" header names make name-based indexing unreliable, so
  # we take column 7 of the Vcells row positionally.
  mb <- tryCatch(as.numeric(gc_after["Vcells", 7L]),
                 error = function(e) NA_real_)

  append_result(n, stage_label, wall, mb, status, error_msg)

  log_line(sprintf("[n=%5d] %-10s %-7s wall=%7.1fs  maxR=%6.0f MB%s",
                   n, stage_label, status, wall, mb,
                   if (!is.na(error_msg)) paste0("  -- ", substr(error_msg, 1, 120)) else ""))

  list(status = status, value = value, wall = wall)
}

# -------------------------------------------------------------------------
# Mixed-type dataset simulator
# -------------------------------------------------------------------------

simulate_dataset <- function(n) {
  set.seed(n)
  tree <- ape::rcoal(n)

  # 3 continuous BM traits plus 2 latent BM traits for binary / ordinal.
  x1 <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  x2 <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  x3 <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  z1 <- ape::rTraitCont(tree, model = "BM", sigma = 1)
  z2 <- ape::rTraitCont(tree, model = "BM", sigma = 1)

  # Binary via sign of z1 centered on zero; ordinal via 3-quantile binning.
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

  # 15% MCAR on observed cells.
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

  # --- simulate ---
  sim <- tryCatch(simulate_dataset(n), error = function(e) {
    append_result(n, "simulate", NA_real_, NA_real_, "error", conditionMessage(e))
    NULL
  })
  if (is.null(sim)) return(invisible(NULL))
  traits <- sim$traits
  tree   <- sim$tree

  # --- preprocess ---
  st <- time_stage(n, "preprocess", function() {
    preprocess_traits(traits, tree, log_transform = FALSE)
  })
  if (st$status != "ok") return(invisible(NULL))
  pd <- st$value

  # --- graph ---
  st <- time_stage(n, "graph", function() {
    build_phylo_graph(tree, k_eigen = 8L, sigma_mult = 0.5)
  })
  if (st$status != "ok") return(invisible(NULL))
  graph <- st$value

  # --- splits ---
  st <- time_stage(n, "splits", function() {
    make_missing_splits(pd$X_scaled, missing_frac = 0.25,
                        val_frac = 0.25, seed = 42L,
                        trait_map = pd$trait_map)
  })
  if (st$status != "ok") return(invisible(NULL))
  splits <- st$value

  # --- baseline ---
  st <- time_stage(n, "baseline", function() {
    fit_baseline(pd, tree, splits = splits, model = "BM")
  })
  if (st$status != "ok") return(invisible(NULL))
  baseline <- st$value

  # --- train ---
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

  # --- predict ---
  st <- time_stage(n, "predict", function() {
    stats::predict(fit, return_se = TRUE, n_imputations = 5L)
  })
  invisible(NULL)
}

# -------------------------------------------------------------------------
# Main loop
# -------------------------------------------------------------------------

log_line(sprintf("bench_scaling starting at %s", format(Sys.time())))
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
        title    = "pigauto v0.3.0 scaling diagnostic",
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
# Markdown summary + verdict
# -------------------------------------------------------------------------

machine <- tryCatch(
  sprintf("%s %s (%s), R %s",
          Sys.info()[["sysname"]],
          Sys.info()[["release"]],
          Sys.info()[["machine"]],
          paste(R.version$major, R.version$minor, sep = ".")),
  error = function(e) "machine info unavailable"
)

# Identify the first n at which *any* stage took more than 5 minutes OR
# errored / timed out. That n is the scaling wall.
results$blocking <- (results$status != "ok") |
  (results$wall_sec > 300 & !is.na(results$wall_sec))
blockers <- results[results$blocking, , drop = FALSE]

verdict <- if (nrow(blockers) == 0) {
  "All stages completed within 5 minutes at every N -- no scaling wall observed in this run."
} else {
  first_n    <- min(blockers$n)
  first_rows <- blockers[blockers$n == first_n, , drop = FALSE]
  stages     <- paste(first_rows$stage, collapse = ", ")
  sprintf("First scaling wall at n = %d in stage(s): %s.", first_n, stages)
}

# Fix A vs Fix D verdict heuristic: if the `graph` stage is the
# earliest-blocking or slowest, Fix A (sparse eigensolver) is the primary
# win. If `baseline` (Rphylopars) is also blocking at the same or an
# earlier n, Fix D (fast Felsenstein) becomes mandatory too.
graph_block    <- any(results$blocking & results$stage == "graph")
baseline_block <- any(results$blocking & results$stage == "baseline")

fix_verdict <- if (graph_block && baseline_block) {
  "Fix A (sparse eigensolver) AND Fix D (fast Felsenstein baseline) are mandatory; both graph and baseline stages block."
} else if (graph_block) {
  "Fix A (sparse eigensolver) is sufficient to unblock the current wall; Fix D can be deferred unless baseline becomes the next bottleneck."
} else if (baseline_block) {
  "Fix D (fast Felsenstein baseline) is the primary bottleneck; Fix A may still be needed at larger N but is secondary here."
} else if (nrow(blockers) == 0) {
  "No fixes strictly required from this run; pigauto already reaches N = 10,000 within acceptable budgets."
} else {
  "Bottleneck is in neither graph nor baseline -- inspect train/predict stages before committing to Fix A or Fix D."
}

md_lines <- c(
  "# pigauto scaling diagnostic (Part 3 baseline, tag v0.3.0)",
  "",
  sprintf("Run on: %s", format(Sys.time())),
  sprintf("Machine: %s", machine),
  sprintf("Commit:  %s (tag v0.3.0)", system("git rev-parse HEAD", intern = TRUE)),
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

md_lines <- c(md_lines, "",
              "## Wall of first failure", "", verdict, "",
              "## Verdict on Part 4 fixes", "", fix_verdict, "")

writeLines(md_lines, out_md)
log_line(sprintf("summary written to %s", out_md))
log_line("done")
