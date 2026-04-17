#!/usr/bin/env Rscript
#
# script/bench_scaling_v090.R
#
# Scaling curve for pigauto v0.9.0 (Level-C baseline + B1 soft E-step +
# B2 rate-aware bandwidth + B3 ordinal threshold + Graph Transformer).
# Measures per-stage wall time, peak R heap, and held-out accuracy as a
# function of tree size, using a mixed-type dataset (2 continuous + 1
# binary + 1 K=4 categorical).
#
# Design
#   * Tree sizes: n in {100, 300, 1000, 3000, 5000, 7500, 10000}
#   * Traits per n:
#       - cont1, cont2  : Brownian motion, sigma = 1
#       - bin           : 2-level factor, thresholded BM
#       - cat           : K = 4 factor, argmax of 4 independent BM liabilities
#   * Missingness: 25% MCAR on observed cells
#   * Epochs: 300 for n<=300, 200 for n<=3000, 100 for n<=10000
#   * Single rep per n (the curve is about runtime, not sampling precision)
#   * Stage timing: preprocess, graph, splits, baseline, train, predict
#   * Stage budget: 45 min; overall budget: 6 h
#   * Checkpoint after every stage so a timeout still leaves a usable RDS.
#
# Output
#   script/bench_scaling_v090.rds    per-stage data.frame
#   script/bench_scaling_v090.png    wall-clock vs n (log-log)
#   script/bench_scaling_v090.md     machine info + table + commentary
#
# Run with
#   cd pigauto && Rscript script/bench_scaling_v090.R

options(warn = 1, stringsAsFactors = FALSE)

here <- Sys.getenv("PIGAUTO_ROOT", unset = NA_character_)
if (is.na(here) || !dir.exists(here)) {
  cands <- c("/Users/z3437171/Dropbox/Github Local/pigauto",
             "/workspaces/pigauto",
             ".")
  here <- Find(function(p) file.exists(file.path(p, "DESCRIPTION")), cands)
  if (is.null(here)) stop("Cannot find pigauto project root")
}

suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(here, quiet = TRUE)
})

out_rds <- file.path(here, "script", "bench_scaling_v090.rds")
out_png <- file.path(here, "script", "bench_scaling_v090.png")
out_md  <- file.path(here, "script", "bench_scaling_v090.md")

n_grid             <- c(100L, 300L, 1000L, 3000L, 5000L, 7500L, 10000L)
stage_budget_sec   <- 45 * 60
overall_budget_sec <- 6 * 3600

epochs_for <- function(n) {
  if (n <= 300L) return(300L)
  if (n <= 3000L) return(200L)
  100L
}

has_rspectra <- requireNamespace("RSpectra", quietly = TRUE)

script_start <- proc.time()[["elapsed"]]

log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start), ..., "\n",
      sep = "")
  flush.console()
}

# -------------------------------------------------------------------------
# Results accumulator with checkpointing
# -------------------------------------------------------------------------

results <- data.frame(
  n         = integer(),
  stage     = character(),
  wall_sec  = numeric(),
  max_mb    = numeric(),
  status    = character(),
  error_msg = character(),
  stringsAsFactors = FALSE
)

accuracy <- data.frame(
  n         = integer(),
  trait     = character(),
  metric    = character(),
  value     = numeric(),
  stringsAsFactors = FALSE
)

checkpoint <- function() {
  saveRDS(list(
    results  = results,
    accuracy = accuracy,
    n_grid   = n_grid,
    commit   = tryCatch(system("git rev-parse HEAD", intern = TRUE),
                        error = function(e) "unknown"),
    session  = utils::sessionInfo()
  ), out_rds)
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
  checkpoint()
}

append_accuracy <- function(n, trait, metric, value) {
  row <- data.frame(
    n      = n, trait = trait, metric = metric, value = value,
    stringsAsFactors = FALSE
  )
  accuracy <<- rbind(accuracy, row)
  checkpoint()
}

# -------------------------------------------------------------------------
# Per-stage timer
# -------------------------------------------------------------------------

time_stage <- function(n, stage_label, expr_fun) {
  gc(reset = TRUE, full = TRUE)
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
# Mixed-type data simulator
# -------------------------------------------------------------------------

simulate_mixed <- function(n) {
  set.seed(n)
  tree <- ape::rcoal(n)

  # Two continuous BM traits
  x1 <- as.numeric(ape::rTraitCont(tree, model = "BM", sigma = 1))
  x2 <- as.numeric(ape::rTraitCont(tree, model = "BM", sigma = 1))

  # Binary via threshold on latent BM
  z_bin <- as.numeric(ape::rTraitCont(tree, model = "BM", sigma = 1))
  bin <- factor(ifelse(z_bin > stats::median(z_bin), "high", "low"),
                levels = c("low", "high"))

  # K = 4 categorical via argmax of 4 independent BM liabilities
  K <- 4L
  liab <- matrix(NA_real_, nrow = n, ncol = K)
  for (k in seq_len(K)) {
    liab[, k] <- as.numeric(ape::rTraitCont(tree, model = "BM", sigma = 1))
  }
  rownames(liab) <- tree$tip.label
  cat_idx <- apply(liab, 1, which.max)
  cat_lev <- c("A", "B", "C", "D")
  cat <- factor(cat_lev[cat_idx], levels = cat_lev)

  traits <- data.frame(
    cont1 = x1,
    cont2 = x2,
    bin   = bin,
    cat   = cat,
    row.names = tree$tip.label,
    stringsAsFactors = FALSE
  )

  # Save the complete ground truth for later accuracy scoring
  truth <- traits

  # Apply 25% MCAR missingness
  set.seed(1000L + n)
  mask <- matrix(stats::runif(nrow(traits) * ncol(traits)) < 0.25,
                 nrow = nrow(traits), ncol = ncol(traits))
  for (j in seq_len(ncol(traits))) {
    traits[mask[, j], j] <- NA
  }

  list(traits = traits, tree = tree, truth = truth, mask = mask)
}

# -------------------------------------------------------------------------
# Accuracy scoring on held-out MCAR cells
# -------------------------------------------------------------------------

score_predictions <- function(n, completed, truth, mask) {
  for (col in colnames(truth)) {
    j <- match(col, colnames(truth))
    held <- which(mask[, j])
    if (!length(held)) next
    pred_vals <- completed[[col]][held]
    true_vals <- truth[[col]][held]
    if (is.factor(true_vals)) {
      # binary or categorical: accuracy
      ok <- sum(as.character(pred_vals) == as.character(true_vals), na.rm = TRUE)
      acc <- ok / length(held)
      append_accuracy(n, col, "accuracy", acc)
    } else {
      resid <- true_vals - pred_vals
      rmse  <- sqrt(mean(resid^2, na.rm = TRUE))
      append_accuracy(n, col, "rmse", rmse)
    }
  }
}

# -------------------------------------------------------------------------
# Per-size runner
# -------------------------------------------------------------------------

run_one <- function(n) {
  log_line(sprintf("===== n = %d =====", n))

  sim <- tryCatch(simulate_mixed(n), error = function(e) {
    append_result(n, "simulate", NA_real_, NA_real_, "error",
                  conditionMessage(e))
    NULL
  })
  if (is.null(sim)) return(invisible(NULL))
  traits <- sim$traits
  tree   <- sim$tree
  truth  <- sim$truth
  mask   <- sim$mask

  st <- time_stage(n, "preprocess", function() {
    preprocess_traits(traits, tree, log_transform = FALSE)
  })
  if (st$status != "ok") return(invisible(NULL))
  pd <- st$value

  st <- time_stage(n, "graph", function() {
    build_phylo_graph(tree, k_eigen = "auto")
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

  st <- time_stage(n, "baseline", function() {
    fit_baseline(pd, tree, splits = splits, graph = graph)
  })
  if (st$status != "ok") return(invisible(NULL))
  baseline <- st$value

  # Drop D to save memory during training, as in bench_scaling_v031
  graph$D <- NULL
  invisible(gc(full = TRUE, verbose = FALSE))

  this_epochs <- epochs_for(n)
  log_line(sprintf("  epochs for this n: %d", this_epochs))

  st <- time_stage(n, "train", function() {
    fit_pigauto(
      data     = pd,
      tree     = tree,
      splits   = splits,
      graph    = graph,
      baseline = baseline,
      epochs   = this_epochs,
      verbose  = FALSE,
      seed     = 1L
    )
  })
  if (st$status != "ok") return(invisible(NULL))
  fit <- st$value

  st <- time_stage(n, "predict", function() {
    stats::predict(fit, return_se = TRUE, n_imputations = 1L)
  })
  if (st$status != "ok") return(invisible(NULL))
  pred <- st$value

  # Score on held-out MCAR cells. We already have `pred` from the predict
  # stage; reconstruct the completed frame via the internal build_completed
  # helper (same function impute() uses under the hood).
  completed <- tryCatch({
    info <- pigauto:::build_completed(traits, pred$imputed, NULL)
    info$completed
  }, error = function(e) {
    log_line(sprintf("  build_completed failed: %s", conditionMessage(e)))
    NULL
  })
  if (!is.null(completed)) {
    score_predictions(n, completed, truth, mask)
  }

  invisible(NULL)
}

# -------------------------------------------------------------------------
# Main loop
# -------------------------------------------------------------------------

log_line(sprintf("bench_scaling_v090 starting at %s", format(Sys.time())))
log_line(sprintf("RSpectra available: %s", has_rspectra))
log_line(sprintf("n_grid = %s", paste(n_grid, collapse = ", ")))

for (n in n_grid) {
  elapsed <- proc.time()[["elapsed"]] - script_start
  if (elapsed > overall_budget_sec) {
    log_line(sprintf("overall budget exhausted at %.0fs, skipping n >= %d",
                     elapsed, n))
    append_result(n, "budget", elapsed, NA_real_, "budget-exhausted")
    break
  }
  run_one(n)
}

log_line(sprintf("loop finished at %s (total %.0fs)",
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
        title    = "pigauto v0.9.0 scaling (Level-C + B1/B2/B3 + GT)",
        subtitle = "Wall time per stage vs tree size (mixed types)",
        x        = "Number of tips (log)",
        y        = "Wall-clock seconds (log)",
        colour   = "Stage"
      ) +
      ggplot2::theme_minimal(base_size = 12)
    ggplot2::ggsave(out_png, plot = p, width = 10, height = 6, dpi = 120)
    log_line(sprintf("plot written to %s", out_png))
  }
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
  "# pigauto v0.9.0 scaling curve",
  "",
  sprintf("Run on: %s", format(Sys.time())),
  sprintf("Machine: %s", machine),
  sprintf("Commit:  %s", tryCatch(system("git rev-parse HEAD", intern = TRUE),
                                  error = function(e) "unknown")),
  sprintf("RSpectra available: %s", has_rspectra),
  "",
  "## Workload per n",
  "",
  "- 2 continuous BM traits",
  "- 1 binary trait (BM threshold)",
  "- 1 categorical trait (K = 4, argmax of independent BM liabilities)",
  "- 25% MCAR missingness",
  sprintf("- Epochs: 300 for n<=300, 200 for n<=3000, 100 for n<=10000"),
  "",
  "## Per-stage wall time",
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

md_lines <- c(md_lines, "", "## Accuracy on held-out MCAR cells",
              "",
              "| N | trait | metric | value |",
              "|---|---|---|---|")
if (nrow(accuracy) > 0) {
  for (i in seq_len(nrow(accuracy))) {
    a <- accuracy[i, ]
    md_lines <- c(md_lines,
      sprintf("| %d | %s | %s | %s |",
              a$n, a$trait, a$metric,
              if (is.na(a$value)) "-" else sprintf("%.4f", a$value)))
  }
}

writeLines(md_lines, out_md)
log_line(sprintf("summary written to %s", out_md))
log_line("done")
