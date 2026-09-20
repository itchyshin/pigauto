#!/usr/bin/env Rscript
# Summarise a pre-run: per-arm walls by n, BACE convergence, failures, and the re-derived budget.
#
#   Rscript script/campaign_sim_prerun_summary.R --dirs prerun_fast,prerun_bace [--out summary.md]
#
# The budget it prints is the number Shinichi approves at G0: measured walls x the design's cell
# and replicate counts, not the planning estimate.

args <- commandArgs(trailingOnly = TRUE)
get <- function(f, d = NULL) { i <- match(f, args); if (is.na(i)) d else args[i + 1L] }
dirs <- strsplit(get("--dirs", "prerun_fast,prerun_bace"), ",")[[1]]
out_md <- get("--out", NULL)
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))

`%||%` <- function(a, b) if (is.null(a)) b else a

fs <- unlist(lapply(dirs, function(d) list.files(d, pattern = "[.]rds$", full.names = TRUE)))
if (!length(fs)) stop("no rds under: ", paste(dirs, collapse = ", "))

rows <- list(); conv <- list(); fails <- list()
for (f in fs) {
  x <- readRDS(f)
  w <- x$walls
  if (length(w)) rows[[length(rows) + 1L]] <- data.frame(
    cell = x$tag, n = x$n, lambda = x$lambda, rho = x$rho, evo = x$evo, miss = x$miss,
    frac = x$miss_frac, arm = names(w), wall = as.numeric(w), stringsAsFactors = FALSE)
  if (length(x$failed)) fails[[length(fails) + 1L]] <- data.frame(
    cell = x$tag, arm = names(x$failed), error = vapply(names(x$failed),
      function(a) substr(x$errors[[a]] %||% "", 1, 90), ""), stringsAsFactors = FALSE)
  d <- x$diag$bace
  if (!is.null(d)) conv[[length(conv) + 1L]] <- data.frame(
    cell = x$tag, n = x$n, converged = isTRUE(d$converged), attempts = d$n_attempts %||% NA,
    drift = d$drift %||% NA_real_, ess_med = d$ess_med %||% NA_real_, ess_min = d$ess_min %||% NA_real_,
    stringsAsFactors = FALSE)
}
W <- do.call(rbind, rows); FL <- do.call(rbind, fails); C <- do.call(rbind, conv)

cat("## Walls per arm by n (seconds, 4 threads, one replicate)\n\n")
tab <- tapply(W$wall, list(W$arm, W$n), function(z) round(mean(z, na.rm = TRUE), 1))
print(tab)

if (!is.null(C)) { cat("\n## BACE convergence\n\n"); print(C, row.names = FALSE) }
cat("\n## Failures\n\n")
if (is.null(FL)) cat("none\n") else print(FL, row.names = FALSE)

# ---- re-derived budget ------------------------------------------------------------------------
# Mean wall per arm per n from this pre-run, applied to the design's cells and replicates.
design <- function(stage) {
  p <- file.path(script_dir, "campaign_sim_design.R")
  read.csv(text = paste(system2("Rscript", c(p, "--stage", stage), stdout = TRUE), collapse = "\n"),
           stringsAsFactors = FALSE)
}
arm_wall <- function(a, n) {
  z <- W$wall[W$arm == a & W$n == n]
  if (length(z)) return(mean(z, na.rm = TRUE))
  z2 <- W$wall[W$arm == a]                     # fall back: scale the nearest measured n
  if (!length(z2)) return(NA_real_)
  ns <- unique(W$n[W$arm == a]); nn <- ns[which.min(abs(ns - n))]
  mean(W$wall[W$arm == a & W$n == nn], na.rm = TRUE) * (n / nn)^2
}
arms_fast <- c("gnn_on", "gnn_off", "gnn_off_rphylopars", "freq", "floor")
cat("\n## Re-derived budget (4-thread slot-hours; measured walls x design)\n\n")
total <- 0
for (stage in c("core", "factorial")) {
  d <- design(stage)
  s <- 0
  for (i in seq_len(nrow(d))) {
    per_fast <- sum(vapply(arms_fast, arm_wall, 0, n = d$n[i]), na.rm = TRUE)
    per_bace <- arm_wall("bace", d$n[i])
    s <- s + d$reps[i] * per_fast + d$bace_reps[i] * ifelse(is.na(per_bace), 0, per_bace)
  }
  cat(sprintf("  %-10s %3d cells -> %8.0f slot-hours\n", stage, nrow(d), s / 3600))
  total <- total + s / 3600
}
cat(sprintf("  %-10s %19s %8.0f slot-hours\n", "TOTAL", "", total))
cat(sprintf("\n  Totoro at 36 slots: %.1f h if it all ran there.\n", total / 36))
cat(sprintf("  Core on Totoro (36 slots) + factorial split over nibi and rorqual (400 slots each):\n"))
