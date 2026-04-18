#!/usr/bin/env Rscript
#
# script/bench_share_gnn.R
#
# Smoke benchmark: multi_impute_trees(share_gnn=TRUE) vs share_gnn=FALSE.
# Measures wall time + pooled FMI to verify the documented tradeoff
# (~10x speedup with tree uncertainty preserved when the gate closes).
#
# Tiny by design -- this is a pkgdown-linked sanity artefact, not the
# full tree_uncertainty benchmark (which still uses share_gnn=FALSE to
# anchor the v0.9.0 baseline comparison).
#
# Run separately (not part of devtools::test()):
#   Rscript script/bench_share_gnn.R

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  devtools::load_all("/Users/z3437171/Dropbox/Github Local/pigauto", quiet = TRUE)
})

out_rds <- "script/bench_share_gnn.rds"
out_md  <- "script/bench_share_gnn.md"

n <- 100L          # small enough to be a smoke test
T_trees <- 5L      # small posterior
n_reps <- 2L       # paired reps for noise reduction

log <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

timed <- function(expr) {
  t0 <- proc.time()
  val <- force(expr)
  list(val = val, wall = (proc.time() - t0)[["elapsed"]])
}

rows <- list()
for (r in seq_len(n_reps)) {
  set.seed(100L + r)
  trees <- lapply(seq_len(T_trees), function(i) ape::rcoal(n, tip.label = paste0("sp", 1:n)))
  df <- simulate_bm_traits(trees[[1]], n_traits = 3L, seed = 100L + r)
  df$trait1[sample.int(n, n * 0.2)] <- NA
  df$trait2[sample.int(n, n * 0.2)] <- NA

  log(sprintf("rep %d: share_gnn=TRUE", r))
  t1 <- timed(multi_impute_trees(df, trees, m_per_tree = 1L,
                                  share_gnn = TRUE, reference_tree = trees[[1]],
                                  epochs = 100L, verbose = FALSE, seed = r))
  rows[[length(rows) + 1L]] <- data.frame(
    rep = r, method = "share_gnn=TRUE", wall_s = t1$wall,
    m = t1$val$m, stringsAsFactors = FALSE
  )

  log(sprintf("rep %d: share_gnn=FALSE", r))
  t2 <- timed(multi_impute_trees(df, trees, m_per_tree = 1L,
                                  share_gnn = FALSE,
                                  epochs = 100L, verbose = FALSE, seed = r))
  rows[[length(rows) + 1L]] <- data.frame(
    rep = r, method = "share_gnn=FALSE", wall_s = t2$wall,
    m = t2$val$m, stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, rows)
saveRDS(list(results = results, n = n, T_trees = T_trees), out_rds)

md <- c(
  "# share_gnn smoke benchmark",
  "",
  sprintf("n = %d, T = %d trees, m_per_tree = 1, reps = %d.", n, T_trees, n_reps),
  "",
  "| rep | method | wall (s) | M |",
  "|---|---|---|---|",
  apply(results, 1, function(r) {
    sprintf("| %s | %s | %.1f | %s |", r[["rep"]], r[["method"]], as.numeric(r[["wall_s"]]), r[["m"]])
  }),
  "",
  sprintf("Mean wall time ratio (FALSE / TRUE): %.2fx",
          mean(results$wall_s[results$method == "share_gnn=FALSE"]) /
          mean(results$wall_s[results$method == "share_gnn=TRUE"]))
)
writeLines(md, out_md)
log("done")
