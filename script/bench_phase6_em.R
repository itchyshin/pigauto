#!/usr/bin/env Rscript
#
# script/bench_phase6_em.R
#
# Smoke benchmark: plug-in baseline vs Phase 6 EM (diagonal) vs
# Phase 7 EM (off-diagonal) on a small correlated-binary simulation.
# Reports per-method wall time, Phase 6 iterations, whether Phase 7
# converged, and per-trait held-out accuracy.
#
# Run separately from devtools::test():
#   Rscript script/bench_phase6_em.R

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  devtools::load_all("/Users/z3437171/Dropbox/Github Local/pigauto",
                      quiet = TRUE)
})

out_rds <- "script/bench_phase6_em.rds"
out_md  <- "script/bench_phase6_em.md"

# ---- Simulate correlated binary + continuous data -----------------------
# Shared continuous BM drives two binary traits. Cross-trait liability
# correlation is strong (~0.7); this is the regime where off-diagonal
# conditioning should move the needle.
n <- 150L
set.seed(2026L)
tree <- ape::rcoal(n, tip.label = paste0("sp", seq_len(n)))
z <- ape::rTraitCont(tree, model = "BM")
df <- data.frame(
  b1 = factor(ifelse(z + stats::rnorm(n, sd = 0.3) > 0, "1", "0"),
               levels = c("0", "1")),
  b2 = factor(ifelse(z + stats::rnorm(n, sd = 0.3) > 0, "1", "0"),
               levels = c("0", "1")),
  c1 = z + stats::rnorm(n, sd = 0.5),
  row.names = tree$tip.label
)
# Hold out 30% at random on b1 and b2 for held-out accuracy scoring
held_b1 <- sample.int(n, 0.30 * n)
held_b2 <- sample.int(n, 0.30 * n)
truth_b1 <- df$b1
truth_b2 <- df$b2
df$b1[held_b1] <- NA
df$b2[held_b2] <- NA

pd <- preprocess_traits(df, tree, log_transform = FALSE)
sp <- make_missing_splits(pd$X_scaled, missing_frac = 0.25,
                           seed = 1L, trait_map = pd$trait_map)
gr <- build_phylo_graph(tree, k_eigen = 4L)

score_accuracy <- function(bl, trait_col, held_idx, truth) {
  # bl$mu for binary trait is logit(P(y=1)). Threshold at 0 = 0.5.
  col_idx <- which(colnames(pd$X_scaled) == trait_col)
  mu_col  <- bl$mu[, col_idx]
  pred    <- factor(ifelse(mu_col > 0, "1", "0"), levels = c("0", "1"))
  mean(pred[held_idx] == truth[held_idx], na.rm = TRUE)
}

timed <- function(expr) {
  t0 <- proc.time()
  val <- force(expr)
  list(val = val, wall = (proc.time() - t0)[["elapsed"]])
}

cat("Running plug-in baseline (em_iterations = 0) ...\n")
bl0 <- timed(fit_baseline(pd, tree, splits = sp, graph = gr,
                            em_iterations = 0L))
cat("Running Phase 6 EM (em_iterations = 5, em_offdiag = FALSE) ...\n")
bl6 <- timed(fit_baseline(pd, tree, splits = sp, graph = gr,
                            em_iterations = 5L, em_offdiag = FALSE))
cat("Running Phase 7 EM (em_iterations = 5, em_offdiag = TRUE) ...\n")
bl7 <- timed(fit_baseline(pd, tree, splits = sp, graph = gr,
                            em_iterations = 5L, em_offdiag = TRUE))

results <- data.frame(
  method = c("plug_in", "phase_6_diag", "phase_7_offdiag"),
  wall_s = c(bl0$wall, bl6$wall, bl7$wall),
  acc_b1 = c(score_accuracy(bl0$val, "b1", held_b1, truth_b1),
             score_accuracy(bl6$val, "b1", held_b1, truth_b1),
             score_accuracy(bl7$val, "b1", held_b1, truth_b1)),
  acc_b2 = c(score_accuracy(bl0$val, "b2", held_b2, truth_b2),
             score_accuracy(bl6$val, "b2", held_b2, truth_b2),
             score_accuracy(bl7$val, "b2", held_b2, truth_b2))
)

saveRDS(list(results = results, n = n), out_rds)

md <- c(
  "# Phase 6 + Phase 7 EM smoke benchmark",
  "",
  sprintf("n = %d, correlated binary (shared BM liability) + 1 continuous.", n),
  "Held out 30% of each binary trait for accuracy scoring.",
  "",
  "| method | wall (s) | acc b1 | acc b2 |",
  "|---|---:|---:|---:|",
  apply(results, 1, function(r) sprintf(
    "| %s | %.2f | %.3f | %.3f |",
    r[["method"]], as.numeric(r[["wall_s"]]),
    as.numeric(r[["acc_b1"]]), as.numeric(r[["acc_b2"]])
  )),
  ""
)
writeLines(md, out_md)
cat("Wrote", out_md, "and", out_rds, "\n")
print(results)
