#!/usr/bin/env Rscript
#
# script/bench_binary_joint.R
#
# A/B benchmark for Phase 3 (threshold-joint baseline on binary traits).
# Scenarios sweep phylogenetic signal on a correlated continuous+binary
# pair simulated by joint BM on the tree.

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/level-c-phase-3",
    quiet = TRUE
  )
})

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/level-c-phase-3"
out_md  <- file.path(here, "script", "bench_binary_joint.md")

n_species <- 200L
reps      <- 3L
signals   <- c(0.3, 0.6, 0.9)

simulate_correlated_bm <- function(tree, rho, seed) {
  set.seed(seed)
  n <- length(tree$tip.label)
  # Joint BM on the tree for 2 correlated traits
  V   <- ape::vcv(tree)
  L   <- chol(V)
  z1  <- as.numeric(t(L) %*% rnorm(n))
  z2  <- as.numeric(t(L) %*% rnorm(n))
  x1  <- z1
  x2  <- rho * z1 + sqrt(1 - rho^2) * z2
  data.frame(
    x     = x1,
    y     = factor(ifelse(x2 > 0, "B", "A"), levels = c("A", "B")),
    row.names = tree$tip.label,
    stringsAsFactors = FALSE
  )
}

run_one <- function(signal, rep_id) {
  set.seed(rep_id * 1000L + round(signal * 100))
  tree <- ape::rtree(n_species)
  # signal acts as the correlation rho between continuous and the
  # binary-generating latent: low rho = LP should tie joint; high rho =
  # joint should win.
  df   <- simulate_correlated_bm(tree, rho = signal,
                                  seed = rep_id * 1000L + round(signal * 100))
  df$y[sample(n_species, n_species %/% 3)] <- NA  # 33% binary missing

  pd     <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.25,
                                 seed = rep_id, trait_map = pd$trait_map)

  # Joint path (default when Rphylopars available)
  bl_joint <- fit_baseline(pd, tree, splits = splits)

  # LP path (force fallback)
  orig <- pigauto:::joint_mvn_available
  assignInNamespace("joint_mvn_available", function() FALSE, ns = "pigauto")
  on.exit(assignInNamespace("joint_mvn_available", orig, ns = "pigauto"),
          add = TRUE)
  bl_lp <- fit_baseline(pd, tree, splits = splits)
  assignInNamespace("joint_mvn_available", orig, ns = "pigauto")

  # Evaluate binary accuracy on TEST split (cells held out by splits)
  test_idx <- splits$test_idx
  y_col    <- which(colnames(pd$X_scaled) == "y")
  test_y   <- pd$X_scaled[, y_col]
  # test_idx is linear over (n, p); decompose
  n <- nrow(pd$X_scaled); p <- ncol(pd$X_scaled)
  row_i <- ((test_idx - 1L) %% n) + 1L
  col_j <- ((test_idx - 1L) %/% n) + 1L
  y_rows <- row_i[col_j == y_col]
  if (length(y_rows) == 0L) return(NULL)
  truth <- pd$X_scaled[y_rows, y_col]

  acc_joint <- mean(as.integer(plogis(bl_joint$mu[y_rows, y_col]) > 0.5) == truth)
  acc_lp    <- mean(as.integer(plogis(bl_lp$mu[y_rows, y_col])    > 0.5) == truth)

  data.frame(signal = signal, rep = rep_id,
             acc_joint = acc_joint, acc_lp = acc_lp,
             lift = acc_joint - acc_lp)
}

all_results <- do.call(rbind, lapply(signals, function(s) {
  do.call(rbind, lapply(seq_len(reps), function(r) run_one(s, r)))
}))

agg <- aggregate(cbind(acc_joint, acc_lp, lift) ~ signal,
                  data = all_results, FUN = mean)

md <- c("# Phase 3 A/B: threshold-joint vs LP baseline on binary traits",
        "",
        sprintf("Run on: %s", format(Sys.time())),
        sprintf("Species per scenario: %d, reps: %d", n_species, reps),
        "",
        "Accuracy on held-out binary test cells:",
        "",
        "```",
        capture.output(print(agg, row.names = FALSE)),
        "```",
        "",
        sprintf("Scenarios with lift >= 2pp: %d / %d",
                sum(agg$lift >= 0.02), nrow(agg)))

writeLines(md, out_md)
cat("\n--- Summary ---\n")
cat(paste(md, collapse = "\n"), "\n")
