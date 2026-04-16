#!/usr/bin/env Rscript
#
# script/bench_categorical_joint.R
#
# Phase 4 A/B benchmark: joint_K vs OVR vs LP baseline for categorical
# traits. Simulates continuous anchor + binary anchor + K-class categorical,
# all jointly Brownian-evolved with cross-trait correlation rho. Binary
# anchor is required to fire the threshold-joint path (Phase 4 scope
# restriction).

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/level-c-phase-4",
    quiet = TRUE
  )
})

here   <- "/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/level-c-phase-4"
out_md <- file.path(here, "script", "bench_categorical_joint.md")

n_species <- 100L
reps      <- 3L
Ks        <- c(3L, 5L)                       # drop K=8 (multi-cat pushes phylopars past numerical limit)
rho       <- 0.4                              # moderate correlation keeps covariance stable

simulate_correlated_cat <- function(tree, K, rho, seed) {
  set.seed(seed)
  n <- length(tree$tip.label)
  V <- ape::vcv(tree)
  L <- chol(V)
  z_anchor_x <- as.numeric(t(L) %*% rnorm(n))      # continuous anchor
  z_anchor_y <- as.numeric(t(L) %*% rnorm(n))      # binary-generating latent
  liab <- matrix(NA_real_, nrow = n, ncol = K)
  for (k in seq_len(K)) {
    z_k <- as.numeric(t(L) %*% rnorm(n))
    # Correlate the categorical liability with BOTH anchors
    liab[, k] <- rho * z_anchor_x + sqrt(1 - rho^2) * z_k
  }
  class_vec <- apply(liab, 1, which.max)
  data.frame(
    x = z_anchor_x,
    y = factor(ifelse(z_anchor_y > 0, "B", "A"), levels = c("A", "B")),
    z = factor(LETTERS[class_vec], levels = LETTERS[seq_len(K)]),
    row.names = tree$tip.label,
    stringsAsFactors = FALSE
  )
}

force_lp <- function(fn) {
  orig <- pigauto:::joint_mvn_available
  assignInNamespace("joint_mvn_available", function() FALSE, ns = "pigauto")
  on.exit(assignInNamespace("joint_mvn_available", orig, ns = "pigauto"),
          add = TRUE)
  res <- fn()
  assignInNamespace("joint_mvn_available", orig, ns = "pigauto")
  res
}

run_one <- function(K, rep_id) {
  seed <- rep_id * 1000L + K
  set.seed(seed)
  tree <- ape::rtree(n_species)
  df   <- simulate_correlated_cat(tree, K = K, rho = rho, seed = seed)
  df$z[sample(n_species, n_species %/% 3)] <- NA

  pd     <- preprocess_traits(df, tree)
  splits <- make_missing_splits(pd$X_scaled, missing_frac = 0.25,
                                 seed = rep_id, trait_map = pd$trait_map)

  bl_joint_K <- fit_baseline(pd, tree, splits = splits, cat_encoding = "joint_K")
  bl_ovr     <- fit_baseline(pd, tree, splits = splits, cat_encoding = "ovr")
  bl_lp      <- force_lp(function() fit_baseline(pd, tree, splits = splits,
                                                   cat_encoding = "joint_K"))

  cat_col_names <- grep("^z=", colnames(pd$X_scaled), value = TRUE)
  cat_col_idx   <- which(colnames(pd$X_scaled) %in% cat_col_names)
  n <- nrow(pd$X_scaled); p <- ncol(pd$X_scaled)
  hold_row_i <- ((splits$test_idx - 1L) %% n) + 1L
  hold_col_j <- ((splits$test_idx - 1L) %/% n) + 1L
  hold_rows <- unique(hold_row_i[hold_col_j %in% cat_col_idx])
  if (length(hold_rows) == 0L) return(NULL)

  truth <- apply(pd$X_scaled[hold_rows, cat_col_names, drop = FALSE], 1, which.max)
  acc <- function(bl) mean(apply(bl$mu[hold_rows, cat_col_names, drop = FALSE],
                                   1, which.max) == truth)

  data.frame(
    K          = K,
    rep        = rep_id,
    n_test     = length(hold_rows),
    acc_joint_K = acc(bl_joint_K),
    acc_ovr    = acc(bl_ovr),
    acc_lp     = acc(bl_lp),
    lift_joint_K = acc(bl_joint_K) - acc(bl_lp),
    lift_ovr    = acc(bl_ovr)     - acc(bl_lp)
  )
}

all_results <- do.call(rbind, lapply(Ks, function(K) {
  do.call(rbind, lapply(seq_len(reps), function(r) run_one(K, r)))
}))

agg <- aggregate(cbind(acc_joint_K, acc_ovr, acc_lp,
                         lift_joint_K, lift_ovr) ~ K,
                  data = all_results, FUN = mean)

winner_per_K <- ifelse(agg$lift_joint_K >= agg$lift_ovr, "joint_K", "ovr")
overall_winner <- if (mean(agg$lift_joint_K) > mean(agg$lift_ovr) + 0.01) {
  "joint_K"
} else if (mean(agg$lift_ovr) > mean(agg$lift_joint_K) + 0.01) {
  "ovr"
} else {
  "ovr (tie-break: fewer assumptions)"
}

md <- c("# Phase 4 A/B: categorical joint_K vs OVR vs LP baseline",
        "",
        sprintf("Run on: %s", format(Sys.time())),
        sprintf("Species per scenario: %d, reps: %d, rho: %.2f",
                n_species, reps, rho),
        "Data: continuous + binary + K-class categorical (all BM with rho-correlation).",
        "",
        "Accuracy on held-out categorical test rows (argmax):",
        "",
        "```",
        capture.output(print(agg, row.names = FALSE)),
        "```",
        "",
        sprintf("Per-K winners: %s",
                paste(Ks, "->", winner_per_K, collapse = ", ")),
        sprintf("Overall winner: **%s**", overall_winner),
        "",
        sprintf("joint_K lifts >= 2pp: %d / %d",
                sum(agg$lift_joint_K >= 0.02), nrow(agg)),
        sprintf("OVR lifts     >= 2pp: %d / %d",
                sum(agg$lift_ovr     >= 0.02), nrow(agg)))

writeLines(md, out_md)
cat("\n--- Summary ---\n")
cat(paste(md, collapse = "\n"), "\n")
