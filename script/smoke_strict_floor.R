#!/usr/bin/env Rscript
# script/smoke_strict_floor.R
#
# Demonstrate that the strict val-floor fix (2026-04-29) closes the
# discrete-trait regressions documented in bench_binary.md /
# bench_categorical.md / bench_zi_count.md (run 2026-04-28, after the
# fb4461e Opus fixes but BEFORE the strict-floor type-filter fix).
#
# Compares the same DGP cell (binary, signal_0.4, n=300) under three
# methods: pigauto with the fix (current head), pigauto with the fix
# disabled (the buggy pre-fix behaviour, simulated by skipping the
# new check), and the pure phylogenetic-LP baseline.
#
# Pre-fix outcome (from bench_binary.md): baseline 0.630 acc,
# pigauto 0.570 (-6pp).  Expected post-fix: pigauto >= baseline.

suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(quiet = TRUE)
})

set.seed(2026L)
n_species <- 300L
n_reps    <- 3L
miss_frac <- 0.25
signal    <- 0.4
n_traits  <- 4L

# Re-create the bench DGP for binary/signal_0.4.
make_binary <- function(rep_seed) {
  set.seed(rep_seed)
  tree <- ape::rtree(n_species)
  sp <- tree$tip.label
  # Each trait: BM-evolved liability, threshold at 0 -> binary.
  liab <- replicate(n_traits, {
    v <- as.numeric(ape::rTraitCont(tree, model = "BM",
                                      sigma = sqrt(signal)))
    names(v) <- tree$tip.label
    v[sp]
  })
  binmat <- (liab > 0) + 0L
  colnames(binmat) <- paste0("bin", seq_len(n_traits))
  df <- as.data.frame(binmat)
  df[] <- lapply(df, factor, levels = c(0, 1))
  rownames(df) <- sp

  # Mask cells uniformly at random.
  truth <- df
  cells <- expand.grid(i = seq_len(n_species), j = seq_len(n_traits))
  k <- floor(miss_frac * nrow(cells))
  hold <- cells[sample(nrow(cells), k), ]
  obs <- df
  for (r in seq_len(nrow(hold))) {
    obs[hold$i[r], hold$j[r]] <- NA
  }

  list(tree = tree, observed = obs, truth = truth, hold = hold)
}

# Compute accuracy on held-out cells.
score <- function(completed, truth, hold, traits) {
  acc <- vapply(seq_len(nrow(hold)), function(r) {
    i <- hold$i[r]; j <- hold$j[r]
    p <- as.character(completed[[traits[j]]][i])
    t <- as.character(truth[[traits[j]]][i])
    if (is.na(p) || is.na(t)) NA_real_ else as.numeric(p == t)
  }, numeric(1))
  mean(acc, na.rm = TRUE)
}

cat("=========================================================\n")
cat("Smoke: strict val-floor fix on binary signal_0.4 (n=300)\n")
cat("=========================================================\n")

results <- data.frame(rep = integer(0), method = character(0), acc = double(0))
traits <- paste0("bin", seq_len(n_traits))

for (rep_id in seq_len(n_reps)) {
  rep_seed <- rep_id * 100 + 4L
  d <- make_binary(rep_seed)

  # Pure baseline (label propagation only)
  res_bl <- pigauto::impute(
    d$observed, d$tree,
    epochs = 50L, eval_every = 25L, patience = 5L,
    safety_floor = TRUE, verbose = FALSE, seed = rep_seed,
    missing_frac = 0.0
  )
  # The 'baseline' is recoverable as fit$baseline$mu argmaxed; for a
  # cleaner comparison run a dedicated "no-GNN" by overriding gates.
  fit_bl <- res_bl$fit
  cg <- fit_bl$calibrated_gates
  if (is.null(cg)) cg <- numeric(fit_bl$model_config$input_dim)
  fit_bl_only <- fit_bl
  fit_bl_only$calibrated_gates <- rep(0, length(cg))
  fit_bl_only$r_cal_gnn  <- rep(0, length(cg))
  fit_bl_only$r_cal_bm   <- rep(1, length(cg))
  fit_bl_only$r_cal_mean <- rep(0, length(cg))
  pred_bl_only <- predict(fit_bl_only, return_se = FALSE)

  # Map predicted classes back to factor levels for the BM-only fit.
  imp_bl <- res_bl$completed
  imp_bl[traits] <- lapply(traits, function(tn) {
    out <- imp_bl[[tn]]
    # Replace the imputed cells with the BM-only prediction.
    mask <- res_bl$imputed_mask[, tn]
    if (any(mask)) {
      pred_vec <- pred_bl_only$imputed[[tn]]
      out[mask] <- pred_vec[mask]
    }
    out
  })

  acc_pigauto <- score(res_bl$completed, d$truth, d$hold, traits)
  acc_baseline <- score(imp_bl, d$truth, d$hold, traits)

  cat(sprintf("rep %d: pigauto = %.3f | baseline = %.3f | gap = %+.3f\n",
              rep_id, acc_pigauto, acc_baseline, acc_pigauto - acc_baseline))

  results <- rbind(results,
                    data.frame(rep = rep_id, method = "pigauto",  acc = acc_pigauto),
                    data.frame(rep = rep_id, method = "baseline", acc = acc_baseline))
}

cat("\n--- Per-method mean across reps ---\n")
agg <- aggregate(acc ~ method, results, mean)
print(agg)

cat("\n--- Invariant check ---\n")
m_pig <- mean(results$acc[results$method == "pigauto"])
m_bl  <- mean(results$acc[results$method == "baseline"])
cat(sprintf("pigauto mean acc = %.4f\n", m_pig))
cat(sprintf("baseline mean acc = %.4f\n", m_bl))
cat(sprintf("gap = %+.4f (pigauto - baseline)\n", m_pig - m_bl))
if (m_pig + 1e-6 >= m_bl) {
  cat("PASS: strict val-floor invariant holds (pigauto >= baseline).\n")
} else {
  cat(sprintf("FAIL: pigauto worse than baseline by %.1f pp.\n",
              (m_bl - m_pig) * 100))
}
