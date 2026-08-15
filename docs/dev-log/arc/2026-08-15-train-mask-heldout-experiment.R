# train_mask_heldout experiment ------------------------------------------------
#
# Question: does P0 fix #4's held-out masking explain the strict val-floor
# failure (test-safety-floor.R:712)?
#
# Design: replicate that test's fixture exactly, then fit under two conditions
# that differ ONLY in whether held-out truth is visible as DAE context:
#
#   FIXED  — mask_heldout_with_baseline() as merged (P0). Held-out val/test
#            cells replaced by the baseline prediction before training.
#   LEAKED — mask_heldout_with_baseline() forced to identity, restoring the
#            pre-fix behaviour where held-out truth IS visible as context.
#
# Both conditions share the fixture, the seed, and the device, so the only
# difference is the leak. R reps each, because the GNN blend is nondeterministic
# (observed: 0.3434 / 0.3442 across runs) while the BM baseline is not.
#
# Prediction if the hypothesis is right:
#   LEAKED passes the 5% floor; FIXED does not; loss_bm identical in both.

suppressMessages(devtools::load_all(".", quiet = TRUE))

REPS <- 4L

make_fixture <- function() {
  set.seed(20260429L)
  n <- 80L
  tree <- ape::rtree(n)
  sp <- tree$tip.label
  bm_draw <- function(seed) {
    withr::with_seed(seed, {
      v <- as.numeric(ape::rTraitCont(tree, model = "BM", sigma = 1))
      names(v) <- tree$tip.label
      v[sp]
    })
  }
  v1 <- bm_draw(11L); v2 <- bm_draw(12L); v3 <- bm_draw(13L); v4 <- bm_draw(14L)
  qs <- stats::quantile(v3, c(0, 1/3, 2/3, 1), na.rm = TRUE); qs[1] <- qs[1] - 1e-9
  cat3 <- factor(c("a", "b", "c")[as.integer(cut(v3, qs, include.lowest = TRUE))],
                 levels = c("a", "b", "c"))
  df <- data.frame(row.names = sp,
                   x_cont = v1,
                   x_bin  = factor(ifelse(v2 > 0, "yes", "no"), levels = c("no", "yes")),
                   x_cat  = cat3,
                   x_cnt  = as.integer(round(pmax(v4 + 5, 0))))
  list(df = df, tree = tree)
}

# Per-trait val loss for the blend and for pure BM — the same surface the test
# uses (0-1 loss for binary/categorical, MSE on latent otherwise).
losses_one_fit <- function(fx) {
  res <- pigauto::impute(fx$df, fx$tree,
                         epochs = 30L, eval_every = 10L, patience = 5L,
                         missing_frac = 0.30, verbose = FALSE,
                         seed = 20260429L)
  fit <- res$fit; data <- res$data; splits <- res$splits
  X_scaled <- data$X_scaled
  val_mask <- matrix(FALSE, nrow(X_scaled), ncol(X_scaled))
  val_mask[splits$val_idx] <- TRUE

  blend_pred <- predict(fit, return_se = FALSE)$imputed_latent
  bm_pred <- fit$baseline$mu
  if (isTRUE(fit$multi_obs)) bm_pred <- bm_pred[fit$obs_to_species, , drop = FALSE]

  out <- list()
  for (tm in fit$trait_map) {
    lc <- tm$latent_cols
    val_rows <- which(val_mask[, lc[1]])
    if (length(val_rows) == 0L) next
    truth <- X_scaled[val_rows, lc, drop = FALSE]
    bm_v <- bm_pred[val_rows, lc, drop = FALSE]
    bl_v <- blend_pred[val_rows, lc, drop = FALSE]
    disc <- tm$type %in% c("binary", "categorical")
    lb <- if (tm$type == "binary") {
      mean(as.numeric(bl_v[, 1] > 0) != truth[, 1], na.rm = TRUE)
    } else if (tm$type == "categorical") {
      mean(max.col(bl_v, ties.method = "first") !=
             max.col(truth, ties.method = "first"), na.rm = TRUE)
    } else mean((bl_v - truth)^2, na.rm = TRUE)
    lm_ <- if (tm$type == "binary") {
      mean(as.numeric(bm_v[, 1] > 0) != truth[, 1], na.rm = TRUE)
    } else if (tm$type == "categorical") {
      mean(max.col(bm_v, ties.method = "first") !=
             max.col(truth, ties.method = "first"), na.rm = TRUE)
    } else mean((bm_v - truth)^2, na.rm = TRUE)
    out[[tm$name]] <- c(blend = lb, bm = lm_)
  }
  out
}

fx <- make_fixture()

run_condition <- function(label, leaked) {
  cat(sprintf("\n===== %s =====\n", label))
  for (r in seq_len(REPS)) {
    res <- if (leaked) {
      testthat::with_mocked_bindings(
        losses_one_fit(fx),
        mask_heldout_with_baseline = function(X_fill, MU, hold_idx) X_fill,
        .package = "pigauto"
      )
    } else {
      losses_one_fit(fx)
    }
    for (nm in names(res)) {
      lb <- res[[nm]]["blend"]; lm_ <- res[[nm]]["bm"]
      thr <- lm_ * 1.05 + 1e-10
      cat(sprintf("  rep%d  %-7s blend=%.4f  bm=%.4f  thr=%.4f  %s\n",
                  r, nm, lb, lm_, thr, if (lb <= thr) "PASS" else "**FAIL**"))
    }
  }
}

run_condition("FIXED  (P0 as merged — held-out truth masked)", leaked = FALSE)
run_condition("LEAKED (pre-fix behaviour restored)",           leaked = TRUE)
cat("\ndone\n")
