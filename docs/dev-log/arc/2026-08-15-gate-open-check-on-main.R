# Gate-open check on pristine main ---------------------------------------------
#
# Question: does the calibration<->predict surface asymmetry PREDATE P0?
#
# On main the calibrated gate closes (r_cal_gnn = 0), so blend == pure BM and
# the input context provably cannot affect the output. That makes the asymmetry
# unobservable — not absent. This forces the gate open to the value P0's own
# calibration chose (0.85 / 0.15) on the SAME fitted model, then measures the
# two surfaces:
#
#   A = predict(fit)                               val truth PINNED   (what users get)
#   B = predict(fit, .mask_observed_idx = val_idx) val truth HIDDEN   (what calibration scores)
#
# No refit between A and B: the model and the (forced) gate are identical, so
# any difference is caused by input context alone.
#
# If A != B here, the asymmetry exists on main too and P0 is only the messenger.
# If A == B here, P0 introduced it.

suppressMessages(devtools::load_all(".", quiet = TRUE))

set.seed(20260429L)
n <- 80L; tree <- ape::rtree(n); sp <- tree$tip.label
bm_draw <- function(seed) withr::with_seed(seed, {
  v <- as.numeric(ape::rTraitCont(tree, model = "BM", sigma = 1))
  names(v) <- tree$tip.label; v[sp] })
v1 <- bm_draw(11L); v2 <- bm_draw(12L); v3 <- bm_draw(13L); v4 <- bm_draw(14L)
qs <- stats::quantile(v3, c(0, 1/3, 2/3, 1), na.rm = TRUE); qs[1] <- qs[1] - 1e-9
cat3 <- factor(c("a", "b", "c")[as.integer(cut(v3, qs, include.lowest = TRUE))],
               levels = c("a", "b", "c"))
df <- data.frame(row.names = sp, x_cont = v1,
                 x_bin = factor(ifelse(v2 > 0, "yes", "no"), levels = c("no", "yes")),
                 x_cat = cat3, x_cnt = as.integer(round(pmax(v4 + 5, 0))))

res <- pigauto::impute(df, tree, epochs = 30L, eval_every = 10L, patience = 5L,
                       missing_frac = 0.30, verbose = FALSE, seed = 20260429L)
fit <- res$fit; splits <- res$splits; X <- res$data$X_scaled
lc <- fit$trait_map[[1]]$latent_cols[1]
rows <- which(matrix(replace(logical(length(X)), splits$val_idx, TRUE),
                     nrow(X), ncol(X))[, lc])
truth <- X[rows, lc]; bm <- fit$baseline$mu[rows, lc]
mse <- function(x) mean((x - truth)^2)

cat(sprintf("\n--- pristine main %s ---\n", substr(system("git rev-parse HEAD", intern = TRUE), 1, 7)))
cat(sprintf("AS-CALIBRATED gate: r_cal_bm=%.2f r_cal_gnn=%.2f  (expect the gate CLOSED)\n",
            fit$r_cal_bm[lc], fit$r_cal_gnn[lc]))
A0 <- predict(fit, return_se = FALSE)$imputed_latent[rows, lc]
B0 <- predict(fit, return_se = FALSE,
              .mask_observed_idx = splits$val_idx)$imputed_latent[rows, lc]
cat(sprintf("  gate as-is : A=%.4f  B=%.4f  max|A-B|=%.4g   (BM=%.4f)\n",
            mse(A0), mse(B0), max(abs(A0 - B0)), mse(bm)))

# Force the gate open to exactly what P0's calibration selected. Same model.
fit2 <- fit
fit2$r_cal_bm[lc]   <- 0.85
fit2$r_cal_gnn[lc]  <- 0.15
if (!is.null(fit2$r_cal)) fit2$r_cal[lc] <- 0.15

A <- predict(fit2, return_se = FALSE)$imputed_latent[rows, lc]
B <- predict(fit2, return_se = FALSE,
             .mask_observed_idx = splits$val_idx)$imputed_latent[rows, lc]

cat(sprintf("\nFORCED gate 0.85/0.15 on the SAME model (%d val cells)\n", length(rows)))
cat(sprintf("  pure BM                              : %.4f\n", mse(bm)))
cat(sprintf("  A  predict surface (truth PINNED)    : %.4f  %s\n", mse(A),
            if (mse(A) <= mse(bm) * 1.05) "within floor" else "**BREACHES**"))
cat(sprintf("  B  calib-like      (truth HIDDEN)    : %.4f  %s\n", mse(B),
            if (mse(B) <= mse(bm) * 1.05) "within floor" else "**BREACHES**"))
cat(sprintf("  max|A-B| = %.4f   corr = %.4f   mean(A-B) = %+.4f\n",
            max(abs(A - B)), stats::cor(A, B), mean(A - B)))
cat("\nA != B here  => asymmetry PREDATES P0 (P0 only opened the gate that reveals it)\n")
cat("A == B here  => P0 introduced it\n")
