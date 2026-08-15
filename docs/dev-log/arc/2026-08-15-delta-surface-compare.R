# Calibration vs predict delta-surface comparison --------------------------------
#
# Question: does gate calibration score a DIFFERENT delta surface than the one
# predict() produces for the same val cells? If so, a gate chosen as helpful at
# calibration can be harmful at predict — which is what the strict val-floor
# failure looks like.
#
# The structural asymmetry, read from the code:
#
#   calibration (fit_pigauto.R:819)  X_init = t_X_eval  (val cells START AT BM)
#                                    pin_mask = t_M_obs (val cells NOT pinned:
#                                    M_obs_mat[val_idx, test_idx] <- FALSE)
#     -> the GNN never sees val truth as context.
#
#   predict  (predict_pigauto.R:339) observed_mask = !is.na(X_scaled)
#     -> val cells ARE non-NA (held out only for training), so they are PINNED
#        TO THEIR OWN TRUTH as DAE context throughout the refine loop.
#
# `.mask_observed_idx` un-pins chosen cells at predict, reproducing the
# calibration-time context. Comparing the two predicts isolates the gap.
#
# Measured per val cell of x_cont, on the same fit (no refit, so the model and
# the calibrated gate are held fixed and ONLY the input context differs):
#   A = predict(fit)                              -> "predict surface"  (val truth pinned)
#   B = predict(fit, .mask_observed_idx = val_idx) -> "calibration-like" (val truth hidden)

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
lc <- fit$trait_map[[1]]$latent_cols[1]          # x_cont

val_mask <- matrix(FALSE, nrow(X), ncol(X)); val_mask[splits$val_idx] <- TRUE
rows <- which(val_mask[, lc])
truth <- X[rows, lc]
bm <- fit$baseline$mu[rows, lc]

# A: the surface predict() actually produces (val truth pinned as context)
A <- predict(fit, return_se = FALSE)$imputed_latent[rows, lc]
# B: the surface calibration effectively scored (val truth hidden)
B <- predict(fit, return_se = FALSE,
             .mask_observed_idx = splits$val_idx)$imputed_latent[rows, lc]

mse <- function(x) mean((x - truth)^2)
cat(sprintf("\ngate: r_cal_bm=%.2f  r_cal_gnn=%.2f\n",
            fit$r_cal_bm[lc], fit$r_cal_gnn[lc]))
cat(sprintf("\nval MSE on x_cont (%d cells)\n", length(rows)))
cat(sprintf("  pure BM baseline                      : %.4f\n", mse(bm)))
cat(sprintf("  A  predict surface  (val truth PINNED) : %.4f   %s\n", mse(A),
            if (mse(A) <= mse(bm) * 1.05) "within floor" else "**BREACHES FLOOR**"))
cat(sprintf("  B  calib-like       (val truth HIDDEN) : %.4f   %s\n", mse(B),
            if (mse(B) <= mse(bm) * 1.05) "within floor" else "**BREACHES FLOOR**"))
cat(sprintf("\n  max |A - B| across val cells : %.4f\n", max(abs(A - B))))
cat(sprintf("  corr(A, B)                   : %.4f\n", stats::cor(A, B)))
cat(sprintf("  mean(A - B)                  : %+.4f\n", mean(A - B)))
cat("\nIf B is inside the floor and A is not, calibration scored a surface it\n")
cat("does not get at predict time, and the gate it chose is not safe there.\n")
