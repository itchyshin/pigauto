#!/usr/bin/env Rscript
# script/compare_se_methods.R
#
# Empirical comparison of three uncertainty sources for continuous traits:
#   1. BM-only SE  -- analytic conditional-MVN (single forward pass)
#   2. MC dropout SD -- between-pass variation across M GNN forward passes
#   3. Combined SE  -- sqrt(BM^2 + MC^2), what pred$se returns with n_imp > 1
#   4. Conformal half-width / 1.96 -- empirical from held-out residuals
#
# Question: how similar / different are these, and which is most conservative?
#
# Run from package root:
#   Rscript script/compare_se_methods.R

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(ggplot2)
})

set.seed(42)
data(avonet300, package = "pigauto")
data(tree300,   package = "pigauto")

# --- prep: create ~20% missingness in three continuous traits -----------------
df <- avonet300
rownames(df) <- df$Species_Key
df$Species_Key <- NULL
traits_keep <- c("Mass", "Beak.Length_Culmen", "Wing.Length",
                 "Trophic.Level", "Primary.Lifestyle")
df <- df[, traits_keep]

set.seed(7)
n  <- nrow(df)
for (col in c("Mass", "Beak.Length_Culmen", "Wing.Length")) {
  df[[col]][sample(n, round(0.20 * n))] <- NA
}

cat("Missing cells per trait:\n")
print(colSums(is.na(df)))
cat("\n")

# --- 1. Single forward pass: BM SE only ---------------------------------------
cat("Fitting model (single pass)...\n")
res1 <- impute(df, tree300, epochs = 500L, verbose = FALSE, seed = 1L)
pred1 <- res1$prediction   # n_imputations = 1L by default

bm_se <- pred1$se          # matrix: BM SE (original scale)

# --- 2. MC dropout: M=30 stochastic forward passes ---------------------------
cat("Running MC dropout (M=30)...\n")
pred_mc <- predict(res1$fit, return_se = TRUE, n_imputations = 30L)
mc_se   <- pred_mc$se      # matrix: sqrt(BM^2 + MC^2) combined SE

# Extract pure MC SD by back-computing from combined (approximate):
# combined^2 = BM^2 + MC^2  =>  MC = sqrt(max(combined^2 - BM^2, 0))
mc_only_se <- sqrt(pmax(mc_se^2 - bm_se^2, 0))

# --- 3. Conformal half-width / 1.96 ------------------------------------------
conf_scores <- pred1$conformal_scores   # named vector, latent scale
conf_se <- matrix(NA_real_, nrow = n, ncol = 3,
                  dimnames = list(rownames(df),
                                  c("Mass", "Beak.Length_Culmen", "Wing.Length")))
for (nm in colnames(conf_se)) {
  tm <- res1$fit$trait_map[[nm]]
  if (!is.null(tm) && !is.null(conf_scores) && nm %in% names(conf_scores) &&
      is.finite(conf_scores[nm])) {
    # conformal score is in latent (z-score) scale; convert to original
    cs_orig <- conf_scores[nm] * tm$sd
    if (isTRUE(tm$log_transform)) {
      # delta method: se_orig = mu * se_log
      mu_hat <- pred1$imputed[[nm]]
      cs_orig <- mu_hat * (conf_scores[nm] / 1.96)   # already approximate
    }
    conf_se[, nm] <- cs_orig / 1.96
  }
}

# --- 4. Summary table: mean SE for missing cells only ------------------------
imask <- res1$imputed_mask

cat("=============================================================\n")
cat("Mean SE at originally-missing cells (continuous traits only)\n")
cat("=============================================================\n\n")

for (nm in c("Mass", "Beak.Length_Culmen", "Wing.Length")) {
  miss_rows <- which(imask[, nm])
  if (length(miss_rows) == 0) next

  bm   <- mean(bm_se[miss_rows, nm],      na.rm = TRUE)
  mc   <- mean(mc_only_se[miss_rows, nm], na.rm = TRUE)
  comb <- mean(mc_se[miss_rows, nm],      na.rm = TRUE)
  conf <- mean(conf_se[miss_rows, nm],    na.rm = TRUE)

  cat(sprintf("Trait: %s  (n_missing = %d)\n", nm, length(miss_rows)))
  cat(sprintf("  BM-only SE              : %8.3f\n", bm))
  cat(sprintf("  MC dropout SD (GNN only): %8.3f\n", mc))
  cat(sprintf("  Combined SE (BM+MC)     : %8.3f\n", comb))
  cat(sprintf("  Conformal/1.96          : %8.3f\n", conf))
  cat(sprintf("  Ratio MC/BM             : %8.3f\n", mc / bm))
  cat(sprintf("  Ratio Conformal/BM      : %8.3f\n", conf / bm))
  cat("\n")
}

# --- 5. Scatter: BM SE vs Combined SE (log scale) ----------------------------
plot_df <- data.frame(
  bm_se    = as.vector(bm_se[, c("Mass","Beak.Length_Culmen","Wing.Length")]),
  comb_se  = as.vector(mc_se[, c("Mass","Beak.Length_Culmen","Wing.Length")]),
  trait    = rep(c("Mass","Beak.Length_Culmen","Wing.Length"), each = n),
  missing  = as.vector(imask[, c("Mass","Beak.Length_Culmen","Wing.Length")])
)
plot_df <- plot_df[plot_df$missing & is.finite(plot_df$bm_se) &
                   is.finite(plot_df$comb_se), ]

p <- ggplot(plot_df, aes(x = bm_se, y = comb_se, colour = trait)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "BM-only SE vs Combined SE (BM + MC dropout)",
       subtitle = "Dashed line = equality; points above = MC adds uncertainty",
       x = "BM-only SE (log scale)", y = "Combined SE (log scale)",
       colour = "Trait") +
  theme_bw()

out_png <- "script/compare_se_methods.png"
ggsave(out_png, p, width = 6, height = 5, dpi = 150)
cat(sprintf("Scatter plot saved to %s\n", out_png))

cat("\nConclusion:\n")
cat("  - BM SE captures phylogenetic imputation uncertainty (model-dependent).\n")
cat("  - MC dropout SD captures GNN architectural uncertainty (model-free).\n")
cat("  - Combined = sqrt(BM^2 + MC^2) is always >= BM alone.\n")
cat("  - Conformal/1.96 is empirically calibrated on held-out residuals:\n")
cat("    it will be wider than BM alone when GNN adds prediction error,\n")
cat("    and is the most defensible 95% CI without Gaussian assumption.\n")
