# AVONET300 regression smoke for the lambda default (spec 5.2): mask 20% of the continuous cells MCAR,
# impute with gnn = FALSE under lambda_mode = "fixed_1" and "estimate", compare z-scored RMSE on the
# masked cells. AVONET has strong phylogenetic signal, so lambda should sit near 1 and RMSE should
# change by less than 5%. Prints AVONET_OK or AVONET_FAIL plus the estimated lambdas.
suppressMessages(devtools::load_all(quiet = TRUE))
data(avonet300); data(tree300)
df <- avonet300; rownames(df) <- df$Species_Key; df$Species_Key <- NULL
cont <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")
set.seed(20260922)
truth <- df[cont]
mask <- matrix(runif(nrow(df) * length(cont)) < 0.20, nrow(df))
mask[is.na(as.matrix(truth))] <- FALSE
df_masked <- df; for (j in seq_along(cont)) df_masked[mask[, j], cont[j]] <- NA
zrmse <- function(res) {
  out <- numeric(length(cont))
  for (j in seq_along(cont)) {
    y <- log(truth[[cont[j]]]); yhat <- log(res$completed[[cont[j]]])
    sel <- mask[, j]; out[j] <- sqrt(mean((yhat[sel] - y[sel])^2)) / sd(y, na.rm = TRUE)
  }
  setNames(out, cont)
}
r1 <- impute(df_masked, tree300, gnn = FALSE, lambda_mode = "fixed_1")
r2 <- impute(df_masked, tree300, gnn = FALSE, lambda_mode = "estimate")
z1 <- zrmse(r1); z2 <- zrmse(r2)
cat("fixed_1 :", round(z1, 4), "\n"); cat("estimate:", round(z2, 4), "\n")
lam <- r2$fit$model_config$lambda_per_trait
cat("lambda_per_trait:", if (is.null(lam)) "NULL" else paste(names(lam), round(lam, 3), collapse = " "), "\n")
rel <- (mean(z2) - mean(z1)) / mean(z1)
cat(sprintf("relative change in mean z-RMSE: %+.2f%%\n", 100 * rel))
cat(if (is.finite(rel) && abs(rel) < 0.05) "AVONET_OK\n" else "AVONET_FAIL\n")
