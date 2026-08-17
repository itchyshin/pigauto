# Verify per-type lambda dispatch on the AVONET lambda-modes mask (seed 2026):
# bayes should now keep the continuous lift AND the joint/OVR categorical accuracy.
# Reference (script/bench_avonet_lambda_modes.md, arc/bace-comparators):
#   OLD bayes: Trophic 0.600 Primary 0.578 | Mass 2301.7 Beak 21.27 Tarsus 24.42
#   fixed_1  : Trophic 0.789 Primary 0.667 | Mass 2625.3
devtools::load_all(".", quiet = TRUE)
suppressPackageStartupMessages(library(ape))
SEED <- 2026L; MISS <- 0.30
e <- new.env(); utils::data("avonet300", package = "pigauto", envir = e)
utils::data("tree300", package = "pigauto", envir = e)
df <- e$avonet300; tree <- e$tree300
rownames(df) <- df$Species_Key; df$Species_Key <- NULL
set.seed(SEED)
mask <- matrix(FALSE, nrow(df), ncol(df), dimnames = list(rownames(df), names(df)))
for (v in names(df)) { obs <- which(!is.na(df[[v]]))
  mask[sample(obs, ceiling(length(obs) * MISS)), v] <- TRUE }
dm <- df; for (v in names(df)) dm[[v]][mask[, v]] <- NA
t0 <- proc.time()[[3]]
res <- pigauto::impute(dm, tree, lambda_mode = "bayes", seed = SEED,
                       verbose = FALSE, n_imputations = 20L)
cat(sprintf("wall %.0fs\n", proc.time()[[3]] - t0))
for (v in names(df)) {
  idx <- which(mask[, v]); tv <- df[[v]][idx]; cv <- res$completed[[v]][idx]
  if (is.factor(tv) || is.character(tv))
    cat(sprintf("%-22s acc  %.3f\n", v, mean(as.character(cv) == as.character(tv), na.rm=TRUE)))
  else cat(sprintf("%-22s rmse %.2f\n", v, sqrt(mean((as.numeric(cv)-as.numeric(tv))^2, na.rm=TRUE))))
}
