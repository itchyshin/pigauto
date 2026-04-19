# =============================================================================
# Phase 3 probe: compare conformal_method = "split" vs "bootstrap"
# =============================================================================
# Same seed-sensitivity probe as probe_seed_sensitivity.R, but runs each
# seed twice: once with the default "split" conformal, once with the new
# "bootstrap" smoothed variant (B = 500 resamples of val residuals).
#
# Expected: bootstrap SHOULD reduce the conformal-score variance across
# seeds substantially (from sd = 2.94 seen in Phase 1c). If it also
# reduces coverage variance without lowering the mean, we ship it.
#
# Run: cd pigauto && Rscript dev/coverage_sim/probe_bootstrap_conformal.R
# =============================================================================

library(BACE)
suppressMessages(devtools::load_all(".", quiet = TRUE))
library(ape); library(MASS)

set.seed(42)
source_funcs <- function(path, fn_names) {
  src <- readLines(path)
  for (fn in fn_names) {
    s <- grep(paste0("^", fn, " <- function"), src)[1]
    d <- 0; e <- NA
    for (i in s:length(src)) {
      d <- d + sum(gregexpr("\\{", src[i])[[1]] > 0) -
               sum(gregexpr("\\}", src[i])[[1]] > 0)
      if (i > s && d == 0) { e <- i; break }
    }
    eval(parse(text = src[s:e]), envir = globalenv())
  }
}
SIGNAL <- list(high=0.90, moderate=0.60, low=0.20); DEP_STRENGTH <- 1.5
source_funcs("dev/coverage_sim/sim_bace_source.R",
             c(".trait_to_numeric", ".calibrate_intercept", "make_phylo_signal"))

inject_singleobs <- function(cd, tree, mech, rate, vars=c("y","x1","x2","x3","x4")) {
  n <- nrow(cd)
  mm <- as.data.frame(matrix(FALSE, n, length(vars), dimnames=list(NULL, vars)))
  md <- cd
  Sp <- if (mech == "phylo_MAR") ape::vcv(tree, corr=TRUE) else NULL
  tmc <- c(y="x3", x1="y", x2="y", x3="x1", x4="y")
  for (v in vars) {
    lp <- switch(mech,
      phylo_MAR = {
        z <- MASS::mvrnorm(1, mu=rep(0, nrow(Sp)), Sigma=Sp)
        names(z) <- rownames(Sp); -DEP_STRENGTH * z[rownames(cd)]
      },
      trait_MAR = -DEP_STRENGTH * as.numeric(scale(.trait_to_numeric(cd[[tmc[[v]]]]))),
      trait_MNAR = -DEP_STRENGTH * as.numeric(scale(.trait_to_numeric(cd[[v]])))
    )
    c_hat <- .calibrate_intercept(lp, rate)
    mm[[v]] <- rbinom(n, 1, plogis(c_hat + lp)) == 1
    md[[v]][mm[[v]]] <- NA
  }
  list(miss_data=md, miss_mask=mm)
}

N <- 150L
sim <- sim_bace(response_type="gaussian",
                predictor_types=c("binary","multinomial3","poisson","threshold3"),
                var_names=c("y","x1","x2","x3","x4"),
                phylo_signal=make_phylo_signal("all_high"),
                n_cases=N, n_species=N, beta_sparsity=0.3,
                missingness=c(0,0,0,0,0))
cd <- sim$data
tree <- sim$tree
rownames(cd) <- cd$species; cd$species <- NULL
cd$x1 <- factor(as.character(cd$x1), levels=c("0","1"))

miss <- inject_singleobs(cd, tree, "trait_MAR", 0.35)
miss_data <- miss$miss_data; miss_mask <- miss$miss_mask

# ----- probe loop ------------------------------------------------------------
SEEDS  <- 1:10
EPOCHS <- 500L

run_seed <- function(seed_i, method) {
  res <- tryCatch(
    impute(miss_data, tree, species_col=NULL, log_transform=FALSE,
           missing_frac=0.25, n_imputations=1L, epochs=EPOCHS, verbose=FALSE,
           seed=seed_i, conformal_method=method),
    error = function(e) e)
  if (inherits(res, "error")) return(NULL)

  fit <- res$fit; tm <- fit$trait_map
  y_lc <- tm[["y"]]$latent_cols[1]
  gate_y <- fit$calibrated_gates[y_lc]
  cs_y   <- fit$conformal_scores[y_lc]

  pred_y <- res$completed$y
  sd_y   <- sd(cd$y)
  nrmse_y <- sqrt(mean((pred_y[miss_mask$y] - cd$y[miss_mask$y])^2)) / sd_y

  n_draws <- 500
  pred_on_miss <- pred_y[miss_mask$y]
  halfwidth_orig <- (cs_y / 1.96) * tm[["y"]]$sd
  dr <- matrix(rnorm(length(pred_on_miss) * n_draws, 0, halfwidth_orig),
               ncol = n_draws) + pred_on_miss
  lo <- apply(dr, 1, quantile, 0.025)
  hi <- apply(dr, 1, quantile, 0.975)
  tv <- cd$y[miss_mask$y]
  cov_y <- mean(tv >= lo & tv <= hi)

  data.frame(seed = seed_i, method = method, gate_y = gate_y,
             cs_y = cs_y, nrmse_y = nrmse_y, cov_y = cov_y,
             stringsAsFactors = FALSE)
}

cat("=== split conformal (backward-compat baseline) ===\n")
split_res <- list()
for (s in SEEDS) {
  t0 <- Sys.time()
  r <- run_seed(s, "split")
  if (!is.null(r)) {
    split_res[[length(split_res) + 1]] <- r
    cat(sprintf("seed=%2d split  gate=%.2f cs=%5.2f NRMSE=%5.2f cov=%.2f (%ds)\n",
                s, r$gate_y, r$cs_y, r$nrmse_y, r$cov_y,
                round(as.numeric(difftime(Sys.time(), t0, units="secs")))))
  }
}

cat("\n=== bootstrap conformal (B = 500) ===\n")
boot_res <- list()
for (s in SEEDS) {
  t0 <- Sys.time()
  r <- run_seed(s, "bootstrap")
  if (!is.null(r)) {
    boot_res[[length(boot_res) + 1]] <- r
    cat(sprintf("seed=%2d boot   gate=%.2f cs=%5.2f NRMSE=%5.2f cov=%.2f (%ds)\n",
                s, r$gate_y, r$cs_y, r$nrmse_y, r$cov_y,
                round(as.numeric(difftime(Sys.time(), t0, units="secs")))))
  }
}

split_df <- do.call(rbind, split_res)
boot_df  <- do.call(rbind, boot_res)

cat("\n=== COMPARISON ===\n")
summarise_df <- function(d, label) {
  cat(sprintf("--- %s ---\n", label))
  for (col in c("gate_y", "cs_y", "nrmse_y", "cov_y")) {
    v <- d[[col]]
    cat(sprintf("  %-8s mean=%.3f median=%.3f sd=%.3f range=[%.2f, %.2f]\n",
                col, mean(v), median(v), sd(v), min(v), max(v)))
  }
}
summarise_df(split_df, "split conformal")
summarise_df(boot_df,  "bootstrap conformal (B=500)")

# Side-by-side table
cat("\nSide-by-side per seed:\n")
comp <- data.frame(
  seed = split_df$seed,
  gate_split = split_df$gate_y, gate_boot = boot_df$gate_y,
  cs_split = split_df$cs_y,     cs_boot   = boot_df$cs_y,
  cov_split = split_df$cov_y,   cov_boot  = boot_df$cov_y
)
print(comp, row.names=FALSE, digits=3)

saveRDS(list(split = split_df, bootstrap = boot_df, comparison = comp),
        "dev/coverage_sim/results/bootstrap_conformal_probe.rds")
cat("\nSaved to dev/coverage_sim/results/bootstrap_conformal_probe.rds\n")
