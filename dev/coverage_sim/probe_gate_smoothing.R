# =============================================================================
# Phase 3b probe: does median-splits gate calibration break the bimodal gate?
# =============================================================================
# Compares four combinations on the SAME sim dataset:
#   (conformal_method x gate_method)
#     split     x single_split   (legacy, backward compat)
#     split     x median_splits  (gate fix only)
#     bootstrap x single_split   (conformal fix only)
#     bootstrap x median_splits  (both fixes)
#
# 10 seeds per combination. Reports gate SD, conformal score SD, NRMSE SD,
# coverage SD. Goal: gate SD drops substantially when median_splits is on;
# coverage SD ideally drops too.
# =============================================================================

library(BACE)
suppressMessages(devtools::load_all(".", quiet = TRUE))
library(ape); library(MASS)

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
      phylo_MAR = { z <- MASS::mvrnorm(1, mu=rep(0, nrow(Sp)), Sigma=Sp)
                    names(z) <- rownames(Sp); -DEP_STRENGTH * z[rownames(cd)] },
      trait_MAR  = -DEP_STRENGTH * as.numeric(scale(.trait_to_numeric(cd[[tmc[[v]]]]))),
      trait_MNAR = -DEP_STRENGTH * as.numeric(scale(.trait_to_numeric(cd[[v]])))
    )
    c_hat <- .calibrate_intercept(lp, rate)
    mm[[v]] <- rbinom(n, 1, plogis(c_hat + lp)) == 1
    md[[v]][mm[[v]]] <- NA
  }
  list(miss_data=md, miss_mask=mm)
}

set.seed(42)
N <- 150L
sim <- sim_bace(response_type="gaussian",
                predictor_types=c("binary","multinomial3","poisson","threshold3"),
                var_names=c("y","x1","x2","x3","x4"),
                phylo_signal=make_phylo_signal("all_high"),
                n_cases=N, n_species=N, beta_sparsity=0.3,
                missingness=c(0,0,0,0,0))
cd <- sim$data; tree <- sim$tree
rownames(cd) <- cd$species; cd$species <- NULL
cd$x1 <- factor(as.character(cd$x1), levels=c("0","1"))

miss <- inject_singleobs(cd, tree, "trait_MAR", 0.35)
miss_data <- miss$miss_data; miss_mask <- miss$miss_mask

SEEDS <- 1:10
EPOCHS <- 500L

run_one <- function(seed_i, cm, gm) {
  res <- tryCatch(
    impute(miss_data, tree, species_col=NULL, log_transform=FALSE,
           missing_frac=0.25, n_imputations=1L, epochs=EPOCHS, verbose=FALSE,
           seed=seed_i,
           conformal_method=cm, gate_method=gm),
    error = function(e) e)
  if (inherits(res, "error")) return(NULL)

  fit <- res$fit; tm <- fit$trait_map
  y_lc <- tm[["y"]]$latent_cols[1]
  gate_y <- fit$calibrated_gates[y_lc]
  cs_y   <- fit$conformal_scores[y_lc]

  pred_y <- res$completed$y
  sd_y <- sd(cd$y)
  nrmse_y <- sqrt(mean((pred_y[miss_mask$y] - cd$y[miss_mask$y])^2)) / sd_y

  n_draws <- 500
  po <- pred_y[miss_mask$y]
  hw <- (cs_y / 1.96) * tm[["y"]]$sd
  dr <- matrix(rnorm(length(po) * n_draws, 0, hw), ncol=n_draws) + po
  lo <- apply(dr, 1, quantile, 0.025); hi <- apply(dr, 1, quantile, 0.975)
  tv <- cd$y[miss_mask$y]
  cov_y <- mean(tv >= lo & tv <= hi)

  data.frame(seed=seed_i, cm=cm, gm=gm, gate_y=gate_y, cs_y=cs_y,
             nrmse_y=nrmse_y, cov_y=cov_y, stringsAsFactors=FALSE)
}

COMBOS <- list(
  list(cm="split",     gm="single_split"),
  list(cm="split",     gm="median_splits"),
  list(cm="bootstrap", gm="single_split"),
  list(cm="bootstrap", gm="median_splits")
)

rows <- list()
for (combo in COMBOS) {
  cat(sprintf("\n=== cm=%s, gm=%s ===\n", combo$cm, combo$gm))
  for (s in SEEDS) {
    t0 <- Sys.time()
    r <- run_one(s, combo$cm, combo$gm)
    if (!is.null(r)) {
      rows[[length(rows) + 1]] <- r
      cat(sprintf("seed=%2d gate=%.2f cs=%5.2f NRMSE=%5.2f cov=%.2f (%ds)\n",
                  s, r$gate_y, r$cs_y, r$nrmse_y, r$cov_y,
                  round(as.numeric(difftime(Sys.time(), t0, units="secs")))))
    } else {
      cat(sprintf("seed=%2d FAILED\n", s))
    }
  }
}

df <- do.call(rbind, rows)
cat("\n\n=== AGGREGATES ===\n")
for (combo in COMBOS) {
  sub <- df[df$cm == combo$cm & df$gm == combo$gm, ]
  cat(sprintf("\n-- cm=%s gm=%s --\n", combo$cm, combo$gm))
  for (col in c("gate_y","cs_y","nrmse_y","cov_y")) {
    v <- sub[[col]]
    cat(sprintf("  %-8s mean=%.3f median=%.3f sd=%.3f range=[%.2f, %.2f]\n",
                col, mean(v), median(v), sd(v), min(v), max(v)))
  }
}

saveRDS(df, "dev/coverage_sim/results/gate_smoothing_probe.rds")
cat("\nSaved to dev/coverage_sim/results/gate_smoothing_probe.rds\n")
