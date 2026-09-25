# Pre-run (lambda-bias lane, 2026-09-24): does predict_method = "exact" (cross-trait covariance)
# close the gap to the frequentist lambda arm? One job = one (lambda, n, seed). Reuses the
# imputation-sim lane's make_cell/run_arms/score_arm unchanged (sourced read-only).
args <- commandArgs(trailingOnly = TRUE)
lam <- as.numeric(args[1]); n <- as.integer(args[2]); seed <- as.integer(args[3]); out <- args[4]
f <- file.path(out, sprintf("l%s_n%d_s%d.rds", lam, n, seed)); if (file.exists(f)) quit(save = "no")
suppressMessages(library(pigauto))
source(file.path(Sys.getenv("SIM_ROOT"), "script", "campaign_gnn_off_lib.R"))
cell <- make_cell("types_mixed", n, seed, miss_frac = 0.3, miss = "mcar", lambda = lam, rho = 0,
                  evo = "BM", thresholds = "fixed", driver = TRUE)
truth <- cell$truth; mask <- cell$mask
cell_obj <- list(truth = truth, tree = cell$tree, mask = mask, df_miss = cell$df_miss,
                 cont_traits = cell$cont_traits, trait_types = cell$trait_types, covs = cell$covs)
t0 <- proc.time()[[3]]
ra <- run_arms(cell_obj, c("gnn_off", "freq_lambda"), list(seed = seed))
t1 <- proc.time()[[3]]
res <- tryCatch(pigauto::impute(cell$df_miss, cell$tree, verbose = FALSE, seed = seed,
                                trait_types = cell$trait_types, gnn = FALSE, predict_method = "exact"),
                error = function(e) e)
t2 <- proc.time()[[3]]
assign("truth", truth, envir = .GlobalEnv); assign("mask", mask, envir = .GlobalEnv)
ex <- if (inherits(res, "error")) NULL else {
  comp <- res$completed[rownames(truth), names(truth)]
  score_arm("gnn_off_exact", comp, res$prediction$conformal_lower, res$prediction$conformal_upper,
            prob = res$prediction$probabilities)
}
tab <- rbind(ra$results, ex)
tab$lambda <- lam; tab$n <- n; tab$seed <- seed
saveRDS(list(results = tab, error = if (inherits(res, "error")) conditionMessage(res) else NULL,
             secs = c(arms = t1 - t0, exact = t2 - t1),
             exact_path = if (inherits(res, "error")) NULL else res$fit$model_config$predict_method), f)
