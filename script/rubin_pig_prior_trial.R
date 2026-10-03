# script/rubin_pig_prior_trial.R
#
# Trial of Sigma_E prior variants for pig_post's lambda = 1 slope bias (cause in
# docs/dev-log/arc/2026-10-02-rubin-pigauto-campaign.md). Patches pigauto:::.mip_fit IN MEMORY for this R session only
# (assignInNamespace); no package file is edited. Refits one campaign dataset with the patched prior and keeps the
# completed datasets for scoring against the oracle (script/rubin_pig_prior_trial_summary.R).
#
#   variant base: S_E = diag(0.01 * obs_var)                (pigauto b565cad, unchanged)
#   variant A   : S_E = diag(1e-4 * obs_var)                (smaller scale)
#   variant B   : S_E = 0.01 * pairwise covariance of the observed latents (residual correlation follows the data)
#   variant C   : S_E = diag((0.1 / n) * obs_var)            (shrinks with n: 1e-4 at n = 1000, 1e-3 at n = 100)
#
#   Rscript script/rubin_pig_prior_trial.R --variant A --n 1000 --lambda 1 --rho 0.5 --seed 1 --out <dir>

RNGkind("L'Ecuyer-CMRG")
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[i + 1L] }
variant <- get_arg("--variant"); n <- as.integer(get_arg("--n")); lambda <- as.numeric(get_arg("--lambda"))
rho <- as.numeric(get_arg("--rho")); seed <- as.integer(get_arg("--seed")); out <- get_arg("--out")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
f_out <- file.path(out, sprintf("prior_%s_n%d_l%s_r%s_s%d.rds", variant, n, format(lambda), format(rho), seed))
if (file.exists(f_out)) quit(save = "no", status = 0)

# variant SEP: pigauto's opt-in residual_prior = "sep" (research branch research/mi-posterior-sep-prior); no patch
old <- "S_E = diag(0.01 * prob$obs_var, K)"
if (identical(variant, "SEP")) variant_control <- list(residual_prior = "sep") else variant_control <- list()
new <- switch(variant,
  base = old, SEP = old,
  A = "S_E = diag(1e-04 * prob$obs_var, K)",
  C = "S_E = diag((0.1 / prob$n) * prob$obs_var, K)",
  B = paste0("S_E = { C0 <- stats::cov(prob$Y, use = \"pairwise.complete.obs\"); C0 <- (C0 + t(C0)) / 2; ",
             "ev <- eigen(C0, symmetric = TRUE); 0.01 * (ev$vectors %*% diag(pmax(ev$values, 1e-3 * max(ev$values)), K) %*% t(ev$vectors)) }"),
  stop("unknown variant"))
fit_src <- paste(deparse(pigauto:::.mip_fit, width.cutoff = 500L), collapse = "\n")
if (!grepl(old, fit_src, fixed = TRUE)) stop("installed .mip_fit does not contain the expected prior line")
f_new <- eval(parse(text = sub(old, new, fit_src, fixed = TRUE)))
environment(f_new) <- asNamespace("pigauto")
utils::assignInNamespace(".mip_fit", f_new, ns = "pigauto")

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
suppressMessages({ for (f in c("campaign_gnn_off_lib.R", "rubin_lib.R", "rubin_freq.R", "rubin_pigauto.R")) source(file.path(here, f)) })
cell <- make_cell("types_mixed", n, seed, miss_frac = 0.30, miss = "mcar", lambda = lambda, rho = rho,
                  thresholds = "fixed", driver = TRUE)
bt <- default_block_traits(cell)
set.seed(seed + 909L); t0 <- proc.time()[["elapsed"]]
pp <- mi_pig_post(cell, 20, seed = seed + 909L, control = variant_control)
saveRDS(list(variant = variant, n = n, lambda = lambda, rho = rho, seed = seed, block_traits = bt,
             truth = cell$truth[, bt], mask = cell$mask[, bt], tree = cell$tree,
             pig_post = lapply(pp$datasets, function(d) d[, bt, drop = FALSE]),
             pig_diag = pp$diag[c("converged", "rhat_max", "ess_min", "n_extensions")],
             wall_s = proc.time()[["elapsed"]] - t0), f_out)
cat("done", basename(f_out), "\n")
