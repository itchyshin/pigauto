# script/rubin_pig_lambda_diag.R
#
# Diagnostic for the pig_post campaign finding (docs/dev-log/arc/2026-10-02-rubin-pigauto-campaign.md): at n = 1000 and
# lambda = 1, pig_post's downstream coverage was 0.85-0.86. Refits one campaign cell EXACTLY (same dataset, same seed
# offset 909, same call) and saves the posterior draws the campaign files did not keep: lambda, Sigma_P, Sigma_E, mu.
#
#   Rscript script/rubin_pig_lambda_diag.R --n 1000 --lambda 1 --rho 0.5 --seed 1 --out <dir>

RNGkind("L'Ecuyer-CMRG")
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[i + 1L] }
n <- as.integer(get_arg("--n")); lambda <- as.numeric(get_arg("--lambda")); rho <- as.numeric(get_arg("--rho"))
seed <- as.integer(get_arg("--seed")); out <- get_arg("--out")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
f_out <- file.path(out, sprintf("lambda_diag_n%d_l%s_r%s_s%d.rds", n, format(lambda), format(rho), seed))
if (file.exists(f_out)) quit(save = "no", status = 0)

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
suppressMessages({ for (f in c("campaign_gnn_off_lib.R", "rubin_lib.R", "rubin_freq.R", "rubin_pigauto.R")) source(file.path(here, f)) })

cell <- make_cell("types_mixed", n, seed, miss_frac = 0.30, miss = "mcar", lambda = lambda, rho = rho,
                  thresholds = "fixed", driver = TRUE)
block_traits <- default_block_traits(cell)
is_prp <- block_traits %in% names(cell$trait_types)[cell$trait_types == "proportion"]
Yt <- as.data.frame(transform_block(cell$df_miss, block_traits, is_prp)); rownames(Yt) <- rownames(cell$df_miss)
set.seed(seed + 909L)                                          # rubin_cell.R: arm_seed(909L) before the arm
t0 <- proc.time()[["elapsed"]]
mi <- pigauto::multi_impute(Yt, cell$tree, m = 20, draws_method = "posterior", log_transform = FALSE,
                            posterior_control = list(seed = seed + 909L))
# complete-data Sigma on the same latent scale (sample covariance of the fully observed truth) for comparison
Ytrue <- transform_block(cell$truth, block_traits, is_prp)
saveRDS(list(n = n, lambda = lambda, rho = rho, seed = seed, block_traits = block_traits,
             params = mi$posterior$params, start = mi$posterior$start, hyper = mi$posterior$hyper,
             diagnostics = mi$posterior$diagnostics, converged = mi$posterior$converged,
             truth_cov = stats::cov(Ytrue), truth_cor = stats::cor(Ytrue),
             wall_s = proc.time()[["elapsed"]] - t0,
             pigauto_sha = pigauto_install_sha()), f_out)
cat("done", basename(f_out), "\n")
