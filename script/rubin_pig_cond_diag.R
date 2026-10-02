# script/rubin_pig_cond_diag.R
#
# Diagnostic for pig_post's slope bias at n = 1000, lambda = 1, rho = 0.5 (docs/dev-log/arc/2026-10-02-rubin-pigauto-
# campaign.md). Refits one campaign cell exactly for pig_post (seed offset 909) and freqA (offset 101) and keeps their
# M = 20 completed datasets, which the campaign files did not, plus the truth, mask and tree. The exact conditional
# under the true model is computed afterwards (script/rubin_pig_cond_summary.R).
#
#   Rscript script/rubin_pig_cond_diag.R --n 1000 --lambda 1 --rho 0.5 --seed 1 --out <dir>

RNGkind("L'Ecuyer-CMRG")
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[i + 1L] }
n <- as.integer(get_arg("--n")); lambda <- as.numeric(get_arg("--lambda")); rho <- as.numeric(get_arg("--rho"))
seed <- as.integer(get_arg("--seed")); out <- get_arg("--out")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
f_out <- file.path(out, sprintf("cond_diag_n%d_l%s_r%s_s%d.rds", n, format(lambda), format(rho), seed))
if (file.exists(f_out)) quit(save = "no", status = 0)

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))
suppressMessages({ for (f in c("campaign_gnn_off_lib.R", "rubin_lib.R", "rubin_freq.R", "rubin_pigauto.R")) source(file.path(here, f)) })

cell <- make_cell("types_mixed", n, seed, miss_frac = 0.30, miss = "mcar", lambda = lambda, rho = rho,
                  thresholds = "fixed", driver = TRUE)
block_traits <- default_block_traits(cell)
keep <- function(sets) lapply(sets, function(d) d[, block_traits, drop = FALSE])

set.seed(seed + 101L); fa <- mi_freq_A(cell, 20)
set.seed(seed + 909L); pp <- mi_pig_post(cell, 20, seed = seed + 909L)
saveRDS(list(n = n, lambda = lambda, rho = rho, seed = seed, block_traits = block_traits,
             truth = cell$truth[, block_traits], mask = cell$mask[, block_traits], tree = cell$tree,
             freqA = keep(Filter(Negate(is.null), fa$datasets)), pig_post = keep(pp$datasets),
             pig_diag = pp$diag[c("converged", "rhat_max", "ess_min")]), f_out)
cat("done", basename(f_out), "\n")
