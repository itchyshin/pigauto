# script/campaign_gnn_off_cell.R
#
# One cell of the with/without-GNN campaign (plan: ~/.claude/plans/valiant-forging-gray.md,
# section ARC B). Runs every arm on ONE (dgp, n, seed) cell with the SAME user-level 30% MCAR
# mask, scores on the masked cells, and writes one rds. A driver fans cells out over cores.
#
# Usage:
#   Rscript script/campaign_gnn_off_cell.R --dgp bm_mixed --n 100 --seed 1 --out results/ \
#           [--arms gnn_on,gnn_off,gnn_off_pure,rphylopars,bace,floor] \
#           [--epochs 2000] [--bace_nitt 50000 --bace_burnin 10000 --bace_thin 25] [--smoke]
#
# Arm "freq" = frequentist stack: Rphylopars (continuous family) + castor Mk (discrete traits).
# Arm "mf_phylo" = missForest on traits + phylogenetic eigenvectors (Gendre et al. 2024 hybrid).
# DGPs: bm_mixed (4 continuous BM + 1 binary + 1 categorical(3) on ape::rcoal(n)), types_mixed (every
#       pigauto trait type: 2 continuous, count, proportion, binary, ordinal, categorical),
#       ou_mixed (simulate_non_bm OU continuous + same discrete), bace_dgp (BACE::sim_bace: y + 2 gaussian
#       + 1 binary predictors, all imputed), avonet (bundled avonet300/tree300, n ignored). Arms and metrics as locked in the plan. `--smoke` shrinks epochs/chains so the
#       whole cell runs in seconds; it is for checking the invocation, never for results.

suppressPackageStartupMessages({
  library(pigauto)
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1L], "--")) return(TRUE)
  args[i + 1L]
}
dgp    <- get_arg("--dgp", "bm_mixed")
n      <- as.integer(get_arg("--n", 100L))
seed   <- as.integer(get_arg("--seed", 1L))
out    <- get_arg("--out", "results")
arms   <- strsplit(get_arg("--arms", "gnn_on,gnn_off,gnn_off_pure,rphylopars,bace,floor"), ",")[[1]]
# gnn_on_full (GNN-on predicting from the tax-free baseline) is derived from the gnn_on fit, not refit.
epochs <- as.integer(get_arg("--epochs", 2000L))
smoke  <- isTRUE(get_arg("--smoke", FALSE))
bace_nitt   <- as.integer(get_arg("--bace_nitt", 50000L))
bace_burnin <- as.integer(get_arg("--bace_burnin", 10000L))
bace_thin   <- as.integer(get_arg("--bace_thin", 25L))
miss_frac <- 0.30
if (smoke) { epochs <- 20L; bace_nitt <- 600L; bace_burnin <- 100L; bace_thin <- 5L }
dir.create(out, showWarnings = FALSE, recursive = TRUE)
tag <- sprintf("%s_n%d_s%d", dgp, n, seed)
log_line <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

# ---- data, mask, scoring: shared with the solver runner --------------------------------------
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "campaign_gnn_off_lib.R"))
cell_data <- make_cell(dgp, n, seed, miss_frac)
truth <- cell_data$truth; tree <- cell_data$tree; mask <- cell_data$mask
df_miss <- cell_data$df_miss; cont_traits <- cell_data$cont_traits; trait_types <- cell_data$trait_types
log_line("cell %s: n=%d traits=%d (continuous %d) masked=%d", tag, nrow(truth), ncol(truth),
         length(cont_traits), sum(mask))

# ---- arms ---------------------------------------------------------------------------------
# Arm-dispatch loop factored into run_arms() (script/campaign_gnn_off_lib.R, Section G) so this
# script and script/campaign_sim_cell.R share one implementation. Behaviour is unchanged: same
# arms, same pigauto calls, same scoring inputs; run_arms() additionally scores a failed arm at the
# floor (never dropped) and returns the calibration long-frame, both new and additive.
cell_obj <- list(truth = truth, tree = tree, mask = mask, df_miss = df_miss,
                  cont_traits = cont_traits, trait_types = trait_types)
opts <- list(seed = seed, epochs = epochs, bace_nitt = bace_nitt, bace_burnin = bace_burnin,
             bace_thin = bace_thin, log_line = log_line)
out_arms <- run_arms(cell_obj, arms, opts)
tab <- out_arms$results
if (!is.null(tab)) tab$dgp <- dgp else tab <- tab
if (!is.null(tab)) { tab$n <- nrow(truth); tab$seed <- seed }
walls <- out_arms$walls; errors <- out_arms$errors; paths <- out_arms$paths
cell <- list(tag = tag, dgp = dgp, n = nrow(truth), seed = seed, arms = arms, smoke = smoke,
             epochs = epochs, bace = c(nitt = bace_nitt, burnin = bace_burnin, thin = bace_thin),
             results = tab, calib = out_arms$calib, walls = unlist(walls), errors = errors, paths = paths,
             pigauto_version = as.character(utils::packageVersion("pigauto")),
             host = Sys.info()[["nodename"]], time = Sys.time())
saveRDS(cell, file.path(out, paste0(tag, if (smoke) "_smoke" else "", ".rds")))
print(tab, digits = 3); print(round(unlist(walls), 1))
if (length(errors)) { cat("ERRORS:\n"); print(errors) }
