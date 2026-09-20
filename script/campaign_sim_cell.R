# script/campaign_sim_cell.R
#
# One cell of the corrected-design four-arm imputation simulation (docs/dev-log/arc/
# 2026-09-19-imputation-simulation-plan.md, sections D, E, M, P, superseded by the exact design in
# this lane's task brief). Builds the cell via make_cell() with lambda/rho/evo/thresholds/driver/miss,
# runs every requested arm via run_arms() (script/campaign_gnn_off_lib.R, shared with
# script/campaign_gnn_off_cell.R so both scripts dispatch arms identically), and writes ONE rds.
#
# Usage:
#   Rscript script/campaign_sim_cell.R --dgp types_mixed --n 300 --seed 1 --out results/ \
#           [--arms gnn_on,gnn_off,gnn_off_rphylopars,freq,bace,floor] \
#           [--lambda 1 --rho 0 --evo BM --miss mcar --frac 0.30 --thresholds sample --driver] \
#           [--epochs 2000 --bace_nitt 50000 --bace_burnin 10000 --bace_thin 25] [--smoke]
#
# Output filename: <dgp>_<evo>_l<lambda>_r<rho>_<miss><frac>_n<n>_s<seed>.rds
# Resumable: a cell whose rds already exists is skipped.

suppressPackageStartupMessages({
  library(pigauto)
  library(ape)
})

RNGkind("L'Ecuyer-CMRG")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1L], "--")) return(TRUE)
  args[i + 1L]
}
dgp        <- get_arg("--dgp", "types_mixed")
n          <- as.integer(get_arg("--n", 100L))
seed       <- as.integer(get_arg("--seed", 1L))
out        <- get_arg("--out", "results")
arms       <- strsplit(get_arg("--arms", "gnn_on,gnn_off,gnn_off_rphylopars,freq,bace,floor"), ",")[[1]]
lambda     <- as.numeric(get_arg("--lambda", 1))
rho        <- as.numeric(get_arg("--rho", 0))
evo        <- get_arg("--evo", NULL)   # NULL -> make_dgp() default per dgp (BM, or OU for ou_mixed)
miss       <- get_arg("--miss", "mcar")
frac       <- as.numeric(get_arg("--frac", 0.30))
thresholds <- get_arg("--thresholds", "sample")
driver     <- isTRUE(get_arg("--driver", FALSE))
epochs     <- as.integer(get_arg("--epochs", 2000L))
smoke      <- isTRUE(get_arg("--smoke", FALSE))
bace_nitt   <- as.integer(get_arg("--bace_nitt", 50000L))
bace_burnin <- as.integer(get_arg("--bace_burnin", 10000L))
bace_thin   <- as.integer(get_arg("--bace_thin", 25L))
if (smoke) { epochs <- 20L; bace_nitt <- 600L; bace_burnin <- 100L; bace_thin <- 5L }

# mar requires an always-observed driver trait
if (miss == "mar" && !driver) { driver <- TRUE }

dir.create(out, showWarnings = FALSE, recursive = TRUE)
default_evo <- if (dgp == "ou_mixed") "OU" else "BM"
evo_tag <- if (is.null(evo)) default_evo else evo
tag <- sprintf("%s_%s_l%s_r%s_%s%s_n%d_s%d", dgp, evo_tag, format(lambda), format(rho), miss, format(frac), n, seed)
out_path <- file.path(out, paste0(tag, if (smoke) "_smoke" else "", ".rds"))
log_line <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

if (file.exists(out_path)) { log_line("cell %s: rds exists, skipping (resume)", tag); quit(save = "no", status = 0) }

# Cap torch threads per D-143/D-249-style sharing rules; harmless when torch is not loaded.
if ("torch" %in% loadedNamespaces() || requireNamespace("torch", quietly = TRUE)) {
  try({
    torch::torch_set_num_threads(as.integer(Sys.getenv("PIG_TORCH_THREADS", "4")))
    torch::torch_set_num_interop_threads(1L)
  }, silent = TRUE)
}

# ---- data, mask, scoring, arm dispatch: shared with campaign_gnn_off_cell.R -------------------
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "campaign_gnn_off_lib.R"))

git_hash <- tryCatch({
  h <- system2("git", c("-C", dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "rev-parse", "HEAD"),
               stdout = TRUE, stderr = FALSE)
  if (length(h) == 1 && nzchar(h)) h else NA_character_
}, error = function(e) NA_character_)

cell_data <- make_cell(dgp, n, seed, miss_frac = frac, miss = miss, lambda = lambda, rho = rho,
                        evo = evo, thresholds = thresholds, driver = driver)
truth <- cell_data$truth; tree <- cell_data$tree; mask <- cell_data$mask
df_miss <- cell_data$df_miss; cont_traits <- cell_data$cont_traits; trait_types <- cell_data$trait_types
log_line("cell %s: n=%d traits=%d (continuous %d) masked=%d realised_frac=%.3f", tag, nrow(truth),
         ncol(truth), length(cont_traits), sum(mask), cell_data$realised_frac)

cell_obj <- list(truth = truth, tree = tree, mask = mask, df_miss = df_miss,
                  cont_traits = cont_traits, trait_types = trait_types)
opts <- list(seed = seed, epochs = epochs, bace_nitt = bace_nitt, bace_burnin = bace_burnin,
             bace_thin = bace_thin, log_line = log_line)
out_arms <- run_arms(cell_obj, arms, opts)
tab <- out_arms$results
if (!is.null(tab)) { tab$dgp <- dgp; tab$n <- nrow(truth); tab$seed <- seed }

cell <- list(tag = tag, dgp = dgp, n = nrow(truth), seed = seed, arms = arms, smoke = smoke,
             lambda = lambda, rho = rho, evo = evo_tag, thresholds = thresholds, driver = driver,
             miss = miss, miss_frac = frac, realised_frac = cell_data$realised_frac,
             epochs = epochs, bace = c(nitt = bace_nitt, burnin = bace_burnin, thin = bace_thin),
             results = tab, calib = out_arms$calib, walls = unlist(out_arms$walls),
             failed = out_arms$failed, errors = out_arms$errors, paths = out_arms$paths,
             L = cell_data$L, truth = truth, mask = mask,
             RNGkind = RNGkind(), sessionInfo = utils::sessionInfo(),
             pigauto_version = as.character(utils::packageVersion("pigauto")),
             git_hash = git_hash, host = Sys.info()[["nodename"]], time = Sys.time())
saveRDS(cell, out_path)
print(tab, digits = 3); print(round(unlist(out_arms$walls), 1))
if (length(out_arms$errors)) { cat("ERRORS:\n"); print(out_arms$errors) }
