#!/usr/bin/env Rscript
# Design table for the four-arm imputation simulation (arc/imputation-sim).
# One row per cell. Drivers read this; gates G10/G11 count against it.
#
#   Rscript script/campaign_sim_design.R --stage core     > cells_core.csv
#   Rscript script/campaign_sim_design.R --stage factorial > cells_factorial.csv
#   Rscript script/campaign_sim_design.R --stage prerun   > cells_prerun.csv
#
# Stages
#   core      BM, MCAR 0.30, lambda {0.3, 0.7, 1}, rho {0, 0.5}, n {100, 300, 1000}     18 cells
#   factorial evo {BM, OU} x lambda {0.3, 1} x rho {0, 0.5} x n {100, 1000}
#             x mechanism {mcar 0.1, mcar 0.3, mar 0.3, clade 0.3}, minus core overlap  56 cells
#             (Q7 trim: n = 300 and lambda = 0.7 appear only in the core)
#   prerun    the 12 core cells at n {100, 1000} + 4 sentinels at n = 1000
#             (OU lambda 0.3, MAR, clade, MCAR 0.1)                                     16 cells
#   avonet    the bundled AVONET300 data, every arm                                      1 cell
# Replicates: 200 for arms freq, gnn_off, gnn_off_rphylopars, gnn_on, floor; BACE on seeds 1..100.

args <- commandArgs(trailingOnly = TRUE)
get <- function(flag, default) { i <- match(flag, args); if (is.na(i)) default else args[i + 1L] }
stage <- get("--stage", "core")

cell_key <- function(d) sprintf("%s_%s_l%s_r%s_%s%s_n%d", d$dgp, d$evo, d$lambda, d$rho, d$miss, d$frac, d$n)

core <- expand.grid(dgp = "types_mixed", evo = "BM", lambda = c(0.3, 0.7, 1), rho = c(0, 0.5),
                    miss = "mcar", frac = 0.3, n = c(100L, 300L, 1000L), stringsAsFactors = FALSE)

fact <- expand.grid(dgp = "types_mixed", evo = c("BM", "OU"), lambda = c(0.3, 1), rho = c(0, 0.5),
                    mech = c("mcar0.1", "mcar0.3", "mar0.3", "clade0.3"), n = c(100L, 1000L),
                    stringsAsFactors = FALSE)
fact$miss <- sub("[0-9.]+$", "", fact$mech)
fact$frac <- as.numeric(sub("^[a-z]+", "", fact$mech))
fact$mech <- NULL
fact <- fact[, names(core)]
fact <- fact[!(cell_key(fact) %in% cell_key(core)), ]

sent <- data.frame(dgp = "types_mixed",
                   evo    = c("OU",  "BM",  "BM",    "BM"),
                   lambda = c(0.3,   1,     1,       1),
                   rho    = c(0.5,   0.5,   0.5,     0.5),
                   miss   = c("mcar", "mar", "clade", "mcar"),
                   frac   = c(0.3,   0.3,   0.3,     0.1),
                   n = 1000L, stringsAsFactors = FALSE)
prerun <- rbind(core[core$n %in% c(100L, 1000L), ], sent)

avonet <- data.frame(dgp = "avonet", evo = "BM", lambda = 1, rho = 0, miss = "mcar", frac = 0.3, n = 300L,
                     stringsAsFactors = FALSE)

tab <- switch(stage, core = core, factorial = fact, prerun = prerun, avonet = avonet,
              stop("unknown --stage ", stage))
tab$stage <- stage
tab$cell <- cell_key(tab)
tab$reps <- if (stage == "prerun") 1L else if (stage == "avonet") 20L else 200L
tab$bace_reps <- if (stage == "prerun") 1L else if (stage == "avonet") 20L else 100L
tab <- tab[, c("stage", "cell", "dgp", "evo", "lambda", "rho", "miss", "frac", "n", "reps", "bace_reps")]
write.csv(tab, stdout(), row.names = FALSE, quote = FALSE)
