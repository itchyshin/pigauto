# script/campaign_gnn_off_aggregate.R
# Aggregate campaign cells (rds from campaign_gnn_off_cell.R) into per (dgp, n, arm, trait) tables with
# mean and Monte Carlo SE over seeds, wall-time summaries, coverage, dispatch paths, and error counts.
# Usage: Rscript script/campaign_gnn_off_aggregate.R <results_dir> <out_prefix>
args <- commandArgs(trailingOnly = TRUE)
res_dir <- args[1]; out <- args[2]
fs <- list.files(res_dir, pattern = "\\.rds$", full.names = TRUE)
fs <- fs[!grepl("_smoke", fs)]
cells <- lapply(fs, readRDS)
cat(length(cells), "cells read\n")
tab <- do.call(rbind, lapply(cells, function(c) { r <- c$results; if (is.null(r) || !nrow(r)) return(NULL); r$dgp <- c$dgp; r$n <- c$n; r$seed <- c$seed; r }))
walls <- do.call(rbind, lapply(cells, function(c) data.frame(dgp = c$dgp, n = c$n, seed = c$seed, arm = names(c$walls), wall = as.numeric(unlist(c$walls)))))
errs <- do.call(rbind, lapply(cells, function(c) if (length(c$errors)) data.frame(dgp = c$dgp, n = c$n, seed = c$seed, arm = names(c$errors), error = unlist(c$errors)) else NULL))
paths <- do.call(rbind, lapply(cells, function(c) { p <- c$paths$gnn_off; if (is.null(p)) return(NULL); data.frame(dgp = c$dgp, n = c$n, seed = c$seed, trait = names(p), path = as.character(p)) }))
mcse <- function(x) stats::sd(x) / sqrt(length(x))
# per trait
per_trait <- do.call(rbind, lapply(split(tab, list(tab$dgp, tab$n, tab$arm, tab$trait, tab$metric), drop = TRUE), function(d)
  data.frame(dgp = d$dgp[1], n = d$n[1], arm = d$arm[1], trait = d$trait[1], metric = d$metric[1], n_seeds = nrow(d),
             mean = mean(d$value, na.rm = TRUE), mcse = mcse(d$value[!is.na(d$value)]),
             coverage = if (all(is.na(d$coverage))) NA_real_ else mean(d$coverage, na.rm = TRUE),
             coverage_mcse = if (all(is.na(d$coverage))) NA_real_ else mcse(d$coverage[!is.na(d$coverage)]))))
# per cell-arm summary: mean over continuous traits (zRMSE) and discrete traits (accuracy), then over seeds
cell_arm <- do.call(rbind, lapply(split(tab, list(tab$dgp, tab$n, tab$seed, tab$arm), drop = TRUE), function(d)
  data.frame(dgp = d$dgp[1], n = d$n[1], seed = d$seed[1], arm = d$arm[1],
             zrmse = mean(d$value[d$metric == "zRMSE"], na.rm = TRUE),
             acc = if (any(d$metric == "accuracy" & !is.na(d$value))) mean(d$value[d$metric == "accuracy"], na.rm = TRUE) else NA_real_,
             coverage = if (all(is.na(d$coverage))) NA_real_ else mean(d$coverage, na.rm = TRUE))))
summary_arm <- do.call(rbind, lapply(split(cell_arm, list(cell_arm$dgp, cell_arm$n, cell_arm$arm), drop = TRUE), function(d)
  data.frame(dgp = d$dgp[1], n = d$n[1], arm = d$arm[1], n_seeds = nrow(d),
             zrmse = mean(d$zrmse), zrmse_mcse = mcse(d$zrmse),
             acc = mean(d$acc, na.rm = TRUE), acc_mcse = if (all(is.na(d$acc))) NA_real_ else mcse(d$acc[!is.na(d$acc)]),
             coverage = mean(d$coverage, na.rm = TRUE))))
summary_arm <- summary_arm[order(summary_arm$dgp, summary_arm$n, summary_arm$zrmse), ]
wall_sum <- aggregate(wall ~ dgp + n + arm, walls, function(x) c(mean = mean(x), max = max(x)))
path_tab <- if (!is.null(paths)) as.data.frame(table(paths$dgp, paths$n, paths$trait, paths$path)) else NULL
saveRDS(list(tab = tab, per_trait = per_trait, cell_arm = cell_arm, summary_arm = summary_arm, walls = walls,
             wall_sum = wall_sum, errors = errs, paths = paths, n_cells = length(cells)), paste0(out, ".rds"))
write.csv(summary_arm, paste0(out, "_summary_arm.csv"), row.names = FALSE)
write.csv(per_trait, paste0(out, "_per_trait.csv"), row.names = FALSE)
cat("errors:", if (is.null(errs)) 0 else nrow(errs), "\n")
print(summary_arm, digits = 3, row.names = FALSE)
