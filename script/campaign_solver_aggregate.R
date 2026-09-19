# script/campaign_solver_aggregate.R
# Aggregate the arc-C solver diagnostic cells. Usage: Rscript script/campaign_solver_aggregate.R <results_dir> <out_prefix>
args <- commandArgs(trailingOnly = TRUE)
res_dir <- args[1]; out <- args[2]
fs <- list.files(res_dir, pattern = "\\.rds$", full.names = TRUE)
cells <- lapply(fs, readRDS); cat(length(cells), "cells read\n")
tab <- do.call(rbind, lapply(cells, function(c) c$results))
walls <- do.call(rbind, lapply(cells, function(c) data.frame(dgp = c$dgp, n = c$n, seed = c$seed, arm = names(c$walls), wall = as.numeric(unlist(c$walls)))))
errs <- do.call(rbind, lapply(cells, function(c) if (length(c$errors)) data.frame(dgp = c$dgp, n = c$n, seed = c$seed, arm = names(c$errors), error = unlist(c$errors)) else NULL))
mcse <- function(x) stats::sd(x) / sqrt(length(x))
# fisher_ml fallback detection: identical continuous results to inhouse in the same cell
fb <- do.call(rbind, lapply(cells, function(c) {
  r <- c$results[c$results$metric == "zRMSE", ]
  a <- r[r$arm == "inhouse_pure", ]; b <- r[r$arm == "fisher_pure", ]
  if (!nrow(a) || !nrow(b)) return(NULL)
  m <- merge(a, b, by = "trait")
  data.frame(dgp = c$dgp, n = c$n, seed = c$seed, fisher_fell_back = all(abs(m$value.x - m$value.y) < 1e-10))
}))
cell_arm <- do.call(rbind, lapply(split(tab, list(tab$dgp, tab$n, tab$seed, tab$arm), drop = TRUE), function(d)
  data.frame(dgp = d$dgp[1], n = d$n[1], seed = d$seed[1], arm = d$arm[1],
             zrmse = mean(d$value[d$metric == "zRMSE"], na.rm = TRUE),
             acc = if (any(d$metric == "accuracy" & !is.na(d$value))) mean(d$value[d$metric == "accuracy"], na.rm = TRUE) else NA_real_,
             coverage = if (all(is.na(d$coverage))) NA_real_ else mean(d$coverage, na.rm = TRUE))))
summary_arm <- do.call(rbind, lapply(split(cell_arm, list(cell_arm$dgp, cell_arm$n, cell_arm$arm), drop = TRUE), function(d)
  data.frame(dgp = d$dgp[1], n = d$n[1], arm = d$arm[1], n_seeds = nrow(d), zrmse = mean(d$zrmse), zrmse_mcse = mcse(d$zrmse),
             acc = mean(d$acc, na.rm = TRUE), acc_mcse = if (all(is.na(d$acc))) NA_real_ else mcse(d$acc[!is.na(d$acc)]),
             coverage = mean(d$coverage, na.rm = TRUE))))
summary_arm <- summary_arm[order(summary_arm$dgp, summary_arm$n, summary_arm$zrmse), ]
# paired differences vs raw_rphylopars per cell (same mask): negative = pigauto arm better
paired <- do.call(rbind, lapply(split(cell_arm, list(cell_arm$dgp, cell_arm$n, cell_arm$seed), drop = TRUE), function(d) {
  ref <- d$zrmse[d$arm == "raw_rphylopars"]; if (!length(ref)) return(NULL)
  data.frame(dgp = d$dgp, n = d$n, seed = d$seed, arm = d$arm, diff = d$zrmse - ref) }))
paired_sum <- do.call(rbind, lapply(split(paired, list(paired$dgp, paired$n, paired$arm), drop = TRUE), function(d)
  data.frame(dgp = d$dgp[1], n = d$n[1], arm = d$arm[1], diff_vs_raw = mean(d$diff), diff_mcse = mcse(d$diff))))
per_trait <- do.call(rbind, lapply(split(tab, list(tab$dgp, tab$n, tab$arm, tab$trait, tab$metric), drop = TRUE), function(d)
  data.frame(dgp = d$dgp[1], n = d$n[1], arm = d$arm[1], trait = d$trait[1], metric = d$metric[1], n_seeds = nrow(d),
             mean = mean(d$value, na.rm = TRUE), mcse = mcse(d$value[!is.na(d$value)]))))
wall_sum <- aggregate(wall ~ dgp + n + arm, walls, mean)
saveRDS(list(tab = tab, cell_arm = cell_arm, summary_arm = summary_arm, paired = paired, paired_sum = paired_sum,
             per_trait = per_trait, walls = walls, wall_sum = wall_sum, errors = errs, fisher_fallback = fb, n_cells = length(cells)),
        paste0(out, ".rds"))
write.csv(summary_arm, paste0(out, "_summary_arm.csv"), row.names = FALSE)
cat("errors:", if (is.null(errs)) 0 else nrow(errs), "\n")
if (!is.null(fb)) { cat("fisher_ml fell back (identical to single_pass) in", sum(fb$fisher_fell_back), "of", nrow(fb), "cells\n") }
print(summary_arm[summary_arm$arm != "floor", ], digits = 3, row.names = FALSE)
cat("\npaired difference vs raw_rphylopars (negative = better):\n")
print(paired_sum[paired_sum$arm %in% c("inhouse_pure","rphylo_pure","rphylo_def","fisher_pure","bayes_pure","cont_only_pure"), ], digits = 3, row.names = FALSE)
