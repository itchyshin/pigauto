# script/campaign_gnn_off_aggregate.R
# Aggregate campaign cells (rds from campaign_gnn_off_cell.R OR campaign_sim_cell.R) into per
# (dgp, n, arm, trait) tables with mean and Monte Carlo SE over seeds, wall-time summaries, coverage,
# dispatch paths, and error counts. Section I additions (only fire when the new-shape fields are
# present -- i.e. rds written by campaign_sim_cell.R): per (cell, arm, trait, metric) MCSE, a paired
# difference vs `--reference` with MCSE of the per-replicate difference (common seeds only), pooled
# ECE (10 equal-mass bins over all masked cells of a cell) with a bootstrap-over-replicates MCSE,
# failure rate, and summary.csv / paired.csv / failures.csv / raw_index.csv + sessionInfo().
# Usage: Rscript script/campaign_gnn_off_aggregate.R <results_dir> <out_prefix> [--reference ARM]
args <- commandArgs(trailingOnly = TRUE)
get_flag <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L] }
reference_arm <- get_flag("--reference", "floor")
pos <- args[!grepl("^--", args)]; pos <- pos[!(pos %in% c(reference_arm))]
res_dir <- pos[1]; out <- pos[2]
# recursive: a pooled directory holds one subdirectory per machine, and the same (cell, seed) can
# appear in two of them carrying DIFFERENT arms (fast arms on Totoro, BACE on the clusters). Both
# are kept; their `results` rows carry distinct arms, so rbind below merges rather than duplicates.
fs <- list.files(res_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
fs <- fs[!grepl("_smoke", fs)]
cells <- lapply(fs, function(f) tryCatch(readRDS(f), error = function(e) NULL))
bad <- sum(vapply(cells, is.null, TRUE))
keep <- !vapply(cells, is.null, TRUE)
fs <- fs[keep]; cells <- cells[keep]
# Two machines can also hold the SAME (cell, seed) carrying the SAME arm, which rbind would count
# twice: measured 2026-09-20, 650 core cell-seeds existed on both nibi and fir, both arm = bace,
# with different values (a BACE re-run, not a copy). Dedupe on (filename, arm set) keeping the
# first host in listing order, and say how many were dropped. The legitimate cross-host case --
# same filename, DIFFERENT arms -- has a different key and is untouched.
dup_key <- vapply(seq_along(fs), function(i)
  paste(basename(fs[i]), paste(sort(as.character(cells[[i]]$arms)), collapse = ","), sep = "|"), "")
dup <- duplicated(dup_key)
if (any(dup)) {
  cat(sprintf("dropped %d duplicate (cell, seed, arm-set) files present on more than one machine\n", sum(dup)))
  print(utils::head(sort(table(sub("^[^|]*[|]", "", dup_key[dup])), decreasing = TRUE), 5))
  fs <- fs[!dup]; cells <- cells[!dup]
}
cat(length(cells), "cells read", if (bad) sprintf("(%d unreadable, skipped)", bad) else "", "\n")
if (!length(cells)) stop("no readable cells under ", res_dir)
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

# ---- Section I: campaign_sim_cell.R rds shape (lambda/rho/evo/miss/calib/failed present) ------
is_sim_shape <- length(cells) > 0 && !is.null(cells[[1]]$lambda)
if (is_sim_shape) {
  # `frac` belongs in the key. Without it MCAR 0.10 and MCAR 0.30 share the key "mcar" and are
  # averaged into one reported cell across the whole factorial, and the paired merge below forms a
  # Cartesian product across the two mechanisms. Found by the harness audit, 2026-09-21.
  regime_of <- function(c) data.frame(lambda = c$lambda %||% NA, rho = c$rho %||% NA, evo = c$evo %||% NA,
                                       miss = c$miss %||% NA, frac = c$miss_frac %||% NA,
                                       thresholds = c$thresholds %||% NA)
  `%||%` <- function(a, b) if (is.null(a)) b else a
  tab2 <- do.call(rbind, lapply(cells, function(c) {
    r <- c$results; if (is.null(r) || !nrow(r)) return(NULL)
    cbind(r, regime_of(c))
  }))
  mcse_boot <- function(x, R = 500L) {
    x <- x[is.finite(x)]; if (length(x) < 2) return(NA_real_)
    b <- replicate(R, mean(sample(x, length(x), replace = TRUE)))
    stats::sd(b)
  }
  cell_key <- c("dgp", "n", "lambda", "rho", "evo", "miss", "frac")
  summary_i <- do.call(rbind, lapply(split(tab2, list(tab2$dgp, tab2$n, tab2$lambda, tab2$rho, tab2$evo, tab2$miss, tab2$frac, tab2$arm, tab2$trait, tab2$metric), drop = TRUE), function(d)
    data.frame(d[1, cell_key], arm = d$arm[1], trait = d$trait[1], metric = d$metric[1], reps = nrow(d),
               mean = mean(d$value, na.rm = TRUE), mcse = mcse(d$value[!is.na(d$value)]))))
  # Coverage is NOT a metric row: the runner carries it as a side column on each zRMSE row, so the
  # split above never emitted it and _summary.csv held no coverage at all (measured 2026-09-20 --
  # the results board's coverage panel, half the pre-registered primary contrast, therefore rendered
  # its empty state on every version). Emit it as its own metric so consumers that filter on
  # `metric` see it, with the same mean-and-MCSE-over-seeds treatment as every other number.
  cov_src <- tab2[tab2$metric == "zRMSE" & !is.na(tab2$coverage), ]
  if (nrow(cov_src)) {
    cov_rows <- do.call(rbind, lapply(split(cov_src, list(cov_src$dgp, cov_src$n, cov_src$lambda, cov_src$rho, cov_src$evo, cov_src$miss, cov_src$frac, cov_src$arm, cov_src$trait), drop = TRUE), function(d)
      data.frame(d[1, cell_key], arm = d$arm[1], trait = d$trait[1], metric = "coverage", reps = nrow(d),
                 mean = mean(d$coverage, na.rm = TRUE), mcse = mcse(d$coverage[!is.na(d$coverage)]))))
    summary_i <- rbind(summary_i, cov_rows)
    cat(sprintf("coverage emitted as a metric: %d rows from %d zRMSE rows carrying it\n",
                nrow(cov_rows), nrow(cov_src)))
  }
  write.csv(summary_i, paste0(out, "_summary.csv"), row.names = FALSE)

  # paired difference vs reference_arm, on common seeds only
  paired_rows <- list()
  for (grp in split(tab2, list(tab2$dgp, tab2$n, tab2$lambda, tab2$rho, tab2$evo, tab2$miss, tab2$frac, tab2$trait, tab2$metric), drop = TRUE)) {
    ref <- grp[grp$arm == reference_arm, ]
    if (!nrow(ref)) next
    for (a in setdiff(unique(grp$arm), reference_arm)) {
      cur <- grp[grp$arm == a, ]
      m <- merge(ref[c("seed", "value")], cur[c("seed", "value")], by = "seed", suffixes = c(".ref", ".arm"))
      if (!nrow(m)) next
      d <- m$value.arm - m$value.ref
      paired_rows[[length(paired_rows) + 1L]] <- data.frame(grp[1, cell_key], trait = grp$trait[1], metric = grp$metric[1],
                                                              arm = a, reference = reference_arm, n_common_seeds = nrow(m),
                                                              mean_diff = mean(d), mcse_diff = mcse(d))
    }
  }
  paired <- do.call(rbind, paired_rows)
  if (!is.null(paired)) write.csv(paired, paste0(out, "_paired.csv"), row.names = FALSE)

  # pooled ECE per (cell, arm), 10 equal-mass bins over all masked cells, bootstrap MCSE over replicates
  calib_all <- do.call(rbind, lapply(cells, function(c) { cb <- c$calib; if (is.null(cb) || !nrow(cb)) return(NULL); cbind(cb, regime_of(c), dgp = c$dgp, n = c$n, seed = c$seed) }))
  ece <- NULL
  if (!is.null(calib_all)) {
    ece_one <- function(conf, correct, n_bins = 10L) {
      ord <- order(conf); conf <- conf[ord]; correct <- correct[ord]
      bins <- cut(seq_along(conf), breaks = n_bins, labels = FALSE)
      sum(vapply(seq_len(n_bins), function(b) {
        idx <- bins == b; if (!any(idx)) return(0)
        (sum(idx) / length(conf)) * abs(mean(correct[idx]) - mean(conf[idx]))
      }, numeric(1)))
    }
    # POOLED, as the methods note claims: bin every masked cell of the whole (cell, arm, trait)
    # together, then bin. The previous line computed ece_one() per replicate and averaged, which is
    # the estimator the note explicitly rejects: on a few dozen cells per replicate each per-seed ECE
    # carries a positive small-sample bias of roughly 0.29 that does not cancel under averaging, and
    # it grows with how many bins an arm's probabilities occupy, so it penalises the sharper arm.
    # Found by the harness audit, 2026-09-21, which measured a perfectly calibrated arm scoring
    # 0.108 under the old estimator against 0.007 under this one.
    #
    # `trait` joins the key as well: binary confidences floor at 0.5 and three-class at 0.33, so
    # pooling them into one number compares nothing.
    ece_key <- c(cell_key, "trait")
    ece <- do.call(rbind, lapply(split(calib_all, list(calib_all$dgp, calib_all$n, calib_all$lambda, calib_all$rho,
                                                        calib_all$evo, calib_all$miss, calib_all$frac,
                                                        calib_all$arm, calib_all$trait), drop = TRUE), function(d) {
      pooled <- ece_one(d$confidence, d$correct)
      # MCSE by bootstrapping REPLICATES and re-pooling, so the resample unit stays the replicate.
      seeds <- unique(d$seed)
      boot <- if (length(seeds) < 2) NA_real_ else stats::sd(replicate(200L, {
        pick <- sample(seeds, length(seeds), replace = TRUE)
        ece_one(unlist(lapply(pick, function(s) d$confidence[d$seed == s])),
                unlist(lapply(pick, function(s) d$correct[d$seed == s])))
      }))
      data.frame(d[1, ece_key], arm = d$arm[1], n_reps = length(seeds), n_cells = nrow(d),
                 ece = pooled, ece_mcse = boot)
    }))
    write.csv(ece, paste0(out, "_ece.csv"), row.names = FALSE)
  }

  # failure rate per (cell, arm)
  fail_rows <- do.call(rbind, lapply(cells, function(c) {
    if (!length(c$failed)) return(NULL)
    data.frame(regime_of(c), dgp = c$dgp, n = c$n, seed = c$seed, arm = names(c$failed))
  }))
  failure_rate <- if (!is.null(fail_rows)) do.call(rbind, lapply(split(fail_rows, list(fail_rows$dgp, fail_rows$n, fail_rows$arm), drop = TRUE), function(d)
    data.frame(dgp = d$dgp[1], n = d$n[1], arm = d$arm[1], n_failed = nrow(d)))) else NULL
  write.csv(if (is.null(failure_rate)) data.frame() else failure_rate, paste0(out, "_failures.csv"), row.names = FALSE)

  raw_index <- do.call(rbind, lapply(seq_along(cells), function(i) {
    c <- cells[[i]]
    data.frame(file = basename(fs[i]), dgp = c$dgp, n = c$n, seed = c$seed, lambda = c$lambda %||% NA,
               rho = c$rho %||% NA, evo = c$evo %||% NA, miss = c$miss %||% NA, host = c$host %||% NA,
               arms = paste(c$arms, collapse = ";"), n_failed = length(c$failed))
  }))
  write.csv(raw_index, paste0(out, "_raw_index.csv"), row.names = FALSE)
  saveRDS(utils::sessionInfo(), paste0(out, "_sessionInfo.rds"))
  cat("Section I outputs written:", paste0(out, c("_summary.csv", "_paired.csv", "_ece.csv", "_failures.csv", "_raw_index.csv")), "\n")
}
