# script/campaign_gnn_off_tables.R
# Emit the markdown tables for the campaign results note from the aggregate rds.
# Usage: Rscript script/campaign_gnn_off_tables.R <agg_rds> <out_md>
args <- commandArgs(trailingOnly = TRUE)
a <- readRDS(args[1]); out <- args[2]
arm_label <- c(gnn_off_pure = "pigauto, GNN off, pure baseline", gnn_off = "pigauto, GNN off",
               rphylopars = "raw Rphylopars (continuous only)", gnn_on_full = "pigauto, GNN on, full baseline",
               gnn_on = "pigauto, GNN on (as shipped)", bace = "BACE (50k iterations, OVR)", floor = "mean / mode floor")
arm_order <- names(arm_label)
fmt <- function(m, s) ifelse(is.na(m), "", sprintf("%.3f (%.3f)", m, s))
s <- a$summary_arm
s$arm <- factor(s$arm, levels = arm_order); s <- s[order(s$dgp, s$n, s$arm), ]
lines <- character()
for (d in unique(s$dgp)) {
  sd <- s[s$dgp == d, ]
  ns <- sort(unique(sd$n))
  lines <- c(lines, "", sprintf("### %s: mean z-RMSE over the continuous traits, mean (MCSE) over seeds", d), "",
             paste0("| arm | ", paste(sprintf("n = %d (seeds %d)", ns, sapply(ns, function(nn) max(sd$n_seeds[sd$n == nn]))), collapse = " | "), " |"),
             paste0("|---|", paste(rep("---:", length(ns)), collapse = "|"), "|"))
  for (arm in arm_order) {
    row <- sapply(ns, function(nn) { r <- sd[sd$n == nn & sd$arm == arm, ]; if (nrow(r)) fmt(r$zrmse, r$zrmse_mcse) else "" })
    lines <- c(lines, paste0("| ", arm_label[[arm]], " | ", paste(row, collapse = " | "), " |"))
  }
  lines <- c(lines, "", sprintf("### %s: mean accuracy over the discrete traits, mean (MCSE); and production-interval coverage (nominal 0.95)", d), "",
             paste0("| arm | ", paste(sprintf("acc n = %d", ns), collapse = " | "), " | ", paste(sprintf("cov n = %d", ns), collapse = " | "), " |"),
             paste0("|---|", paste(rep("---:", 2 * length(ns)), collapse = "|"), "|"))
  for (arm in arm_order) {
    acc <- sapply(ns, function(nn) { r <- sd[sd$n == nn & sd$arm == arm, ]; if (nrow(r) && is.finite(r$acc)) fmt(r$acc, r$acc_mcse) else "" })
    cov <- sapply(ns, function(nn) { r <- sd[sd$n == nn & sd$arm == arm, ]; if (nrow(r) && is.finite(r$coverage)) sprintf("%.3f", r$coverage) else "" })
    if (all(acc == "") && all(cov == "")) next
    lines <- c(lines, paste0("| ", arm_label[[arm]], " | ", paste(acc, collapse = " | "), " | ", paste(cov, collapse = " | "), " |"))
  }
}
# per-trait table for avonet (the real data), continuous traits
pt <- a$per_trait[a$per_trait$dgp == "avonet" & a$per_trait$metric == "zRMSE" & a$per_trait$arm != "floor", ]
if (nrow(pt)) {
  lines <- c(lines, "", "### AVONET300 per continuous trait: z-RMSE, mean (MCSE)", "",
             paste0("| arm | ", paste(unique(pt$trait), collapse = " | "), " |"), paste0("|---|", paste(rep("---:", length(unique(pt$trait))), collapse = "|"), "|"))
  for (arm in arm_order[arm_order != "floor"]) {
    row <- sapply(unique(pt$trait), function(tr) { r <- pt[pt$arm == arm & pt$trait == tr, ]; if (nrow(r)) fmt(r$mean, r$mcse) else "" })
    lines <- c(lines, paste0("| ", arm_label[[arm]], " | ", paste(row, collapse = " | "), " |"))
  }
  pa <- a$per_trait[a$per_trait$dgp == "avonet" & a$per_trait$metric == "accuracy" & a$per_trait$arm %in% c("gnn_off_pure","gnn_off","gnn_on_full","gnn_on","bace"), ]
  lines <- c(lines, "", "### AVONET300 per discrete trait: accuracy, mean (MCSE)", "",
             paste0("| arm | ", paste(unique(pa$trait), collapse = " | "), " |"), paste0("|---|", paste(rep("---:", length(unique(pa$trait))), collapse = "|"), "|"))
  for (arm in c("gnn_off_pure","gnn_off","gnn_on_full","gnn_on","bace")) {
    row <- sapply(unique(pa$trait), function(tr) { r <- pa[pa$arm == arm & pa$trait == tr, ]; if (nrow(r)) fmt(r$mean, r$mcse) else "" })
    lines <- c(lines, paste0("| ", arm_label[[arm]], " | ", paste(row, collapse = " | "), " |"))
  }
}
# wall time
w <- a$walls[a$walls$arm != "gnn_on_full", ]
ws <- aggregate(wall ~ n + arm, w, mean); ws$arm <- factor(ws$arm, levels = arm_order); ws <- ws[order(ws$n, ws$arm), ]
ns <- sort(unique(ws$n))
lines <- c(lines, "", "### Wall time per arm, seconds, mean over DGPs and seeds (4 threads per cell)", "",
           paste0("| arm | ", paste(sprintf("n = %d", ns), collapse = " | "), " |"), paste0("|---|", paste(rep("---:", length(ns)), collapse = "|"), "|"))
for (arm in arm_order[arm_order != "gnn_on_full"]) {
  row <- sapply(ns, function(nn) { r <- ws[ws$n == nn & ws$arm == arm, ]; if (nrow(r)) sprintf("%.1f", r$wall) else "" })
  lines <- c(lines, paste0("| ", arm_label[[arm]], " | ", paste(row, collapse = " | "), " |"))
}
# dispatch paths
if (!is.null(a$paths)) {
  pt2 <- as.data.frame(table(dgp = a$paths$dgp, path = a$paths$path)); pt2 <- pt2[pt2$Freq > 0, ]
  lines <- c(lines, "", "### Baseline dispatch recorded by `fit$baseline$path` (GNN-off arm; trait-cells over all n and seeds)", "",
             "| dgp | path | count |", "|---|---|---:|", sprintf("| %s | %s | %d |", pt2$dgp, pt2$path, pt2$Freq))
}
lines <- c(lines, "", sprintf("Cells aggregated: %d. Errors: %d.", a$n_cells, if (is.null(a$errors)) 0L else nrow(a$errors)))
writeLines(lines, out); cat("wrote", out, length(lines), "lines\n")
