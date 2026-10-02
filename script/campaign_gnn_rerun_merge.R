# script/campaign_gnn_rerun_merge.R
#
# 2026-10-02 rerun of the with/without-GNN campaign: the GNN-off arms and the GNN-on arm ran as
# two passes of campaign_gnn_off_cell.R (the off pass needs no torch), writing one rds per cell to
# two directories. This merges them cell by cell so campaign_gnn_off_aggregate.R can read one
# directory. A cell missing from either pass is reported and skipped.
#
# Usage: Rscript script/campaign_gnn_rerun_merge.R <off_dir> <on_dir> <out_dir>
args <- commandArgs(trailingOnly = TRUE)
off_dir <- args[1]; on_dir <- args[2]; out_dir <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
keep <- function(d) { f <- list.files(d, pattern = "\\.rds$"); f[!grepl("_smoke", f)] }
f_off <- keep(off_dir); f_on <- keep(on_dir)
both <- intersect(f_off, f_on)
cat(length(f_off), "off cells,", length(f_on), "on cells,", length(both), "in both\n")
miss <- setdiff(union(f_off, f_on), both)
if (length(miss)) cat("missing from one pass:", paste(miss, collapse = ", "), "\n")
for (f in both) {
  a <- readRDS(file.path(off_dir, f)); b <- readRDS(file.path(on_dir, f))
  m <- a
  m$arms    <- c(a$arms, b$arms)
  m$results <- rbind(a$results, b$results)
  m$walls   <- c(a$walls, b$walls)
  m$errors  <- c(a$errors, b$errors)
  m$paths   <- c(a$paths, b$paths)
  m$epochs  <- b$epochs
  saveRDS(m, file.path(out_dir, f))
}
cat("merged", length(both), "cells into", out_dir, "\n")
