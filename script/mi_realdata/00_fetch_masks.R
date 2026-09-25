#!/usr/bin/env Rscript
# Extract mask receipts from the Mondrian real-data harness on branch
# arc/mondrian-realdata into
# script/mi_realdata/inputs/<dataset>-<arm>-m<seed>/mask_receipt.rds, and,
# with --with-conformal, a small derived conformal_metrics.rds beside it.
#
# mask_receipt.rds is list(dataset, seed, arm, truth, masked, mask, tree,
# structured) -- see script/mondrian_confirmation/01_run_masked_confirmation.R
# on arc/mondrian-realdata. `truth` is the pre-masking observed data (used
# for the reference PGLS slope); `masked` is what multi_impute() sees;
# `mask` is a logical matrix, TRUE where an originally-observed cell was
# masked out for this arm/seed.
#
# conformal_metrics.rds (--with-conformal, added 2026-09-24) holds, for the
# split and the Mondrian receipt of the same cell, only the per-trait
# summary that 01_run.R needs (lib.R::derive_conformal_metrics(): trait,
# n_masked, n_interval, coverage, mean_width, median_width), plus the
# source ref, commit and blob ids. The raw split.rds / mondrian.rds objects
# are never copied. With these files in place, 01_run.R runs without git
# (MI_REALDATA_OFFLINE=1), which Totoro and DRAC compute nodes need.
#
# Usage (on a machine with the git history, e.g. the Mac):
#   Rscript 00_fetch_masks.R --all [--with-conformal]
#   Rscript 00_fetch_masks.R dataset arm seed [--with-conformal]
#
# --all fetches all 10 planned cells (script/mi_realdata/lib.R::planned_cells()).
# Existing files are skipped, never overwritten.

# Rscript passes a script path containing spaces as `~+~` in --file=
# (e.g. ".../Github~+~Local/..."); undo that before normalizePath().
here <- function() {
  a <- commandArgs(FALSE)
  f <- gsub("~+~", " ", sub("^--file=", "", a[grepl("^--file=", a)]), fixed = TRUE)
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
source(file.path(here(), "lib.R"))

args <- commandArgs(trailingOnly = TRUE)
with_conformal <- "--with-conformal" %in% args
args <- setdiff(args, "--with-conformal")
if (length(args) == 1L && identical(args[[1L]], "--all")) {
  todo <- planned_cells()
} else if (length(args) == 3L) {
  todo <- data.frame(dataset = args[[1L]], arm = args[[2L]], seed = as.integer(args[[3L]]))
} else {
  stop("expected: --all | dataset arm seed, optionally with --with-conformal", call. = FALSE)
}

inputs_dir <- file.path(here(), "inputs")

fetch_one <- function(dataset, arm, seed) {
  name <- cell_name(dataset, arm, seed)
  out_file <- file.path(inputs_dir, name, "mask_receipt.rds")
  if (file.exists(out_file)) {
    cat(sprintf("[skip] %s already present at %s\n", name, out_file))
    return(invisible(TRUE))
  }
  rel <- file.path(mondrian_prefix, name, "mask_receipt.rds")
  git_show_to_file(rel, out_file)
  chk <- readRDS(out_file)
  stopifnot(is.list(chk), is.data.frame(chk$truth), is.data.frame(chk$masked),
            is.matrix(chk$mask), inherits(chk$tree, "phylo"))
  cat(sprintf("[fetched] %s -> %s (%d tips x %d traits)\n",
              name, out_file, nrow(chk$truth), ncol(chk$truth)))
  invisible(TRUE)
}

# A git or read failure stops before anything is written. A source receipt
# whose own status is not "ok" is written as it is (01_run.R records it and
# G8 then fails naming it) and flagged here with [WARN].
fetch_conformal_one <- function(dataset, arm, seed) {
  name <- cell_name(dataset, arm, seed)
  out_file <- conformal_metrics_path(inputs_dir, name)
  if (file.exists(out_file)) {
    cat(sprintf("[skip] %s conformal metrics already present at %s\n", name, out_file))
    return(invisible(TRUE))
  }
  cm <- fetch_conformal_metrics(dataset, arm, seed)
  if (is.na(cm$source$split_blob) || is.na(cm$source$mondrian_blob)) {
    stop("could not resolve the split/mondrian blob ids on ", mondrian_ref, call. = FALSE)
  }
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(cm, out_file)
  tag <- if (identical(cm$split$status, "ok") && identical(cm$mondrian$status, "ok")) "fetched" else "WARN"
  cat(sprintf("[%s] %s conformal metrics -> %s (split=%s, mondrian=%s, %d traits, %d bytes)\n",
              tag, name, out_file, cm$split$status, cm$mondrian$status,
              if (is.null(cm$split$metrics)) 0L else nrow(cm$split$metrics), file.size(out_file)))
  invisible(TRUE)
}

ok <- TRUE
for (i in seq_len(nrow(todo))) {
  nm <- cell_name(todo$dataset[[i]], todo$arm[[i]], todo$seed[[i]])
  steps <- list(masks = fetch_one)
  if (with_conformal) steps$conformal <- fetch_conformal_one
  for (st in names(steps)) {
    tryCatch(
      steps[[st]](todo$dataset[[i]], todo$arm[[i]], todo$seed[[i]]),
      error = function(e) {
        ok <<- FALSE
        cat(sprintf("[FAIL] %s (%s): %s\n", nm, st, conditionMessage(e)))
      }
    )
  }
}
if (!ok) stop("one or more inputs failed to fetch (see [FAIL] lines above)", call. = FALSE)
cat("FETCH_MASKS_OK\n")
