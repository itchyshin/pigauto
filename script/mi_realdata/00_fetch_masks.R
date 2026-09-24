#!/usr/bin/env Rscript
# Extract mask receipts (mask_receipt.rds ONLY -- not split.rds/mondrian.rds)
# from the Mondrian real-data harness on branch arc/mondrian-realdata into
# script/mi_realdata/inputs/<dataset>-<arm>-m<seed>/mask_receipt.rds.
#
# mask_receipt.rds is list(dataset, seed, arm, truth, masked, mask, tree,
# structured) -- see script/mondrian_confirmation/01_run_masked_confirmation.R
# on arc/mondrian-realdata. `truth` is the pre-masking observed data (used
# for the reference PGLS slope); `masked` is what multi_impute() sees;
# `mask` is a logical matrix, TRUE where an originally-observed cell was
# masked out for this arm/seed.
#
# Usage:
#   Rscript 00_fetch_masks.R --all
#   Rscript 00_fetch_masks.R dataset arm seed
#
# --all fetches all 10 planned cells (script/mi_realdata/lib.R::planned_cells()).

here <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
source(file.path(here(), "lib.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1L && identical(args[[1L]], "--all")) {
  todo <- planned_cells()
} else if (length(args) == 3L) {
  todo <- data.frame(dataset = args[[1L]], arm = args[[2L]], seed = as.integer(args[[3L]]))
} else {
  stop("expected: --all | dataset arm seed", call. = FALSE)
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

ok <- TRUE
for (i in seq_len(nrow(todo))) {
  res <- tryCatch(
    fetch_one(todo$dataset[[i]], todo$arm[[i]], todo$seed[[i]]),
    error = function(e) {
      ok <<- FALSE
      cat(sprintf("[FAIL] %s: %s\n", cell_name(todo$dataset[[i]], todo$arm[[i]], todo$seed[[i]]),
                  conditionMessage(e)))
      FALSE
    }
  )
}
if (!ok) stop("one or more mask receipts failed to fetch (see [FAIL] lines above)", call. = FALSE)
cat("FETCH_MASKS_OK\n")
