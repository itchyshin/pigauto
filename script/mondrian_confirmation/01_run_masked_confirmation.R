#!/usr/bin/env Rscript
# Usage: Rscript 01_run_masked_confirmation.R input.rds output_dir [seed] [epochs]
# input.rds is list(data = data.frame, tree = ape::phylo, dataset = character).
# It deliberately masks only originally observed cells and uses the identical
# mask for split and Mondrian fits.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("expected: input.rds output_dir [seed] [epochs]", call. = FALSE)
input <- readRDS(args[[1L]])
if (!is.list(input) || !is.data.frame(input$data) || !inherits(input$tree, "phylo")) {
  stop("input must be list(data = data.frame, tree = phylo, dataset = character)", call. = FALSE)
}
seed <- if (length(args) >= 3L) as.integer(args[[3L]]) else 20260818L
epochs <- if (length(args) >= 4L) as.integer(args[[4L]]) else 500L
`%||%` <- function(x, y) if (is.null(x)) y else x
if (is.null(rownames(input$data)) || !identical(rownames(input$data), input$tree$tip.label)) {
  stop("input data rows must exactly match tree tip order", call. = FALSE)
}
dir.create(args[[2L]], recursive = TRUE, showWarnings = FALSE)
set.seed(seed)
truth <- input$data
mask <- matrix(FALSE, nrow(truth), ncol(truth), dimnames = dimnames(truth))
for (nm in names(truth)) {
  observed <- which(!is.na(truth[[nm]]))
  if (length(observed) < 20L) next
  mask[sample(observed, ceiling(0.2 * length(observed))), nm] <- TRUE
}
masked <- truth
for (nm in names(masked)) masked[[nm]][mask[, nm]] <- NA
saveRDS(list(dataset = input$dataset %||% basename(args[[1L]]), seed = seed,
             truth = truth, masked = masked, mask = mask, tree = input$tree),
        file.path(args[[2L]], "mask_receipt.rds"))
one_method <- function(method) tryCatch({
  t0 <- proc.time()[["elapsed"]]
  fit <- pigauto::impute(masked, input$tree, seed = seed, epochs = epochs,
                         n_imputations = 1L, verbose = FALSE,
                         conformal_method = method)
  lo <- fit$prediction$conformal_lower; hi <- fit$prediction$conformal_upper
  rows <- lapply(names(truth), function(nm) {
    i <- which(mask[, nm]); if (!length(i) || !is.numeric(truth[[nm]]) || is.null(lo) || !(nm %in% colnames(lo))) return(NULL)
    valid <- i[is.finite(as.numeric(truth[[nm]][i])) & is.finite(lo[i, nm]) & is.finite(hi[i, nm])]
    if (!length(valid)) return(NULL)
    data.frame(trait = nm, n_masked = length(i), n_interval = length(valid),
      coverage = mean(truth[[nm]][valid] >= lo[valid, nm] & truth[[nm]][valid] <= hi[valid, nm]),
      width = mean(hi[valid, nm] - lo[valid, nm]), stringsAsFactors = FALSE)
  })
  mondrian <- if (identical(method, "mondrian")) fit$fit$conformal_mondrian else NULL
  list(status = "ok", method = method, elapsed_s = proc.time()[["elapsed"]] - t0,
       metrics = do.call(rbind, Filter(Negate(is.null), rows)), mondrian = mondrian)
}, error = function(e) list(status = "error", method = method, error = conditionMessage(e)))
for (method in c("split", "mondrian")) saveRDS(one_method(method), file.path(args[[2L]], paste0(method, ".rds")))
