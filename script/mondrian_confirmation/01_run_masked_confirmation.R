#!/usr/bin/env Rscript
# Usage: Rscript 01_run_masked_confirmation.R input.rds output_dir arm [seed] [epochs] [methods]
# input.rds is list(data = data.frame, tree = ape::phylo, dataset = character).
# It deliberately masks only originally observed cells and uses the identical
# mask for split and Mondrian fits. `methods` is a comma-separated subset of
# `split,mondrian`, so expensive methods can be run as independent receipts.
# `arm` is one of:
#   mcar       - 20% of observed cells per trait masked uniformly at random.
#   structured - per trait, fit glm(missing_in_real_data ~ first k phylo
#                eigenvectors, binomial), k = min(10, n/20); 20% of observed
#                cells masked with probability proportional to the fitted
#                propensity (see docs/dev-log/mondrian-realdata/00-preregistration.md).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("expected: input.rds output_dir arm [seed] [epochs] [methods]", call. = FALSE)
if (requireNamespace("torch", quietly = TRUE)) {
  try(torch::torch_set_num_threads(1L), silent = TRUE)
  try(torch::torch_set_num_interop_threads(1L), silent = TRUE)
}
input <- readRDS(args[[1L]])
if (!is.list(input) || !is.data.frame(input$data) || !inherits(input$tree, "phylo")) {
  stop("input must be list(data = data.frame, tree = phylo, dataset = character)", call. = FALSE)
}
arm <- match.arg(args[[3L]], c("mcar", "structured"))
seed <- if (length(args) >= 4L) as.integer(args[[4L]]) else 20260818L
epochs <- if (length(args) >= 5L) as.integer(args[[5L]]) else 500L
methods <- if (length(args) >= 6L) trimws(strsplit(args[[6L]], ",", fixed = TRUE)[[1L]]) else c("split", "mondrian")
if (!length(methods) || anyDuplicated(methods) || !all(methods %in% c("split", "mondrian"))) {
  stop("methods must be a non-empty, comma-separated subset of split,mondrian", call. = FALSE)
}
`%||%` <- function(x, y) if (is.null(x)) y else x
if (is.null(rownames(input$data)) || !identical(rownames(input$data), input$tree$tip.label)) {
  stop("input data rows must exactly match tree tip order", call. = FALSE)
}
dir.create(args[[2L]], recursive = TRUE, showWarnings = FALSE)
set.seed(seed)
truth <- input$data

# ---- Structured-arm mask: propensity-weighted on real missingness -------
build_structured_mask <- function(truth, tree) {
  R <- cov2cor(ape::vcv(tree))
  eig <- eigen(R, symmetric = TRUE)
  n <- nrow(truth)
  mask <- matrix(FALSE, nrow(truth), ncol(truth), dimnames = dimnames(truth))
  k_used <- stats::setNames(integer(0), character(0))
  prop_summary <- list()
  for (nm in names(truth)) {
    observed <- which(!is.na(truth[[nm]]))
    if (length(observed) < 20L) next
    k <- max(1L, min(10L, floor(n / 20)))
    ev <- eig$vectors[, seq_len(k), drop = FALSE]
    colnames(ev) <- paste0("ev", seq_len(k))
    miss <- as.integer(is.na(truth[[nm]]))
    df <- as.data.frame(ev)
    df$miss <- miss
    fit <- tryCatch(
      stats::glm(stats::as.formula(paste("miss ~", paste(colnames(ev), collapse = " + "))),
                 data = df, family = stats::binomial()),
      error = function(e) NULL
    )
    prob <- if (is.null(fit)) {
      rep(1, length(observed))
    } else {
      p <- suppressWarnings(stats::predict(fit, newdata = df[observed, , drop = FALSE], type = "response"))
      p[!is.finite(p) | p <= 0] <- 1e-6
      p
    }
    n_mask <- ceiling(0.2 * length(observed))
    sel <- sample(observed, min(n_mask, length(observed)), prob = prob)
    mask[sel, nm] <- TRUE
    k_used[nm] <- k
    prop_summary[[nm]] <- as.list(summary(prob))
  }
  list(mask = mask, k = k_used, propensity_summary = prop_summary)
}

mask_file <- file.path(args[[2L]], "mask_receipt.rds")
if (file.exists(mask_file)) {
  prior <- readRDS(mask_file)
  valid_prior <- is.list(prior) &&
    identical(prior$seed, seed) &&
    identical(prior$arm, arm) &&
    identical(prior$dataset, input$dataset %||% basename(args[[1L]])) &&
    identical(rownames(prior$truth), rownames(truth)) &&
    identical(names(prior$truth), names(truth)) &&
    identical(prior$tree$tip.label, input$tree$tip.label)
  if (!valid_prior) {
    stop("existing mask receipt does not match this input, seed, arm, or dataset", call. = FALSE)
  }
  truth <- prior$truth
  masked <- prior$masked
  mask <- prior$mask
} else {
  if (identical(arm, "mcar")) {
    mask <- matrix(FALSE, nrow(truth), ncol(truth), dimnames = dimnames(truth))
    for (nm in names(truth)) {
      observed <- which(!is.na(truth[[nm]]))
      if (length(observed) < 20L) next
      mask[sample(observed, ceiling(0.2 * length(observed))), nm] <- TRUE
    }
    struct_receipt <- NULL
  } else {
    sm <- build_structured_mask(truth, input$tree)
    mask <- sm$mask
    struct_receipt <- list(k = sm$k, propensity_summary = sm$propensity_summary)
  }
  masked <- truth
  for (nm in names(masked)) masked[[nm]][mask[, nm]] <- NA
  saveRDS(list(dataset = input$dataset %||% basename(args[[1L]]), seed = seed, arm = arm,
               truth = truth, masked = masked, mask = mask, tree = input$tree,
               structured = if (identical(arm, "structured")) struct_receipt else NULL),
          mask_file)
}

# ---- Mirror of R/predict_pigauto.R:820-866 mondrian stratum logic -------
# Reading only: calls pigauto's internal mondrian_locality() against the
# fit's own training-observed mask, exactly as predict.pigauto_fit() does.
stratify_rows <- function(fit_obj, mondrian_info, D_sq, trait_map, row_idx_by_trait) {
  out <- list()
  if (is.null(mondrian_info) || is.null(D_sq) || is.null(fit_obj$X_scaled)) {
    for (nm in names(row_idx_by_trait)) out[[nm]] <- rep(NA_character_, length(row_idx_by_trait[[nm]]))
    return(out)
  }
  obs_mask_train <- !is.na(fit_obj$X_scaled)
  hold <- c(fit_obj$splits$val_idx, fit_obj$splits$test_idx)
  if (length(hold)) obs_mask_train[hold] <- FALSE
  for (tm in trait_map) {
    nm <- tm$name
    idx <- row_idx_by_trait[[nm]]
    if (is.null(idx) || !length(idx)) next
    mo <- mondrian_info[[nm]]
    if (is.null(mo) || isTRUE(mo$fallback)) {
      out[[nm]] <- rep(NA_character_, length(idx)); next
    }
    lc <- tm$latent_cols
    score_col <- if (identical(tm$type, "zi_count")) lc[2L] else lc[1L]
    obs_idx <- which(obs_mask_train[, score_col])
    if (!length(obs_idx)) { out[[nm]] <- rep(NA_character_, length(idx)); next }
    locality <- pigauto:::mondrian_locality(D_sq, obs_idx, idx, k = 5L)
    out[[nm]] <- ifelse(!is.finite(locality), NA_character_,
                         ifelse(locality > mo$threshold, "far", "near"))
  }
  out
}

one_method <- function(method, mondrian_info_for_stratum) tryCatch({
  t0 <- proc.time()[["elapsed"]]
  fit <- pigauto::impute(masked, input$tree, seed = seed, epochs = epochs,
                         n_imputations = 1L, verbose = FALSE,
                         conformal_method = method)
  lo <- fit$prediction$conformal_lower; hi <- fit$prediction$conformal_upper
  row_idx_by_trait <- list()
  cells <- lapply(names(truth), function(nm) {
    i <- which(mask[, nm]); if (!length(i) || !is.numeric(truth[[nm]]) || is.null(lo) || !(nm %in% colnames(lo))) return(NULL)
    valid <- i[is.finite(as.numeric(truth[[nm]][i])) & is.finite(lo[i, nm]) & is.finite(hi[i, nm])]
    if (!length(valid)) return(NULL)
    row_idx_by_trait[[nm]] <<- valid
    data.frame(trait = nm, row = valid, truth = truth[[nm]][valid],
               lo = lo[valid, nm], hi = hi[valid, nm], stringsAsFactors = FALSE)
  })
  cells <- do.call(rbind, Filter(Negate(is.null), cells))
  rows <- lapply(names(truth), function(nm) {
    i <- which(mask[, nm]); if (!length(i) || !is.numeric(truth[[nm]]) || is.null(lo) || !(nm %in% colnames(lo))) return(NULL)
    valid <- i[is.finite(as.numeric(truth[[nm]][i])) & is.finite(lo[i, nm]) & is.finite(hi[i, nm])]
    if (!length(valid)) return(NULL)
    data.frame(trait = nm, n_masked = length(i), n_interval = length(valid),
      coverage = mean(truth[[nm]][valid] >= lo[valid, nm] & truth[[nm]][valid] <= hi[valid, nm]),
      width = mean(hi[valid, nm] - lo[valid, nm]), stringsAsFactors = FALSE)
  })
  mondrian <- if (identical(method, "mondrian")) fit$fit$conformal_mondrian else NULL
  mi <- mondrian_info_for_stratum %||% mondrian
  if (!is.null(cells) && !is.null(mi)) {
    strat <- stratify_rows(fit$fit, mi, fit$fit$graph$D_sq, fit$fit$trait_map, row_idx_by_trait)
    cells$stratum <- unlist(lapply(seq_len(nrow(cells)), function(r) {
      s <- strat[[cells$trait[r]]]
      pos <- match(cells$row[r], row_idx_by_trait[[cells$trait[r]]])
      if (is.null(s) || is.na(pos)) NA_character_ else s[pos]
    }))
  } else if (!is.null(cells)) {
    cells$stratum <- NA_character_
  }
  list(status = "ok", method = method, arm = arm, elapsed_s = proc.time()[["elapsed"]] - t0,
       metrics = do.call(rbind, Filter(Negate(is.null), rows)), mondrian = mondrian,
       cells = cells)
}, error = function(e) list(status = "error", method = method, arm = arm, error = conditionMessage(e)))

`%||%` <- function(x, y) if (is.null(x)) y else x

# Mondrian's own conformal_mondrian is the stratification reference for BOTH
# methods (comparability): fit mondrian first when requested/available so
# split's cells can be labelled the same way.
mondrian_ref <- NULL
mondrian_result_file <- file.path(args[[2L]], "mondrian.rds")
if (file.exists(mondrian_result_file)) {
  prior_mond <- readRDS(mondrian_result_file)
  if (identical(prior_mond$status, "ok")) mondrian_ref <- prior_mond$mondrian
}
ordered_methods <- methods[order(match(methods, c("mondrian", "split")))]
for (method in ordered_methods) {
  result_file <- file.path(args[[2L]], paste0(method, ".rds"))
  if (!file.exists(result_file)) {
    res <- one_method(method, mondrian_ref)
    if (identical(method, "mondrian") && identical(res$status, "ok")) mondrian_ref <- res$mondrian
    saveRDS(res, result_file)
  } else if (identical(method, "mondrian")) {
    prior_mond <- readRDS(result_file)
    if (identical(prior_mond$status, "ok")) mondrian_ref <- prior_mond$mondrian
  }
}
