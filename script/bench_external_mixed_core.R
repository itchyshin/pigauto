#!/usr/bin/env Rscript
# Retained Stage-B mixed AVONET comparator.  This is a runner, not evidence:
# it requires PIGAUTO_OUT_DIR and is never invoked by package tests or CI.
options(warn = 1, stringsAsFactors = FALSE)
`%||%` <- function(x, y) if (is.null(x)) y else x
suppressPackageStartupMessages({
  library(ape)
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH", unset = "")
  if (nzchar(pkg_path) && file.exists(file.path(pkg_path, "DESCRIPTION"))) {
    devtools::load_all(pkg_path, quiet = TRUE)
  } else library(pigauto)
})
if (requireNamespace("torch", quietly = TRUE)) {
  try(torch::torch_set_num_threads(1L), silent = TRUE)
  try(torch::torch_set_num_interop_threads(1L), silent = TRUE)
}

MISS_FRAC <- 0.30
MASK_SEEDS <- 20260901L:20260905L
CONT <- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")
DISC <- c("Trophic.Level", "Primary.Lifestyle", "Migration")
TRAITS <- c(CONT, DISC)
out_dir <- Sys.getenv("PIGAUTO_OUT_DIR", unset = "")
if (!nzchar(out_dir)) stop("Set PIGAUTO_OUT_DIR to a dedicated retained result directory", call. = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
receipt_dir <- file.path(out_dir, "receipts")
dir.create(receipt_dir, recursive = TRUE, showWarnings = FALSE)

load_data <- function() {
  e <- new.env(parent = emptyenv())
  utils::data("avonet300", package = "pigauto", envir = e)
  utils::data("tree300", package = "pigauto", envir = e)
  d <- e$avonet300; rownames(d) <- d$Species_Key; d$Species_Key <- NULL
  list(df = d[, TRAITS, drop = FALSE], tree = e$tree300)
}
make_mask <- function(d, seed) {
  set.seed(seed); z <- matrix(FALSE, nrow(d), ncol(d), dimnames = list(rownames(d), names(d)))
  for (v in names(d)) { o <- which(!is.na(d[[v]])); z[sample(o, ceiling(length(o) * MISS_FRAC)), v] <- TRUE }
  z
}
fit_safe <- function(expr) {
  t0 <- proc.time()[["elapsed"]]; x <- tryCatch(expr, error = function(e) e)
  list(value = if (inherits(x, "error")) NULL else x,
       error = if (inherits(x, "error")) conditionMessage(x) else NULL,
       wall_s = proc.time()[["elapsed"]] - t0)
}
fit_pigauto <- function(d, tree, seed, method) fit_safe(pigauto::impute(d, tree, seed = seed, verbose = FALSE, predict_method = method))
fit_missforest <- function(d) {
  if (!requireNamespace("missForest", quietly = TRUE)) return(list(value = NULL, error = "missForest not installed", wall_s = 0))
  fit_safe(missForest::missForest(d, verbose = FALSE)$ximp)
}
fit_rphylopars <- function(d, tree) {
  if (!requireNamespace("Rphylopars", quietly = TRUE)) return(list(value = NULL, error = "Rphylopars not installed", wall_s = 0))
  ans <- fit_safe(Rphylopars::phylopars(data.frame(species = rownames(d), d, check.names = FALSE), tree = tree, model = "BM", phylo_correlated = TRUE, pheno_correlated = TRUE, REML = TRUE))
  if (!is.null(ans$value)) ans$value <- tryCatch(as.data.frame(ans$value$anc_recon[rownames(d), names(d), drop = FALSE]), error = function(e) e)
  if (inherits(ans$value, "error")) { ans$error <- conditionMessage(ans$value); ans$value <- NULL }
  ans
}
fit_phylolm <- function(d, tree) {
  if (!requireNamespace("phylolm", quietly = TRUE)) return(list(value = NULL, error = "phylolm not installed", wall_s = 0))
  fit_safe({
    out <- d
    for (v in names(d)) {
      y <- d[[v]]; names(y) <- rownames(d); obs <- names(y)[!is.na(y)]
      if (length(obs) < 5L) stop(sprintf("%s has fewer than five observed species", v))
      tr <- ape::keep.tip(tree, obs)
      dat <- data.frame(y = unname(y[obs]), row.names = obs)
      mod <- phylolm::phylolm(y ~ 1, data = dat, phy = tr, model = "lambda")
      b <- unname(stats::coef(mod)[["(Intercept)"]]); R <- stats::cov2cor(ape::vcv(tree))[names(y), names(y)]
      oi <- which(!is.na(y)); mi <- which(is.na(y)); pred <- y
      if (length(mi)) pred[mi] <- b + as.numeric(R[mi, oi, drop = FALSE] %*% solve(R[oi, oi, drop = FALSE] + diag(1e-6, length(oi)), y[oi] - b))
      out[[v]] <- unname(pred)
    }; out
  })
}
fit_bace <- function(d, tree) {
  if (!requireNamespace("BACE", quietly = TRUE)) return(list(value = NULL, error = "BACE not installed", wall_s = 0))
  fit_safe({
    tb <- tree; tb$edge.length[tb$edge.length == 0] <- 1e-8
    x <- d; x$Species <- rownames(x); f <- lapply(names(d), function(v) paste0(v, " ~ ", paste(setdiff(names(d), v), collapse = " + ")))
    BACE::bace(fixformula = f, ran_phylo_form = "~ 1 |Species", phylo = tb, data = x, nitt = 2000L, burnin = 500L, thin = 5L, runs = 2L, n_final = 2L, verbose = FALSE, skip_conv = TRUE)
  })
}
bace_completed <- function(fit, d, mask) {
  draws <- fit$imputed_datasets %||% fit$imputed_data %||% if (!is.null(fit$data)) list(fit$data) else NULL
  if (is.null(draws) || !length(draws)) stop("BACE output shape not recognised")
  out <- d
  for (v in names(d)) for (i in which(mask[, v])) {
    vals <- vapply(draws, function(z) as.character(z[[v]][i]), character(1)); vals <- vals[!is.na(vals)]
    if (!length(vals)) next
    if (is.factor(out[[v]])) out[[v]][i] <- factor(names(sort(table(vals), decreasing = TRUE))[1], levels = levels(out[[v]]), ordered = is.ordered(out[[v]])) else out[[v]][i] <- stats::median(as.numeric(vals))
  }; out
}
as_completed <- function(name, fit, d, mask) {
  if (is.null(fit$value)) return(list(completed = NULL, error = fit$error, wall_s = fit$wall_s))
  z <- tryCatch(if (name %in% c("pigauto_default", "pigauto_exact")) fit$value$completed else if (name == "bace") bace_completed(fit$value, d, mask) else fit$value, error = function(e) e)
  if (inherits(z, "error")) return(list(completed = NULL, error = conditionMessage(z), wall_s = fit$wall_s))
  list(completed = z, error = NULL, wall_s = fit$wall_s)
}
score <- function(pred, truth, hidden, is_cont, train) {
  i <- which(hidden & !is.na(pred) & !is.na(truth)); n <- length(i)
  if (!n) return(c(n_cells = 0, rmse_raw = NA, rmse_norm = NA, pearson_r = NA, accuracy = NA, brier = NA))
  if (is_cont) {
    err <- as.numeric(pred[i]) - truth[i]; s <- train[["sd"]]
    c(n_cells = n, rmse_raw = sqrt(mean(err^2)), rmse_norm = if (is.finite(s) && s > 0) sqrt(mean((err / s)^2)) else NA, pearson_r = if (n > 2 && stats::sd(pred[i]) > 0) stats::cor(pred[i], truth[i]) else NA, accuracy = NA, brier = NA)
  } else {
    p <- as.character(pred[i]); y <- as.character(truth[i]); c(n_cells = n, rmse_raw = NA, rmse_norm = NA, pearson_r = NA, accuracy = mean(p == y), brier = NA)
  }
}

dat <- load_data(); truth <- dat$df; rows <- list()
for (seed in MASK_SEEDS) {
  mask <- make_mask(truth, seed); masked <- truth
  for (v in names(masked)) masked[[v]][mask[, v]] <- NA
  fits <- list(
    pigauto_default = fit_pigauto(masked, dat$tree, seed, "per_column"),
    pigauto_exact = fit_pigauto(masked, dat$tree, seed, "exact"),
    rphylopars = fit_rphylopars(masked[, CONT, drop = FALSE], dat$tree),
    phylolm_lambda = fit_phylolm(masked[, CONT, drop = FALSE], dat$tree),
    missforest = fit_missforest(masked), bace = fit_bace(masked, dat$tree))
  methods <- lapply(names(fits), function(n) as_completed(n, fits[[n]], masked, mask)); names(methods) <- names(fits)
  for (mn in names(methods)) for (v in TRAITS) {
    applicable <- !(mn %in% c("rphylopars", "phylolm_lambda") && v %in% DISC)
    mo <- methods[[mn]]; st <- if (!applicable) "not_applicable" else if (is.null(mo$completed) || is.null(mo$completed[[v]])) "failed" else "ok"
    tr <- truth[[v]][!mask[, v]]; sc <- if (st == "ok") score(mo$completed[[v]], truth[[v]], mask[, v], v %in% CONT, c(mean = mean(tr, na.rm = TRUE), sd = stats::sd(tr, na.rm = TRUE))) else score(rep(NA, nrow(truth)), truth[[v]], mask[, v], v %in% CONT, c(sd = NA))
    rows[[length(rows) + 1L]] <- cbind(
      data.frame(regime = "B-mixed", mask_seed = seed, method = mn, trait = v,
                 status = st, error = if (st == "failed") mo$error else NA_character_,
                 wall_s = mo$wall_s, stringsAsFactors = FALSE),
      as.data.frame(as.list(sc), stringsAsFactors = FALSE)
    )
  }
  last <- rows[(length(rows) - length(methods) * length(TRAITS) + 1L):length(rows)]
  saveRDS(list(protocol = "stage_b_mixed_core_v1", mask_seed = seed, miss_frac = MISS_FRAC, source_sha = Sys.getenv("PIGAUTO_SOURCE_SHA", unset = NA_character_), mask = mask, truth = truth, masked = masked, methods = methods, rows = last), file.path(receipt_dir, sprintf("stage_b_mixed_core_mask-%d.rds", seed)))
}
results <- do.call(rbind, rows); rownames(results) <- NULL
summary <- do.call(rbind, lapply(split(results, interaction(results$method, results$trait, drop = TRUE)), function(x) data.frame(method = x$method[1], trait = x$trait[1], n_ok = sum(x$status == "ok"), n_failed = sum(x$status == "failed"), n_not_applicable = sum(x$status == "not_applicable"), rmse_raw_mean = mean(x$rmse_raw, na.rm = TRUE), rmse_norm_mean = mean(x$rmse_norm, na.rm = TRUE), accuracy_mean = mean(x$accuracy, na.rm = TRUE), stringsAsFactors = FALSE)))
saveRDS(list(results = results, summary = summary, protocol = "stage_b_mixed_core_v1", mask_seeds = MASK_SEEDS), file.path(out_dir, "stage_b_mixed_core.rds"))
writeLines(c("# Stage B mixed core: retained comparator output", "", "This is a pre-registered, five-mask runner output. It reports per-method/trait status; continuous RMSE and discrete accuracy are never pooled.", "", "## Summary", "```", capture.output(print(summary, row.names = FALSE)), "```", "", "## Non-claims", "No parity, default-change, or general mixed-type superiority claim follows from this descriptive protocol."), file.path(out_dir, "stage_b_mixed_core.md"))
