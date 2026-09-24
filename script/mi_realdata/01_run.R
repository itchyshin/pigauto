#!/usr/bin/env Rscript
# Run multi_impute(draws_method = "posterior") on one real-data mask cell
# and compare against the Mondrian split/mondrian conformal receipts for
# the SAME masked cells (read only -- conformal is never rerun here).
#
# Usage: Rscript 01_run.R dataset arm seed outdir
#
# Requires script/mi_realdata/inputs/<dataset>-<arm>-m<seed>/mask_receipt.rds
# to already exist (run 00_fetch_masks.R first).
#
# Writes <outdir>/<dataset>-<arm>-m<seed>/mi_posterior.rds, a list with
# `status` "ok" or "error" (see below), that 02_summarise.R / 03_acceptance.R
# consume.

here <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
here_dir <- here()
source(file.path(here_dir, "lib.R"))
source(file.path(here_dir, "pairs.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("expected: dataset arm seed outdir", call. = FALSE)
dataset <- args[[1L]]; arm <- args[[2L]]; seed <- as.integer(args[[3L]]); outdir <- args[[4L]]
name <- cell_name(dataset, arm, seed)

receipt_dir <- file.path(outdir, name)
dir.create(receipt_dir, recursive = TRUE, showWarnings = FALSE)
receipt_file <- file.path(receipt_dir, "mi_posterior.rds")

write_receipt <- function(x) {
  saveRDS(x, receipt_file)
  invisible(x)
}

# ---- 1. load the mask receipt (truth / masked / mask / tree) --------------
mask_file <- file.path(here_dir, "inputs", name, "mask_receipt.rds")
if (!file.exists(mask_file)) {
  stop("missing ", mask_file, " -- run `Rscript 00_fetch_masks.R ", dataset, " ", arm, " ", seed,
       "` (or --all) first.", call. = FALSE)
}
mr <- readRDS(mask_file)
truth <- mr$truth; masked <- mr$masked; mask <- mr$mask; tree <- mr$tree
if (!identical(rownames(masked), tree$tip.label)) {
  stop("masked data rownames do not match tree tip order in ", mask_file, call. = FALSE)
}

# ---- 2. restrict to continuous-family traits -------------------------------
cls <- classify_traits(masked)
kept <- cls$trait[cls$keep]
dropped <- cls[!cls$keep, , drop = FALSE]
cat(sprintf("[%s] %d traits total, %d kept (continuous), %d dropped\n",
            name, nrow(cls), length(kept), nrow(dropped)))
if (nrow(dropped)) {
  for (i in seq_len(nrow(dropped))) {
    cat(sprintf("  dropped %-25s %-12s %s\n", dropped$trait[i], dropped$class[i], dropped$reason[i]))
  }
}
if (length(kept) == 0L) {
  write_receipt(list(status = "error", stage = "trait_selection", dataset = dataset, arm = arm,
                      seed = seed, name = name, dropped_traits = dropped,
                      error = "no continuous-family traits remain after filtering"))
  stop("[", name, "] no continuous-family traits remain; posterior method needs >= 1.", call. = FALSE)
}
traits_sub <- masked[, kept, drop = FALSE]

# ---- 3. run multi_impute(draws_method = "posterior") ----------------------
suppressPackageStartupMessages(library(pigauto))

t0 <- proc.time()[["elapsed"]]
mi <- tryCatch(
  pigauto::multi_impute(
    traits_sub, tree, m = 20L, draws_method = "posterior",
    posterior_control = list(keep_draws = 1000L),
    verbose = FALSE, seed = seed
  ),
  error = function(e) e
)
wall_time_s <- proc.time()[["elapsed"]] - t0

if (inherits(mi, "condition")) {
  msg <- conditionMessage(mi)
  write_receipt(list(status = "error", stage = "multi_impute", dataset = dataset, arm = arm,
                      seed = seed, name = name, kept_traits = kept, dropped_traits = dropped,
                      wall_time_s = wall_time_s, error = msg))
  stop(
    "[", name, "] multi_impute(draws_method = 'posterior') is not available or failed:\n  ", msg, "\n",
    "This harness is coded exactly against the frozen API in ",
    "docs/dev-log/mi-posterior/design.md section 4 (mi$datasets, mi$posterior$cell_interval, ",
    "mi$posterior$diagnostics, mi$posterior$params, mi$draws_method == \"posterior\"). ",
    "If R/mi_posterior.R implementing draws_method = \"posterior\" has not landed on this ",
    "branch (arc/mi-posterior) yet, this failure is expected -- re-run this script once it has.",
    call. = FALSE
  )
}

stopifnot(
  identical(mi$draws_method, "posterior"),
  is.list(mi$posterior),
  is.data.frame(mi$posterior$cell_interval),
  is.list(mi$datasets), length(mi$datasets) == 20L
)
cat(sprintf("[%s] multi_impute(posterior) OK in %.1fs, m = %d\n", name, wall_time_s, length(mi$datasets)))

# ---- 3b. convergence: max split R-hat / min bulk ESS over Sigma_P, ---------
# Sigma_E, lambda (design-review addition, 2026-09-24), plus the
# diagnostics object's own attr(,"converged"). Frozen API unchanged --
# this only reads mi$posterior$diagnostics as documented.
convergence <- summarize_convergence(mi$posterior$diagnostics)
cat(sprintf("[%s] convergence: max_rhat=%.3f min_ess=%.0f converged=%s\n",
            name, convergence$max_rhat, convergence$min_ess, convergence$converged))

# ---- 4. per-cell model-based coverage / width, from mi$posterior ----------
# design.md section 4 documents mi$posterior$cell_interval as
# data.frame(row, trait, lower, upper, median) but doesn't pin down
# whether `row` is a tip index into tree$tip.label / rownames(masked)
# (pigauto's usual internal convention, e.g. splits$val_idx) or the
# species name itself. Handle both without guessing wrong silently.
ci <- mi$posterior$cell_interval
resolve_row_idx <- function(row_col, species_names) {
  if (is.numeric(row_col) || is.integer(row_col)) return(as.integer(row_col))
  match(as.character(row_col), species_names)
}
ci$.row_idx <- resolve_row_idx(ci$row, rownames(masked))

model_coverage <- do.call(rbind, lapply(kept, function(nm) {
  midx <- which(mask[, nm])
  sub <- ci[ci$trait == nm & ci$.row_idx %in% midx, , drop = FALSE]
  n_matched <- nrow(sub)
  if (n_matched == 0L) {
    return(data.frame(trait = nm, n_masked = length(midx), n_matched = 0L,
                       coverage = NA_real_, median_width = NA_real_))
  }
  tv <- truth[[nm]][sub$.row_idx]
  cov <- mean(tv >= sub$lower & tv <= sub$upper, na.rm = TRUE)
  wid <- stats::median(sub$upper - sub$lower, na.rm = TRUE)
  data.frame(trait = nm, n_masked = length(midx), n_matched = n_matched,
             coverage = cov, median_width = wid)
}))
cat("[", name, "] model-based per-cell coverage:\n"); print(model_coverage)

# ---- 5. conformal comparison -- READ ONLY, never rerun ---------------------
read_conformal <- function(leaf) {
  res <- tryCatch(read_mondrian_rds(dataset, arm, seed, leaf), error = function(e) e)
  if (inherits(res, "condition")) {
    return(list(status = "error", error = conditionMessage(res), metrics = NULL))
  }
  if (!identical(res$status, "ok") || is.null(res$metrics)) {
    return(list(status = res$status %||% "error", error = res$error %||% "no metrics in receipt",
                metrics = NULL))
  }
  list(status = "ok", error = NULL, metrics = res$metrics[res$metrics$trait %in% kept, , drop = FALSE])
}
split_conformal <- read_conformal("split.rds")
mondrian_conformal <- read_conformal("mondrian.rds")

# ---- 6. downstream slope check for pre-registered pairs on this dataset ---
trait_map <- mi$data$trait_map

# Paired design: the reference and the MI legs use the SAME species (rows
# where both traits were observed before the Mondrian mask) and the SAME
# analysis model (phylolm, model = "lambda", on the pruned tree), so the only
# difference between them is that masked cells are imputed in the MI leg.
# phylolm is O(n) per likelihood evaluation; dense gls(corPagel) is not
# feasible at these n. MI pooling is Rubin's rules with the Barnard-Rubin df
# (complete-data df = n - 2), written out here because pool_mi() has no
# phylolm adapter.
pair_frame <- function(src, sp, y_nm, x_nm, y_log, x_log) {
  d <- data.frame(row.names = sp)
  d$y <- if (y_log) log(src[sp, y_nm]) else src[sp, y_nm]
  d$x <- if (x_log) log(src[sp, x_nm]) else src[sp, x_nm]
  d
}
fit_phylolm <- function(d, tr) {
  fit <- tryCatch(suppressWarnings(phylolm::phylolm(y ~ x, data = d, phy = tr, model = "lambda")),
                  error = function(e) e)
  if (inherits(fit, "condition")) return(NULL)
  co <- summary(fit)$coefficients
  c(slope = unname(co["x", "Estimate"]), se = unname(co["x", "StdErr"]), lambda = unname(fit$optpar))
}

pair_setup <- function(pair) {
  y_nm <- pair$response; x_nm <- pair$predictor
  if (!(y_nm %in% kept) || !(x_nm %in% kept)) {
    return(list(status = "skipped", error = "response/predictor not in continuous-family kept traits"))
  }
  y_log <- isTRUE(trait_map[[y_nm]]$log_transform)
  x_log <- isTRUE(trait_map[[x_nm]]$log_transform)
  sp <- rownames(truth)[!is.na(truth[[y_nm]]) & !is.na(truth[[x_nm]])]
  d <- pair_frame(truth, sp, y_nm, x_nm, y_log, x_log)
  sp <- rownames(d)[is.finite(d$y) & is.finite(d$x)]
  if (length(sp) < 10L) return(list(status = "error", error = "fewer than 10 originally-observed finite rows for this pair"))
  tr <- tryCatch(ape::keep.tip(tree, sp), error = function(e) e)
  if (inherits(tr, "condition")) return(list(status = "error", error = conditionMessage(tr)))
  list(status = "ok", y_nm = y_nm, x_nm = x_nm, y_log = y_log, x_log = x_log, sp = tr$tip.label, tr = tr)
}

ref_slope_one <- function(ps) {
  if (!identical(ps$status, "ok")) return(ps[c("status", "error")])
  f <- fit_phylolm(pair_frame(truth, ps$sp, ps$y_nm, ps$x_nm, ps$y_log, ps$x_log), ps$tr)
  if (is.null(f)) return(list(status = "error", error = "phylolm reference fit failed"))
  list(status = "ok", n = length(ps$sp), slope = f[["slope"]], se = f[["se"]], lambda = f[["lambda"]],
       y_log = ps$y_log, x_log = ps$x_log)
}

mi_slope_one <- function(ps) {
  if (!identical(ps$status, "ok")) return(ps[c("status", "error")])
  fits <- lapply(mi$datasets, function(dat)
    fit_phylolm(pair_frame(dat, ps$sp, ps$y_nm, ps$x_nm, ps$y_log, ps$x_log), ps$tr))
  ok <- !vapply(fits, is.null, logical(1))
  if (sum(ok) < 2L) {
    return(list(status = "error", error = sprintf("only %d/%d phylolm fits succeeded", sum(ok), length(fits))))
  }
  est <- vapply(fits[ok], `[[`, numeric(1), "slope")
  se  <- vapply(fits[ok], `[[`, numeric(1), "se")
  m <- length(est); qbar <- mean(est); ubar <- mean(se^2); b <- stats::var(est)
  tvar <- ubar + (1 + 1 / m) * b
  lam <- (1 + 1 / m) * b / tvar
  nu_com <- length(ps$sp) - 2
  nu_old <- (m - 1) / max(lam, 1e-12)^2
  nu_obs <- (nu_com + 1) / (nu_com + 3) * nu_com * (1 - lam)
  list(status = "ok", m_used = m, m_total = length(fits), n = length(ps$sp),
       slope = qbar, se = sqrt(tvar), df = 1 / (1 / nu_old + 1 / nu_obs), fmi = lam,
       y_log = ps$y_log, x_log = ps$x_log)
}

pairs_here <- pairs_for_dataset(mi_realdata_pairs, dataset)
pair_results <- lapply(pairs_here, function(p) {
  ps <- pair_setup(p)
  ref <- ref_slope_one(ps)
  mip <- mi_slope_one(ps)
  rel_diff <- NA_real_; se_ratio <- NA_real_
  if (identical(ref$status, "ok") && identical(mip$status, "ok") && is.finite(ref$slope) && ref$slope != 0) {
    rel_diff <- (mip$slope - ref$slope) / ref$slope
    se_ratio <- mip$se / ref$se
  }
  list(dataset = p$dataset, response = p$response, predictor = p$predictor,
       rationale = p$rationale, reference = ref, mi = mip,
       rel_diff = rel_diff, se_ratio = se_ratio)
})
cat(sprintf("[%s] pair slopes: %d pre-registered pairs for dataset '%s'\n", name, length(pair_results), dataset))
for (pr in pair_results) {
  cat(sprintf("  %s ~ %s: ref=%s mi=%s\n", pr$response, pr$predictor,
              if (identical(pr$reference$status, "ok")) sprintf("%.4g", pr$reference$slope) else pr$reference$status,
              if (identical(pr$mi$status, "ok")) sprintf("%.4g", pr$mi$slope) else pr$mi$status))
}

# ---- 7. write receipt -------------------------------------------------------
receipt <- list(
  status = "ok", dataset = dataset, arm = arm, seed = seed, name = name,
  n_species = nrow(masked), kept_traits = kept, dropped_traits = dropped,
  wall_time_s = wall_time_s, model_coverage = model_coverage,
  split_conformal = split_conformal, mondrian_conformal = mondrian_conformal,
  diagnostics = mi$posterior$diagnostics, convergence = convergence, pairs = pair_results
)
write_receipt(receipt)
cat("[", name, "] wrote ", receipt_file, "\n", sep = "")
cat("MI_REALDATA_CELL_OK\n")
