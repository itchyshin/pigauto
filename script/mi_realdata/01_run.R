#!/usr/bin/env Rscript
# Run multi_impute(draws_method = "posterior") on one real-data mask cell
# and compare against the Mondrian split/mondrian conformal receipts for
# the SAME masked cells (read only -- conformal is never rerun here).
#
# Usage: Rscript 01_run.R dataset arm seed outdir
#
# Requires script/mi_realdata/inputs/<dataset>-<arm>-m<seed>/mask_receipt.rds
# (run 00_fetch_masks.R first). Conformal results come from
# inputs/<cell>/conformal_metrics.rds when present
# (`00_fetch_masks.R --with-conformal`), otherwise from `git show` on
# arc/mondrian-realdata.
#
# Environment (all optional):
#   MI_REALDATA_OFFLINE=1  never call git: conformal_metrics.rds must exist
#                          (Totoro / DRAC compute nodes).
#   PIGAUTO_PKG_PATH       package root to devtools::load_all(); otherwise
#                          library(pigauto).
#   MI_POST_SHA            code commit recorded in the receipt; otherwise
#                          `git rev-parse HEAD` when available. Required
#                          with MI_REALDATA_OFFLINE=1 (G8 rejects a receipt
#                          without a SHA, so the run refuses up front).
#   MI_POST_NITER, MI_POST_BURNIN
#                          override posterior_control n_iter / burnin for
#                          quick local smokes only. Unset means pigauto's
#                          defaults. An override is recorded in the receipt
#                          and makes the cell fail G8 (03_acceptance.R).
# There is no n_cores control: the chains of one cell run sequentially in
# one R process, so parallelise across cells (12_totoro_run.sh, or one
# 1-CPU array task per cell in 10_fir.sbatch).
#
# Writes <outdir>/<dataset>-<arm>-m<seed>/mi_posterior.rds, a list with
# `status` "ok" or "error" (see below), that 02_summarise.R / 03_acceptance.R
# consume. The first thing it does is overwrite that file with a status
# "running" placeholder, so a crash or kill (e.g. out of memory) can never
# leave an older "ok" receipt from a previous run or SHA in place.

# Rscript passes a script path containing spaces as `~+~` in --file=
# (e.g. ".../Github~+~Local/..."); undo that before normalizePath().
here <- function() {
  a <- commandArgs(FALSE)
  f <- gsub("~+~", " ", sub("^--file=", "", a[grepl("^--file=", a)]), fixed = TRUE)
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
saveRDS(list(status = "running", stage = "started", dataset = dataset, arm = arm, seed = seed,
             name = name, started_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
             error = "01_run.R started but did not write a final receipt (crashed, killed or still running)"),
        receipt_file)

offline <- identical(Sys.getenv("MI_REALDATA_OFFLINE"), "1")

# Code provenance: MI_POST_SHA (set by 12_totoro_run.sh from code_dir/SHA),
# else git HEAD of this checkout, flagged when R/ or DESCRIPTION has
# uncommitted changes (load_all() would then run code that no SHA names).
code_sha <- Sys.getenv("MI_POST_SHA")
code_sha_source <- "MI_POST_SHA"
if (!nzchar(code_sha)) {
  code_sha <- if (offline) NA_character_ else git_rev_parse("HEAD", dir = here_dir)
  code_sha_source <- if (is.na(code_sha)) "unknown" else "git HEAD"
  if (!is.na(code_sha)) {
    top <- git_rev_parse("--show-toplevel", dir = here_dir)
    dirty <- if (is.na(top)) character(0) else
      tryCatch(suppressWarnings(system2("git", c("-C", top, "status", "--porcelain", "--", "R", "DESCRIPTION"),
                                        stdout = TRUE, stderr = FALSE)), error = function(e) character(0))
    if (length(dirty)) code_sha_source <- "git HEAD, R/ or DESCRIPTION has uncommitted changes"
  }
}

write_receipt <- function(x) {
  x$receipt_schema <- receipt_schema_current
  x$code_sha <- code_sha
  x$code_sha_source <- code_sha_source
  saveRDS(x, receipt_file)
  invisible(x)
}

if (offline && is.na(code_sha)) {
  msg <- paste0("MI_REALDATA_OFFLINE=1 needs MI_POST_SHA (the commit of this code): G8 rejects a ",
                "receipt without a code SHA, so the run would be wasted.")
  write_receipt(list(status = "error", stage = "provenance", dataset = dataset, arm = arm,
                      seed = seed, name = name, error = msg))
  stop("[", name, "] ", msg, call. = FALSE)
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

# ---- 2b. conformal results: READ ONLY, never rerun ------------------------
# Read before the sampler so a missing input fails in seconds, not after
# the run. A receipt whose own status is not "ok" is recorded as it is
# (G8 then fails and names it); an unreadable source stops here.
cm_file <- conformal_metrics_path(file.path(here_dir, "inputs"), name)
if (file.exists(cm_file)) {
  cm <- readRDS(cm_file)
  conformal_source <- list(kind = "inputs file", file = cm_file, source = cm$source)
} else if (offline) {
  msg <- paste0("MI_REALDATA_OFFLINE=1 and ", cm_file, " does not exist; run ",
                "`Rscript 00_fetch_masks.R --all --with-conformal` on a machine with the git history ",
                "and copy inputs/ over.")
  write_receipt(list(status = "error", stage = "conformal_inputs", dataset = dataset, arm = arm,
                      seed = seed, name = name, error = msg))
  stop("[", name, "] ", msg, call. = FALSE)
} else {
  cm <- tryCatch(fetch_conformal_metrics(dataset, arm, seed), error = function(e) e)
  if (inherits(cm, "condition")) {
    write_receipt(list(status = "error", stage = "conformal_inputs", dataset = dataset, arm = arm,
                        seed = seed, name = name, error = conditionMessage(cm)))
    stop("[", name, "] could not read conformal results: ", conditionMessage(cm), call. = FALSE)
  }
  conformal_source <- list(kind = "git show", source = cm$source)
}
keep_kept <- function(co) {
  if (identical(co$status, "ok")) co$metrics <- co$metrics[co$metrics$trait %in% kept, , drop = FALSE]
  co
}
split_conformal <- keep_kept(cm$split)
mondrian_conformal <- keep_kept(cm$mondrian)
cat(sprintf("[%s] conformal (%s): split=%s mondrian=%s\n", name, conformal_source$kind,
            split_conformal$status, mondrian_conformal$status))

# ---- 3. run multi_impute(draws_method = "posterior") ----------------------
pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH")
if (nzchar(pkg_path)) {
  suppressPackageStartupMessages(devtools::load_all(pkg_path, quiet = TRUE))
} else {
  suppressPackageStartupMessages(library(pigauto))
}

# Sampler settings: pigauto defaults (4 chains, run one after another),
# keep_draws = 1000 for the per-cell intervals. MI_POST_NITER /
# MI_POST_BURNIN are for quick local smokes only and are recorded.
env_int <- function(var) {
  v <- Sys.getenv(var)
  if (!nzchar(v)) return(NULL)
  out <- suppressWarnings(as.integer(v))
  if (is.na(out)) stop(var, " must be an integer, got '", v, "'", call. = FALSE)
  out
}
overrides <- Filter(Negate(is.null), list(n_iter = env_int("MI_POST_NITER"),
                                          burnin = env_int("MI_POST_BURNIN")))
posterior_control <- c(list(keep_draws = 1000L), overrides)
if (length(overrides)) {
  cat(sprintf("[%s] SAMPLER OVERRIDE (smoke only, fails G8): %s\n", name,
              paste(names(overrides), unlist(overrides), sep = " = ", collapse = ", ")))
}

t0 <- proc.time()[["elapsed"]]
mi <- tryCatch(
  pigauto::multi_impute(
    traits_sub, tree, m = 20L, draws_method = "posterior",
    posterior_control = posterior_control,
    verbose = FALSE, seed = seed
  ),
  error = function(e) e
)
wall_time_s <- proc.time()[["elapsed"]] - t0

if (inherits(mi, "condition")) {
  msg <- conditionMessage(mi)
  write_receipt(list(status = "error", stage = "multi_impute", dataset = dataset, arm = arm,
                      seed = seed, name = name, kept_traits = kept, dropped_traits = dropped,
                      wall_time_s = wall_time_s, sampler = list(overrides = overrides), error = msg))
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
# Widths are summarised like for like with the conformal receipts: mean
# AND median of (upper - lower) over the masked cells, on the original
# trait scale (cell_interval is decoded, as are conformal lo/hi).
# n_matched counts only masked cells with a FINITE lower and upper bound,
# and nothing below drops NA: one missing interval or truth makes
# n_matched < n_masked or coverage NA, and G8 names the trait (schema 3,
# 2026-09-24 repair).
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
  n_matched <- sum(is.finite(sub$lower) & is.finite(sub$upper))
  if (nrow(sub) == 0L) {
    return(data.frame(trait = nm, n_masked = length(midx), n_matched = 0L,
                       coverage = NA_real_, mean_width = NA_real_, median_width = NA_real_))
  }
  tv <- truth[[nm]][sub$.row_idx]
  cov <- mean(tv >= sub$lower & tv <= sub$upper)
  w <- sub$upper - sub$lower
  data.frame(trait = nm, n_masked = length(midx), n_matched = n_matched,
             coverage = cov, mean_width = mean(w), median_width = stats::median(w))
}))
cat("[", name, "] model-based per-cell coverage:\n"); print(model_coverage)

# ---- 5. conformal comparison: read in step 2b (split_conformal,
# mondrian_conformal), restricted to the kept traits.

# ---- 6. downstream slope check for pre-registered pairs on this dataset ---
# Paired design: the reference and the MI legs use the SAME species (rows
# where both traits were observed before the Mondrian mask) and the SAME
# analysis model (phylolm, model = "lambda", on the pruned tree), so the only
# difference between them is that masked cells are imputed in the MI leg.
# phylolm is O(n) per likelihood evaluation; dense gls(corPagel) is not
# feasible at these n. MI pooling is Rubin's rules with the Barnard-Rubin df
# (complete-data df = n - 2), written out here because pool_mi() has no
# phylolm adapter. The analysis scale (y_log, x_log) is pre-registered per
# pair in pairs.R and used for both legs; pigauto's trait_map$log_transform
# is not consulted (it would re-log the already-logged PanTHERIA traits).
pair_frame <- function(src, sp, y_nm, x_nm, y_log, x_log) {
  d <- data.frame(row.names = sp)
  d$y <- if (y_log) log(src[sp, y_nm]) else src[sp, y_nm]
  d$x <- if (x_log) log(src[sp, x_nm]) else src[sp, x_nm]
  d
}
fit_phylolm <- function(d, tr) {
  # A non-finite value (e.g. log of a non-positive imputed value) would make
  # phylolm drop or reject rows and break the pairing: count it as a failed
  # fit instead.
  if (!all(is.finite(d$y)) || !all(is.finite(d$x))) return(NULL)
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
  y_log <- pair$y_log; x_log <- pair$x_log
  stopifnot(is.logical(y_log), length(y_log) == 1L, !is.na(y_log),
            is.logical(x_log), length(x_log) == 1L, !is.na(x_log))
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
  frames <- lapply(mi$datasets, function(dat) pair_frame(dat, ps$sp, ps$y_nm, ps$x_nm, ps$y_log, ps$x_log))
  n_nonfinite <- sum(vapply(frames, function(d) !all(is.finite(d$y)) || !all(is.finite(d$x)), logical(1)))
  fits <- lapply(frames, fit_phylolm, tr = ps$tr)
  ok <- !vapply(fits, is.null, logical(1))
  if (sum(ok) < 2L) {
    return(list(status = "error", n_nonfinite = n_nonfinite,
                error = sprintf("only %d/%d phylolm fits succeeded (%d completions non-finite on the analysis scale)",
                                sum(ok), length(fits), n_nonfinite)))
  }
  est <- vapply(fits[ok], `[[`, numeric(1), "slope")
  se  <- vapply(fits[ok], `[[`, numeric(1), "se")
  m <- length(est); qbar <- mean(est); ubar <- mean(se^2); b <- stats::var(est)
  tvar <- ubar + (1 + 1 / m) * b
  lam <- (1 + 1 / m) * b / tvar
  nu_com <- length(ps$sp) - 2
  nu_old <- (m - 1) / max(lam, 1e-12)^2
  nu_obs <- (nu_com + 1) / (nu_com + 3) * nu_com * (1 - lam)
  list(status = "ok", m_used = m, m_total = length(fits), n_nonfinite = n_nonfinite, n = length(ps$sp),
       slope = qbar, se = sqrt(tvar), df = 1 / (1 / nu_old + 1 / nu_obs), fmi = lam,
       y_log = ps$y_log, x_log = ps$x_log)
}

pairs_here <- pairs_for_dataset(mi_realdata_pairs, dataset)
pair_results <- lapply(pairs_here, function(p) {
  ps <- pair_setup(p)
  ref <- ref_slope_one(ps)
  mip <- mi_slope_one(ps)
  c(list(dataset = p$dataset, response = p$response, predictor = p$predictor,
         rationale = p$rationale, reference = ref, mi = mip),
    pair_contrast(ref, mip))  # diff, rel_diff, diff_ref_se, se_ratio
})
cat(sprintf("[%s] pair slopes: %d pre-registered pairs for dataset '%s'\n", name, length(pair_results), dataset))
for (pr in pair_results) {
  cat(sprintf("  %s ~ %s: ref=%s mi=%s rel_diff=%s\n", pr$response, pr$predictor,
              if (identical(pr$reference$status, "ok")) sprintf("%.4g", pr$reference$slope) else pr$reference$status,
              if (identical(pr$mi$status, "ok")) sprintf("%.4g", pr$mi$slope) else pr$mi$status,
              format(signif(pr$rel_diff, 3))))
}

# ---- 7. write receipt -------------------------------------------------------
ctl <- mi$posterior$control
receipt <- list(
  status = "ok", dataset = dataset, arm = arm, seed = seed, name = name,
  n_species = nrow(masked), kept_traits = kept, dropped_traits = dropped,
  wall_time_s = wall_time_s, model_coverage = model_coverage,
  split_conformal = split_conformal, mondrian_conformal = mondrian_conformal,
  conformal_source = conformal_source,
  sampler = list(overrides = overrides,
                 control = ctl[intersect(c("n_chains", "n_iter", "burnin", "thin", "keep_draws",
                                           "param_uncertainty"), names(ctl))],
                 sweeps = mi$posterior$sweeps, wall_s = mi$posterior$wall_s),
  pkg_source = if (nzchar(pkg_path)) paste0("devtools::load_all(", pkg_path, ")") else "library(pigauto)",
  diagnostics = mi$posterior$diagnostics, convergence = convergence, pairs = pair_results
)
write_receipt(receipt)
cat("[", name, "] wrote ", receipt_file, "\n", sep = "")
cat("MI_REALDATA_CELL_OK\n")
