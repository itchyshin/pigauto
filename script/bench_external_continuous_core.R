#!/usr/bin/env Rscript
#
# script/bench_external_continuous_core.R
#
# pigauto vs external rival R packages, fit on the SAME masked real
# data (bundled AVONET300). This is the first committed bench that
# compares pigauto against genuinely external methods rather than
# only against BACE (script/bench_avonet_bace.R) or internal-only
# baselines (column mean, phylogenetic BLUP).
#
# Methods compared, all on the SAME 30%-MCAR mask per rep:
#   1. pigauto::impute()            -- production defaults, all 7 traits
#   2. pigauto::impute(exact)       -- opt-in exact conditional route
#   3. Rphylopars::phylopars()      -- standalone, continuous traits only
#   4. phylolm(model = "lambda")    -- per-trait phylogenetic BLUP,
#                                       continuous traits only
#   4. column mean                  -- floor
#
# IMPORTANT CAVEAT (see also the .md report): Rphylopars is ALSO the
# solver pigauto uses internally for its joint-baseline path
# (R/joint_mvn_baseline.R, R/joint_threshold_baseline.R via
# R/joint_mvn_solver.R, both model = "BM"). So the pigauto-vs-phylopars
# comparison here measures pigauto's GNN + gating + conformal layer
# against the raw phylopars solver it can call internally -- it is
# NOT a comparison against an unrelated method.
#
# Scoring: pigauto is scored ONLY on the 4 continuous traits (Mass,
# Beak.Length_Culmen, Tarsus.Length, Wing.Length) for comparability
# with the continuous-only rivals, even though it is fit on all 7
# AVONET traits (its unified mixed-type imputation is the point of
# the package, not something the rivals can do at all).
#
# Metric: z-scored RMSE (scale = training-portion mean/sd per trait
# per rep, i.e. computed from the non-held-out observed cells only --
# no leakage from the held-out truth) + Pearson r on held-out cells.
# This is the continuous-core component of the Stage-B protocol.  It uses the
# locked five masks (seeds 20260901--20260905) and writes a complete receipt
# after every mask.  It does not replace the separately required mixed-type
# regime, and it must not be used to support a general parity or default claim.
#
# A rival that errors on a trait is recorded as NA with the error
# message verbatim in the `error` column -- never silently dropped.
#
# Output:
#   <PIGAUTO_OUT_DIR>/stage_b_continuous_core.rds
#   <PIGAUTO_OUT_DIR>/stage_b_continuous_core.md
#   <PIGAUTO_OUT_DIR>/receipts/stage_b_continuous_core_mask-<seed>.rds
#
# CPU discipline: torch pinned to 1 thread; run under `nice -n 15`.

options(warn = 1, stringsAsFactors = FALSE)
`%||%` <- function(x, y) if (is.null(x)) y else x
suppressPackageStartupMessages({
  library(ape)
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH", unset = "")
  if (nzchar(pkg_path) && dir.exists(pkg_path) &&
      file.exists(file.path(pkg_path, "DESCRIPTION"))) {
    devtools::load_all(pkg_path, quiet = TRUE)
    message("Loaded pigauto from source: ", pkg_path)
  } else {
    library(pigauto)
    message("Loaded pigauto from installed library")
  }
})

if (requireNamespace("torch", quietly = TRUE)) {
  try(torch::torch_set_num_threads(1L), silent = TRUE)
  try(torch::torch_set_num_interop_threads(1L), silent = TRUE)
}

MISS_FRAC <- 0.30
MASK_SEEDS <- 20260901L:20260905L
N_REPS <- length(MASK_SEEDS)

script_start <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start),
      ..., "\n", sep = "")
  flush.console()
}

out_dir <- Sys.getenv("PIGAUTO_OUT_DIR", unset = "")
if (!nzchar(out_dir)) stop("Set PIGAUTO_OUT_DIR to a dedicated retained result directory", call. = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
receipt_dir <- file.path(out_dir, "receipts")
dir.create(receipt_dir, recursive = TRUE, showWarnings = FALSE)
out_rds <- file.path(out_dir, "stage_b_continuous_core.rds")
out_md  <- file.path(out_dir, "stage_b_continuous_core.md")

# -------------------------------------------------------------------------
# 1. Dataset loader(s) -- a function per dataset so a second one can be
#    added later without touching the main sweep.
# -------------------------------------------------------------------------

load_avonet300 <- function() {
  e <- new.env(parent = emptyenv())
  utils::data("avonet300", package = "pigauto", envir = e)
  utils::data("tree300",   package = "pigauto", envir = e)
  df   <- e$avonet300
  tree <- e$tree300
  rownames(df) <- df$Species_Key
  df$Species_Key <- NULL
  stopifnot(all(rownames(df) == tree$tip.label))
  list(
    name        = "avonet300",
    df          = df,
    tree        = tree,
    cont_traits = c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length")
  )
}

DATASETS <- list(avonet300 = load_avonet300)

# -------------------------------------------------------------------------
# 2. MCAR masking -- same mask matrix (all traits) feeds pigauto; the
#    rivals' continuous-only subset inherits the identical held-out
#    cells for the 4 scored traits.
# -------------------------------------------------------------------------

make_mask <- function(df, miss_frac, seed) {
  set.seed(seed)
  mask <- matrix(FALSE, nrow(df), ncol(df),
                  dimnames = list(rownames(df), names(df)))
  for (v in names(df)) {
    obs_idx <- which(!is.na(df[[v]]))
    to_hide <- sample(obs_idx, ceiling(length(obs_idx) * miss_frac))
    mask[to_hide, v] <- TRUE
  }
  mask
}

# -------------------------------------------------------------------------
# 3. Methods
# -------------------------------------------------------------------------

fit_pigauto_method <- function(df_miss_all, tree, seed, predict_method) {
  res <- tryCatch(
    pigauto::impute(df_miss_all, tree, seed = seed, verbose = FALSE,
                    predict_method = predict_method),
    error = function(e) e)
  if (inherits(res, "error")) return(list(completed = NULL, error = conditionMessage(res)))
  list(completed = res$completed, error = NULL)
}

fit_rphylopars_method <- function(df_miss4, tree) {
  df_in <- data.frame(species = rownames(df_miss4), df_miss4,
                        stringsAsFactors = FALSE)
  fit <- tryCatch(
    Rphylopars::phylopars(df_in, tree = tree, model = "BM",
                            phylo_correlated = TRUE, pheno_correlated = TRUE,
                            REML = TRUE),
    error = function(e) e)
  if (inherits(fit, "error")) return(list(completed = NULL, error = conditionMessage(fit)))
  out <- tryCatch(
    as.data.frame(fit$anc_recon[rownames(df_miss4), colnames(df_miss4), drop = FALSE]),
    error = function(e) e)
  if (inherits(out, "error")) return(list(completed = NULL, error = conditionMessage(out)))
  list(completed = out, error = NULL)
}

# Per-trait phylolm(model = "lambda") BLUP -- mirrors the arm in
# script/bench_lambda_sweep.R (m_phylolm_blup, ~line 155) with the
# fixed-effect covariates stripped (no predictors here: y ~ 1) since
# this bench has no exogenous covariates, only the trait itself.
fit_phylolm_trait <- function(y_named, tree) {
  obs_names <- names(y_named)[!is.na(y_named)]
  if (length(obs_names) < 5L) {
    return(list(pred = y_named, error = "fewer than 5 observed species"))
  }
  sp_df <- data.frame(species = names(y_named), y = as.numeric(y_named),
                        stringsAsFactors = FALSE)
  rownames(sp_df) <- sp_df$species
  tree_obs <- tryCatch(ape::keep.tip(tree, obs_names), error = function(e) e)
  if (inherits(tree_obs, "error")) return(list(pred = y_named, error = conditionMessage(tree_obs)))
  fit <- tryCatch(
    phylolm::phylolm(y ~ 1, data = sp_df[obs_names, , drop = FALSE],
                       phy = tree_obs, model = "lambda"),
    error = function(e) e)
  if (inherits(fit, "error")) return(list(pred = y_named, error = conditionMessage(fit)))
  beta_hat <- unname(stats::coef(fit)[["(Intercept)"]])
  R <- stats::cov2cor(ape::vcv(tree))
  R <- R[sp_df$species, sp_df$species]
  oi <- which(sp_df$species %in% obs_names)
  mi <- setdiff(seq_len(nrow(sp_df)), oi)
  pred <- rep(beta_hat, nrow(sp_df))
  pred[oi] <- sp_df$y[oi]
  if (length(mi)) {
    e_obs <- sp_df$y[oi] - beta_hat
    blup <- as.numeric(R[mi, oi, drop = FALSE] %*%
                          solve(R[oi, oi, drop = FALSE] + diag(1e-6, length(oi)), e_obs))
    pred[mi] <- beta_hat + blup
  }
  names(pred) <- sp_df$species
  list(pred = pred[names(y_named)], error = NULL)
}

fit_phylolm_method <- function(df_miss4, tree) {
  out  <- df_miss4
  errs <- list()
  for (v in names(df_miss4)) {
    y   <- stats::setNames(df_miss4[[v]], rownames(df_miss4))
    res <- fit_phylolm_trait(y, tree)
    out[[v]] <- unname(res$pred)
    if (!is.null(res$error)) errs[[v]] <- res$error
  }
  list(completed = out, errors = errs)
}

fit_column_mean_method <- function(df_miss4, mask4) {
  out <- df_miss4
  for (v in names(out)) {
    tr_mean <- mean(out[[v]], na.rm = TRUE)
    out[[v]][mask4[, v]] <- tr_mean
  }
  list(completed = out, error = NULL)
}

fit_missforest_method <- function(df_miss_all) {
  if (!requireNamespace("missForest", quietly = TRUE)) {
    return(list(completed = NULL, error = "missForest not installed"))
  }
  fit <- tryCatch(missForest::missForest(df_miss_all, verbose = FALSE), error = function(e) e)
  if (inherits(fit, "error")) return(list(completed = NULL, error = conditionMessage(fit)))
  list(completed = fit$ximp, error = NULL)
}

fit_bace_method <- function(df_miss_all, tree, mask_all) {
  if (!requireNamespace("BACE", quietly = TRUE)) {
    return(list(completed = NULL, error = "BACE not installed"))
  }
  tree_b <- tree
  tree_b$edge.length[tree_b$edge.length == 0] <- 1e-8
  data_b <- df_miss_all
  data_b$Species <- rownames(data_b)
  traits <- names(df_miss_all)
  formulae <- lapply(traits, function(v) paste0(v, " ~ ", paste(setdiff(traits, v), collapse = " + ")))
  fit <- tryCatch(BACE::bace(fixformula = formulae, ran_phylo_form = "~ 1 |Species",
                             phylo = tree_b, data = data_b, nitt = 2000L,
                             burnin = 500L, thin = 5L, runs = 2L, n_final = 2L,
                             verbose = FALSE, skip_conv = TRUE), error = function(e) e)
  if (inherits(fit, "error")) return(list(completed = NULL, error = conditionMessage(fit)))
  draws <- fit$imputed_datasets %||% fit$imputed_data %||% if (!is.null(fit$data)) list(fit$data) else NULL
  if (is.null(draws) || !length(draws)) return(list(completed = NULL, error = "BACE output shape not recognised"))
  out <- df_miss_all
  for (v in traits) {
    for (i in which(mask_all[, v])) {
      vals <- vapply(draws, function(d) as.character(d[[v]][i]), character(1))
      vals <- vals[!is.na(vals)]
      if (!length(vals)) next
      if (is.factor(out[[v]])) {
        mode <- names(sort(table(vals), decreasing = TRUE))[1L]
        out[[v]][i] <- factor(mode, levels = levels(out[[v]]), ordered = is.ordered(out[[v]]))
      } else {
        out[[v]][i] <- stats::median(as.numeric(vals), na.rm = TRUE)
      }
    }
  }
  list(completed = out, error = NULL)
}

# -------------------------------------------------------------------------
# 4. Scoring -- z-scored RMSE + Pearson r on held-out cells. z-scale
#    (mean/sd) is computed from the TRAINING portion only (non-held-out
#    observed cells for that trait/rep) to avoid leaking the truth
#    scale into the metric.
# -------------------------------------------------------------------------

score_trait <- function(pred_v, truth_v, mask_v, train_mean, train_sd) {
  idx <- which(mask_v & is.finite(truth_v) & is.finite(pred_v))
  n <- length(idx)
  if (n < 5L || !is.finite(train_sd) || train_sd < 1e-10) {
    return(c(rmse_z = NA_real_, pearson_r = NA_real_, n_cells = n))
  }
  z_t <- (truth_v[idx] - train_mean) / train_sd
  z_p <- (pred_v[idx]  - train_mean) / train_sd
  rmse <- sqrt(mean((z_p - z_t)^2))
  r <- if (stats::sd(z_p) > 1e-10) suppressWarnings(stats::cor(z_p, z_t)) else NA_real_
  c(rmse_z = rmse, pearson_r = r, n_cells = n)
}

# -------------------------------------------------------------------------
# 5. Main sweep
# -------------------------------------------------------------------------

all_rows <- list()

for (ds_name in names(DATASETS)) {
  d <- DATASETS[[ds_name]]()
  df <- d$df; tree <- d$tree; cont_traits <- d$cont_traits
  log_line(sprintf("Dataset %s: n=%d species, %d traits (%d scored continuous)",
                    d$name, nrow(df), ncol(df), length(cont_traits)))

  for (rep_id in seq_len(N_REPS)) {
    mask_seed <- MASK_SEEDS[[rep_id]]
    log_line(sprintf("--- rep %d/%d (mask_seed=%d) ---", rep_id, N_REPS, mask_seed))

    mask_test <- make_mask(df, MISS_FRAC, mask_seed)
    df_miss_all <- df
    for (v in names(df)) df_miss_all[[v]][mask_test[, v]] <- NA
    df_miss4 <- df_miss_all[, cont_traits, drop = FALSE]
    mask4    <- mask_test[, cont_traits, drop = FALSE]

    train_stats <- lapply(cont_traits, function(v) {
      obs <- df[[v]][!mask4[, v]]
      c(mean = mean(obs, na.rm = TRUE), sd = stats::sd(obs, na.rm = TRUE))
    })
    names(train_stats) <- cont_traits

    methods <- list()

    t0 <- proc.time()[["elapsed"]]
    m <- fit_pigauto_method(df_miss_all, tree, seed = mask_seed,
                            predict_method = "per_column")
    wall <- proc.time()[["elapsed"]] - t0
    log_line(sprintf("  pigauto done in %.1fs%s", wall,
                      if (!is.null(m$error)) paste0(" -- ERROR: ", m$error) else ""))
    methods$pigauto_default <- list(completed = m$completed, wall = wall,
                              errors = if (!is.null(m$error))
                                stats::setNames(as.list(rep(m$error, length(cont_traits))), cont_traits)
                              else list())

    t0 <- proc.time()[["elapsed"]]
    m <- fit_pigauto_method(df_miss_all, tree, seed = mask_seed,
                            predict_method = "exact")
    wall <- proc.time()[["elapsed"]] - t0
    log_line(sprintf("  pigauto_exact done in %.1fs%s", wall,
                      if (!is.null(m$error)) paste0(" -- ERROR: ", m$error) else ""))
    methods$pigauto_exact <- list(completed = m$completed, wall = wall,
                              errors = if (!is.null(m$error))
                                stats::setNames(as.list(rep(m$error, length(cont_traits))), cont_traits)
                              else list())

    t0 <- proc.time()[["elapsed"]]
    m <- fit_rphylopars_method(df_miss4, tree)
    wall <- proc.time()[["elapsed"]] - t0
    log_line(sprintf("  rphylopars done in %.1fs%s", wall,
                      if (!is.null(m$error)) paste0(" -- ERROR: ", m$error) else ""))
    methods$rphylopars <- list(completed = m$completed, wall = wall,
                                 errors = if (!is.null(m$error))
                                   stats::setNames(as.list(rep(m$error, length(cont_traits))), cont_traits)
                                 else list())

    t0 <- proc.time()[["elapsed"]]
    m <- fit_phylolm_method(df_miss4, tree)
    wall <- proc.time()[["elapsed"]] - t0
    log_line(sprintf("  phylolm_lambda done in %.1fs%s", wall,
                      if (length(m$errors)) paste0(" -- errors on: ", paste(names(m$errors), collapse = ", ")) else ""))
    methods$phylolm_lambda <- list(completed = m$completed, wall = wall, errors = m$errors)

    t0 <- proc.time()[["elapsed"]]
    m <- fit_column_mean_method(df_miss4, mask4)
    wall <- proc.time()[["elapsed"]] - t0
    log_line(sprintf("  column_mean done in %.1fs", wall))
    methods$column_mean <- list(completed = m$completed, wall = wall, errors = list())

    t0 <- proc.time()[["elapsed"]]
    m <- fit_missforest_method(df_miss_all)
    wall <- proc.time()[["elapsed"]] - t0
    methods$missforest <- list(completed = m$completed, wall = wall,
      errors = if (!is.null(m$error)) stats::setNames(as.list(rep(m$error, length(cont_traits))), cont_traits) else list())

    t0 <- proc.time()[["elapsed"]]
    m <- fit_bace_method(df_miss_all, tree, mask_test)
    wall <- proc.time()[["elapsed"]] - t0
    methods$bace <- list(completed = m$completed, wall = wall,
      errors = if (!is.null(m$error)) stats::setNames(as.list(rep(m$error, length(cont_traits))), cont_traits) else list())

    for (m_name in names(methods)) {
      mobj <- methods[[m_name]]
      for (v in cont_traits) {
        err_v <- mobj$errors[[v]]
        if (is.null(mobj$completed) || is.null(mobj$completed[[v]])) {
          if (is.null(err_v)) err_v <- "method produced no completed output"
          all_rows[[length(all_rows) + 1L]] <- data.frame(
            dataset = d$name, method = m_name, trait = v, rep = rep_id,
            mask_seed = mask_seed, rmse_z = NA_real_, pearson_r = NA_real_,
            n_cells = NA_integer_, wall_s = mobj$wall, error = err_v,
            stringsAsFactors = FALSE)
          next
        }
        sc <- score_trait(mobj$completed[[v]], df[[v]], mask4[, v],
                            train_stats[[v]]["mean"], train_stats[[v]]["sd"])
        all_rows[[length(all_rows) + 1L]] <- data.frame(
          dataset = d$name, method = m_name, trait = v, rep = rep_id,
          mask_seed = mask_seed,
          rmse_z = unname(sc["rmse_z"]), pearson_r = unname(sc["pearson_r"]),
          n_cells = unname(sc["n_cells"]), wall_s = mobj$wall,
          error = if (is.null(err_v)) NA_character_ else err_v,
          stringsAsFactors = FALSE)
      }
    }
    saveRDS(list(
      protocol = "stage_b_continuous_core_v1",
      dataset = d$name,
      rep = rep_id,
      mask_seed = mask_seed,
      miss_frac = MISS_FRAC,
      source_sha = Sys.getenv("PIGAUTO_SOURCE_SHA", unset = NA_character_),
      mask = mask_test,
      truth = df,
      masked = df_miss_all,
      methods = methods,
      rows = all_rows[(length(all_rows) - length(methods) * length(cont_traits) + 1L):length(all_rows)]
    ), file.path(receipt_dir, sprintf("stage_b_continuous_core_mask-%d.rds", mask_seed)))
  }
}

results_df <- do.call(rbind, all_rows)
rownames(results_df) <- NULL

# -------------------------------------------------------------------------
# 6. Aggregate across reps: mean +/- MCSE (sd/sqrt(n_reps_ok))
# -------------------------------------------------------------------------

agg_key <- interaction(results_df$dataset, results_df$method, results_df$trait, drop = TRUE)
summary_rows <- lapply(split(results_df, agg_key), function(g) {
  rmse_ok <- g$rmse_z[!is.na(g$rmse_z)]
  r_ok    <- g$pearson_r[!is.na(g$pearson_r)]
  errs    <- unique(g$error[!is.na(g$error)])
  data.frame(
    dataset = g$dataset[1], method = g$method[1], trait = g$trait[1],
    n_reps_ok   = length(rmse_ok),
    rmse_z_mean = if (length(rmse_ok)) mean(rmse_ok) else NA_real_,
    rmse_z_mcse = if (length(rmse_ok) > 1) stats::sd(rmse_ok) / sqrt(length(rmse_ok)) else NA_real_,
    pearson_r_mean = if (length(r_ok)) mean(r_ok) else NA_real_,
    pearson_r_mcse = if (length(r_ok) > 1) stats::sd(r_ok) / sqrt(length(r_ok)) else NA_real_,
    wall_s_mean = mean(g$wall_s, na.rm = TRUE),
    errors = if (length(errs)) paste(errs, collapse = " | ") else NA_character_,
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, summary_rows)
summary_df <- summary_df[order(summary_df$dataset, summary_df$trait, summary_df$method), ]
rownames(summary_df) <- NULL

error_rows <- results_df[!is.na(results_df$error), c("dataset", "method", "trait", "rep", "error")]
rownames(error_rows) <- NULL

total_wall <- proc.time()[["elapsed"]] - script_start

saveRDS(list(
  results   = results_df,
  summary   = summary_df,
  errors    = error_rows,
  miss_frac = MISS_FRAC,
  mask_seeds = MASK_SEEDS,
  n_reps    = N_REPS,
  total_wall_s = total_wall
), out_rds)

# -------------------------------------------------------------------------
# 7. Markdown report
# -------------------------------------------------------------------------

fmt_num <- function(x, digits = 3) ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))

summary_display <- summary_df
summary_display$rmse_z <- sprintf("%s +/- %s", fmt_num(summary_display$rmse_z_mean), fmt_num(summary_display$rmse_z_mcse))
summary_display$pearson_r <- sprintf("%s +/- %s", fmt_num(summary_display$pearson_r_mean), fmt_num(summary_display$pearson_r_mcse))
summary_display <- summary_display[, c("dataset", "trait", "method", "n_reps_ok", "rmse_z", "pearson_r", "wall_s_mean", "errors")]

error_block <- if (nrow(error_rows)) {
  c("## Rival failures (verbatim)", "",
    "Any method that errored on a trait is recorded as NA above; the",
    "verbatim error message is listed here per (method, trait, rep) --",
    "never silently dropped.", "", "```",
    capture.output(print(error_rows, row.names = FALSE, max = 500)),
    "```")
} else {
  c("## Rival failures", "", "None -- all methods produced output for all traits/reps.")
}

md <- c(
  "# Stage B continuous core: pigauto vs external comparators -- AVONET300",
  "",
  "## Regime",
  "",
  sprintf("- Dataset: bundled AVONET300 (n = %d species, single-obs).", nrow(DATASETS[[1]]()$df)),
  sprintf("- Scored traits (4 continuous, all methods): %s.",
          paste(DATASETS[[1]]()$cont_traits, collapse = ", ")),
  "- pigauto is fit on all 7 AVONET traits (its unified mixed-type",
  "  imputation across continuous/categorical/ordinal is the point of",
  "  the package) but scored only on the 4 continuous traits above,",
  "  for comparability with the continuous-only rivals.",
  sprintf("- Missingness: %.0f%% MCAR on observed cells, %d locked masks (seeds %s).",
          MISS_FRAC * 100, N_REPS, paste(MASK_SEEDS, collapse = ", ")),
  "- Metric: z-scored RMSE (scale = training-portion mean/sd per",
  "  trait/rep, no leakage from held-out truth) + Pearson r on",
  "  held-out cells. Reported as mean +/- MCSE (sd/sqrt(n_reps_ok)).",
  "",
  "## Methods",
  "",
  "1. `pigauto::impute(predict_method = \"per_column\")` -- production defaults (all 7 traits fit; scored on 4).",
  "2. `pigauto::impute(predict_method = \"exact\")` -- opt-in exact conditional route; it does not change the default.",
  "3. `Rphylopars::phylopars(model = \"BM\")` standalone -- continuous traits only.",
  "4. `phylolm(model = \"lambda\")` per trait, no covariates (y ~ 1) --",
  "   phylogenetic BLUP, mirrors the arm in `script/bench_lambda_sweep.R`",
  "   (`m_phylolm_blup`, ~line 155) with the fixed-effect covariates dropped.",
  "5. `missForest::missForest()` -- non-phylogenetic mixed-type comparator.",
  "6. `BACE::bace()` -- phylogenetic chained-equation comparator (errors retained).",
  "7. Column mean -- floor.",
  "",
  "## IMPORTANT: what the pigauto-vs-phylopars comparison actually measures",
  "",
  "`Rphylopars::phylopars()` is ALSO the solver pigauto calls internally",
  "for its joint-baseline path (`R/joint_mvn_baseline.R`,",
  "`R/joint_threshold_baseline.R`, via `R/joint_mvn_solver.R`, both using",
  "`model = \"BM\"`). So `pigauto` vs `rphylopars` here is NOT a comparison",
  "against an unrelated package -- it measures what pigauto's GNN + gating",
  "+ conformal layer adds on top of the raw phylopars solver it can",
  "already call as its own baseline.",
  "",
  "This is only the continuous-core regime. It is not the mixed-type Stage-B",
  "regime and cannot be used to erase pigauto's mixed-type capability or to",
  "support a parity, superiority, or default-change claim.",
  "",
  "## Per-trait summary (mean +/- MCSE across reps)",
  "",
  "```",
  capture.output(print(summary_display, row.names = FALSE, max = 500)),
  "```",
  "",
  error_block,
  "",
  sprintf("Total wall time: %.1f min.", total_wall / 60)
)
writeLines(md, out_md)

log_line("=== DONE ===")
log_line("  rds: ", out_rds)
log_line("  md:  ", out_md)
log_line(sprintf("  total wall: %.1f min", total_wall / 60))
