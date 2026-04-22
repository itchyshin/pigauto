#!/usr/bin/env Rscript
#
# script/bench_avonet_bace.R
#
# pigauto vs BACE head-to-head on AVONET data. Two modes:
#
#   * default (bundled AVONET 300 -- ~5-10 min total, local smoke)
#   * PIGAUTO_BENCH_FULL=1   -- uses full avonet_full + tree_full
#                               (9,993 species). BACE may be slow at
#                               this scale; set PIGAUTO_N_SUBSET to
#                               subset first (e.g. 3000).
#
# Traits (4 continuous + 2 categorical + 1 ordinal) mirror AVONET:
#   Mass, Beak.Length_Culmen, Tarsus.Length, Wing.Length,
#   Trophic.Level, Primary.Lifestyle, Migration.
#
# The script fixes the long-standing broken BACE::bace() call in the
# Vulcan GPU bundle: that driver passed `data=, tree=, n_iter=, ovr=`
# while the real BACE API takes `fixformula=, ran_phylo_form=, phylo=,
# nitt=`. Every Vulcan bench that claimed to "include BACE" actually
# silently skipped it with "BACE run failed: argument fixformula is
# missing".
#
# Output:
#   script/bench_avonet_bace.rds
#   script/bench_avonet_bace.md

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  # Prefer source load if PIGAUTO_PKG_PATH is set; otherwise library().
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

SEED      <- 2026L
MISS_FRAC <- 0.30

FULL    <- identical(Sys.getenv("PIGAUTO_BENCH_FULL"), "1")
N_SUB <- {
  v <- Sys.getenv("PIGAUTO_N_SUBSET", unset = "")
  if (nzchar(v)) as.integer(v) else NA_integer_
}
N_IMP <- {
  v <- Sys.getenv("PIGAUTO_N_IMPUTATIONS", unset = "")
  if (nzchar(v)) as.integer(v) else 20L
}

script_start <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start),
      ..., "\n", sep = "")
  flush.console()
}

out_rds <- "script/bench_avonet_bace.rds"
out_md  <- "script/bench_avonet_bace.md"

# -------------------------------------------------------------------------
# 1. Load AVONET + tree
# -------------------------------------------------------------------------

e <- new.env(parent = emptyenv())
if (FULL) {
  log_line("Loading avonet_full + tree_full (9,993 species) ...")
  utils::data("avonet_full", package = "pigauto", envir = e)
  utils::data("tree_full",   package = "pigauto", envir = e)
  df   <- e$avonet_full
  tree <- e$tree_full
  rownames(df) <- df$Species_Key
  df$Species_Key <- NULL
} else {
  log_line("Loading bundled AVONET 300 (avonet300 + tree300) ...")
  utils::data("avonet300", package = "pigauto", envir = e)
  utils::data("tree300",   package = "pigauto", envir = e)
  df   <- e$avonet300
  tree <- e$tree300
  rownames(df) <- df$Species_Key
  df$Species_Key <- NULL
}

stopifnot(all(rownames(df) == tree$tip.label))
log_line(sprintf("Aligned: %d species x %d traits", nrow(df), ncol(df)))

if (!is.na(N_SUB) && N_SUB < nrow(df)) {
  set.seed(SEED)
  keep <- sample(rownames(df), N_SUB)
  df   <- df[keep, , drop = FALSE]
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
  df   <- df[tree$tip.label, , drop = FALSE]
  log_line(sprintf("N_SUBSET = %d: subset -> %d species", N_SUB, nrow(df)))
}

# -------------------------------------------------------------------------
# 2. MCAR 30% held-out mask
# -------------------------------------------------------------------------

set.seed(SEED)
df_truth <- df
mask_test <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df),
                     dimnames = list(rownames(df), names(df)))
for (v in names(df)) {
  obs_idx <- which(!is.na(df[[v]]))
  to_hide <- sample(obs_idx, ceiling(length(obs_idx) * MISS_FRAC))
  mask_test[to_hide, v] <- TRUE
}
df_miss <- df
for (v in names(df)) df_miss[[v]][mask_test[, v]] <- NA
log_line("Held-out test cells per trait:")
print(colSums(mask_test))

# -------------------------------------------------------------------------
# 3. Eval helper
# -------------------------------------------------------------------------

safe_cor <- function(x, y) {
  idx <- which(is.finite(x) & is.finite(y))
  if (length(idx) < 2L) return(NA_real_)
  if (stats::sd(x[idx]) == 0 || stats::sd(y[idx]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[idx], y[idx]))
}

eval_completed <- function(completed, truth, mask, method, wall_s) {
  if (is.null(completed)) return(NULL)
  rows <- list()
  for (v in colnames(mask)) {
    idx <- which(mask[, v])
    if (!length(idx)) next
    t_v <- truth[[v]][idx]
    c_v <- completed[[v]][idx]
    if (is.factor(t_v) || is.ordered(t_v) || is.character(t_v)) {
      acc <- mean(as.character(c_v) == as.character(t_v), na.rm = TRUE)
      rows[[length(rows) + 1L]] <- data.frame(
        method = method, trait = v, metric = "accuracy",
        value = acc, n_cells = length(idx), wall_s = wall_s)
    } else {
      rmse <- sqrt(mean((as.numeric(t_v) - as.numeric(c_v))^2, na.rm = TRUE))
      pear <- safe_cor(as.numeric(t_v), as.numeric(c_v))
      rows[[length(rows) + 1L]] <- data.frame(
        method = method, trait = v, metric = "rmse",
        value = rmse, n_cells = length(idx), wall_s = wall_s)
      rows[[length(rows) + 1L]] <- data.frame(
        method = method, trait = v, metric = "pearson_r",
        value = pear, n_cells = length(idx), wall_s = wall_s)
    }
  }
  if (!length(rows)) NULL else do.call(rbind, rows)
}

# -------------------------------------------------------------------------
# 4. Methods
# -------------------------------------------------------------------------

run_mean <- function() {
  out <- df_miss
  for (v in names(out)) {
    if (is.factor(out[[v]])) {
      mode <- names(sort(table(out[[v]], useNA = "no"), decreasing = TRUE))[1]
      out[[v]][mask_test[, v]] <- factor(mode, levels = levels(out[[v]]),
                                           ordered = is.ordered(out[[v]]))
    } else {
      out[[v]][mask_test[, v]] <- mean(out[[v]], na.rm = TRUE)
    }
  }
  list(completed = out)
}

run_pigauto <- function() {
  res <- pigauto::impute(df_miss, tree,
                           log_transform = TRUE,
                           missing_frac  = 0.20,
                           n_imputations = N_IMP,
                           epochs        = 500L,
                           verbose       = TRUE,
                           seed          = SEED)
  list(completed = res$completed)
}

# Correct BACE::bace() call -- uses fixformula + ran_phylo_form + phylo
# + data (with a Species column). Each trait with missing data gets its
# own imputation formula in the chained-equations scheme.
run_bace <- function() {
  if (!requireNamespace("BACE", quietly = TRUE)) {
    log_line("BACE not installed -- skipping")
    return(list(completed = NULL))
  }
  # BACE (via MCMCglmm) refuses phylogenies with zero-length edges.
  # Replace zero edges with a tiny positive value so MCMCglmm accepts
  # the tree. 1e-8 is below the resolution of any real evolutionary
  # distance so this is a no-op for the model.
  tree_b <- tree
  if (any(tree_b$edge.length == 0, na.rm = TRUE)) {
    n_zero <- sum(tree_b$edge.length == 0, na.rm = TRUE)
    log_line(sprintf("Tree has %d zero-length edges -- patching to 1e-8 for BACE", n_zero))
    tree_b$edge.length[tree_b$edge.length == 0] <- 1e-8
  }

  # Build BACE-shaped data frame with Species column matching phylo tips.
  df_b <- df_miss
  df_b$Species <- rownames(df_miss)

  # One formula per trait: each trait ~ the other traits (chained).
  # BACE uses MCMCglmm under the hood which does not take dots in
  # formulas, so we build them explicitly from the column names.
  all_traits <- setdiff(names(df_b), "Species")
  fixformula <- lapply(all_traits, function(v) {
    others <- setdiff(all_traits, v)
    paste0(v, " ~ ", paste(others, collapse = " + "))
  })

  out <- tryCatch({
    BACE::bace(
      fixformula    = fixformula,
      ran_phylo_form = "~ 1 |Species",
      phylo         = tree_b,
      data          = df_b,
      nitt          = 2000L,
      burnin        = 500L,
      thin          = 5L,
      runs          = 2L,
      n_final       = 2L,
      verbose       = FALSE,
      skip_conv     = TRUE   # small n: skip convergence retry
    )
  }, error = function(e) {
    log_line("BACE run failed: ", conditionMessage(e))
    NULL
  })

  if (is.null(out)) return(list(completed = NULL))

  # BACE returns a list with $imputed_datasets (list of M completed
  # data.frames). Pool by mean / mode across draws for the final
  # completed dataset.
  imputed_sets <- tryCatch({
    if ("imputed_datasets" %in% names(out)) out$imputed_datasets
    else if ("imputed_data" %in% names(out)) out$imputed_data
    else if ("data" %in% names(out)) list(out$data)
    else NULL
  }, error = function(e) NULL)

  if (is.null(imputed_sets) || !length(imputed_sets)) {
    log_line("BACE output shape not recognised (keys: ",
              paste(names(out), collapse = ", "), ")")
    return(list(completed = NULL))
  }

  # Pool: for each cell in mask_test, take mean (continuous) or mode
  # (factor) across the M imputed sets.
  M <- length(imputed_sets)
  completed <- df_miss
  completed$Species <- NULL   # drop the BACE key before returning
  for (v in names(completed)) {
    if (!any(mask_test[, v])) next
    draws <- sapply(imputed_sets, function(d) d[[v]])
    if (!is.matrix(draws)) draws <- matrix(draws, ncol = M)
    if (is.factor(completed[[v]])) {
      for (i in which(mask_test[, v])) {
        vals <- as.character(draws[i, ])
        vals <- vals[!is.na(vals)]
        if (!length(vals)) next
        mode_val <- names(sort(table(vals), decreasing = TRUE))[1]
        completed[[v]][i] <- factor(mode_val, levels = levels(completed[[v]]),
                                     ordered = is.ordered(completed[[v]]))
      }
    } else {
      for (i in which(mask_test[, v])) {
        vals <- as.numeric(draws[i, ])
        completed[[v]][i] <- stats::median(vals, na.rm = TRUE)
      }
    }
  }

  list(completed = completed)
}

# -------------------------------------------------------------------------
# 5. Run all methods
# -------------------------------------------------------------------------

timed <- function(tag, fn) {
  log_line("[STAGE] ", tag, " starting")
  t0 <- proc.time()[["elapsed"]]
  val <- fn()
  wall <- proc.time()[["elapsed"]] - t0
  log_line("[STAGE] ", tag, " done in ", round(wall, 1), " s")
  list(val = val, wall = wall)
}

r_mean <- timed("mean_baseline",   run_mean)
r_pig  <- timed("pigauto_default", run_pigauto)
r_bace <- timed("bace_default",    run_bace)

ev_mean <- eval_completed(r_mean$val$completed, df_truth, mask_test,
                            "mean_baseline",   r_mean$wall)
ev_pig  <- eval_completed(r_pig$val$completed,  df_truth, mask_test,
                            "pigauto_default", r_pig$wall)
ev_bace <- eval_completed(r_bace$val$completed, df_truth, mask_test,
                            "bace_default",    r_bace$wall)

all_rows <- do.call(rbind, Filter(Negate(is.null),
                                     list(ev_mean, ev_pig, ev_bace)))

saveRDS(list(
  results   = all_rows,
  seed      = SEED,
  miss_frac = MISS_FRAC,
  n_species = nrow(df),
  full      = FULL,
  n_subset  = N_SUB,
  n_imputations = N_IMP,
  bace_ran  = !is.null(ev_bace)
), out_rds)

md <- c(
  "# AVONET x pigauto + BACE head-to-head",
  "",
  sprintf("Mode: %s. n = %d species x %d traits.",
          if (FULL) "full" else "bundled AVONET 300",
          nrow(df), ncol(df)),
  sprintf("Seed = %d, miss_frac = %.2f, n_imputations = %d.",
          SEED, MISS_FRAC, N_IMP),
  if (is.null(ev_bace)) "**BACE skipped** (not installed or failed)." else "",
  "",
  "## Per-trait metrics",
  "",
  "```",
  capture.output(print(all_rows, row.names = FALSE, max = 500)),
  "```"
)
writeLines(md, out_md)

log_line("=== DONE ===")
log_line("  rds: ", out_rds)
log_line("  md : ", out_md)
