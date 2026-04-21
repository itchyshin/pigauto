#!/usr/bin/env Rscript
#
# script/bench_fishbase.R
#
# Third real-data benchmark completing the vertebrate breadth triad
# (birds via AVONET, mammals via PanTHERIA, fish via FishBase + fishtree).
#
# Data
#   Phylogeny:  fishtree::fishtree_phylogeny()  -- Rabosky et al. 2018
#                all-fish tree with branch lengths (~11,600 species).
#   Traits:     rfishbase::species()  + rfishbase::ecology()  curated by
#                the FishBase team; merged on Species.
#
# Traits (mixed types, 6):
#   Length           numeric  (max TL, cm)        continuous (auto-log)
#   Weight           numeric  (max wt, g)         continuous (auto-log)
#   DepthRangeDeep   numeric  (m)                 continuous (auto-log)
#   Vulnerability    numeric  (0-100)             continuous
#   Troph            numeric  (trophic level)     continuous
#   BodyShapeI       factor   (shape class)       categorical
#
# Design
#   - Missingness:  30% MCAR held out per-trait (seed 2026).
#   - val_frac:     0.5  (matches the other realised benches).
#   - Scope:        controlled by PIGAUTO_N_SUBSET env var. Default is the
#                   full matched intersection of fishtree + FishBase
#                   (expect ~5,000-8,000 species).  Local testing:
#                   PIGAUTO_N_SUBSET=500 for <5-min iterations.
#   - Methods:      mean/mode baseline vs pigauto_default.  BACE skipped
#                   by default (opt in with PIGAUTO_RUN_BACE=1).
#   - Metrics:      per-trait RMSE/pearson_r (continuous), accuracy
#                   (categorical), plus conformal + MC-dropout coverage.
#
# Output
#   script/bench_fishbase.rds         tidy data.frame of results
#   script/bench_fishbase.md          human-readable summary
#
# Run
#   # Quick local smoke:
#   PIGAUTO_N_SUBSET=500 Rscript script/bench_fishbase.R
#   # Full matched set:
#   Rscript script/bench_fishbase.R
#
# Caching: rfishbase + fishtree cache downloads into ~/.cache/rfishbase
# and ~/Library/Caches/org.R-project.R/R/fishtree so rerunning the
# script does not re-hit the network.

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  library(rfishbase)
  library(fishtree)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto",
    quiet = TRUE
  )
})

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "bench_fishbase.rds")
out_md  <- file.path(here, "script", "bench_fishbase.md")

SEED      <- 2026L
MISS_FRAC <- 0.30
N_SUBSET <- {
  v <- Sys.getenv("PIGAUTO_N_SUBSET", unset = "")
  if (nzchar(v)) as.integer(v) else NA_integer_
}
RUN_BACE <- identical(Sys.getenv("PIGAUTO_RUN_BACE"), "1") &&
            requireNamespace("BACE", quietly = TRUE)

script_start <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start),
      ..., "\n", sep = "")
  flush.console()
}

# -------------------------------------------------------------------------
# 1. Load phylogeny + trait tables
# -------------------------------------------------------------------------

log_line("Loading fishtree phylogeny (Rabosky et al. 2018) ...")
tree <- fishtree::fishtree_phylogeny()
log_line(sprintf("  tree tips: %d  (FIN subset by default)", length(tree$tip.label)))

# fishtree labels are "Genus_species". FishBase uses "Genus species".
tree$tip.label <- gsub("_", " ", tree$tip.label)

log_line("Loading rfishbase taxonomy (SpecCode + Species binomial) ...")
taxa <- rfishbase::load_taxa()
# rfishbase 5.x's load_taxa() already stores the binomial in $Species
# (e.g. "Aapticheilichthys websteri", not just the epithet). Use as-is.
taxa <- taxa[, c("SpecCode", "Species"), drop = FALSE]
log_line(sprintf("  taxa rows: %d", nrow(taxa)))

log_line("Loading rfishbase species() trait fields ...")
sp_tbl <- rfishbase::species()
keep_fields <- c("SpecCode", "Length", "Weight", "BodyShapeI",
                  "DepthRangeDeep", "Vulnerability")
keep_fields <- intersect(keep_fields, names(sp_tbl))
sp_tbl <- sp_tbl[, keep_fields, drop = FALSE]
sp_tbl <- merge(taxa, sp_tbl, by = "SpecCode", all.x = FALSE)
log_line(sprintf("  species() rows after taxa join: %d", nrow(sp_tbl)))

log_line("Loading rfishbase ecology() (Troph) ...")
ec_tbl <- rfishbase::ecology()
if ("SpecCode" %in% names(ec_tbl)) {
  troph_cols <- intersect(c("FoodTroph", "DietTroph"), names(ec_tbl))
  ec_tbl <- ec_tbl[, c("SpecCode", troph_cols), drop = FALSE]
  # prefer DietTroph, fall back to FoodTroph
  if (all(c("DietTroph", "FoodTroph") %in% names(ec_tbl))) {
    ec_tbl$Troph <- ifelse(is.finite(ec_tbl$DietTroph),
                            ec_tbl$DietTroph, ec_tbl$FoodTroph)
  } else if ("DietTroph" %in% names(ec_tbl)) {
    ec_tbl$Troph <- ec_tbl$DietTroph
  } else if ("FoodTroph" %in% names(ec_tbl)) {
    ec_tbl$Troph <- ec_tbl$FoodTroph
  } else {
    ec_tbl$Troph <- NA_real_
  }
  ec_tbl <- ec_tbl[!is.na(ec_tbl$Troph), c("SpecCode", "Troph"),
                    drop = FALSE]
  ec_tbl <- ec_tbl[!duplicated(ec_tbl$SpecCode), , drop = FALSE]
  log_line(sprintf("  ecology rows with Troph: %d", nrow(ec_tbl)))
} else {
  ec_tbl <- data.frame(SpecCode = integer(0), Troph = numeric(0))
  log_line("  ecology() missing SpecCode — Troph will be all NA")
}

# -------------------------------------------------------------------------
# 2. Merge and align to tree
# -------------------------------------------------------------------------

df <- merge(sp_tbl, ec_tbl, by = "SpecCode", all.x = TRUE)
df$SpecCode <- NULL

# Keep only species in the tree
df <- df[df$Species %in% tree$tip.label, , drop = FALSE]

# Clean BodyShapeI: FishBase has casing/spacing inconsistencies
# ("elongated" vs "Elongated"; "fusiform / normal" with spaces).
# Normalise to lowercase + trimmed, then collapse rare levels to NA
# to avoid rank issues downstream.
if ("BodyShapeI" %in% names(df)) {
  df$BodyShapeI[df$BodyShapeI == ""] <- NA
  df$BodyShapeI <- tolower(trimws(df$BodyShapeI))
  df$BodyShapeI <- gsub("\\s+", " ", df$BodyShapeI)
  df$BodyShapeI <- factor(df$BodyShapeI)
  counts <- table(df$BodyShapeI, useNA = "no")
  keep_levels <- names(counts[counts >= 50L])
  df$BodyShapeI[!as.character(df$BodyShapeI) %in% keep_levels] <- NA
  df$BodyShapeI <- droplevels(factor(df$BodyShapeI))
}

# Continuous columns: force numeric
for (v in c("Length", "Weight", "DepthRangeDeep", "Vulnerability", "Troph")) {
  if (v %in% names(df)) df[[v]] <- as.numeric(df[[v]])
}

# Require at least one trait present per species
trait_cols <- intersect(c("Length", "Weight", "DepthRangeDeep",
                           "Vulnerability", "Troph", "BodyShapeI"),
                          names(df))
has_any <- rowSums(!is.na(df[, trait_cols, drop = FALSE])) > 0L
df <- df[has_any, , drop = FALSE]

rownames(df) <- df$Species
df$Species <- NULL

# Prune tree to matching species
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(df)))
df   <- df[tree$tip.label, , drop = FALSE]

log_line(sprintf("Matched fishtree × FishBase: %d species × %d traits",
                  nrow(df), ncol(df)))

if (!is.na(N_SUBSET) && N_SUBSET < nrow(df)) {
  set.seed(SEED)
  keep <- sample(rownames(df), N_SUBSET)
  df   <- df[keep, , drop = FALSE]
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
  df   <- df[tree$tip.label, , drop = FALSE]
  log_line(sprintf("N_SUBSET = %d: random subset -> %d species",
                    N_SUBSET, nrow(df)))
}

# -------------------------------------------------------------------------
# 3. MCAR 30% held-out mask
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
# 4. Eval helper (captures RMSE/accuracy + both coverage types)
# -------------------------------------------------------------------------

safe_cor <- function(x, y) {
  idx <- which(is.finite(x) & is.finite(y))
  if (length(idx) < 2L) return(NA_real_)
  if (stats::sd(x[idx]) == 0 || stats::sd(y[idx]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[idx], y[idx]))
}

eval_completed <- function(completed, truth, mask, method, wall_s,
                            res_obj = NULL) {
  if (is.null(completed)) return(NULL)
  lo      <- if (!is.null(res_obj)) res_obj$prediction$conformal_lower else NULL
  hi      <- if (!is.null(res_obj)) res_obj$prediction$conformal_upper else NULL
  mi_list <- if (!is.null(res_obj)) res_obj$prediction$imputed_datasets else NULL

  rows <- list()
  for (v in colnames(mask)) {
    idx <- which(mask[, v])
    if (!length(idx)) next
    t_v <- truth[[v]][idx]; c_v <- completed[[v]][idx]
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
      t_num <- as.numeric(t_v)
      if (!is.null(lo) && !is.null(hi) && v %in% colnames(lo)) {
        lo_v <- lo[idx, v]; hi_v <- hi[idx, v]
        valid <- is.finite(lo_v) & is.finite(hi_v) & is.finite(t_num)
        if (any(valid)) {
          hits <- t_num[valid] >= lo_v[valid] & t_num[valid] <= hi_v[valid]
          rows[[length(rows) + 1L]] <- data.frame(
            method = method, trait = v, metric = "coverage95_conformal",
            value = mean(hits), n_cells = sum(valid), wall_s = wall_s)
        }
      }
      if (!is.null(mi_list) && length(mi_list) > 1L &&
          v %in% names(mi_list[[1]])) {
        draws_mat <- vapply(mi_list, function(d) as.numeric(d[idx, v]),
                             numeric(length(idx)))
        if (!is.matrix(draws_mat))
          draws_mat <- matrix(draws_mat, ncol = length(mi_list))
        if (nrow(draws_mat) == length(idx) && ncol(draws_mat) > 1L) {
          q_lo <- apply(draws_mat, 1L, stats::quantile, probs = 0.025,
                         na.rm = TRUE)
          q_hi <- apply(draws_mat, 1L, stats::quantile, probs = 0.975,
                         na.rm = TRUE)
          valid <- is.finite(q_lo) & is.finite(q_hi) & is.finite(t_num)
          if (any(valid)) {
            hits <- t_num[valid] >= q_lo[valid] & t_num[valid] <= q_hi[valid]
            rows[[length(rows) + 1L]] <- data.frame(
              method = method, trait = v, metric = "coverage95_mcdropout",
              value = mean(hits), n_cells = sum(valid), wall_s = wall_s)
          }
        }
      }
    }
  }
  if (!length(rows)) NULL else do.call(rbind, rows)
}

# -------------------------------------------------------------------------
# 5. Methods
# -------------------------------------------------------------------------

run_mean <- function() {
  out <- df_miss
  for (v in names(out)) {
    if (is.factor(out[[v]])) {
      mode <- names(sort(table(out[[v]], useNA = "no"),
                           decreasing = TRUE))[1]
      out[[v]][mask_test[, v]] <- factor(mode, levels = levels(out[[v]]),
                                           ordered = is.ordered(out[[v]]))
    } else {
      out[[v]][mask_test[, v]] <- mean(out[[v]], na.rm = TRUE)
    }
  }
  list(completed = out, res = NULL)
}

run_pigauto <- function() {
  res <- pigauto::impute(df_miss, tree,
                           log_transform  = TRUE,
                           missing_frac   = 0.20,
                           n_imputations  = 20L,
                           epochs         = 500L,
                           verbose        = TRUE,
                           seed           = SEED)
  list(completed = res$completed, res = res)
}

run_bace <- function() {
  if (!RUN_BACE) {
    log_line("BACE skipped (PIGAUTO_RUN_BACE != 1 or BACE not installed)")
    return(list(completed = NULL, res = NULL))
  }
  out <- tryCatch({
    BACE::bace(data = df_miss, tree = tree,
                n_iter = 2000L, burnin = 500L, thin = 5L,
                ovr = TRUE, verbose = FALSE)
  }, error = function(e) {
    log_line("BACE run failed: ", conditionMessage(e)); NULL
  })
  if (is.null(out)) return(list(completed = NULL, res = NULL))
  completed <- tryCatch({
    if ("completed" %in% names(out)) out$completed
    else if ("data" %in% names(out)) out$data
    else if ("imputed_data" %in% names(out)) out$imputed_data
    else NULL
  }, error = function(e) NULL)
  list(completed = completed, res = NULL)
}

# -------------------------------------------------------------------------
# 6. Run methods
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
                            "mean_baseline",   r_mean$wall, NULL)
ev_pig  <- eval_completed(r_pig$val$completed,  df_truth, mask_test,
                            "pigauto_default", r_pig$wall, r_pig$val$res)
ev_bace <- eval_completed(r_bace$val$completed, df_truth, mask_test,
                            "bace_default",    r_bace$wall, NULL)

all_rows <- do.call(rbind, Filter(Negate(is.null),
                                     list(ev_mean, ev_pig, ev_bace)))

saveRDS(list(
  results   = all_rows,
  seed      = SEED,
  miss_frac = MISS_FRAC,
  n_species = nrow(df),
  n_subset  = N_SUBSET,
  bace_ran  = !is.null(ev_bace)
), out_rds)

md <- c(
  "# FishBase + fishtree x pigauto + BACE",
  "",
  sprintf("n = %d species x %d traits (fishtree x FishBase matched).",
          nrow(df), ncol(df)),
  sprintf("Seed = %d, miss_frac = %.2f.", SEED, MISS_FRAC),
  if (is.null(ev_bace)) "**BACE skipped** (not installed or opted out)." else "",
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
