#!/usr/bin/env Rscript
#
# script/bench_phase_c_cross_dataset.R
#
# Phase C: cross-dataset bench to test whether AVONET findings
# (Phase G clamp_outliers helps Mass tail; Phase H mode pooling
# helps K=3 ordinal but loses on K=5) generalise to other taxa.
#
# Datasets
#   * PanTHERIA (mammals): n~5,000, 8 mixed-type traits.
#     Log-cont mass (analogue of AVONET Mass), gestation, longevity;
#     count litter_size; ordinal diet_breadth (K=5) + habitat_breadth;
#     categorical terrestriality.
#   * AmphiBIO (amphibians): n~6,000, 6 mixed-type traits.
#     Log-cont body_size + body_mass; binary Diu/Noc; ordinal
#     diet_breadth (K=5); categorical Habitat.
#
# Configs (per dataset)
#   1. default                     -- baseline; no opt-ins
#   2. clamp_outliers = TRUE       -- Phase G tail-safety on log-cont
#   3. pool_method = "mode"        -- Phase H mode pooling for ordinal
#   4. both clamp + mode           -- combined opt-ins
#
# Pre-registered questions
#   Q1. Does clamp_outliers help OR hurt log-cont RMSE on a different
#       taxon (mammals)? AVONET evidence: -74 % on seed-2030 Mass,
#       small change on others.  PanTHERIA result will tell us if
#       the pattern generalises.
#   Q2. Does pool_method = "mode" help OR hurt the AmphiBIO K=5
#       ordinal (diet_breadth)?  AVONET Migration K=3 evidence:
#       mode wins +6.6 pp.  Phase H+ simulated K=5 evidence: mode
#       LOSES.  Real K=5 data tells us which pattern dominates.
#   Q3. Does the Mass-instability finding (one of three AVONET seeds
#       produced ~11x baseline RMSE) appear on PanTHERIA mammals?
#       Mammals span 6 orders of magnitude in mass (shrew to whale).
#       If pigauto's BM joint baseline emits singular-Sigma warnings
#       on mammals too, we expect similar instability.
#
# Single seed (2026) for now -- 4 cells per dataset x 2 datasets = 8
# cells.  Wall: ~5 min PanTHERIA per cell + ~2 min AmphiBIO per cell.
# Total ~30 min.  Multi-seed would multiply wall by N_seeds.
#
# Output:
#   script/bench_phase_c_cross_dataset.{rds,md}

options(warn = 1, stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ape)
  pkg_path <- Sys.getenv("PIGAUTO_PKG_PATH", unset = "")
  if (nzchar(pkg_path) && dir.exists(pkg_path) &&
      file.exists(file.path(pkg_path, "DESCRIPTION"))) {
    devtools::load_all(pkg_path, quiet = TRUE)
  } else {
    library(pigauto)
  }
})

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto"
out_rds <- file.path(here, "script", "bench_phase_c_cross_dataset.rds")
out_md  <- file.path(here, "script", "bench_phase_c_cross_dataset.md")
cache_dir <- file.path(here, "script", "data-cache")

# Phase C v2 (2026-05-01): multi-seed + AmphiBIO loader fix
#   * 3 seeds (2026 / 2027 / 2028) to separate config effects from
#     GNN-training MC noise on continuous traits.
#   * AmphiBIO loader drops Diu/Noc (presence-only encoding => degenerate
#     binary; the v1 bench saw a -39.8 pp artefact on Diu under
#     pool_method = "mode" caused by class-imbalance after the NA -> 0
#     hack).  Keeps body_mass / body_size (continuous), diet_breadth
#     (ordinal K=5), habitat (categorical K=4).
SEEDS     <- c(2026L, 2027L, 2028L)
MISS_FRAC <- 0.30
N_IMP     <- 20L
# Subset to match AVONET-bench scale.  PanTHERIA full (n~4600) and
# AmphiBIO full (n~6800) cause OOM on the 200-epoch / N_IMP=20 path
# at standard memory.  n=1500 mirrors the AVONET multi-seed bench
# (where Phase G/H benches were run) so cross-dataset comparison is
# apples-to-apples.
N_SUB     <- 1500L

t0 <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - t0),
      ..., "\n", sep = "")
  flush.console()
}

# ---- Load PanTHERIA --------------------------------------------------------

load_pantheria <- function() {
  pantheria_local <- file.path(cache_dir, "pantheria.txt")
  tree_local      <- file.path(cache_dir, "mammal_tree.tre")
  pan <- utils::read.table(pantheria_local, header = TRUE, sep = "\t",
                            na.strings = "-999", quote = "",
                            stringsAsFactors = FALSE, comment.char = "")
  pan$species_key <- gsub("[^A-Za-z0-9_]", "",
                          paste(pan$MSW93_Genus, pan$MSW93_Species, sep = "_"))
  tree <- ape::read.tree(tree_local)
  tree$tip.label <- gsub("[^A-Za-z0-9_]", "", tree$tip.label)
  overlap <- intersect(tree$tip.label, pan$species_key)
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, overlap))
  pan  <- pan[pan$species_key %in% overlap, ]
  pan  <- pan[match(tree$tip.label, pan$species_key), ]
  rownames(pan) <- pan$species_key

  # Build the mixed-type trait df (4 log-cont, 1 count, 2 ordinal,
  # 1 categorical) -- same selection as bench_pantheria_full.R
  df <- data.frame(row.names = pan$species_key)
  df$body_mass_g          <- as.numeric(pan$X5.1_AdultBodyMass_g)
  df$head_body_length_mm  <- as.numeric(pan$X13.1_AdultHeadBodyLen_mm)
  df$gestation_d          <- as.numeric(pan$X9.1_GestationLen_d)
  df$max_longevity_m      <- as.numeric(pan$X17.1_MaxLongevity_m)
  df$litter_size          <- as.integer(round(as.numeric(pan$X15.1_LitterSize)))
  df$diet_breadth         <- ordered(as.integer(pan$X6.1_DietBreadth),
                                     levels = 1:5)
  df$habitat_breadth      <- ordered(as.integer(pan$X12.1_HabitatBreadth),
                                     levels = 1:3)
  df$terrestriality       <- factor(as.integer(pan$X12.2_Terrestriality),
                                    levels = 1:2)
  # Replace negative / non-finite values in cont with NA
  for (v in c("body_mass_g", "head_body_length_mm", "gestation_d",
              "max_longevity_m")) {
    x <- df[[v]]
    x[x <= 0 | !is.finite(x)] <- NA
    df[[v]] <- x
  }
  list(df = df, tree = tree)
}

# ---- Load AmphiBIO ---------------------------------------------------------

load_amphibio <- function() {
  csv_local  <- file.path(cache_dir, "AmphiBIO_v1.csv")
  tree_local <- file.path(cache_dir, "amphibio_tree.rda")
  amph <- utils::read.csv(csv_local, stringsAsFactors = FALSE)
  e <- new.env(parent = emptyenv())
  load(tree_local, envir = e)
  tree <- e[[ls(e)[1]]]   # whatever the rda variable is called
  amph$species_key <- gsub("[^A-Za-z0-9_]", "_", amph$Species)
  tree$tip.label <- gsub("[^A-Za-z0-9_]", "_", tree$tip.label)
  overlap <- intersect(tree$tip.label, amph$species_key)
  tree <- ape::drop.tip(tree, setdiff(tree$tip.label, overlap))
  amph <- amph[amph$species_key %in% overlap, ]
  amph <- amph[match(tree$tip.label, amph$species_key), ]
  rownames(amph) <- amph$species_key

  # Build mixed-type trait df: 2 log-cont, 2 binary, 1 ordinal, 1 categorical
  df <- data.frame(row.names = amph$species_key)
  df$body_size_mm <- {
    x <- as.numeric(amph$Body_size_mm)
    x[x <= 0 | !is.finite(x)] <- NA
    x
  }
  df$body_mass_g <- {
    x <- as.numeric(amph$Body_mass_g)
    x[x <= 0 | !is.finite(x)] <- NA
    x
  }
  # AmphiBIO Diu / Noc are presence-only (only 1 and NA values; no
  # explicit 0 = "not diurnal").  Phase C v1 saw a -39.8 pp artefact
  # on Diu under mode pooling that was actually class-imbalance noise
  # (NA -> 0 created a 13 % vs 87 % imbalance with mode-class accuracy
  # = 0.85 by default, easy to flip).  Phase C v2 drops Diu / Noc from
  # the AmphiBIO bench entirely: not a real binary signal.
  # Diet breadth: count of 1's across the 6 diet indicator columns
  # (Leaves, Flowers, Seeds, Fruits, Arthro, Vert).  Range 0-6; clip to
  # 1-5 for ordinal encoding.  Species with all-NA diet -> NA breadth.
  diet_cols <- intersect(c("Leaves", "Flowers", "Seeds", "Fruits",
                           "Arthro", "Vert"), names(amph))
  diet_mat <- sapply(diet_cols, function(cn) {
    as.integer(!is.na(amph[[cn]]) & amph[[cn]] == 1L)
  })
  any_diet_data <- rowSums(!is.na(amph[, diet_cols, drop = FALSE])) > 0L
  breadth <- rowSums(diet_mat, na.rm = TRUE)
  breadth[!any_diet_data] <- NA_integer_
  breadth[breadth == 0L]  <- NA_integer_
  breadth <- pmin(breadth, 5L)
  df$diet_breadth <- ordered(as.integer(breadth), levels = 1:5)
  # Habitat (Fos/Ter/Aqu/Arb) -- pick first non-NA, encode as factor
  hab_cols <- intersect(c("Fos", "Ter", "Aqu", "Arb"), names(amph))
  hab_assignment <- apply(amph[, hab_cols, drop = FALSE], 1, function(row) {
    idx <- which(!is.na(row) & row == 1L)
    if (length(idx) == 0L) NA_character_ else hab_cols[idx[1]]
  })
  df$habitat <- factor(hab_assignment, levels = hab_cols)
  list(df = df, tree = tree)
}

# ---- Per-cell run ----------------------------------------------------------

run_one_cell <- function(dataset_name, df, tree, config, seed) {
  log_line(sprintf("Cell %s + %s seed=%d ...", dataset_name, config, seed))

  # Apply mask
  set.seed(seed)
  mask_test <- matrix(FALSE, nrow = nrow(df), ncol = ncol(df),
                      dimnames = list(rownames(df), names(df)))
  for (v in names(df)) {
    obs_idx <- which(!is.na(df[[v]]))
    if (length(obs_idx) == 0L) next
    n_hide <- ceiling(length(obs_idx) * MISS_FRAC)
    to_hide <- sample(obs_idx, n_hide)
    mask_test[to_hide, v] <- TRUE
  }
  df_truth <- df
  df_miss  <- df
  for (v in names(df)) df_miss[[v]][mask_test[, v]] <- NA

  args <- list(traits = df_miss, tree = tree,
               epochs = 200L,
               n_imputations = N_IMP,
               verbose = FALSE,
               seed = seed)
  if (config == "default") {
    # nothing extra
  } else if (config == "clamp") {
    args$clamp_outliers <- TRUE
    args$clamp_factor   <- 5
  } else if (config == "mode") {
    args$pool_method    <- "mode"
  } else if (config == "both") {
    args$clamp_outliers <- TRUE
    args$clamp_factor   <- 5
    args$pool_method    <- "mode"
  } else {
    stop("Unknown config: ", config)
  }

  res <- tryCatch(do.call(pigauto::impute, args),
                  error = function(e) {
                    log_line(sprintf("  ERROR: %s", conditionMessage(e)))
                    NULL
                  })
  if (is.null(res)) return(NULL)

  out <- list()
  for (tr in names(df)) {
    idx <- which(mask_test[, tr])
    if (length(idx) == 0L) next
    truth <- df_truth[[tr]][idx]
    pred  <- res$completed[[tr]][idx]
    if (is.factor(truth) || is.ordered(truth) || is.character(truth)) {
      acc <- mean(as.character(pred) == as.character(truth), na.rm = TRUE)
      mode_class <- names(sort(table(df_miss[[tr]]), decreasing = TRUE))[1]
      base_acc <- mean(as.character(truth) == mode_class, na.rm = TRUE)
      out[[tr]] <- data.frame(
        dataset = dataset_name, seed = seed, config = config, trait = tr,
        type = if (is.ordered(truth)) "ordinal"
               else if (nlevels(as.factor(truth)) == 2L) "binary"
               else "categorical",
        metric = "accuracy",
        value  = acc, baseline = base_acc, n = length(idx),
        stringsAsFactors = FALSE
      )
    } else {
      truth_n <- as.numeric(truth)
      pred_n  <- as.numeric(pred)
      rmse <- sqrt(mean((truth_n - pred_n)^2, na.rm = TRUE))
      base_rmse <- sqrt(mean(
        (truth_n - mean(df_miss[[tr]], na.rm = TRUE))^2, na.rm = TRUE))
      out[[tr]] <- data.frame(
        dataset = dataset_name, seed = seed, config = config, trait = tr,
        type = "continuous_or_count",
        metric = "rmse",
        value  = rmse, baseline = base_rmse, n = length(idx),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

# ---- Sweep -----------------------------------------------------------------

CONFIGS <- c("default", "clamp", "mode", "both")

log_line("Loading PanTHERIA ...")
pan <- load_pantheria()
log_line(sprintf("  PanTHERIA full: %d species, %d traits",
                 nrow(pan$df), ncol(pan$df)))
# Subset to N_SUB random species (use first SEED for the subset draw
# so the species set is fixed across config x seed cells; only the
# masking + GNN training varies across seeds, isolating the noise we
# want to measure).
if (nrow(pan$df) > N_SUB) {
  set.seed(SEEDS[1])
  keep <- sample(rownames(pan$df), N_SUB)
  pan$tree <- ape::drop.tip(pan$tree, setdiff(pan$tree$tip.label, keep))
  pan$df   <- pan$df[pan$tree$tip.label, , drop = FALSE]
  log_line(sprintf("  PanTHERIA subset: %d species (matches AVONET bench scale)",
                   nrow(pan$df)))
}

log_line("Loading AmphiBIO ...")
amph <- load_amphibio()
log_line(sprintf("  AmphiBIO full: %d species, %d traits",
                 nrow(amph$df), ncol(amph$df)))
if (nrow(amph$df) > N_SUB) {
  set.seed(SEEDS[1] + 1L)
  keep <- sample(rownames(amph$df), N_SUB)
  amph$tree <- ape::drop.tip(amph$tree, setdiff(amph$tree$tip.label, keep))
  amph$df   <- amph$df[amph$tree$tip.label, , drop = FALSE]
  log_line(sprintf("  AmphiBIO subset: %d species (matches AVONET bench scale)",
                   nrow(amph$df)))
}

results <- list()
total_cells <- length(SEEDS) * length(CONFIGS) * 2L
i <- 0L
for (s in SEEDS) {
  for (cfg in CONFIGS) {
    i <- i + 1L
    log_line(sprintf("==== cell %d/%d ====", i, total_cells))
    results[[length(results) + 1L]] <- run_one_cell("pantheria",
                                                     pan$df, pan$tree, cfg, s)
  }
  for (cfg in CONFIGS) {
    i <- i + 1L
    log_line(sprintf("==== cell %d/%d ====", i, total_cells))
    results[[length(results) + 1L]] <- run_one_cell("amphibio",
                                                     amph$df, amph$tree, cfg, s)
  }
}
all_rows <- do.call(rbind, results)
saveRDS(all_rows, out_rds)

# ---- Summary: per-config mean +/- SD across seeds --------------------------

# Per (dataset, trait, config), aggregate mean +/- SD across seeds, then
# compute lift_pct vs default-config mean.
agg <- stats::aggregate(value ~ dataset + trait + type + config + metric,
                        data = all_rows,
                        FUN = function(x) {
                          c(mean = mean(x, na.rm = TRUE),
                            sd   = stats::sd(x, na.rm = TRUE),
                            n    = sum(!is.na(x)))
                        })
agg <- data.frame(
  dataset = agg$dataset, trait = agg$trait, type = agg$type,
  config  = agg$config, metric = agg$metric,
  mean    = as.numeric(agg$value[, "mean"]),
  sd      = as.numeric(agg$value[, "sd"]),
  n_seeds = as.numeric(agg$value[, "n"]),
  stringsAsFactors = FALSE
)
default_agg <- subset(agg, config == "default")
agg$lift_pct_mean <- NA_real_
for (i in seq_len(nrow(agg))) {
  d <- subset(default_agg,
              dataset == agg$dataset[i] & trait == agg$trait[i])
  if (nrow(d) == 1L && d$mean != 0) {
    agg$lift_pct_mean[i] <- if (agg$metric[i] == "rmse") {
      100 * (d$mean - agg$mean[i]) / d$mean
    } else {
      100 * (agg$mean[i] - d$mean)
    }
  }
}
agg <- agg[order(agg$dataset, agg$trait, agg$config), ]

md <- c(
  "# Phase C v2 cross-dataset bench: PanTHERIA + AmphiBIO (multi-seed)",
  "",
  sprintf("Run: %s", format(Sys.time())),
  sprintf("Seeds: %s.  miss_frac = %.2f, N_IMP = %d.  AmphiBIO Diu / Noc dropped (presence-only encoding).",
          paste(SEEDS, collapse = ", "), MISS_FRAC, N_IMP),
  "",
  "## Aggregated results (mean +/- SD across seeds)",
  "",
  "```",
  capture.output(print(agg[, c("dataset", "trait", "type", "config",
                                 "metric", "mean", "sd", "lift_pct_mean")],
                        row.names = FALSE)),
  "```",
  "",
  "## Raw per-seed results",
  "",
  "```",
  capture.output(print(all_rows[order(all_rows$dataset, all_rows$trait,
                                        all_rows$config, all_rows$seed),
                                  c("dataset", "seed", "trait", "config",
                                    "metric", "value", "baseline")],
                        row.names = FALSE)),
  "```",
  ""
)
writeLines(md, out_md)
log_line("DONE -- ", out_rds)
log_line("DONE -- ", out_md)
