#!/usr/bin/env Rscript
#
# script/bench_multi_obs_mixed.R
#
# Multi-observation + mixed-type benchmark: Phase 10 + B1 soft E-step.
#
# Why this exists
#   Phase 10 removed the !multi_obs guards in fit_baseline() so the joint
#   MVN, threshold-joint, and OVR categorical paths now accept multi-obs
#   input via species-level aggregation. But bench_multi_obs.R simulates
#   only one continuous trait (CTmax), so none of the Level-C dispatch
#   conditions fire regardless of the guards. This benchmark simulates
#   mixed types (2 continuous + 1 binary + 1 categorical K=4) with
#   multiple obs per species, forcing the Level-C paths to activate.
#
# Comparison
#   Four methods per (scenario, rep):
#     - species_mean:    observed-species mean propagated (reference)
#     - pigauto_hard:    Phase 10 default (threshold-at-0.5 / argmax)
#     - pigauto_soft:    B1 soft E-step (Rao-Blackwell convex combination)
#     - pigauto_LP:      same call with joint_mvn_available forced FALSE
#   The hard vs LP delta measures Phase 10's lift; soft vs hard measures B1.
#
# Run with
#   cd pigauto && Rscript script/bench_multi_obs_mixed.R

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/b1-soft-liability",
    quiet = TRUE
  )
})

here    <- "/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/b1-soft-liability"
out_rds <- file.path(here, "script", "bench_multi_obs_mixed.rds")
out_md  <- file.path(here, "script", "bench_multi_obs_mixed.md")

script_start <- proc.time()[["elapsed"]]
log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start), ...,
      "\n", sep = "")
  flush.console()
}

# -------------------------------------------------------------------------
# Constants
# -------------------------------------------------------------------------

n_species         <- 150L
reps              <- 2L
epochs            <- 100L
obs_lambda_pois   <- 5L
obs_min           <- 1L
obs_max           <- 15L
sp_missing_frac   <- 0.30
within_miss_frac  <- 0.20
within_noise_cont <- 0.30   # SD of within-species noise on x1, x2
within_flip_bin   <- 0.10   # probability a binary obs disagrees with species
within_flip_cat   <- 0.10   # probability a categorical obs is NOT the species mode

scenarios <- list(
  high_phylo_correlated = list(lambda = 1.0, rho = 0.7),
  mod_phylo_correlated  = list(lambda = 0.5, rho = 0.6),
  low_phylo_correlated  = list(lambda = 0.2, rho = 0.6)
)

# -------------------------------------------------------------------------
# Helpers (reuse bench_discriminative_phase9.R logic)
# -------------------------------------------------------------------------

tree_with_lambda <- function(tree, lambda) {
  if (abs(lambda - 1) < 1e-6) return(tree)
  t2 <- tree
  is_internal <- t2$edge[, 2] > length(t2$tip.label)
  t2$edge.length[is_internal] <- t2$edge.length[is_internal] * lambda
  t2
}

simulate_species_mixed <- function(tree, rho, seed) {
  set.seed(seed)
  n <- length(tree$tip.label)
  V <- ape::vcv(tree)
  L <- chol(V)
  z_x1 <- as.numeric(t(L) %*% rnorm(n))
  z_x2 <- rho * z_x1 + sqrt(1 - rho^2) * as.numeric(t(L) %*% rnorm(n))
  z_y  <- rho * z_x1 + sqrt(1 - rho^2) * as.numeric(t(L) %*% rnorm(n))
  z_liab <- matrix(0, nrow = n, ncol = 4)
  for (k in seq_len(4)) {
    z_k <- as.numeric(t(L) %*% rnorm(n))
    z_liab[, k] <- rho * z_x1 + sqrt(1 - rho^2) * z_k
  }
  z_class <- apply(z_liab, 1, which.max)
  data.frame(
    x1 = z_x1, x2 = z_x2,
    y  = factor(ifelse(z_y > 0, "B", "A"), levels = c("A", "B")),
    z  = factor(LETTERS[z_class], levels = LETTERS[1:4]),
    row.names = tree$tip.label,
    stringsAsFactors = FALSE
  )
}

simulate_multi_obs_mixed <- function(tree, rho, seed) {
  sp_df <- simulate_species_mixed(tree, rho = rho, seed = seed)
  set.seed(seed + 1e4L)
  sp_vec <- rownames(sp_df)
  n_sp   <- length(sp_vec)
  n_obs_per <- pmin(pmax(rpois(n_sp, obs_lambda_pois), obs_min), obs_max)
  names(n_obs_per) <- sp_vec

  species_vec <- rep(sp_vec, n_obs_per)
  n_total     <- length(species_vec)

  # Species-level values repeated per obs
  x1 <- rep(sp_df$x1, n_obs_per) + rnorm(n_total, 0, within_noise_cont)
  x2 <- rep(sp_df$x2, n_obs_per) + rnorm(n_total, 0, within_noise_cont)

  # Binary: flip some obs
  y_sp  <- rep(as.character(sp_df$y), n_obs_per)
  flips <- runif(n_total) < within_flip_bin
  y_obs <- y_sp
  y_obs[flips] <- ifelse(y_sp[flips] == "A", "B", "A")
  y <- factor(y_obs, levels = c("A", "B"))

  # Categorical: some obs disagree with the species modal class
  z_sp_char <- rep(as.character(sp_df$z), n_obs_per)
  flips_z   <- runif(n_total) < within_flip_cat
  z_alt <- sapply(z_sp_char[flips_z], function(curr) {
    others <- setdiff(LETTERS[1:4], curr)
    sample(others, 1L)
  })
  z_obs <- z_sp_char
  z_obs[flips_z] <- z_alt
  z <- factor(z_obs, levels = LETTERS[1:4])

  data.frame(
    species = species_vec,
    x1 = x1, x2 = x2, y = y, z = z,
    stringsAsFactors = FALSE
  )
}

# Mask entire species + within-species obs (returns df with NAs)
mask_traits <- function(df, sp_miss_frac, within_miss_frac, seed) {
  set.seed(seed + 999L)
  sp_vec <- unique(df$species)
  n_sp   <- length(sp_vec)
  n_miss <- floor(sp_miss_frac * n_sp)
  missing_sp <- sample(sp_vec, n_miss)

  df_masked <- df
  trait_cols <- c("x1", "x2", "y", "z")

  # Full species mask (all obs NA for trait cols)
  for (col in trait_cols) {
    df_masked[df_masked$species %in% missing_sp, col] <- NA
  }

  # Within-species random obs mask on OBSERVED species
  obs_rows <- which(!(df$species %in% missing_sp))
  for (col in trait_cols) {
    pick <- sample(obs_rows,
                   size = round(within_miss_frac * length(obs_rows)))
    df_masked[pick, col] <- NA
  }

  df_masked
}

# -------------------------------------------------------------------------
# LP-forced helper (from bench_discriminative_phase9.R)
# -------------------------------------------------------------------------

with_lp_forced <- function(fn) {
  orig <- get("joint_mvn_available", envir = asNamespace("pigauto"))
  assignInNamespace("joint_mvn_available", function() FALSE, ns = "pigauto")
  on.exit(
    assignInNamespace("joint_mvn_available", orig, ns = "pigauto"),
    add = TRUE
  )
  fn()
}

# -------------------------------------------------------------------------
# Evaluation: obs-level metrics on held-out cells
# -------------------------------------------------------------------------

eval_obs_level <- function(df_truth, df_masked, completed, imputed_mask) {
  # completed has the imputed values. imputed_mask is TRUE where cells were
  # filled. We evaluate ONLY on cells that were originally observed in
  # df_truth AND masked in df_masked (the held-out cells).
  metrics <- list()
  trait_cols <- c("x1", "x2", "y", "z")

  for (col in trait_cols) {
    held_out <- !is.na(df_truth[[col]]) & is.na(df_masked[[col]])
    if (!any(held_out)) next
    truth <- df_truth[[col]][held_out]
    pred  <- completed[[col]][held_out]

    if (col %in% c("x1", "x2")) {
      metrics[[col]] <- list(
        metric = "RMSE",
        value  = sqrt(mean((as.numeric(truth) - as.numeric(pred))^2,
                             na.rm = TRUE))
      )
    } else if (col == "y") {
      metrics[[col]] <- list(
        metric = "accuracy",
        value  = mean(as.character(truth) == as.character(pred), na.rm = TRUE)
      )
    } else if (col == "z") {
      metrics[[col]] <- list(
        metric = "accuracy",
        value  = mean(as.character(truth) == as.character(pred), na.rm = TRUE)
      )
    }
  }
  metrics
}

# Species-mean baseline: for each trait, fill missing obs of each species
# with the species mean (continuous) or species mode (discrete) of the
# OBSERVED obs of that species. Species with no observed obs stay NA.
species_mean_impute <- function(df_masked) {
  trait_cols <- c("x1", "x2", "y", "z")
  df_out     <- df_masked
  sp_vec     <- unique(df_masked$species)

  for (col in trait_cols) {
    vals <- df_masked[[col]]
    out  <- vals
    for (s in sp_vec) {
      idx <- which(df_masked$species == s)
      obs <- vals[idx]
      obs_non_na <- obs[!is.na(obs)]
      if (length(obs_non_na) == 0L) next
      if (col %in% c("x1", "x2")) {
        fill <- mean(obs_non_na, na.rm = TRUE)
      } else {
        tab <- table(obs_non_na)
        fill <- names(tab)[which.max(tab)]
        # Preserve factor levels
        if (is.factor(vals)) {
          fill <- factor(fill, levels = levels(vals))
        }
      }
      na_idx <- idx[is.na(out[idx])]
      out[na_idx] <- fill
    }
    df_out[[col]] <- out
  }
  df_out
}

# -------------------------------------------------------------------------
# Per-cell driver
# -------------------------------------------------------------------------

run_one <- function(sc_name, sc, rep_id) {
  seed <- rep_id * 1000L + nchar(sc_name)
  log_line("=== ", sc_name, " (lambda=", sc$lambda, ", rho=", sc$rho,
           ") rep ", rep_id, " (seed=", seed, ") ===")
  tree0 <- ape::rtree(n_species)
  tree  <- tree_with_lambda(tree0, sc$lambda)

  df_truth  <- simulate_multi_obs_mixed(tree, rho = sc$rho, seed = seed)
  df_masked <- mask_traits(df_truth, sp_miss_frac = sp_missing_frac,
                           within_miss_frac = within_miss_frac, seed = seed)
  n_total <- nrow(df_truth)
  log_line("  n_obs = ", n_total, ", n_species = ", length(unique(df_truth$species)))

  results <- list()

  # ----- Method 1: species_mean -----
  t0 <- proc.time()[["elapsed"]]
  df_sm <- species_mean_impute(df_masked)
  wall_sm <- proc.time()[["elapsed"]] - t0
  metrics_sm <- eval_obs_level(df_truth, df_masked, df_sm, NULL)
  log_line("  [1/4] species_mean done in ", round(wall_sm, 1), "s")

  # ----- Method 2: pigauto hard (Phase 10 default) -----
  t0 <- proc.time()[["elapsed"]]
  res_hard <- tryCatch(
    impute(df_masked, tree, species_col = "species",
           epochs = as.integer(epochs), verbose = FALSE, seed = seed,
           multi_obs_aggregation = "hard"),
    error = function(e) { log_line("  hard FAILED: ", conditionMessage(e)); NULL }
  )
  wall_hard <- proc.time()[["elapsed"]] - t0
  metrics_hard <- if (is.null(res_hard)) list() else
    eval_obs_level(df_truth, df_masked, res_hard$completed, res_hard$imputed_mask)
  log_line("  [2/4] pigauto_hard done in ", round(wall_hard, 1), "s")

  # ----- Method 3: pigauto soft (B1) -----
  t0 <- proc.time()[["elapsed"]]
  res_soft <- tryCatch(
    impute(df_masked, tree, species_col = "species",
           epochs = as.integer(epochs), verbose = FALSE, seed = seed,
           multi_obs_aggregation = "soft"),
    error = function(e) { log_line("  soft FAILED: ", conditionMessage(e)); NULL }
  )
  wall_soft <- proc.time()[["elapsed"]] - t0
  metrics_soft <- if (is.null(res_soft)) list() else
    eval_obs_level(df_truth, df_masked, res_soft$completed, res_soft$imputed_mask)
  log_line("  [3/4] pigauto_soft done in ", round(wall_soft, 1), "s")

  # ----- Method 4: pigauto_LP -----
  t0 <- proc.time()[["elapsed"]]
  res_lp <- tryCatch(
    with_lp_forced(function()
      impute(df_masked, tree, species_col = "species",
             epochs = as.integer(epochs), verbose = FALSE, seed = seed)),
    error = function(e) { log_line("  LP FAILED: ", conditionMessage(e)); NULL }
  )
  wall_lp <- proc.time()[["elapsed"]] - t0
  metrics_lp <- if (is.null(res_lp)) list() else
    eval_obs_level(df_truth, df_masked, res_lp$completed, res_lp$imputed_mask)
  log_line("  [4/4] pigauto_LP done in ", round(wall_lp, 1), "s")

  # Pack results
  pack <- function(method, metrics, wall) {
    if (length(metrics) == 0L) {
      return(data.frame(scenario = sc_name, rep = rep_id, method = method,
                         trait = NA_character_, metric = NA_character_,
                         value = NA_real_, wall_s = wall,
                         stringsAsFactors = FALSE))
    }
    do.call(rbind, lapply(names(metrics), function(tr) {
      m <- metrics[[tr]]
      data.frame(scenario = sc_name, rep = rep_id, method = method,
                  trait = tr, metric = m$metric, value = m$value,
                  wall_s = wall, stringsAsFactors = FALSE)
    }))
  }

  rbind(
    pack("species_mean",   metrics_sm,   wall_sm),
    pack("pigauto_hard",   metrics_hard, wall_hard),
    pack("pigauto_soft",   metrics_soft, wall_soft),
    pack("pigauto_LP",     metrics_lp,   wall_lp)
  )
}

# -------------------------------------------------------------------------
# Sweep
# -------------------------------------------------------------------------

all_results <- list()
for (sc_name in names(scenarios)) {
  for (r in seq_len(reps)) {
    res <- tryCatch(run_one(sc_name, scenarios[[sc_name]], r),
                     error = function(e) {
                       log_line("FAILED: ", conditionMessage(e)); NULL
                     })
    if (!is.null(res)) all_results[[length(all_results) + 1L]] <- res
  }
}
all_results <- do.call(rbind, all_results)
total_wall <- proc.time()[["elapsed"]] - script_start

# Aggregate: mean over reps, one row per (scenario, method, trait)
agg <- aggregate(value ~ scenario + method + trait + metric,
                  data = all_results, FUN = mean, na.rm = TRUE)

# Wide pivot: method as columns, per (scenario, trait, metric)
wide_rows <- list()
for (sc in unique(agg$scenario)) {
  for (tr in unique(agg$trait[!is.na(agg$trait)])) {
    sub <- agg[agg$scenario == sc & agg$trait == tr, , drop = FALSE]
    if (nrow(sub) == 0L) next
    row <- list(scenario = sc, trait = tr,
                metric = sub$metric[1])
    for (m in c("species_mean", "pigauto_hard", "pigauto_soft", "pigauto_LP")) {
      v <- sub$value[sub$method == m]
      row[[m]] <- if (length(v) == 0L) NA_real_ else v
    }
    row$hard_vs_LP <- if (row$metric == "RMSE") {
      row$pigauto_LP - row$pigauto_hard   # lower = better
    } else {
      row$pigauto_hard - row$pigauto_LP   # higher = better
    }
    row$soft_vs_hard <- if (row$metric == "RMSE") {
      row$pigauto_hard - row$pigauto_soft   # lower = better
    } else {
      row$pigauto_soft - row$pigauto_hard   # higher = better
    }
    wide_rows[[length(wide_rows) + 1L]] <- as.data.frame(row,
                                                          stringsAsFactors = FALSE)
  }
}
wide <- do.call(rbind, wide_rows)

# -------------------------------------------------------------------------
# Save
# -------------------------------------------------------------------------

saveRDS(list(
  results = all_results,
  wide    = wide,
  agg     = agg,
  meta = list(
    n_species = n_species,
    reps      = reps,
    epochs    = epochs,
    scenarios = scenarios,
    total_wall = total_wall
  )
), out_rds)

# -------------------------------------------------------------------------
# Markdown summary
# -------------------------------------------------------------------------

md <- c(
  "# Multi-obs + mixed-type benchmark (Phase 10 + B1 soft E-step)",
  "",
  sprintf("Run on: %s", format(Sys.time())),
  sprintf("Species per scenario: %d, reps: %d, epochs: %d",
          n_species, reps, epochs),
  sprintf("Total wall time: %.1f min", total_wall / 60),
  "",
  "Four methods, obs-level metrics on held-out cells:",
  "- **species_mean**: observed-species mean (continuous) or mode (discrete); reference.",
  "- **pigauto_hard**: `impute(..., multi_obs_aggregation = 'hard')` — Phase 10 default (threshold-at-0.5 / argmax).",
  "- **pigauto_soft**: `impute(..., multi_obs_aggregation = 'soft')` — B1 Rao-Blackwell soft E-step.",
  "- **pigauto_LP**: same call with `joint_mvn_available()` forced FALSE (legacy LP path).",
  "",
  "hard_vs_LP: `pigauto_LP - pigauto_hard` for RMSE, `pigauto_hard - pigauto_LP` for accuracy (positive = hard better).",
  "soft_vs_hard: `pigauto_hard - pigauto_soft` for RMSE, `pigauto_soft - pigauto_hard` for accuracy (positive = soft better).",
  "",
  "```",
  capture.output(print(wide, row.names = FALSE, digits = 4)),
  "```",
  "",
  sprintf("Phase 10 exit criterion (>=2pp lift on >=1 trait in >=2 scenarios): **%s**",
          {
            lift_thresh <- 0.02
            n_qualifying_scenarios <- length(unique(
              wide$scenario[!is.na(wide$hard_vs_LP) &
                            abs(wide$hard_vs_LP) >= lift_thresh &
                            sign(wide$hard_vs_LP) == 1]
            ))
            if (n_qualifying_scenarios >= 2L) "MET" else "NOT MET"
          }),
  "",
  sprintf("B1 exit criterion (soft >= hard on >=1 trait in >=2 scenarios): **%s**",
          {
            n_soft_better <- length(unique(
              wide$scenario[!is.na(wide$soft_vs_hard) &
                            wide$soft_vs_hard > 0]
            ))
            if (n_soft_better >= 2L) "MET" else "NOT MET"
          })
)
writeLines(md, out_md)

log_line("Total wall: ", round(total_wall / 60, 1), " min")
log_line("Wrote ", out_rds)
log_line("Wrote ", out_md)
log_line("done")
