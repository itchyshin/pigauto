#' Run a simulation benchmark for pigauto
#'
#' Generates trait data under various evolutionary models, introduces
#' missing data, fits both the Brownian motion baseline and the full
#' pigauto GNN, and compares performance.  This is the recommended way
#' to assess pigauto on data with known properties before applying it
#' to real data.
#'
#' @details
#' **Available scenarios:**
#' \describe{
#'   \item{\code{"BM"}}{Pure Brownian motion — the baseline is exact,
#'     so the GNN should tie or slightly improve via inter-trait
#'     correlations.}
#'   \item{\code{"OU"}}{Ornstein-Uhlenbeck — stabilising selection
#'     constrains variation.  BM over-estimates evolutionary variance.}
#'   \item{\code{"regime_shift"}}{Two-regime BM — clade-specific optima
#'     create bimodal distributions that BM cannot capture.}
#'   \item{\code{"nonlinear"}}{Non-linear inter-trait relationships —
#'     the GNN's multi-layer message passing can capture quadratic and
#'     interaction effects that BM's linear covariance misses.}
#'   \item{\code{"mixed"}}{Mixed trait types: 2 continuous + 1 binary +
#'     1 categorical (3 levels).  Tests the full type pipeline.}
#' }
#'
#' @param n_species integer.  Number of tips in the simulated tree
#'   (default 100).
#' @param n_traits integer.  Number of continuous traits (default 4).
#'   Ignored for \code{scenario = "mixed"}, which generates a fixed
#'   trait set.
#' @param scenarios character vector.  Subset of
#'   \code{c("BM", "OU", "regime_shift", "nonlinear", "mixed")}.
#'   Default runs all.
#' @param missing_frac numeric.  Fraction of observed cells held out
#'   (default 0.25).
#' @param n_reps integer.  Number of replicate trees per scenario
#'   (default 3).
#' @param epochs integer.  Maximum GNN training epochs (default 500).
#' @param verbose logical.  Print progress (default \code{TRUE}).
#' @param ... additional arguments passed to \code{\link{fit_pigauto}}.
#' @return An object of class \code{"pigauto_benchmark"} with:
#'   \describe{
#'     \item{results}{data.frame with columns: \code{scenario},
#'       \code{rep}, \code{method}, \code{trait}, \code{type},
#'       \code{metric}, \code{value}, \code{n_test}.}
#'     \item{summary}{data.frame averaged across replicates.}
#'     \item{scenarios}{character vector of scenarios run.}
#'     \item{n_reps}{integer.}
#'     \item{n_species}{integer.}
#'   }
#' @examples
#' \dontrun{
#' bench <- simulate_benchmark(n_species = 50, epochs = 200, n_reps = 2)
#' bench$summary
#' plot(bench)
#' }
#' @export
simulate_benchmark <- function(
    n_species  = 100L,
    n_traits   = 4L,
    scenarios  = c("BM", "OU", "regime_shift", "nonlinear", "mixed"),
    missing_frac = 0.25,
    n_reps     = 3L,
    epochs     = 500L,
    verbose    = TRUE,
    ...
) {
  all_scenarios <- c("BM", "OU", "regime_shift", "nonlinear", "mixed")
  scenarios <- match.arg(scenarios, all_scenarios, several.ok = TRUE)

  all_results <- list()
  idx <- 0L

  for (scen in scenarios) {
    if (verbose) message("\n=== Scenario: ", scen, " ===")

    for (rep in seq_len(n_reps)) {
      rep_seed <- rep * 100L + match(scen, all_scenarios)
      if (verbose) message("  Rep ", rep, "/", n_reps,
                           " (seed ", rep_seed, ")")

      # ---- Generate tree and traits ------------------------------------------
      set.seed(rep_seed)
      tree <- ape::rtree(as.integer(n_species))

      if (scen == "mixed") {
        df <- simulate_mixed_traits(tree, seed = rep_seed)
      } else if (scen == "BM") {
        df <- simulate_bm_traits(tree, n_traits, seed = rep_seed)
      } else {
        df <- simulate_non_bm(tree, n_traits = n_traits,
                              scenario = scen, seed = rep_seed)
      }

      # ---- Pipeline ----------------------------------------------------------
      result <- tryCatch({
        pd  <- preprocess_traits(df, tree, log_transform = FALSE)
        spl <- make_missing_splits(pd$X_scaled, missing_frac = missing_frac,
                                   seed = rep_seed, trait_map = pd$trait_map)
        bl  <- fit_baseline(pd, tree, splits = spl)
        fit <- fit_pigauto(pd, tree, splits = spl, baseline = bl,
                           epochs = as.integer(epochs),
                           verbose = FALSE, seed = rep_seed, ...)
        ev <- evaluate(fit, data = pd, splits = spl)
        ev$scenario <- scen
        ev$rep <- rep
        ev
      }, error = function(e) {
        if (verbose) message("    ERROR: ", conditionMessage(e))
        NULL
      })

      if (!is.null(result)) {
        idx <- idx + 1L
        all_results[[idx]] <- result
      }
    }
  }

  results <- do.call(rbind, all_results)
  rownames(results) <- NULL

  # Summarise across replicates
  summary_df <- do.call(rbind, lapply(
    split(results, paste(results$scenario, results$method,
                         results$trait, results$metric)),
    function(d) {
      data.frame(
        scenario = d$scenario[1],
        method   = d$method[1],
        trait    = d$trait[1],
        type     = d$type[1],
        metric   = d$metric[1],
        mean     = mean(d$value, na.rm = TRUE),
        sd       = stats::sd(d$value, na.rm = TRUE),
        n_reps   = sum(!is.na(d$value)),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(summary_df) <- NULL

  structure(
    list(
      results    = results,
      summary    = summary_df,
      scenarios  = scenarios,
      n_reps     = n_reps,
      n_species  = n_species
    ),
    class = "pigauto_benchmark"
  )
}


# ---- Print / summary methods -----------------------------------------------

#' @export
print.pigauto_benchmark <- function(x, ...) {
  cat("pigauto simulation benchmark\n")
  cat(sprintf("  Scenarios: %s\n", paste(x$scenarios, collapse = ", ")))
  cat(sprintf("  Species: %d | Replicates: %d\n", x$n_species, x$n_reps))
  cat(sprintf("  Total evaluations: %d rows\n", nrow(x$results)))

  # Quick comparison: BM vs pigauto mean RMSE
  rmse_rows <- x$summary[x$summary$metric == "rmse", ]
  if (nrow(rmse_rows) > 0L) {
    cat("\n  Mean RMSE by method:\n")
    for (m in unique(rmse_rows$method)) {
      sub <- rmse_rows[rmse_rows$method == m, ]
      cat(sprintf("    %s: %.4f\n", m, mean(sub$mean, na.rm = TRUE)))
    }
  }
  invisible(x)
}


#' @export
summary.pigauto_benchmark <- function(object, ...) {
  cat("pigauto simulation benchmark summary\n")
  cat(strrep("=", 60), "\n")

  for (scen in object$scenarios) {
    sub <- object$summary[object$summary$scenario == scen, ]
    cat("\nScenario: ", scen, "\n")
    cat(strrep("-", 50), "\n")

    for (m in c("baseline", "pigauto")) {
      ms <- sub[sub$method == m, ]
      if (nrow(ms) == 0L) next

      # RMSE for continuous traits
      rmse_rows <- ms[ms$metric == "rmse" & ms$type == "continuous", ]
      if (nrow(rmse_rows) > 0L) {
        cat(sprintf("  %8s  RMSE: %.4f (+/- %.4f)\n", m,
                    mean(rmse_rows$mean), mean(rmse_rows$sd)))
      }

      # Accuracy for discrete traits
      acc_rows <- ms[ms$metric == "accuracy" &
                       ms$type %in% c("binary", "categorical"), ]
      if (nrow(acc_rows) > 0L) {
        cat(sprintf("  %8s  Accuracy: %.1f%% (+/- %.1f%%)\n", m,
                    100 * mean(acc_rows$mean), 100 * mean(acc_rows$sd)))
      }
    }
  }

  # Relative improvement
  cat("\n", strrep("=", 60), "\n")
  cat("Relative improvement (pigauto vs baseline):\n")
  cat(strrep("-", 50), "\n")

  for (scen in object$scenarios) {
    sub <- object$summary[object$summary$scenario == scen, ]
    bl  <- sub[sub$method == "baseline" & sub$metric == "rmse", ]
    pa  <- sub[sub$method == "pigauto"  & sub$metric == "rmse", ]
    if (nrow(bl) > 0L && nrow(pa) > 0L) {
      bl_mean <- mean(bl$mean, na.rm = TRUE)
      pa_mean <- mean(pa$mean, na.rm = TRUE)
      pct <- 100 * (bl_mean - pa_mean) / bl_mean
      cat(sprintf("  %18s: %+.1f%% RMSE\n", scen, pct))
    }
  }
  cat("\n")
  invisible(object$summary)
}


#' Plot a pigauto benchmark
#'
#' Creates a multi-panel comparison plot showing BM baseline vs pigauto
#' performance across all simulation scenarios.
#'
#' @param x pigauto_benchmark object.
#' @param metric character. Which metric to plot (default \code{"rmse"}).
#' @param ... passed to base plot functions.
#' @return Invisible \code{NULL}.
#' @importFrom graphics par barplot axis text legend mtext
#' @importFrom grDevices adjustcolor
#' @export
plot.pigauto_benchmark <- function(x, metric = "rmse", ...) {
  metric <- match.arg(metric, c("rmse", "pearson_r", "mae", "accuracy",
                                "spearman_rho", "brier"))

  sub <- x$summary[x$summary$metric == metric, ]
  if (nrow(sub) == 0L) {
    message("No results for metric '", metric, "'")
    return(invisible(NULL))
  }

  scenarios <- unique(sub$scenario)
  n_scen <- length(scenarios)
  methods  <- c("baseline", "pigauto")

  opar <- par(mfrow = c(1, n_scen), mar = c(5, 4, 3, 1), oma = c(0, 0, 2, 0))
  on.exit(par(opar), add = TRUE)

  col_bl <- "#636363"
  col_pa <- "#2b8cbe"

  for (scen in scenarios) {
    ss <- sub[sub$scenario == scen, ]
    traits <- unique(ss$trait)

    vals_bl <- vapply(traits, function(t) {
      v <- ss$mean[ss$trait == t & ss$method == "baseline"]
      if (length(v) == 0L) NA_real_ else v[1]
    }, numeric(1))

    vals_pa <- vapply(traits, function(t) {
      v <- ss$mean[ss$trait == t & ss$method == "pigauto"]
      if (length(v) == 0L) NA_real_ else v[1]
    }, numeric(1))

    # Remove traits with only NA
    ok <- !is.na(vals_bl) | !is.na(vals_pa)
    if (!any(ok)) { plot.new(); next }
    traits <- traits[ok]; vals_bl <- vals_bl[ok]; vals_pa <- vals_pa[ok]

    mat <- rbind(vals_bl, vals_pa)
    colnames(mat) <- traits

    barplot(mat, beside = TRUE, col = c(col_bl, col_pa),
            main = scen, ylab = toupper(metric),
            las = 2, cex.names = 0.8, ...)

    if (scen == scenarios[1]) {
      legend("topright", legend = c("Baseline", "pigauto"),
             fill = c(col_bl, col_pa), bty = "n", cex = 0.8)
    }
  }

  mtext("pigauto Simulation Benchmark", outer = TRUE, line = 0.5,
        cex = 1.2, font = 2)
  invisible(NULL)
}


# ---- Internal: simulate pure BM traits ----------------------------------------

simulate_bm_traits <- function(tree, n_traits, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sp <- tree$tip.label
  df <- data.frame(row.names = sp)
  for (j in seq_len(n_traits)) {
    trait <- ape::rTraitCont(tree, model = "BM", sigma = 1.0,
                             root.value = j * 0.5)
    df[[paste0("trait", j)]] <- as.numeric(trait[sp])
  }
  df
}


# ---- Internal: simulate mixed-type traits ------------------------------------

simulate_mixed_traits <- function(tree, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sp <- tree$tip.label
  n  <- length(sp)

  # 2 continuous traits (BM)
  t1 <- ape::rTraitCont(tree, model = "BM", sigma = 1.0, root.value = 0)
  t2 <- ape::rTraitCont(tree, model = "BM", sigma = 1.0, root.value = 1)

  # 1 binary trait (threshold on a latent BM)
  latent_bin <- ape::rTraitCont(tree, model = "BM", sigma = 1.0,
                                root.value = 0)
  bin_vals <- ifelse(as.numeric(latent_bin[sp]) > 0, "yes", "no")

  # 1 categorical trait (3 levels from thresholds on two latent BMs)
  latent_cat1 <- ape::rTraitCont(tree, model = "BM", sigma = 1.0,
                                 root.value = 0)
  latent_cat2 <- ape::rTraitCont(tree, model = "BM", sigma = 1.0,
                                 root.value = 0)
  lc1 <- as.numeric(latent_cat1[sp])
  lc2 <- as.numeric(latent_cat2[sp])
  cat_vals <- ifelse(lc1 > 0.5, "A",
                     ifelse(lc2 > 0, "B", "C"))

  data.frame(
    row.names = sp,
    cont1 = abs(as.numeric(t1[sp])) + 0.1,
    cont2 = abs(as.numeric(t2[sp])) + 0.1,
    binary_trait = factor(bin_vals, levels = c("no", "yes")),
    cat_trait    = factor(cat_vals, levels = c("A", "B", "C"))
  )
}
