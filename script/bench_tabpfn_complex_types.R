#!/usr/bin/env Rscript
#
# script/bench_tabpfn_complex_types.R
#
# Branch-local TabPFN benchmark for pigauto's non-scalar families:
# zero-inflated counts and multi-proportion compositions.
#
# ZI-count is fitted as a hybrid model: TabPFNClassifier for the non-zero gate
# and TabPFNRegressor for the conditional non-zero magnitude.
# Multi-proportion is fitted as K independent CLR regressions and evaluated
# through pigauto's existing compositional metrics.
#
# Smoke:
#   PIGAUTO_TABPFN_OUT_STEM=bench_tabpfn_complex_types_smoke \
#   PIGAUTO_TABPFN_SCALES=100 PIGAUTO_TABPFN_REPS=1 \
#   PIGAUTO_TABPFN_SCENARIOS=zi_moderate,multi_moderate \
#   PIGAUTO_TABPFN_VARIANTS=lappe \
#   PIGAUTO_TABPFN_RUN_PIGAUTO=false \
#     Rscript script/bench_tabpfn_complex_types.R

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(".", quiet = TRUE)
})

script_start <- proc.time()[["elapsed"]]

log_line <- function(...) {
  cat(
    sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start),
    ...,
    "\n",
    sep = ""
  )
  flush.console()
}

env_chr <- function(name, default) {
  x <- Sys.getenv(name, unset = NA_character_)
  if (is.na(x) || !nzchar(x)) default else x
}

env_int <- function(name, default) {
  as.integer(env_chr(name, as.character(default)))
}

env_num <- function(name, default) {
  as.numeric(env_chr(name, as.character(default)))
}

env_lgl <- function(name, default) {
  x <- tolower(env_chr(name, if (isTRUE(default)) "true" else "false"))
  x %in% c("1", "true", "t", "yes", "y")
}

env_vec <- function(name, default) {
  x <- env_chr(name, paste(default, collapse = ","))
  trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
}

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_stem <- env_chr("PIGAUTO_TABPFN_OUT_STEM", "bench_tabpfn_complex_types")
if (!grepl("^[A-Za-z0-9_.-]+$", out_stem)) {
  stop("PIGAUTO_TABPFN_OUT_STEM must be a simple file stem, not a path.")
}
out_rds <- file.path(repo_root, "script", paste0(out_stem, ".rds"))
out_md <- file.path(repo_root, "script", paste0(out_stem, ".md"))
runner_py <- file.path(repo_root, "script", "run_tabpfn_phylo.py")

scales <- as.integer(env_vec("PIGAUTO_TABPFN_SCALES", c("150")))
scenarios <- env_vec(
  "PIGAUTO_TABPFN_SCENARIOS",
  c(
    "zi_moderate",
    "zi_sparse",
    "multi_moderate",
    "multi_high_phylo",
    "multi_K8"
  )
)
variants <- env_vec("PIGAUTO_TABPFN_VARIANTS", c("plain", "lappe", "lappe_nfa"))
n_reps <- env_int("PIGAUTO_TABPFN_REPS", 1L)
seed0 <- env_int("PIGAUTO_TABPFN_SEED", 20260611L)
missing_frac <- env_num("PIGAUTO_TABPFN_MISSING_FRAC", 0.25)
val_frac <- env_num("PIGAUTO_TABPFN_VAL_FRAC", 0.5)
epochs <- env_int("PIGAUTO_TABPFN_EPOCHS", 200L)
nfa_k <- env_int("PIGAUTO_TABPFN_NFA_K", 15L)
dry_run <- env_lgl("PIGAUTO_TABPFN_DRY_RUN", FALSE)
run_pigauto <- env_lgl("PIGAUTO_TABPFN_RUN_PIGAUTO", !dry_run)
python <- env_chr("PIGAUTO_TABPFN_PYTHON", "python3")
device <- env_chr("PIGAUTO_TABPFN_DEVICE", "auto")

valid_variants <- c("plain", "lappe", "lappe_nfa", "knn")
if (!all(variants %in% valid_variants)) {
  stop(
    "Unknown variants: ",
    paste(setdiff(variants, valid_variants), collapse = ", ")
  )
}

valid_scenarios <- c(
  "zi_moderate",
  "zi_sparse",
  "zi_many_zeros",
  "multi_moderate",
  "multi_high_phylo",
  "multi_K8"
)
if (!all(scenarios %in% valid_scenarios)) {
  stop(
    "Unknown scenarios: ",
    paste(setdiff(scenarios, valid_scenarios), collapse = ", ")
  )
}

metadata <- list(
  benchmark = "tabpfn_complex_types",
  scales = scales,
  scenarios = scenarios,
  variants = variants,
  n_reps = n_reps,
  seed0 = seed0,
  missing_frac = missing_frac,
  val_frac = val_frac,
  epochs = epochs,
  nfa_k = nfa_k,
  dry_run = dry_run,
  run_pigauto = run_pigauto,
  python = python,
  device = device,
  out_stem = out_stem,
  commit = tryCatch(
    system("git rev-parse HEAD", intern = TRUE),
    error = function(e) "unknown"
  )
)

logit <- function(p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  log(p / (1 - p))
}

scenario_family <- function(scenario) {
  if (startsWith(scenario, "zi_")) "zi_count" else "multi_proportion"
}

scenario_params <- function(scenario) {
  switch(
    scenario,
    zi_moderate = list(zero_frac = 0.5, mean_nz = 20, overdispersion = NULL),
    zi_sparse = list(zero_frac = 0.5, mean_nz = 5, overdispersion = 2),
    zi_many_zeros = list(zero_frac = 0.8, mean_nz = 20, overdispersion = NULL),
    multi_moderate = list(K = 5L, signal = 0.6),
    multi_high_phylo = list(K = 5L, signal = 0.9),
    multi_K8 = list(K = 8L, signal = 0.6),
    stop("Unknown scenario: ", scenario)
  )
}

generate_data <- function(n, scenario, seed) {
  pars <- scenario_params(scenario)
  set.seed(seed)
  tree <- ape::rtree(n)
  family <- scenario_family(scenario)

  if (family == "zi_count") {
    df <- simulate_zi_count_traits(
      tree,
      n_traits = 2L,
      zero_frac = pars$zero_frac,
      mean_nz = pars$mean_nz,
      overdispersion = pars$overdispersion,
      seed = seed + 1L
    )
    list(
      df = df,
      tree = tree,
      trait_types = setNames(rep("zi_count", ncol(df)), names(df)),
      multi_proportion_groups = NULL,
      family = family
    )
  } else {
    df <- simulate_multi_proportion_traits(
      tree,
      K = pars$K,
      signal = pars$signal,
      seed = seed + 1L
    )
    list(
      df = df,
      tree = tree,
      trait_types = NULL,
      multi_proportion_groups = list(comp = names(df)),
      family = family
    )
  }
}

target_rows_from_split <- function(splits, target_col, n) {
  split_to_rows <- function(idx) {
    row_i <- ((idx - 1L) %% n) + 1L
    col_j <- ((idx - 1L) %/% n) + 1L
    row_i[col_j == target_col]
  }
  list(
    val = split_to_rows(splits$val_idx),
    test = split_to_rows(splits$test_idx)
  )
}

mean_impute <- function(pd, splits) {
  X <- pd$X_scaled
  X_train <- X
  X_train[c(splits$val_idx, splits$test_idx)] <- NA_real_
  means <- colMeans(X_train, na.rm = TRUE)
  means[!is.finite(means)] <- 0
  matrix(
    means,
    nrow = nrow(X),
    ncol = ncol(X),
    byrow = TRUE,
    dimnames = dimnames(X)
  )
}

neighbor_regression_features <- function(y, D, train_rows, k) {
  n <- length(y)
  out <- matrix(NA_real_, nrow = n, ncol = 4L)
  colnames(out) <- c("nfa_mean", "nfa_sd", "nfa_min_dist", "nfa_n")
  train_rows <- train_rows[is.finite(y[train_rows])]
  if (length(train_rows) < 2L) {
    return(out)
  }

  for (i in seq_len(n)) {
    candidates <- setdiff(train_rows, i)
    if (!length(candidates)) {
      next
    }
    d <- D[i, candidates]
    ord <- order(d, decreasing = FALSE, na.last = NA)
    if (!length(ord)) {
      next
    }
    keep <- candidates[ord[seq_len(min(k, length(ord)))]]
    vals <- y[keep]
    out[i, "nfa_mean"] <- mean(vals, na.rm = TRUE)
    out[i, "nfa_sd"] <- if (length(vals) > 1L) {
      stats::sd(vals, na.rm = TRUE)
    } else {
      0
    }
    out[i, "nfa_min_dist"] <- min(D[i, keep], na.rm = TRUE)
    out[i, "nfa_n"] <- length(vals)
  }
  out
}

neighbor_class_features <- function(labels, D, train_rows, k) {
  n <- length(labels)
  props <- matrix(NA_real_, nrow = n, ncol = 2L)
  colnames(props) <- c("nfa_p_zero", "nfa_p_nonzero")
  extra <- matrix(NA_real_, nrow = n, ncol = 3L)
  colnames(extra) <- c("nfa_entropy", "nfa_min_dist", "nfa_n")
  train_rows <- train_rows[!is.na(labels[train_rows])]
  if (length(train_rows) < 2L) {
    return(cbind(props, extra))
  }

  for (i in seq_len(n)) {
    candidates <- setdiff(train_rows, i)
    if (!length(candidates)) {
      next
    }
    d <- D[i, candidates]
    ord <- order(d, decreasing = FALSE, na.last = NA)
    if (!length(ord)) {
      next
    }
    keep <- candidates[ord[seq_len(min(k, length(ord)))]]
    labs <- labels[keep]
    tab <- tabulate(labs + 1L, nbins = 2L)
    p <- tab / sum(tab)
    props[i, ] <- p
    extra[i, "nfa_entropy"] <- -sum(ifelse(p > 0, p * log(p), 0))
    extra[i, "nfa_min_dist"] <- min(D[i, keep], na.rm = TRUE)
    extra[i, "nfa_n"] <- length(keep)
  }
  cbind(props, extra)
}

feature_matrix <- function(
  pd,
  graph,
  splits,
  target_cols,
  variant,
  train_rows,
  target_values = NULL,
  class_labels = NULL
) {
  X <- pd$X_scaled
  X[c(splits$val_idx, splits$test_idx)] <- NA_real_
  X <- X[, setdiff(seq_len(ncol(X)), target_cols), drop = FALSE]
  parts <- list()
  if (ncol(X)) {
    parts$traits <- X
  }
  if (variant %in% c("lappe", "lappe_nfa")) {
    coords <- graph$coords
    colnames(coords) <- paste0("lappe_", seq_len(ncol(coords)))
    parts$lappe <- coords
  }
  if (variant %in% c("lappe_nfa", "knn")) {
    if (!is.null(class_labels)) {
      parts$nfa <- neighbor_class_features(
        class_labels,
        graph$D,
        train_rows,
        nfa_k
      )
    } else {
      parts$nfa <- neighbor_regression_features(
        target_values,
        graph$D,
        train_rows,
        nfa_k
      )
    }
  }
  if (!length(parts)) {
    stop("No predictors available for this target and variant.")
  }
  out <- do.call(cbind, parts)
  storage.mode(out) <- "double"
  out
}

single_trait_map <- function(pd, target) {
  setNames(list(pd$trait_map[[target]]), target)
}

run_tabpfn <- function(train_x, train_y, pred_x, task, seed) {
  tmpdir <- tempfile("tabpfn_complex_types_")
  dir.create(tmpdir, recursive = TRUE)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  train_path <- file.path(tmpdir, "train.csv")
  pred_path <- file.path(tmpdir, "predict.csv")
  out_path <- file.path(tmpdir, "pred.csv")
  meta_path <- file.path(tmpdir, "metadata.json")

  train_frame <- data.frame(.target = train_y, train_x, check.names = FALSE)
  pred_frame <- data.frame(pred_x, check.names = FALSE)
  utils::write.csv(train_frame, train_path, row.names = FALSE, na = "")
  utils::write.csv(pred_frame, pred_path, row.names = FALSE, na = "")

  cmd_args <- c(
    shQuote(runner_py),
    "--task",
    task,
    "--train",
    shQuote(train_path),
    "--predict",
    shQuote(pred_path),
    "--target",
    ".target",
    "--out",
    shQuote(out_path),
    "--metadata",
    shQuote(meta_path),
    "--device",
    device,
    "--seed",
    as.character(seed)
  )
  cmd_out <- system2(python, cmd_args, stdout = TRUE, stderr = TRUE)
  status <- attr(cmd_out, "status")
  if (!is.null(status) && status != 0L) {
    stop(paste(cmd_out, collapse = "\n"))
  }

  out <- utils::read.csv(out_path, check.names = FALSE)
  if (task == "classification") {
    p1 <- if ("prob_class_1" %in% names(out)) {
      out[["prob_class_1"]]
    } else {
      as.numeric(out$prediction == 1L)
    }
    list(prediction = as.integer(out$prediction), p1 = as.numeric(p1))
  } else {
    list(prediction = as.numeric(out$prediction), p1 = NULL)
  }
}

tag_eval <- function(
  ev,
  method,
  scenario,
  scale_n,
  rep_id,
  variant,
  wall_sec,
  status = "ok"
) {
  if (is.null(ev) || !nrow(ev)) {
    return(NULL)
  }
  ev <- ev[ev$split == "test", , drop = FALSE]
  if (!nrow(ev)) {
    return(NULL)
  }
  ev$method <- method
  ev$scenario <- scenario
  ev$scale_n <- scale_n
  ev$rep <- rep_id
  ev$variant <- variant
  ev$wall_sec <- wall_sec
  ev$status <- status
  ev
}

standard_methods <- function(
  pd,
  tree,
  graph,
  splits,
  scenario,
  scale_n,
  rep_id,
  seed
) {
  rows <- list()

  t0 <- proc.time()[["elapsed"]]
  mean_pred <- mean_impute(pd, splits)
  ev_mean <- evaluate_imputation(
    mean_pred,
    pd$X_scaled,
    splits,
    trait_map = pd$trait_map
  )
  rows[[length(rows) + 1L]] <- tag_eval(
    ev_mean,
    "mean",
    scenario,
    scale_n,
    rep_id,
    NA_character_,
    proc.time()[["elapsed"]] - t0
  )

  t1 <- proc.time()[["elapsed"]]
  baseline <- fit_baseline(pd, tree, splits = splits, graph = graph)
  ev_bl <- evaluate_imputation(
    baseline$mu,
    pd$X_scaled,
    splits,
    trait_map = pd$trait_map
  )
  rows[[length(rows) + 1L]] <- tag_eval(
    ev_bl,
    "baseline",
    scenario,
    scale_n,
    rep_id,
    NA_character_,
    proc.time()[["elapsed"]] - t1
  )

  if (run_pigauto) {
    graph_fit <- graph
    graph_fit$D <- NULL
    invisible(gc(full = TRUE, verbose = FALSE))
    t2 <- proc.time()[["elapsed"]]
    fit <- fit_pigauto(
      pd,
      tree,
      splits = splits,
      baseline = baseline,
      graph = graph_fit,
      epochs = epochs,
      verbose = FALSE,
      seed = seed
    )
    pred <- stats::predict(fit, return_se = FALSE, n_imputations = 1L)
    ev_pg <- evaluate_imputation(pred, pd$X_scaled, splits)
    rows[[length(rows) + 1L]] <- tag_eval(
      ev_pg,
      "pigauto",
      scenario,
      scale_n,
      rep_id,
      NA_character_,
      proc.time()[["elapsed"]] - t2
    )
  }

  do.call(rbind, rows)
}

tabpfn_zi_target <- function(
  pd,
  graph,
  splits,
  target,
  variant,
  seed,
  scenario,
  scale_n,
  rep_id
) {
  tm <- pd$trait_map[[target]]
  lc <- tm$latent_cols
  n <- nrow(pd$X_scaled)
  rows <- target_rows_from_split(splits, lc[1L], n)
  gate <- pd$X_scaled[, lc[1L]]
  mag <- pd$X_scaled[, lc[2L]]
  val_rows <- rows$val[is.finite(gate[rows$val])]
  test_rows <- rows$test[is.finite(gate[rows$test])]
  train_rows <- setdiff(which(is.finite(gate)), c(val_rows, test_rows))

  base_pred <- pd$X_scaled
  method <- paste0("tabpfn_", variant)
  status <- "ok"
  wall <- NA_real_

  if (
    length(train_rows) < 5L ||
      length(test_rows) < 2L ||
      length(unique(gate[train_rows])) < 2L
  ) {
    status <- "skipped_too_few_gate_cells"
  } else if (dry_run) {
    status <- "dry_run"
  } else {
    t0 <- proc.time()[["elapsed"]]
    fit_result <- tryCatch(
      {
        x_gate <- feature_matrix(
          pd,
          graph,
          splits,
          lc,
          variant,
          train_rows,
          class_labels = as.integer(gate)
        )
        gate_fit <- run_tabpfn(
          x_gate[train_rows, , drop = FALSE],
          as.integer(gate[train_rows]),
          x_gate[test_rows, , drop = FALSE],
          "classification",
          seed
        )

        mag_train <- train_rows[is.finite(mag[train_rows])]
        if (length(mag_train) < 5L) {
          stop("Too few non-zero magnitude training rows.")
        }
        x_mag <- feature_matrix(
          pd,
          graph,
          splits,
          lc,
          variant,
          mag_train,
          target_values = mag
        )
        mag_fit <- run_tabpfn(
          x_mag[mag_train, , drop = FALSE],
          mag[mag_train],
          x_mag[test_rows, , drop = FALSE],
          "regression",
          seed + 17L
        )
        list(ok = TRUE, p_nz = gate_fit$p1, mag = mag_fit$prediction)
      },
      error = function(e) {
        list(ok = FALSE, message = conditionMessage(e))
      }
    )
    wall <- proc.time()[["elapsed"]] - t0

    if (isTRUE(fit_result$ok)) {
      base_pred[test_rows, lc[1L]] <- logit(fit_result$p_nz)
      base_pred[test_rows, lc[2L]] <- fit_result$mag
    } else {
      status <- paste0("error: ", fit_result$message)
    }
  }

  if (status != "ok") {
    out <- data.frame(
      split = "test",
      trait = target,
      type = tm$type,
      n = length(test_rows),
      rmse = NA_real_,
      pearson_r = NA_real_,
      coverage_95 = NA_real_,
      mae = NA_real_,
      spearman_rho = NA_real_,
      accuracy = NA_real_,
      brier = NA_real_,
      zero_accuracy = NA_real_,
      aitchison = NA_real_,
      rmse_clr = NA_real_,
      simplex_mae = NA_real_,
      method = method,
      scenario = scenario,
      scale_n = scale_n,
      rep = rep_id,
      variant = variant,
      wall_sec = wall,
      status = status,
      stringsAsFactors = FALSE
    )
    return(out)
  }

  ev <- evaluate_imputation(
    base_pred,
    pd$X_scaled,
    splits,
    trait_map = single_trait_map(pd, target)
  )
  tag_eval(ev, method, scenario, scale_n, rep_id, variant, wall, status)
}

tabpfn_multi_target <- function(
  pd,
  graph,
  splits,
  target,
  variant,
  seed,
  scenario,
  scale_n,
  rep_id
) {
  tm <- pd$trait_map[[target]]
  lc <- tm$latent_cols
  n <- nrow(pd$X_scaled)
  rows <- target_rows_from_split(splits, lc[1L], n)
  y <- pd$X_scaled[, lc, drop = FALSE]
  complete <- stats::complete.cases(y)
  val_rows <- rows$val[complete[rows$val]]
  test_rows <- rows$test[complete[rows$test]]
  train_rows <- setdiff(which(complete), c(val_rows, test_rows))

  base_pred <- pd$X_scaled
  method <- paste0("tabpfn_", variant)
  status <- "ok"
  wall <- NA_real_

  if (length(train_rows) < 5L || length(test_rows) < 2L) {
    status <- "skipped_too_few_rows"
  } else if (dry_run) {
    status <- "dry_run"
  } else {
    t0 <- proc.time()[["elapsed"]]
    fit_result <- tryCatch(
      {
        pred <- matrix(NA_real_, nrow = length(test_rows), ncol = length(lc))
        for (k in seq_along(lc)) {
          y_k <- pd$X_scaled[, lc[k]]
          x_k <- feature_matrix(
            pd,
            graph,
            splits,
            lc,
            variant,
            train_rows,
            target_values = y_k
          )
          fit_k <- run_tabpfn(
            x_k[train_rows, , drop = FALSE],
            y_k[train_rows],
            x_k[test_rows, , drop = FALSE],
            "regression",
            seed + k
          )
          pred[, k] <- fit_k$prediction
        }
        list(ok = TRUE, pred = pred)
      },
      error = function(e) {
        list(ok = FALSE, message = conditionMessage(e))
      }
    )
    wall <- proc.time()[["elapsed"]] - t0

    if (isTRUE(fit_result$ok)) {
      base_pred[test_rows, lc] <- fit_result$pred
    } else {
      status <- paste0("error: ", fit_result$message)
    }
  }

  if (status != "ok") {
    out <- data.frame(
      split = "test",
      trait = target,
      type = tm$type,
      n = length(test_rows),
      rmse = NA_real_,
      pearson_r = NA_real_,
      coverage_95 = NA_real_,
      mae = NA_real_,
      spearman_rho = NA_real_,
      accuracy = NA_real_,
      brier = NA_real_,
      zero_accuracy = NA_real_,
      aitchison = NA_real_,
      rmse_clr = NA_real_,
      simplex_mae = NA_real_,
      method = method,
      scenario = scenario,
      scale_n = scale_n,
      rep = rep_id,
      variant = variant,
      wall_sec = wall,
      status = status,
      stringsAsFactors = FALSE
    )
    return(out)
  }

  ev <- evaluate_imputation(
    base_pred,
    pd$X_scaled,
    splits,
    trait_map = single_trait_map(pd, target)
  )
  tag_eval(ev, method, scenario, scale_n, rep_id, variant, wall, status)
}

tabpfn_methods <- function(pd, graph, splits, scenario, scale_n, rep_id, seed) {
  rows <- list()
  for (variant in variants) {
    for (target in names(pd$trait_map)) {
      tm <- pd$trait_map[[target]]
      log_line(
        "TabPFN ",
        variant,
        " target=",
        target,
        " type=",
        tm$type,
        " n=",
        scale_n,
        " scenario=",
        scenario,
        " rep=",
        rep_id
      )
      rows[[length(rows) + 1L]] <- if (tm$type == "zi_count") {
        tabpfn_zi_target(
          pd,
          graph,
          splits,
          target,
          variant,
          seed,
          scenario,
          scale_n,
          rep_id
        )
      } else {
        tabpfn_multi_target(
          pd,
          graph,
          splits,
          target,
          variant,
          seed,
          scenario,
          scale_n,
          rep_id
        )
      }
    }
  }
  do.call(rbind, rows)
}

mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

write_summary <- function(results, metadata, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("# TabPFN Complex-Type Benchmark\n", con)
  writeLines(
    sprintf("- Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    con
  )
  writeLines(sprintf("- Git commit: `%s`", metadata$commit), con)
  writeLines(
    sprintf("- Scales: `%s`", paste(metadata$scales, collapse = ", ")),
    con
  )
  writeLines(
    sprintf("- Scenarios: `%s`", paste(metadata$scenarios, collapse = ", ")),
    con
  )
  writeLines(
    sprintf("- Variants: `%s`", paste(metadata$variants, collapse = ", ")),
    con
  )
  writeLines(sprintf("- Replicates: `%d`", metadata$n_reps), con)
  writeLines("", con)
  writeLines("## Design", con)
  writeLines("", con)
  writeLines(
    paste(
      "ZI-count uses one TabPFN classifier for the non-zero gate and one",
      "TabPFN regressor for the conditional magnitude. Multi-proportion uses",
      "independent TabPFN regressions on the z-scored CLR latent columns and is",
      "scored by pigauto's existing Aitchison, CLR RMSE, simplex MAE, and",
      "dominant-component accuracy metrics."
    ),
    con
  )
  writeLines("", con)
  writeLines("## Status Counts", con)
  writeLines("", con)
  status_counts <- as.data.frame(table(
    method = results$method,
    status = results$status,
    useNA = "ifany"
  ))
  status_counts <- status_counts[status_counts$Freq > 0L, , drop = FALSE]
  writeLines(capture.output(print(status_counts, row.names = FALSE)), con)
  writeLines("", con)
  writeLines("## Test Summary", con)
  writeLines("", con)
  scored <- results[results$status == "ok", , drop = FALSE]
  if (nrow(scored)) {
    groups <- split(
      scored,
      interaction(
        scored$method,
        scored$scenario,
        scored$scale_n,
        scored$trait,
        drop = TRUE
      )
    )
    agg <- do.call(
      rbind,
      lapply(groups, function(x) {
        data.frame(
          method = x$method[1L],
          scenario = x$scenario[1L],
          scale_n = x$scale_n[1L],
          trait = x$trait[1L],
          type = x$type[1L],
          rmse = mean_or_na(x$rmse),
          mae = mean_or_na(x$mae),
          zero_accuracy = mean_or_na(x$zero_accuracy),
          brier = mean_or_na(x$brier),
          aitchison = mean_or_na(x$aitchison),
          rmse_clr = mean_or_na(x$rmse_clr),
          simplex_mae = mean_or_na(x$simplex_mae),
          accuracy = mean_or_na(x$accuracy),
          stringsAsFactors = FALSE
        )
      })
    )
    agg <- agg[
      order(agg$scenario, agg$scale_n, agg$type, agg$trait, agg$method),
      ,
      drop = FALSE
    ]
    writeLines(capture.output(print(agg, row.names = FALSE, digits = 4)), con)
  } else {
    writeLines("No scored rows were produced. Inspect `Status Counts`.", con)
  }
}

all_results <- list()

for (scale_n in scales) {
  for (scenario in scenarios) {
    for (rep_id in seq_len(n_reps)) {
      seed <- seed0 +
        scale_n * 1000L +
        match(scenario, valid_scenarios) * 100L +
        rep_id
      log_line(
        "Preparing complex type n=",
        scale_n,
        " scenario=",
        scenario,
        " rep=",
        rep_id
      )
      dat <- generate_data(scale_n, scenario, seed)
      pd <- preprocess_traits(
        dat$df,
        dat$tree,
        log_transform = FALSE,
        trait_types = dat$trait_types,
        multi_proportion_groups = dat$multi_proportion_groups
      )
      graph <- build_phylo_graph(dat$tree, k_eigen = "auto")
      splits <- make_missing_splits(
        pd$X_scaled,
        missing_frac = missing_frac,
        val_frac = val_frac,
        seed = seed,
        trait_map = pd$trait_map
      )

      all_results[[length(all_results) + 1L]] <-
        standard_methods(
          pd,
          dat$tree,
          graph,
          splits,
          scenario,
          scale_n,
          rep_id,
          seed
        )
      all_results[[length(all_results) + 1L]] <-
        tabpfn_methods(pd, graph, splits, scenario, scale_n, rep_id, seed)

      results <- do.call(rbind, all_results)
      saveRDS(list(results = results, metadata = metadata), out_rds)
      write_summary(results, metadata, out_md)
    }
  }
}

results <- do.call(rbind, all_results)
saveRDS(list(results = results, metadata = metadata), out_rds)
write_summary(results, metadata, out_md)
log_line("Wrote ", out_rds)
log_line("Wrote ", out_md)
