#!/usr/bin/env Rscript
#
# script/bench_tabpfn_sim_types.R
#
# Branch-local TabPFN benchmark on simulated mixed scalar trait types.
# This extends the AVONET-only TabPFN checks to pigauto's scalar families:
# continuous, count, proportion, binary, ordinal, and categorical.
#
# Smoke:
#   PIGAUTO_TABPFN_SCALES=50 PIGAUTO_TABPFN_REPS=1 \
#   PIGAUTO_TABPFN_SCENARIOS=mixed_moderate \
#   PIGAUTO_TABPFN_VARIANTS=lappe \
#   PIGAUTO_TABPFN_RUN_PIGAUTO=false \
#     Rscript script/bench_tabpfn_sim_types.R
#
# Local grid:
#   PIGAUTO_TABPFN_OUT_STEM=bench_tabpfn_sim_types_local \
#   PIGAUTO_TABPFN_SCALES=75,150 PIGAUTO_TABPFN_REPS=2 \
#   PIGAUTO_TABPFN_SCENARIOS=mixed_moderate,mixed_high_phylo,mixed_sparse_imbalanced \
#   PIGAUTO_TABPFN_VARIANTS=plain,lappe,lappe_nfa \
#   PIGAUTO_TABPFN_RUN_PIGAUTO=true \
#     Rscript script/bench_tabpfn_sim_types.R

options(warn = 1, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(".", quiet = TRUE)
})

script_start <- proc.time()[["elapsed"]]

log_line <- function(...) {
  cat(sprintf("[%6.1fs] ", proc.time()[["elapsed"]] - script_start), ...,
      "\n", sep = "")
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
out_stem <- env_chr("PIGAUTO_TABPFN_OUT_STEM", "bench_tabpfn_sim_types")
if (!grepl("^[A-Za-z0-9_.-]+$", out_stem)) {
  stop("PIGAUTO_TABPFN_OUT_STEM must be a simple file stem, not a path.")
}
out_rds <- file.path(repo_root, "script", paste0(out_stem, ".rds"))
out_md <- file.path(repo_root, "script", paste0(out_stem, ".md"))
runner_py <- file.path(repo_root, "script", "run_tabpfn_phylo.py")

scales <- as.integer(env_vec("PIGAUTO_TABPFN_SCALES", c("75", "150")))
scenarios <- env_vec("PIGAUTO_TABPFN_SCENARIOS", c(
  "mixed_moderate",
  "mixed_high_phylo",
  "mixed_sparse_imbalanced"
))
variants <- env_vec("PIGAUTO_TABPFN_VARIANTS",
                    c("plain", "lappe", "lappe_nfa", "knn"))
n_reps <- env_int("PIGAUTO_TABPFN_REPS", 1L)
seed0 <- env_int("PIGAUTO_TABPFN_SEED", 20260611L)
missing_frac <- env_num("PIGAUTO_TABPFN_MISSING_FRAC", 0.25)
val_frac <- env_num("PIGAUTO_TABPFN_VAL_FRAC", 0.5)
epochs <- env_int("PIGAUTO_TABPFN_EPOCHS", 300L)
nfa_k <- env_int("PIGAUTO_TABPFN_NFA_K", 15L)
dry_run <- env_lgl("PIGAUTO_TABPFN_DRY_RUN", FALSE)
run_pigauto <- env_lgl("PIGAUTO_TABPFN_RUN_PIGAUTO", !dry_run)
python <- env_chr("PIGAUTO_TABPFN_PYTHON", "python3")
device <- env_chr("PIGAUTO_TABPFN_DEVICE", "auto")
target_env <- env_vec("PIGAUTO_TABPFN_TARGETS", character())

valid_variants <- c("plain", "lappe", "lappe_nfa", "knn")
if (!all(variants %in% valid_variants)) {
  stop("Unknown variants: ",
       paste(setdiff(variants, valid_variants), collapse = ", "))
}

valid_scenarios <- c("mixed_moderate", "mixed_high_phylo",
                     "mixed_sparse_imbalanced")
if (!all(scenarios %in% valid_scenarios)) {
  stop("Unknown scenarios: ",
       paste(setdiff(scenarios, valid_scenarios), collapse = ", "))
}

metadata <- list(
  benchmark = "tabpfn_sim_types",
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
  commit = tryCatch(system("git rev-parse HEAD", intern = TRUE),
                    error = function(e) "unknown")
)

regression_types <- c("continuous", "count", "proportion")
classification_types <- c("binary", "ordinal", "categorical")

rmse_vec <- function(truth, pred) {
  ok <- is.finite(truth) & is.finite(pred)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((truth[ok] - pred[ok])^2))
}

cor_vec <- function(truth, pred) {
  ok <- is.finite(truth) & is.finite(pred)
  if (sum(ok) < 2L) return(NA_real_)
  if (stats::sd(truth[ok]) == 0 || stats::sd(pred[ok]) == 0) return(NA_real_)
  stats::cor(truth[ok], pred[ok])
}

mode_int <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(0L)
  as.integer(names(sort(table(x), decreasing = TRUE))[1L])
}

accuracy_vec <- function(truth, pred) {
  ok <- !is.na(truth) & !is.na(pred)
  if (!any(ok)) return(NA_real_)
  mean(truth[ok] == pred[ok])
}

balanced_accuracy_vec <- function(truth, pred) {
  ok <- !is.na(truth) & !is.na(pred)
  if (!any(ok)) return(NA_real_)
  per_class <- tapply(truth[ok] == pred[ok], truth[ok], mean)
  mean(per_class, na.rm = TRUE)
}

scenario_params <- function(scenario) {
  switch(
    scenario,
    mixed_moderate = list(
      signal = 0.6, threshold_quantile = 0.5, mean_count = 20L,
      overdispersion = NULL, boundary_frac = 0.0
    ),
    mixed_high_phylo = list(
      signal = 0.9, threshold_quantile = 0.5, mean_count = 20L,
      overdispersion = NULL, boundary_frac = 0.0
    ),
    mixed_sparse_imbalanced = list(
      signal = 0.6, threshold_quantile = 0.8, mean_count = 5L,
      overdispersion = 2, boundary_frac = 0.2
    ),
    stop("Unknown scenario: ", scenario)
  )
}

rename_cols <- function(df, names_out) {
  names(df) <- names_out
  df
}

generate_mixed_scalar_data <- function(n, scenario, seed) {
  pars <- scenario_params(scenario)
  set.seed(seed)
  tree <- ape::rtree(n)

  cont <- rename_cols(
    simulate_bm_traits(tree, n_traits = 2L, seed = seed + 1L),
    c("continuous_1", "continuous_2")
  )
  count <- rename_cols(
    simulate_count_traits(tree, n_traits = 1L,
                          mean_count = pars$mean_count,
                          overdispersion = pars$overdispersion,
                          seed = seed + 2L),
    "count_1"
  )
  prop <- rename_cols(
    simulate_proportion_traits(tree, n_traits = 1L, signal = pars$signal,
                               boundary_frac = pars$boundary_frac,
                               seed = seed + 3L),
    "proportion_1"
  )
  binary <- rename_cols(
    simulate_binary_traits(tree, n_traits = 1L, signal = pars$signal,
                           threshold_quantile = pars$threshold_quantile,
                           seed = seed + 4L),
    "binary_1"
  )
  ordinal <- rename_cols(
    simulate_ordinal_traits(tree, n_traits = 1L, n_levels = 5L,
                            signal = pars$signal, seed = seed + 5L),
    "ordinal_1"
  )
  categorical <- rename_cols(
    simulate_categorical_traits(tree, n_traits = 1L, n_levels = 3L,
                                signal = pars$signal, seed = seed + 6L),
    "categorical_1"
  )

  df <- data.frame(cont, count, prop, binary, ordinal, categorical,
                   check.names = FALSE)
  rownames(df) <- tree$tip.label
  list(
    df = df,
    tree = tree,
    trait_types = c(proportion_1 = "proportion")
  )
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

class_labels_from_truth <- function(X, tm) {
  n <- nrow(X)
  labels <- rep(NA_integer_, n)
  if (tm$type == "categorical") {
    mat <- X[, tm$latent_cols, drop = FALSE]
    ok <- stats::complete.cases(mat)
    labels[ok] <- max.col(mat[ok, , drop = FALSE], ties.method = "first") - 1L
  } else if (tm$type == "ordinal") {
    z <- X[, tm$latent_cols[1L]]
    ok <- is.finite(z)
    cls <- round(z[ok] * tm$sd + tm$mean)
    cls <- pmax(pmin(cls, length(tm$levels) - 1L), 0L)
    labels[ok] <- as.integer(cls)
  } else if (tm$type == "binary") {
    z <- X[, tm$latent_cols[1L]]
    ok <- is.finite(z)
    labels[ok] <- as.integer(round(z[ok]))
  } else {
    stop("Unsupported classification type: ", tm$type)
  }
  labels
}

class_labels_from_pred <- function(pred_mat, tm) {
  n <- nrow(pred_mat)
  labels <- rep(NA_integer_, n)
  if (tm$type == "categorical") {
    mat <- pred_mat[, tm$latent_cols, drop = FALSE]
    ok <- stats::complete.cases(mat)
    labels[ok] <- max.col(mat[ok, , drop = FALSE], ties.method = "first") - 1L
  } else if (tm$type == "ordinal") {
    z <- pred_mat[, tm$latent_cols[1L]]
    ok <- is.finite(z)
    cls <- round(z[ok] * tm$sd + tm$mean)
    cls <- pmax(pmin(cls, length(tm$levels) - 1L), 0L)
    labels[ok] <- as.integer(cls)
  } else if (tm$type == "binary") {
    z <- pred_mat[, tm$latent_cols[1L]]
    ok <- is.finite(z)
    labels[ok] <- as.integer(1 / (1 + exp(-z[ok])) >= 0.5)
  } else {
    stop("Unsupported classification type: ", tm$type)
  }
  labels
}

naive_impute <- function(pd, splits) {
  X <- pd$X_scaled
  X_train <- X
  X_train[c(splits$val_idx, splits$test_idx)] <- NA_real_
  out <- matrix(NA_real_, nrow = nrow(X), ncol = ncol(X),
                dimnames = dimnames(X))

  for (tm in pd$trait_map) {
    cols <- tm$latent_cols
    if (tm$type == "categorical") {
      labels <- class_labels_from_truth(X_train, tm)
      maj <- mode_int(labels)
      out[, cols] <- 0
      out[, cols[maj + 1L]] <- 1
    } else if (tm$type == "binary") {
      labels <- class_labels_from_truth(X_train, tm)
      out[, cols[1L]] <- mode_int(labels)
    } else if (tm$type == "ordinal") {
      labels <- class_labels_from_truth(X_train, tm)
      maj <- mode_int(labels)
      out[, cols[1L]] <- (maj - tm$mean) / tm$sd
    } else {
      means <- colMeans(X_train[, cols, drop = FALSE], na.rm = TRUE)
      means[!is.finite(means)] <- 0
      out[, cols] <- matrix(means, nrow = nrow(X), ncol = length(cols),
                            byrow = TRUE)
    }
  }

  out
}

neighbor_regression_features <- function(y, D, train_rows, k) {
  n <- length(y)
  out <- matrix(NA_real_, nrow = n, ncol = 4L)
  colnames(out) <- c("nfa_mean", "nfa_sd", "nfa_min_dist", "nfa_n")

  train_rows <- train_rows[is.finite(y[train_rows])]
  if (length(train_rows) < 2L) return(out)

  for (i in seq_len(n)) {
    candidates <- setdiff(train_rows, i)
    if (!length(candidates)) next
    d <- D[i, candidates]
    ord <- order(d, decreasing = FALSE, na.last = NA)
    if (!length(ord)) next
    keep <- candidates[ord[seq_len(min(k, length(ord)))]]
    vals <- y[keep]
    out[i, "nfa_mean"] <- mean(vals, na.rm = TRUE)
    out[i, "nfa_sd"] <- if (length(vals) > 1L) stats::sd(vals, na.rm = TRUE) else 0
    out[i, "nfa_min_dist"] <- min(D[i, keep], na.rm = TRUE)
    out[i, "nfa_n"] <- length(vals)
  }

  out
}

neighbor_class_features <- function(labels, D, train_rows, k, n_classes) {
  n <- length(labels)
  props <- matrix(NA_real_, nrow = n, ncol = n_classes)
  colnames(props) <- paste0("nfa_p_class_", seq_len(n_classes) - 1L)
  extra <- matrix(NA_real_, nrow = n, ncol = 3L)
  colnames(extra) <- c("nfa_entropy", "nfa_min_dist", "nfa_n")

  train_rows <- train_rows[!is.na(labels[train_rows])]
  if (length(train_rows) < 2L) return(cbind(props, extra))

  for (i in seq_len(n)) {
    candidates <- setdiff(train_rows, i)
    if (!length(candidates)) next
    d <- D[i, candidates]
    ord <- order(d, decreasing = FALSE, na.last = NA)
    if (!length(ord)) next
    keep <- candidates[ord[seq_len(min(k, length(ord)))]]
    labs <- labels[keep]
    tab <- tabulate(labs + 1L, nbins = n_classes)
    p <- tab / sum(tab)
    props[i, ] <- p
    extra[i, "nfa_entropy"] <- -sum(ifelse(p > 0, p * log(p), 0))
    extra[i, "nfa_min_dist"] <- min(D[i, keep], na.rm = TRUE)
    extra[i, "nfa_n"] <- length(keep)
  }

  cbind(props, extra)
}

feature_matrix <- function(pd, graph, splits, target_cols, variant,
                           train_rows, target_values = NULL,
                           class_labels = NULL, n_classes = NULL) {
  X <- pd$X_scaled
  X[c(splits$val_idx, splits$test_idx)] <- NA_real_
  X <- X[, setdiff(seq_len(ncol(X)), target_cols), drop = FALSE]

  parts <- list(traits = X)
  if (variant %in% c("lappe", "lappe_nfa")) {
    coords <- graph$coords
    colnames(coords) <- paste0("lappe_", seq_len(ncol(coords)))
    parts$lappe <- coords
  }
  if (variant %in% c("lappe_nfa", "knn")) {
    if (!is.null(class_labels)) {
      parts$nfa <- neighbor_class_features(class_labels, graph$D, train_rows,
                                           nfa_k, n_classes)
    } else {
      parts$nfa <- neighbor_regression_features(target_values, graph$D,
                                                train_rows, nfa_k)
    }
  }

  out <- do.call(cbind, parts)
  storage.mode(out) <- "double"
  out
}

run_tabpfn <- function(train_x, train_y, pred_x, task, seed) {
  tmpdir <- tempfile("tabpfn_sim_types_")
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
    "--task", task,
    "--train", shQuote(train_path),
    "--predict", shQuote(pred_path),
    "--target", ".target",
    "--out", shQuote(out_path),
    "--metadata", shQuote(meta_path),
    "--device", device,
    "--seed", as.character(seed)
  )

  cmd_out <- system2(python, cmd_args, stdout = TRUE, stderr = TRUE)
  status <- attr(cmd_out, "status")
  if (!is.null(status) && status != 0L) {
    stop(paste(cmd_out, collapse = "\n"))
  }

  pred <- utils::read.csv(out_path)$prediction
  if (task == "classification") as.integer(pred) else as.numeric(pred)
}

result_row <- function(method, scenario, scale_n, rep_id, target, type, task,
                       n_train, n_val, n_test, n_features, truth, pred,
                       wall_sec, status) {
  if (task == "classification") {
    rmse <- NA_real_
    pearson_r <- NA_real_
    accuracy <- accuracy_vec(truth, pred)
    balanced_accuracy <- balanced_accuracy_vec(truth, pred)
  } else {
    rmse <- rmse_vec(truth, pred)
    pearson_r <- cor_vec(truth, pred)
    accuracy <- NA_real_
    balanced_accuracy <- NA_real_
  }

  data.frame(
    method = method,
    scenario = scenario,
    scale_n = scale_n,
    rep = rep_id,
    trait = target,
    type = type,
    task = task,
    split = "test",
    n_train = n_train,
    n_val = n_val,
    n_test = n_test,
    n_features = n_features,
    rmse = rmse,
    pearson_r = pearson_r,
    accuracy = accuracy,
    balanced_accuracy = balanced_accuracy,
    wall_sec = wall_sec,
    status = status,
    stringsAsFactors = FALSE
  )
}

target_truth_and_pred <- function(pd, splits, target, pred_mat) {
  tm <- pd$trait_map[[target]]
  n <- nrow(pd$X_scaled)
  rows <- target_rows_from_split(splits, tm$latent_cols[1L], n)

  if (tm$type %in% regression_types) {
    y <- pd$X_scaled[, tm$latent_cols[1L]]
    test_rows <- rows$test[is.finite(y[rows$test])]
    list(
      task = "regression",
      truth = y[test_rows],
      pred = pred_mat[test_rows, tm$latent_cols[1L]],
      test_rows = test_rows
    )
  } else {
    labels <- class_labels_from_truth(pd$X_scaled, tm)
    pred_labels <- class_labels_from_pred(pred_mat, tm)
    test_rows <- rows$test[!is.na(labels[rows$test])]
    list(
      task = "classification",
      truth = labels[test_rows],
      pred = pred_labels[test_rows],
      test_rows = test_rows
    )
  }
}

append_standard_rows <- function(rows, pd, splits, targets, pred_mat, method,
                                 scenario, scale_n, rep_id, wall_sec) {
  n <- nrow(pd$X_scaled)
  for (target in targets) {
    tm <- pd$trait_map[[target]]
    split_rows <- target_rows_from_split(splits, tm$latent_cols[1L], n)
    train_rows <- setdiff(
      seq_len(n),
      unique(c(split_rows$val, split_rows$test))
    )
    vals <- target_truth_and_pred(pd, splits, target, pred_mat)
    rows[[length(rows) + 1L]] <- result_row(
      method = method,
      scenario = scenario,
      scale_n = scale_n,
      rep_id = rep_id,
      target = target,
      type = tm$type,
      task = vals$task,
      n_train = length(train_rows),
      n_val = length(split_rows$val),
      n_test = length(vals$test_rows),
      n_features = NA_integer_,
      truth = vals$truth,
      pred = vals$pred,
      wall_sec = wall_sec,
      status = "ok"
    )
  }
  rows
}

standard_methods <- function(pd, tree, graph, splits, targets, seed,
                             scenario, scale_n, rep_id) {
  rows <- list()

  t0 <- proc.time()[["elapsed"]]
  naive <- naive_impute(pd, splits)
  naive_wall <- proc.time()[["elapsed"]] - t0
  rows <- append_standard_rows(rows, pd, splits, targets, naive, "mean_mode",
                               scenario, scale_n, rep_id, naive_wall)

  t1 <- proc.time()[["elapsed"]]
  baseline <- fit_baseline(pd, tree, splits = splits, graph = graph)
  baseline_wall <- proc.time()[["elapsed"]] - t1
  rows <- append_standard_rows(rows, pd, splits, targets, baseline$mu,
                               "baseline", scenario, scale_n, rep_id,
                               baseline_wall)

  if (run_pigauto) {
    graph_fit <- graph
    graph_fit$D <- NULL
    invisible(gc(full = TRUE, verbose = FALSE))
    t2 <- proc.time()[["elapsed"]]
    fit <- fit_pigauto(pd, tree, splits = splits, baseline = baseline,
                       graph = graph_fit, epochs = epochs, verbose = FALSE,
                       seed = seed)
    pred <- stats::predict(fit, return_se = TRUE, n_imputations = 1L)
    pigauto_wall <- proc.time()[["elapsed"]] - t2
    rows <- append_standard_rows(rows, pd, splits, targets, pred$imputed_latent,
                                 "pigauto", scenario, scale_n, rep_id,
                                 pigauto_wall)
  }

  do.call(rbind, rows)
}

tabpfn_one_target <- function(pd, graph, splits, target, variant, seed,
                              scenario, scale_n, rep_id) {
  tm <- pd$trait_map[[target]]
  n <- nrow(pd$X_scaled)
  split_rows <- target_rows_from_split(splits, tm$latent_cols[1L], n)

  if (tm$type %in% regression_types) {
    task <- "regression"
    y <- pd$X_scaled[, tm$latent_cols[1L]]
    val_rows <- split_rows$val[is.finite(y[split_rows$val])]
    test_rows <- split_rows$test[is.finite(y[split_rows$test])]
    train_rows <- setdiff(which(is.finite(y)), c(val_rows, test_rows))
    x <- feature_matrix(pd, graph, splits, tm$latent_cols, variant,
                        train_rows, target_values = y)
    truth <- y[test_rows]
    train_y <- y[train_rows]
  } else {
    task <- "classification"
    labels <- class_labels_from_truth(pd$X_scaled, tm)
    val_rows <- split_rows$val[!is.na(labels[split_rows$val])]
    test_rows <- split_rows$test[!is.na(labels[split_rows$test])]
    train_rows <- setdiff(which(!is.na(labels)), c(val_rows, test_rows))
    n_classes <- if (is.null(tm$levels)) 2L else length(tm$levels)
    x <- feature_matrix(pd, graph, splits, tm$latent_cols, variant,
                        train_rows, class_labels = labels,
                        n_classes = n_classes)
    truth <- labels[test_rows]
    train_y <- labels[train_rows]
  }

  method <- paste0("tabpfn_", variant)
  n_features <- ncol(x)
  if (length(train_rows) < 5L || length(test_rows) < 2L) {
    return(result_row(
      method, scenario, scale_n, rep_id, target, tm$type, task,
      length(train_rows), length(val_rows), length(test_rows), n_features,
      truth, rep(NA, length(truth)), NA_real_, "skipped_too_few_cells"
    ))
  }
  if (task == "classification" && length(unique(train_y)) < 2L) {
    return(result_row(
      method, scenario, scale_n, rep_id, target, tm$type, task,
      length(train_rows), length(val_rows), length(test_rows), n_features,
      truth, rep(NA_integer_, length(truth)), NA_real_,
      "skipped_too_few_classes"
    ))
  }
  if (dry_run) {
    return(result_row(
      method, scenario, scale_n, rep_id, target, tm$type, task,
      length(train_rows), length(val_rows), length(test_rows), n_features,
      truth, rep(NA, length(truth)), NA_real_, "dry_run"
    ))
  }

  t0 <- proc.time()[["elapsed"]]
  fit_result <- tryCatch({
    pred <- run_tabpfn(x[train_rows, , drop = FALSE], train_y,
                       x[test_rows, , drop = FALSE], task, seed)
    list(ok = TRUE, pred = pred)
  }, error = function(e) {
    list(ok = FALSE, message = conditionMessage(e))
  })
  wall <- proc.time()[["elapsed"]] - t0

  if (!isTRUE(fit_result$ok)) {
    return(result_row(
      method, scenario, scale_n, rep_id, target, tm$type, task,
      length(train_rows), length(val_rows), length(test_rows), n_features,
      truth, rep(NA, length(truth)), wall,
      paste0("error: ", fit_result$message)
    ))
  }

  result_row(
    method, scenario, scale_n, rep_id, target, tm$type, task,
    length(train_rows), length(val_rows), length(test_rows), n_features,
    truth, fit_result$pred, wall, "ok"
  )
}

mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

write_summary <- function(results, metadata, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("# TabPFN Simulated Mixed Scalar-Type Benchmark\n", con)
  writeLines(sprintf("- Generated: %s",
                     format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), con)
  writeLines(sprintf("- Git commit: `%s`", metadata$commit), con)
  writeLines(sprintf("- Scales: `%s`",
                     paste(metadata$scales, collapse = ", ")), con)
  writeLines(sprintf("- Scenarios: `%s`",
                     paste(metadata$scenarios, collapse = ", ")), con)
  writeLines(sprintf("- Variants: `%s`",
                     paste(metadata$variants, collapse = ", ")), con)
  writeLines(sprintf("- Replicates: `%d`", metadata$n_reps), con)
  writeLines(sprintf("- Dry run: `%s`", metadata$dry_run), con)
  writeLines("", con)

  writeLines("## Design", con)
  writeLines("", con)
  writeLines(paste(
    "Each cell simulates one mixed scalar dataset with two continuous traits",
    "and one trait each for count, proportion, binary, ordinal, and",
    "categorical data. Missing cells are held out with pigauto's",
    "`make_missing_splits()`. TabPFN uses regression for continuous, count,",
    "and proportion latent targets and classification for binary, ordinal,",
    "and categorical targets. This branch-local check does not cover `zi_count`,",
    "`multi_proportion`, classification prediction sets, multi-tree Rubin",
    "pooling, multi-observation covariates, or active-imputation guidance."
  ), con)
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

  scored <- results[results$status == "ok", , drop = FALSE]
  writeLines("## Test Summary", con)
  writeLines("", con)
  if (nrow(scored)) {
    groups <- split(scored, interaction(scored$method, scored$scenario,
                                        scored$scale_n, scored$trait,
                                        drop = TRUE))
    agg <- do.call(rbind, lapply(groups, function(x) {
      data.frame(
        method = x$method[1L],
        scenario = x$scenario[1L],
        scale_n = x$scale_n[1L],
        trait = x$trait[1L],
        type = x$type[1L],
        task = x$task[1L],
        rmse = mean_or_na(x$rmse),
        pearson_r = mean_or_na(x$pearson_r),
        accuracy = mean_or_na(x$accuracy),
        balanced_accuracy = mean_or_na(x$balanced_accuracy),
        stringsAsFactors = FALSE
      )
    }))
    agg <- agg[order(agg$scenario, agg$scale_n, agg$type, agg$trait,
                     agg$method), , drop = FALSE]
    writeLines(capture.output(print(agg, row.names = FALSE, digits = 4)), con)
  } else if (isTRUE(metadata$dry_run)) {
    writeLines("No scored rows yet. This is expected for dry runs.", con)
  } else {
    writeLines("No scored rows were produced. Inspect `Status Counts`.", con)
  }
}

all_results <- list()

for (scale_n in scales) {
  for (scenario in scenarios) {
    for (rep_id in seq_len(n_reps)) {
      seed <- seed0 + scale_n * 1000L + match(scenario, valid_scenarios) * 100L + rep_id
      log_line("Preparing simulated mixed scalar n=", scale_n,
               " scenario=", scenario, " rep=", rep_id)
      dat <- generate_mixed_scalar_data(scale_n, scenario, seed)
      pd <- preprocess_traits(dat$df, dat$tree, log_transform = FALSE,
                              trait_types = dat$trait_types)
      graph <- build_phylo_graph(dat$tree, k_eigen = "auto")
      splits <- make_missing_splits(pd$X_scaled, missing_frac = missing_frac,
                                    val_frac = val_frac, seed = seed,
                                    trait_map = pd$trait_map)

      targets <- names(Filter(function(tm) {
        tm$type %in% c(regression_types, classification_types)
      }, pd$trait_map))
      if (length(target_env)) {
        targets <- intersect(targets, target_env)
      }
      if (!length(targets)) {
        stop("No targets selected.")
      }

      log_line("Targets: ", paste(targets, collapse = ", "))
      all_results[[length(all_results) + 1L]] <-
        standard_methods(pd, dat$tree, graph, splits, targets, seed,
                         scenario, scale_n, rep_id)

      for (variant in variants) {
        for (target in targets) {
          log_line("TabPFN ", variant, " target=", target,
                   " n=", scale_n, " scenario=", scenario,
                   " rep=", rep_id)
          all_results[[length(all_results) + 1L]] <- tabpfn_one_target(
            pd, graph, splits, target, variant, seed, scenario, scale_n, rep_id
          )
        }
      }

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
