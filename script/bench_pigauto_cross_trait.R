#!/usr/bin/env Rscript
#
# script/bench_pigauto_cross_trait.R
#
# Branch-local pigauto-side diagnostic following the TabPFN reconciliation.
# Tests whether pigauto can close the same-row TabPFN gap by enabling stronger
# within-row cross-trait machinery, and whether a ridge-stabilized same-row
# model already explains much of the continuous AVONET gain.
#
# Smoke:
#   PIGAUTO_XTRAIT_OUT_STEM=bench_pigauto_cross_trait_smoke \
#   PIGAUTO_XTRAIT_SCALES=50 PIGAUTO_XTRAIT_REPS=1 \
#   PIGAUTO_XTRAIT_TARGETS=Mass \
#   PIGAUTO_XTRAIT_CONFIGS=default,trait_attention \
#     Rscript script/bench_pigauto_cross_trait.R

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

env_vec <- function(name, default) {
  x <- env_chr(name, paste(default, collapse = ","))
  trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
}

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_stem <- env_chr("PIGAUTO_XTRAIT_OUT_STEM", "bench_pigauto_cross_trait")
if (!grepl("^[A-Za-z0-9_.-]+$", out_stem)) {
  stop("PIGAUTO_XTRAIT_OUT_STEM must be a simple file stem, not a path.")
}
out_rds <- file.path(repo_root, "script", paste0(out_stem, ".rds"))
out_md <- file.path(repo_root, "script", paste0(out_stem, ".md"))

scales <- as.integer(env_vec("PIGAUTO_XTRAIT_SCALES", c("300")))
n_reps <- env_int("PIGAUTO_XTRAIT_REPS", 3L)
seed0 <- env_int("PIGAUTO_XTRAIT_SEED", 20260611L)
missing_frac <- env_num("PIGAUTO_XTRAIT_MISSING_FRAC", 0.25)
val_frac <- env_num("PIGAUTO_XTRAIT_VAL_FRAC", 0.5)
epochs <- env_int("PIGAUTO_XTRAIT_EPOCHS", 500L)
run_pigauto <- tolower(env_chr("PIGAUTO_XTRAIT_RUN_PIGAUTO", "true")) %in%
  c("true", "1", "yes", "y")
ridge_lambda <- env_num("PIGAUTO_XTRAIT_RIDGE_LAMBDA", 1)
ridge_min_sd <- env_num("PIGAUTO_XTRAIT_RIDGE_MIN_SD", 1e-6)
target_env <- env_vec("PIGAUTO_XTRAIT_TARGETS", character())
config_names <- env_vec(
  "PIGAUTO_XTRAIT_CONFIGS",
  c("default", "trait_attention", "relaxed_gate")
)

valid_configs <- c(
  "default",
  "trait_attention",
  "trait_attention_wide",
  "relaxed_gate"
)
if (!all(config_names %in% valid_configs)) {
  stop(
    "Unknown configs: ",
    paste(setdiff(config_names, valid_configs), collapse = ", ")
  )
}

metadata <- list(
  benchmark = "pigauto_cross_trait",
  scales = scales,
  n_reps = n_reps,
  seed0 = seed0,
  missing_frac = missing_frac,
  val_frac = val_frac,
  epochs = epochs,
  run_pigauto = run_pigauto,
  ridge_lambda = ridge_lambda,
  ridge_min_sd = ridge_min_sd,
  configs = config_names,
  out_stem = out_stem,
  commit = tryCatch(
    system("git rev-parse HEAD", intern = TRUE),
    error = function(e) "unknown"
  )
)

rmse_vec <- function(truth, pred) {
  sqrt(mean((truth - pred)^2, na.rm = TRUE))
}

cor_vec <- function(truth, pred) {
  ok <- is.finite(truth) & is.finite(pred)
  if (sum(ok) < 2L) {
    return(NA_real_)
  }
  stats::cor(truth[ok], pred[ok])
}

load_avonet_subset <- function(n, seed) {
  env <- new.env(parent = emptyenv())
  if (n <= 300L) {
    data("avonet300", package = "pigauto", envir = env)
    data("tree300", package = "pigauto", envir = env)
    df <- env$avonet300
    tree <- env$tree300
  } else {
    data("avonet_full", package = "pigauto", envir = env)
    data("tree_full", package = "pigauto", envir = env)
    df <- env$avonet_full
    tree <- env$tree_full
  }

  if (n > nrow(df)) {
    stop(
      "Requested n = ",
      n,
      " but AVONET source has only ",
      nrow(df),
      " rows."
    )
  }

  if (n < nrow(df)) {
    set.seed(seed)
    species <- sample(df$Species_Key, n)
    tree <- ape::keep.tip(tree, species)
    df <- df[match(tree$tip.label, df$Species_Key), , drop = FALSE]
  } else {
    df <- df[match(tree$tip.label, df$Species_Key), , drop = FALSE]
  }

  rownames(df) <- df$Species_Key
  df$Species_Key <- NULL
  list(df = df, tree = tree)
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

same_row_features <- function(
  pd,
  graph,
  splits,
  target_col,
  include_lappe = FALSE
) {
  X <- pd$X_scaled
  X[c(splits$val_idx, splits$test_idx)] <- NA_real_
  X <- X[, setdiff(seq_len(ncol(X)), target_col), drop = FALSE]
  parts <- list(traits = X)
  if (include_lappe) {
    coords <- graph$coords
    colnames(coords) <- paste0("lappe_", seq_len(ncol(coords)))
    parts$lappe <- coords
  }
  out <- do.call(cbind, parts)
  storage.mode(out) <- "double"
  out
}

impute_feature_means <- function(x_train, x_pred) {
  means <- colMeans(x_train, na.rm = TRUE)
  means[!is.finite(means)] <- 0
  fill <- function(x) {
    for (j in seq_len(ncol(x))) {
      bad <- !is.finite(x[, j])
      if (any(bad)) {
        x[bad, j] <- means[[j]]
      }
    }
    x
  }
  list(train = fill(x_train), pred = fill(x_pred))
}

ridge_predict <- function(x_train, y_train, x_pred, lambda) {
  scales <- apply(x_train, 2L, stats::sd)
  keep <- is.finite(scales) & scales > ridge_min_sd
  if (!any(keep)) {
    return(rep(mean(y_train), nrow(x_pred)))
  }
  x_train <- x_train[, keep, drop = FALSE]
  x_pred <- x_pred[, keep, drop = FALSE]
  centers <- colMeans(x_train)
  scales <- scales[keep]
  train_scaled <- sweep(sweep(x_train, 2L, centers, "-"), 2L, scales, "/")
  pred_scaled <- sweep(sweep(x_pred, 2L, centers, "-"), 2L, scales, "/")

  design_train <- cbind(.intercept = 1, train_scaled)
  design_pred <- cbind(.intercept = 1, pred_scaled)
  penalty <- diag(c(0, rep(lambda, ncol(train_scaled))))

  beta <- solve(
    crossprod(design_train) + penalty,
    crossprod(design_train, y_train)
  )
  as.numeric(design_pred %*% beta)
}

ridge_one_target <- function(
  pd,
  graph,
  splits,
  target,
  method,
  include_lappe,
  scale_n,
  rep_id
) {
  tm <- pd$trait_map[[target]]
  target_col <- tm$latent_cols[1L]
  n <- nrow(pd$X_scaled)
  y <- pd$X_scaled[, target_col]
  rows <- target_rows_from_split(splits, target_col, n)
  test_rows <- rows$test[is.finite(y[rows$test])]
  val_rows <- rows$val[is.finite(y[rows$val])]
  train_rows <- setdiff(which(is.finite(y)), c(val_rows, test_rows))

  base <- empty_row(method, scale_n, rep_id, target, "pending")
  base$n_train <- length(train_rows)
  base$n_val <- length(val_rows)
  base$n_test <- length(test_rows)

  if (length(train_rows) < 5L || length(test_rows) < 2L) {
    base$status <- "skipped_too_few_cells"
    return(base)
  }

  x <- same_row_features(
    pd,
    graph,
    splits,
    target_col,
    include_lappe = include_lappe
  )
  base$n_features <- ncol(x)
  filled <- impute_feature_means(
    x[train_rows, , drop = FALSE],
    x[test_rows, , drop = FALSE]
  )

  t0 <- proc.time()[["elapsed"]]
  pred <- tryCatch(
    {
      ridge_predict(
        filled$train,
        y[train_rows],
        filled$pred,
        lambda = ridge_lambda
      )
    },
    error = function(e) {
      structure(rep(NA_real_, length(test_rows)), error = conditionMessage(e))
    }
  )
  base$wall_sec <- proc.time()[["elapsed"]] - t0
  err <- attr(pred, "error", exact = TRUE)
  if (!is.null(err)) {
    base$status <- paste0("error: ", err)
    return(base)
  }

  base$rmse <- rmse_vec(y[test_rows], pred)
  base$pearson_r <- cor_vec(y[test_rows], pred)
  base$status <- "ok"
  base
}

empty_row <- function(method, scale_n, rep_id, target, status) {
  data.frame(
    method = method,
    scale_n = scale_n,
    rep = rep_id,
    trait = target,
    split = "test",
    n_train = NA_integer_,
    n_val = NA_integer_,
    n_test = NA_integer_,
    n_features = NA_integer_,
    rmse = NA_real_,
    pearson_r = NA_real_,
    coverage_95 = NA_real_,
    r_cal_bm = NA_real_,
    r_cal_gnn = NA_real_,
    r_cal_mean = NA_real_,
    wall_sec = NA_real_,
    status = status,
    stringsAsFactors = FALSE
  )
}

gate_map <- function(fit, pd, targets, field) {
  vals <- fit[[field]]
  if (is.null(vals)) {
    return(setNames(rep(NA_real_, length(targets)), targets))
  }
  out <- vapply(
    targets,
    function(target) {
      col <- pd$trait_map[[target]]$latent_cols[1L]
      if (length(vals) < col) NA_real_ else as.numeric(vals[[col]])
    },
    numeric(1L)
  )
  out
}

fit_config_args <- function(config) {
  switch(
    config,
    default = list(),
    trait_attention = list(
      use_trait_attention = TRUE,
      n_trait_heads = 2L,
      trait_embed_dim = 32L
    ),
    trait_attention_wide = list(
      use_trait_attention = TRUE,
      n_trait_heads = 4L,
      trait_embed_dim = 64L,
      hidden_dim = 128L
    ),
    relaxed_gate = list(
      gate_cap = 1.0,
      lambda_gate = 0,
      phylo_signal_gate = FALSE
    ),
    stop("Unknown config: ", config)
  )
}

standardize_eval <- function(
  ev,
  method,
  scale_n,
  rep_id,
  targets,
  wall_sec,
  gate_bm = NULL,
  gate_gnn = NULL,
  gate_mean = NULL
) {
  if (is.null(ev) || !nrow(ev)) {
    return(NULL)
  }
  keep <- ev$split == "test" & ev$trait %in% targets
  ev <- ev[keep, , drop = FALSE]
  if (!nrow(ev)) {
    return(NULL)
  }
  gates <- function(x) {
    if (is.null(x)) rep(NA_real_, nrow(ev)) else as.numeric(x[ev$trait])
  }
  data.frame(
    method = method,
    scale_n = scale_n,
    rep = rep_id,
    trait = ev$trait,
    split = ev$split,
    n_train = NA_integer_,
    n_val = NA_integer_,
    n_test = ev$n,
    n_features = NA_integer_,
    rmse = ev$rmse,
    pearson_r = ev$pearson_r,
    coverage_95 = ev$coverage_95,
    r_cal_bm = gates(gate_bm),
    r_cal_gnn = gates(gate_gnn),
    r_cal_mean = gates(gate_mean),
    wall_sec = wall_sec,
    status = "ok",
    stringsAsFactors = FALSE
  )
}

pigauto_config <- function(
  pd,
  tree,
  graph,
  splits,
  targets,
  config,
  seed,
  scale_n,
  rep_id,
  baseline = NULL
) {
  graph_fit <- graph
  graph_fit$D <- NULL
  invisible(gc(full = TRUE, verbose = FALSE))

  method <- paste0("pigauto_", config)
  t0 <- proc.time()[["elapsed"]]
  fit_result <- tryCatch(
    {
      args <- c(
        list(
          data = pd,
          tree = tree,
          splits = splits,
          baseline = baseline,
          graph = graph_fit,
          epochs = epochs,
          verbose = FALSE,
          seed = seed
        ),
        fit_config_args(config)
      )
      fit <- do.call(fit_pigauto, args)
      pred <- stats::predict(fit, return_se = TRUE, n_imputations = 1L)
      ev <- evaluate_imputation(pred, pd$X_scaled, splits)
      list(ok = TRUE, fit = fit, ev = ev)
    },
    error = function(e) {
      list(ok = FALSE, message = conditionMessage(e))
    }
  )
  wall <- proc.time()[["elapsed"]] - t0

  if (!isTRUE(fit_result$ok)) {
    return(do.call(
      rbind,
      lapply(targets, function(target) {
        row <- empty_row(
          method,
          scale_n,
          rep_id,
          target,
          paste0("error: ", fit_result$message)
        )
        row$wall_sec <- wall
        row
      })
    ))
  }

  standardize_eval(
    fit_result$ev,
    method,
    scale_n,
    rep_id,
    targets,
    wall,
    gate_bm = gate_map(fit_result$fit, pd, targets, "r_cal_bm"),
    gate_gnn = gate_map(fit_result$fit, pd, targets, "r_cal_gnn"),
    gate_mean = gate_map(fit_result$fit, pd, targets, "r_cal_mean")
  )
}

mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

write_summary <- function(results, metadata, path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("# Pigauto Cross-Trait Diagnostic\n", con)
  writeLines(
    sprintf("- Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    con
  )
  writeLines(sprintf("- Git commit: `%s`", metadata$commit), con)
  writeLines(
    sprintf("- Scales: `%s`", paste(metadata$scales, collapse = ", ")),
    con
  )
  writeLines(sprintf("- Replicates: `%d`", metadata$n_reps), con)
  writeLines(sprintf("- Run pigauto configs: `%s`", metadata$run_pigauto), con)
  writeLines(sprintf("- Ridge lambda: `%s`", metadata$ridge_lambda), con)
  writeLines(sprintf("- Ridge min SD: `%s`", metadata$ridge_min_sd), con)
  writeLines(
    sprintf(
      "- Pigauto configs: `%s`",
      paste(metadata$configs, collapse = ", ")
    ),
    con
  )
  writeLines("", con)

  writeLines("## Aim", con)
  writeLines("", con)
  writeLines(
    paste(
      "This benchmark tests whether pigauto-side cross-trait settings can close",
      "the continuous AVONET same-row TabPFN gap, and whether a",
      "ridge-stabilized same-row model already captures much of the signal."
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

  scored <- results[results$status == "ok", , drop = FALSE]
  if (!nrow(scored)) {
    writeLines("No scored rows were produced.", con)
    return(invisible(NULL))
  }

  writeLines("## Mean RMSE By Method", con)
  writeLines("", con)
  groups <- split(
    scored,
    interaction(scored$method, scored$scale_n, drop = TRUE)
  )
  agg <- do.call(
    rbind,
    lapply(groups, function(x) {
      data.frame(
        method = x$method[1L],
        scale_n = x$scale_n[1L],
        rmse = mean_or_na(x$rmse),
        pearson_r = mean_or_na(x$pearson_r),
        coverage_95 = mean_or_na(x$coverage_95),
        wall_sec = mean_or_na(x$wall_sec),
        stringsAsFactors = FALSE
      )
    })
  )
  agg <- agg[order(agg$scale_n, agg$rmse, agg$method), , drop = FALSE]
  writeLines(capture.output(print(agg, row.names = FALSE, digits = 4)), con)
  writeLines("", con)

  writeLines("## Mean RMSE By Trait", con)
  writeLines("", con)
  groups <- split(
    scored,
    interaction(scored$method, scored$scale_n, scored$trait, drop = TRUE)
  )
  trait_agg <- do.call(
    rbind,
    lapply(groups, function(x) {
      data.frame(
        method = x$method[1L],
        scale_n = x$scale_n[1L],
        trait = x$trait[1L],
        rmse = mean_or_na(x$rmse),
        pearson_r = mean_or_na(x$pearson_r),
        coverage_95 = mean_or_na(x$coverage_95),
        stringsAsFactors = FALSE
      )
    })
  )
  trait_agg <- trait_agg[
    order(trait_agg$scale_n, trait_agg$trait, trait_agg$rmse, trait_agg$method),
    ,
    drop = FALSE
  ]
  writeLines(
    capture.output(print(trait_agg, row.names = FALSE, digits = 4)),
    con
  )
  writeLines("", con)

  gate_rows <- scored[startsWith(scored$method, "pigauto_"), , drop = FALSE]
  if (nrow(gate_rows)) {
    writeLines("## Gate Audit", con)
    writeLines("", con)
    gate_groups <- split(
      gate_rows,
      interaction(
        gate_rows$method,
        gate_rows$scale_n,
        gate_rows$trait,
        drop = TRUE
      )
    )
    gate_agg <- do.call(
      rbind,
      lapply(gate_groups, function(x) {
        data.frame(
          method = x$method[1L],
          scale_n = x$scale_n[1L],
          trait = x$trait[1L],
          rmse = mean_or_na(x$rmse),
          r_cal_bm = mean_or_na(x$r_cal_bm),
          r_cal_gnn = mean_or_na(x$r_cal_gnn),
          r_cal_mean = mean_or_na(x$r_cal_mean),
          stringsAsFactors = FALSE
        )
      })
    )
    gate_agg <- gate_agg[
      order(gate_agg$scale_n, gate_agg$trait, gate_agg$method),
      ,
      drop = FALSE
    ]
    writeLines(
      capture.output(print(gate_agg, row.names = FALSE, digits = 4)),
      con
    )
  }
}

all_results <- list()

for (scale_n in scales) {
  for (rep_id in seq_len(n_reps)) {
    seed <- seed0 + scale_n * 100L + rep_id
    log_line("Preparing AVONET n=", scale_n, " rep=", rep_id)
    dat <- load_avonet_subset(scale_n, seed)
    pd <- preprocess_traits(dat$df, dat$tree, log_transform = TRUE)
    graph <- build_phylo_graph(dat$tree, k_eigen = "auto")
    splits <- make_missing_splits(
      pd$X_scaled,
      missing_frac = missing_frac,
      val_frac = val_frac,
      seed = seed,
      trait_map = pd$trait_map
    )

    continuous_targets <- names(Filter(
      function(tm) {
        identical(tm$type, "continuous")
      },
      pd$trait_map
    ))
    if (length(target_env)) {
      continuous_targets <- intersect(continuous_targets, target_env)
    }
    if (!length(continuous_targets)) {
      stop("No continuous targets selected.")
    }
    log_line("Targets: ", paste(continuous_targets, collapse = ", "))

    t0 <- proc.time()[["elapsed"]]
    baseline <- fit_baseline(pd, dat$tree, splits = splits, graph = graph)
    ev_bl <- evaluate_imputation(
      baseline$mu,
      pd$X_scaled,
      splits,
      trait_map = pd$trait_map
    )
    all_results[[length(all_results) + 1L]] <- standardize_eval(
      ev_bl,
      "baseline",
      scale_n,
      rep_id,
      continuous_targets,
      proc.time()[["elapsed"]] - t0
    )

    for (target in continuous_targets) {
      log_line(
        "Ridge same-row target=",
        target,
        " n=",
        scale_n,
        " rep=",
        rep_id
      )
      all_results[[length(all_results) + 1L]] <- ridge_one_target(
        pd,
        graph,
        splits,
        target,
        "ridge_same_row",
        FALSE,
        scale_n,
        rep_id
      )
      log_line(
        "Ridge same-row+LapPE target=",
        target,
        " n=",
        scale_n,
        " rep=",
        rep_id
      )
      all_results[[length(all_results) + 1L]] <- ridge_one_target(
        pd,
        graph,
        splits,
        target,
        "ridge_same_row_lappe",
        TRUE,
        scale_n,
        rep_id
      )
    }

    if (run_pigauto) {
      for (config in config_names) {
        log_line("Pigauto config=", config, " n=", scale_n, " rep=", rep_id)
        all_results[[length(all_results) + 1L]] <- pigauto_config(
          pd,
          dat$tree,
          graph,
          splits,
          continuous_targets,
          config,
          seed,
          scale_n,
          rep_id,
          baseline = baseline
        )
      }
    }

    results <- do.call(rbind, all_results)
    saveRDS(list(results = results, metadata = metadata), out_rds)
    write_summary(results, metadata, out_md)
  }
}

results <- do.call(rbind, all_results)
saveRDS(list(results = results, metadata = metadata), out_rds)
write_summary(results, metadata, out_md)
log_line("Wrote ", out_rds)
log_line("Wrote ", out_md)
