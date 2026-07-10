# pigauto v0.10.0 known-DGP MI validation: configuration and data generation.
#
# Run directly to create a deterministic manifest, or source from 1_run.R.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

.arg_value <- function(args, key, default = NULL) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^--", key, "="), "", hit[[length(hit)]])
}

.env_value <- function(name, default = NULL) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

.as_int <- function(x, name, lower = 1L) {
  value <- suppressWarnings(as.integer(x))
  if (length(value) != 1L || is.na(value) || value < lower) {
    stop(name, " must be an integer >= ", lower, "; got ", x, ".",
         call. = FALSE)
  }
  value
}

.csv_values <- function(x, allowed, name) {
  values <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  bad <- setdiff(values, allowed)
  if (length(bad) > 0L) {
    stop(name, " contains unsupported values: ", paste(bad, collapse = ", "),
         call. = FALSE)
  }
  unique(values)
}

mi_validation_config <- function(args = commandArgs(trailingOnly = TRUE)) {
  profile <- .arg_value(
    args, "profile", .env_value("PIGAUTO_MI_PROFILE", "full")
  )
  if (!(profile %in% c("smoke", "pilot", "full"))) {
    stop("profile must be one of smoke, pilot, or full.", call. = FALSE)
  }

  defaults <- switch(
    profile,
    smoke = list(reps = 1L, m = 2L, epochs = 2L),
    pilot = list(reps = 10L, m = 50L, epochs = 500L),
    full  = list(reps = 500L, m = 50L, epochs = 500L)
  )
  root <- normalizePath(getwd(), mustWork = TRUE)
  default_output <- file.path(root, "script", "mi-validation-v010", "results",
                              profile)

  get_setting <- function(arg, env, default) {
    .arg_value(args, arg, .env_value(env, as.character(default)))
  }

  list(
    profile = profile,
    base_seed = .as_int(
      get_setting("base-seed", "PIGAUTO_MI_BASE_SEED", 20260710L),
      "base_seed", 1L
    ),
    reps = .as_int(
      get_setting("reps", "PIGAUTO_MI_REPS", defaults$reps), "reps", 1L
    ),
    m = .as_int(
      get_setting("m", "PIGAUTO_MI_M", defaults$m), "m", 2L
    ),
    epochs = .as_int(
      get_setting("epochs", "PIGAUTO_MI_EPOCHS", defaults$epochs),
      "epochs", 1L
    ),
    missing_fraction = 0.30,
    dgps = .csv_values(
      get_setting("dgps", "PIGAUTO_MI_DGPS", "lm,glm,lmer"),
      c("lm", "glm", "lmer"), "dgps"
    ),
    regimes = .csv_values(
      get_setting("regimes", "PIGAUTO_MI_REGIMES", "phylogeny,auxiliary"),
      c("phylogeny", "auxiliary"), "regimes"
    ),
    output_dir = normalizePath(
      get_setting("output", "PIGAUTO_MI_OUTPUT", default_output),
      mustWork = FALSE
    ),
    eval_every = .as_int(
      get_setting("eval-every", "PIGAUTO_MI_EVAL_EVERY", 50L),
      "eval_every", 1L
    ),
    patience = .as_int(
      get_setting("patience", "PIGAUTO_MI_PATIENCE", 50L),
      "patience", 1L
    )
  )
}

atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = ".partial-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, version = 3)
  if (!file.rename(tmp, path)) {
    stop("Could not atomically move ", tmp, " to ", path, ".", call. = FALSE)
  }
  invisible(path)
}

make_manifest <- function(config) {
  grid <- expand.grid(
    dgp = config$dgps,
    regime = config$regimes,
    replicate = seq_len(config$reps),
    stringsAsFactors = FALSE
  )
  grid <- grid[order(match(grid$dgp, config$dgps),
                     match(grid$regime, config$regimes),
                     grid$replicate), , drop = FALSE]
  rownames(grid) <- NULL
  grid$task_id <- seq_len(nrow(grid))
  dgp_code <- match(grid$dgp, c("lm", "glm", "lmer"))
  regime_code <- match(grid$regime, c("phylogeny", "auxiliary"))
  grid$seed <- as.integer(
    config$base_seed + dgp_code * 1000000L + regime_code * 100000L +
      grid$replicate
  )
  grid$n_species <- c(lm = 250L, glm = 400L, lmer = 100L)[grid$dgp]
  grid$n_observations <- ifelse(grid$dgp == "lmer", 4L * grid$n_species,
                                grid$n_species)
  grid$missing_fraction <- config$missing_fraction
  grid$m <- config$m
  grid$epochs <- config$epochs
  grid$profile <- config$profile
  grid[, c("task_id", "dgp", "regime", "replicate", "seed", "n_species",
           "n_observations", "missing_fraction", "m", "epochs", "profile")]
}

write_manifest <- function(config, overwrite = FALSE) {
  path <- file.path(config$output_dir, "manifest.rds")
  manifest <- make_manifest(config)
  if (file.exists(path) && !overwrite) {
    existing <- readRDS(path)
    if (!identical(existing, manifest)) {
      stop("Existing manifest does not match the requested configuration: ",
           path, ". Use a different --output path or rerun 0_prepare.R to ",
           "replace it explicitly.", call. = FALSE)
    }
    return(existing)
  }
  atomic_save_rds(manifest, path)
  manifest
}

.zscore <- function(x) {
  out <- as.numeric(scale(x))
  if (any(!is.finite(out))) stop("Cannot standardize a constant vector.")
  out
}

.simulate_x <- function(phylo_value, z, regime, noise_sd = 0.25) {
  if (regime == "phylogeny") {
    raw <- 1.50 * phylo_value + 0.20 * z + stats::rnorm(length(z), sd = noise_sd)
  } else {
    raw <- 0.25 * phylo_value + 1.00 * (z^2 - 1) +
      stats::rnorm(length(z), sd = noise_sd)
  }
  .zscore(raw)
}

.mask_mar_x <- function(data, y_numeric, fraction, seed) {
  n_missing <- as.integer(round(nrow(data) * fraction))
  score <- 0.85 * .zscore(y_numeric) + 0.55 * .zscore(data$z)
  weights <- stats::plogis(score)
  set.seed(seed + 70000L)
  missing_rows <- sample.int(nrow(data), n_missing, replace = FALSE,
                             prob = weights)
  observed <- data
  observed$x[missing_rows] <- NA_real_
  list(data = observed, rows = missing_rows)
}

simulate_validation_dgp <- function(dgp, regime, seed,
                                    missing_fraction = 0.30) {
  stopifnot(dgp %in% c("lm", "glm", "lmer"),
            regime %in% c("phylogeny", "auxiliary"))
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("The validation harness requires package 'ape'.", call. = FALSE)
  }
  set.seed(seed)

  if (dgp %in% c("lm", "glm")) {
    n <- if (dgp == "lm") 250L else 400L
    tree <- ape::rtree(n)
    species <- tree$tip.label
    # Index by tip name rather than relying on an incidental vector order.
    set.seed(seed + 1L)
    phylo_named <- ape::rTraitCont(tree, model = "BM")
    phylo_value <- .zscore(phylo_named[species])
    z <- stats::rnorm(n)
    x <- .simulate_x(phylo_value, z, regime)

    if (dgp == "lm") {
      y <- 1 + 0.60 * x - 0.40 * z + stats::rnorm(n, sd = 1)
      full <- data.frame(x = x, y = y, z = z, row.names = species)
      y_numeric <- y
      truth <- c(x = 0.60, z = -0.40)
    } else {
      probability <- stats::plogis(-0.20 + 0.80 * x - 0.50 * z)
      y_numeric <- stats::rbinom(n, 1L, probability)
      y <- factor(ifelse(y_numeric == 1L, "yes", "no"),
                  levels = c("no", "yes"))
      full <- data.frame(x = x, y = y, z = z, row.names = species)
      truth <- c(x = 0.80, z = -0.50)
    }
    masked <- .mask_mar_x(full, y_numeric, missing_fraction, seed)
    observed <- masked$data
    traits <- observed[, c("x", "y"), drop = FALSE]
    species_col <- NULL
    truth_vc <- numeric(0)
  } else {
    n_species <- 100L
    n_per_species <- 4L
    tree <- ape::rtree(n_species)
    species_names <- tree$tip.label
    set.seed(seed + 1L)
    phylo_named <- ape::rTraitCont(tree, model = "BM")
    phylo_species <- .zscore(phylo_named[species_names])
    species <- rep(species_names, each = n_per_species)
    z <- stats::rnorm(length(species))
    x <- .simulate_x(rep(phylo_species, each = n_per_species), z, regime,
                     noise_sd = 0.30)
    random_intercept <- stats::rnorm(n_species, sd = 0.80)
    names(random_intercept) <- species_names
    y <- 1 + 0.60 * x - 0.40 * z + random_intercept[species] +
      stats::rnorm(length(species), sd = 1)
    full <- data.frame(species = species, x = x, y = y, z = z,
                       stringsAsFactors = FALSE)
    masked <- .mask_mar_x(full, y, missing_fraction, seed)
    observed <- masked$data
    traits <- observed[, c("species", "x", "y"), drop = FALSE]
    species_col <- "species"
    truth <- c(x = 0.60, z = -0.40)
    truth_vc <- c(tau2 = 0.80^2, sigma2 = 1)
  }

  covariates <- data.frame(z = full$z)
  list(
    dgp = dgp,
    regime = regime,
    seed = seed,
    tree = tree,
    full = full,
    observed = observed,
    traits = traits,
    covariates = covariates,
    species_col = species_col,
    missing_rows = masked$rows,
    truth = truth,
    truth_vc = truth_vc
  )
}

fit_downstream <- function(data, dgp) {
  if (dgp == "lm") {
    stats::lm(y ~ x + z, data = data)
  } else if (dgp == "glm") {
    stats::glm(y ~ x + z, data = data, family = stats::binomial())
  } else {
    if (!requireNamespace("lme4", quietly = TRUE)) {
      stop("The lmer DGP requires package 'lme4'.", call. = FALSE)
    }
    lme4::lmer(
      y ~ x + z + (1 | species), data = data, REML = FALSE,
      control = lme4::lmerControl(
        optimizer = "bobyqa",
        check.conv.singular = "ignore"
      )
    )
  }
}

complete_imputed_dataset <- function(data, dgp_object) {
  data$z <- dgp_object$full$z
  data
}

if (sys.nframe() == 0L) {
  config <- mi_validation_config()
  manifest <- write_manifest(config, overwrite = TRUE)
  cat("Wrote", nrow(manifest), "tasks to",
      file.path(config$output_dir, "manifest.rds"), "\n")
}
