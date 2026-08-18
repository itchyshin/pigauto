#!/usr/bin/env Rscript
# Stage A active-recovery campaign helpers.  This file is deliberately
# independent of the public API: it is evidence machinery, not a new method.

active_recovery_entropy <- function(p) {
  p <- pmin(pmax(as.numeric(p), 1e-12), 1 - 1e-12)
  -p * log(p) - (1 - p) * log(1 - p)
}

active_recovery_source_sha <- function() {
  out <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
                  error = function(e) NA_character_)
  if (length(out) == 1L) out else NA_character_
}

active_recovery_draw <- function(n, lambda, family, seed) {
  stopifnot(n >= 20L, lambda > 0, lambda <= 1,
            family %in% c("continuous", "binary"))
  set.seed(seed)
  tree <- ape::rtree(n)
  R <- pigauto:::phylo_cor_matrix(tree)
  Sigma <- lambda * R + (1 - lambda) * diag(n)
  L <- chol(Sigma + diag(1e-10, n))
  z <- as.numeric(t(L) %*% stats::rnorm(n))
  names(z) <- tree$tip.label
  truth <- if (identical(family, "continuous")) {
    data.frame(trait = z, row.names = tree$tip.label)
  } else {
    factor(ifelse(z > stats::median(z), "yes", "no"), levels = c("no", "yes"))
  }
  if (identical(family, "binary")) {
    truth <- data.frame(trait = truth, row.names = tree$tip.label)
  }
  list(tree = tree, truth = truth, latent = z, R = R)
}

active_recovery_split <- function(species, seed) {
  n <- length(species)
  stopifnot(n >= 20L)
  set.seed(seed)
  perm <- sample(species, n)
  n_initial <- floor(0.60 * n)
  n_candidate <- floor(0.20 * n)
  list(initial = perm[seq_len(n_initial)],
       candidate = perm[n_initial + seq_len(n_candidate)],
       test = perm[(n_initial + n_candidate + 1L):n])
}

active_recovery_mask <- function(truth, initial) {
  out <- truth
  out[-match(initial, rownames(out)), "trait"] <- NA
  out
}

active_recovery_data_hash <- function(data) {
  raw <- serialize(data, connection = NULL, version = 2L)
  # A deterministic two-part content receipt, sufficient to detect a wrong
  # treatment application without adding a package dependency to a script.
  x <- as.numeric(as.integer(raw))
  mod <- 2147483647
  sprintf("%010.0f-%010.0f", sum(x) %% mod,
          sum((seq_along(x) %% 65521) * x) %% mod)
}

active_recovery_policy_scores <- function(result, family, candidates) {
  spp <- result$data$species_names
  idx <- match(candidates, spp)
  if (anyNA(idx)) stop("candidate species not found in fitted data", call. = FALSE)
  if (identical(family, "continuous")) {
    x <- result$prediction$se[idx, 1L]
  } else {
    p <- result$prediction$probabilities$trait[idx]
    x <- active_recovery_entropy(p)
  }
  stats::setNames(as.numeric(x), candidates)
}

active_recovery_select <- function(result, family, candidates, policy,
                                   random_seed = NULL) {
  stopifnot(policy %in% c("active", "random", "uncertainty"))
  if (identical(policy, "random")) {
    if (!is.null(random_seed)) set.seed(random_seed)
    return(sample(candidates, 1L))
  }
  if (identical(policy, "uncertainty")) {
    scores <- active_recovery_policy_scores(result, family, candidates)
    return(names(scores)[which.max(scores)])
  }
  suggested <- pigauto::suggest_next_observation(
    result, top_n = length(result$data$species_names), by = "cell",
    types = family
  )
  suggested <- suggested[suggested$trait == "trait" &
                           suggested$species %in% candidates, , drop = FALSE]
  if (!nrow(suggested)) stop("no active recommendation in candidate set", call. = FALSE)
  suggested$species[which.max(suggested$delta)]
}

active_recovery_restore_one <- function(masked, truth, species) {
  before <- sum(!is.na(masked$trait))
  out <- masked
  i <- match(species, rownames(out))
  if (is.na(i) || !is.na(out$trait[i])) stop("restoration target must be missing", call. = FALSE)
  out$trait[i] <- truth$trait[i]
  if (sum(!is.na(out$trait)) != before + 1L) stop("restoration was not exactly one cell", call. = FALSE)
  out
}

active_recovery_score <- function(result, truth, test, family) {
  idx <- match(test, result$data$species_names)
  y <- truth$trait[match(test, rownames(truth))]
  if (identical(family, "continuous")) {
    pred <- as.numeric(result$prediction$imputed$trait[idx])
    mse <- mean((pred - y)^2)
    list(primary = mse / stats::var(y), mse = mse, rmse = sqrt(mse),
         log_loss = NA_real_, accuracy = NA_real_)
  } else {
    p <- as.numeric(result$prediction$probabilities$trait[idx])
    y01 <- as.integer(y == levels(truth$trait)[2L])
    brier <- mean((p - y01)^2)
    list(primary = brier, mse = NA_real_, rmse = NA_real_,
         log_loss = -mean(y01 * log(pmax(p, 1e-12)) +
                           (1 - y01) * log(pmax(1 - p, 1e-12))),
         accuracy = mean((p >= 0.5) == as.logical(y01)))
  }
}

active_recovery_fit <- function(data, tree, seed, epochs) {
  pigauto::impute(data, tree, epochs = as.integer(epochs), seed = as.integer(seed),
                  n_imputations = 1L, verbose = FALSE)
}

active_recovery_one <- function(n, lambda, family, replicate, master_seed,
                                epochs = 100L) {
  dgp_seed <- master_seed + replicate * 1000L
  draw <- active_recovery_draw(n, lambda, family, dgp_seed)
  split <- active_recovery_split(draw$tree$tip.label, dgp_seed + 1L)
  initial <- active_recovery_mask(draw$truth, split$initial)
  initial_hash <- active_recovery_data_hash(initial)
  fit_seed <- dgp_seed + 2L
  t0 <- proc.time()[["elapsed"]]
  initial_fit <- active_recovery_fit(initial, draw$tree, fit_seed, epochs)
  initial_score <- active_recovery_score(initial_fit, draw$truth, split$test, family)
  baseline_row <- c(list(policy = "no_acquisition", selected_species = NA_character_,
                         restored_value = NA_character_, initial_hash = initial_hash,
                         changed_hash = initial_hash, n_initial = length(split$initial),
                         n_candidate = length(split$candidate), n_test = length(split$test),
                         elapsed_s = proc.time()[["elapsed"]] - t0), initial_score)
  rows <- c(list(baseline_row), lapply(c("active", "random", "uncertainty"), function(policy) {
    selected <- active_recovery_select(initial_fit, family, split$candidate, policy,
                                       random_seed = dgp_seed + match(policy, c("active", "random", "uncertainty")))
    changed <- active_recovery_restore_one(initial, draw$truth, selected)
    changed_hash <- active_recovery_data_hash(changed)
    refit <- active_recovery_fit(changed, draw$tree, fit_seed, epochs)
    score <- active_recovery_score(refit, draw$truth, split$test, family)
    c(list(policy = policy, selected_species = selected, restored_value = as.character(draw$truth[selected, "trait"]),
           initial_hash = initial_hash, changed_hash = changed_hash,
           n_initial = length(split$initial), n_candidate = length(split$candidate),
           n_test = length(split$test), elapsed_s = proc.time()[["elapsed"]] - t0), score)
  }))
  data.frame(n = n, lambda = lambda, family = family, replicate = replicate,
             source_sha = active_recovery_source_sha(),
             status = "ok", error = NA_character_,
             do.call(rbind, lapply(rows, as.data.frame)), row.names = NULL,
             stringsAsFactors = FALSE)
}
