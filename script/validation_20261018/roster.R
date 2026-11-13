# Contract-based expected result roster and independent payload audit.
# This file contains no fit or simulation calls. Expected keys come from the
# planned job fields and frozen generator/adapter contracts, never observed
# successful result rows.

.vr_roster_value <- function(job, name, default = NULL) {
  if (is.data.frame(job)) {
    if (!(name %in% names(job))) return(default)
    value <- job[[name]][[1L]]
  } else value <- job[[name]]
  if (is.null(value) || !length(value) ||
      (is.atomic(value) && all(is.na(value)))) default else value
}

.vr_roster_bind_fill <- function(xs) {
  if (!length(xs)) return(data.frame())
  nms <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(x) {
    for (nm in setdiff(nms, names(x))) x[[nm]] <- NA
    x[, nms, drop = FALSE]
  })
  do.call(rbind, xs)
}

.vr_roster_traits <- function(scenario) {
  scenario <- tolower(as.character(scenario))
  if (grepl("^types_mixed", scenario)) {
    list(numeric = c("c1", "c2", "cnt", "prp"),
         discrete = c("bin", "ord", "cat3"))
  } else if (grepl("^(bm|ou)_mixed", scenario)) {
    list(numeric = paste0("c", 1:4), discrete = c("bin", "cat3"))
  } else NULL
}

.vr_roster_bien_traits <- function(traits_path, tree_path, subset_n = 300L) {
  if (!nzchar(traits_path) || !nzchar(tree_path) ||
      !file.exists(traits_path) || !file.exists(tree_path))
    return(list(available = FALSE, reason = "explicit BIEN cache paths are absent or do not exist"))
  trait_means <- tryCatch(readRDS(traits_path), error = identity)
  tree_pkg <- tryCatch(readRDS(tree_path), error = identity)
  if (inherits(trait_means, "error") || inherits(tree_pkg, "error") ||
      !is.list(trait_means) || is.null(names(trait_means)) ||
      !is.list(tree_pkg) || is.null(tree_pkg$scenario.3) ||
      !inherits(tree_pkg$scenario.3, "phylo"))
    return(list(available = FALSE, reason = "BIEN cache does not match the registered cache schema"))
  norm <- function(x) gsub("_", " ", gsub("[[:space:]]+", " ", trimws(as.character(x))))
  species_lists <- lapply(trait_means, function(d) {
    if (is.data.frame(d) && "species" %in% names(d)) norm(d$species) else character()
  })
  all_species <- Reduce(union, species_lists, init = character())
  observed <- matrix(FALSE, nrow = length(all_species), ncol = length(trait_means),
                     dimnames = list(NULL, names(trait_means)))
  for (nm in names(trait_means)) {
    d <- trait_means[[nm]]
    if (!is.data.frame(d) || !("species" %in% names(d))) next
    value_col <- intersect(c("mean_value", "trait_value", "value"), names(d))[1L]
    if (is.na(value_col)) next
    idx <- match(all_species, norm(d$species))
    observed[, nm] <- is.finite(suppressWarnings(as.numeric(d[[value_col]][idx])))
  }
  eligible <- all_species[rowSums(observed) >= 2L]
  matched <- intersect(norm(tree_pkg$scenario.3$tip.label), eligible)
  if (length(matched) < as.integer(subset_n))
    return(list(available = FALSE, reason = sprintf(
      "BIEN cache has %d matched species with at least two traits; %d required",
      length(matched), as.integer(subset_n))))
  list(available = TRUE, traits = names(trait_means))
}

.vr_is_nonnegative_integer <- function(x) {
  is.numeric(x) && length(x) > 0L && !anyNA(x) && all(is.finite(x)) &&
    all(x >= 0) && all(x == floor(x))
}

.vr_calibration_payload_status <- function(hit, n_bins, ece) {
  if (!nrow(hit)) return("missing")
  required <- c("bin", "bin_lower", "bin_upper", "n", "mean_confidence",
                "mean_correct", "absolute_gap", "status")
  if (!all(required %in% names(hit)) || nrow(hit) != n_bins ||
      anyNA(hit$bin) || anyDuplicated(hit$bin) || !setequal(as.integer(hit$bin), seq_len(n_bins)) ||
      !.vr_is_nonnegative_integer(hit$n) ||
      anyNA(hit$status) ||
      any(!hit$status %in% c("ok", "failed", "unavailable")) ||
      length(unique(hit$status)) != 1L) return("failed")
  ord <- order(as.integer(hit$bin)); hit <- hit[ord, , drop = FALSE]
  edges <- seq(0, 1, length.out = n_bins + 1L)
  if (any(!is.finite(hit$bin_lower)) || any(!is.finite(hit$bin_upper)) ||
      any(abs(hit$bin_lower - edges[-length(edges)]) > 1e-10) ||
      any(abs(hit$bin_upper - edges[-1L]) > 1e-10)) return("failed")
  positive <- hit$n > 0L
  if (any(positive & (!is.finite(hit$mean_confidence) | hit$mean_confidence < 0 |
                      hit$mean_confidence > 1 | !is.finite(hit$mean_correct) |
                      hit$mean_correct < 0 | hit$mean_correct > 1 |
                      !is.finite(hit$absolute_gap) | hit$absolute_gap < 0 |
                      hit$absolute_gap > 1 |
                      abs(hit$absolute_gap - abs(hit$mean_confidence - hit$mean_correct)) > 1e-10)) ||
      any(!positive & (is.finite(hit$mean_confidence) | is.finite(hit$mean_correct) |
                       is.finite(hit$absolute_gap)))) return("failed")
  if (!is.data.frame(ece) || nrow(ece) != 1L || !"status" %in% names(ece) ||
      is.na(ece$status[[1L]]) || !ece$status[[1L]] %in% c("ok", "failed", "unavailable")) return("failed")
  if (ece$status[[1L]] != hit$status[[1L]]) return("failed")
  if (ece$status[[1L]] == "ok") {
    if (!"estimate" %in% names(ece) || !is.finite(ece$estimate[[1L]]) ||
        ece$estimate[[1L]] < 0 || ece$estimate[[1L]] > 1) return("failed")
  } else if ("estimate" %in% names(ece) && is.finite(ece$estimate[[1L]])) {
    return("failed")
  }
  if (hit$status[[1L]] == "failed") return("failed")
  if (hit$status[[1L]] == "unavailable") return("unavailable")
  "present"
}

.vr_posterior_interval_status <- function(hit) {
  if (!nrow(hit)) return("missing")
  required <- c("interval_status", "n_expected", "n_finite_scored",
                "n_missing_or_nonfinite", "coverage", "mean_width")
  if (!all(required %in% names(hit)) || nrow(hit) != 1L ||
      is.na(hit$interval_status[[1L]]) ||
      !hit$interval_status[[1L]] %in% c("ok", "failed", "unavailable") ||
      !.vr_is_nonnegative_integer(hit$n_expected) ||
      !.vr_is_nonnegative_integer(hit$n_finite_scored) ||
      !.vr_is_nonnegative_integer(hit$n_missing_or_nonfinite) ||
      hit$n_finite_scored[[1L]] + hit$n_missing_or_nonfinite[[1L]] != hit$n_expected[[1L]])
    return("failed")
  status <- hit$interval_status[[1L]]
  if (status == "unavailable") {
    valid <- hit$n_expected[[1L]] == 0L && hit$n_finite_scored[[1L]] == 0L &&
      hit$n_missing_or_nonfinite[[1L]] == 0L && !is.finite(hit$coverage[[1L]]) &&
      !is.finite(hit$mean_width[[1L]])
    return(if (valid) "unavailable" else "failed")
  }
  if (status == "ok") {
    valid <- hit$n_expected[[1L]] > 0L && hit$n_finite_scored[[1L]] == hit$n_expected[[1L]] &&
      hit$n_missing_or_nonfinite[[1L]] == 0L && is.finite(hit$coverage[[1L]]) &&
      hit$coverage[[1L]] >= 0 && hit$coverage[[1L]] <= 1 &&
      is.finite(hit$mean_width[[1L]]) && hit$mean_width[[1L]] >= 0
    return(if (valid) "present" else "failed")
  }
  "failed"
}

validation_compare_posterior_denominators <- function(intervals, expected) {
  required <- c("job_id", "trait", "n_expected")
  if (!is.data.frame(intervals) || !all(required %in% names(intervals)) ||
      !all(c("job_id", "trait", "n_expected_truth") %in% names(expected)))
    stop("intervals and expected denominator tables have invalid schemas")
  key <- function(x) paste(x$job_id, x$trait, sep = "\034")
  if (anyDuplicated(key(intervals)) || anyDuplicated(key(expected)))
    stop("denominator tables contain duplicate job/trait keys")
  i <- match(key(expected), key(intervals))
  observed <- intervals$n_expected[i]
  data.frame(job_id = expected$job_id, trait = expected$trait,
    n_expected_truth = expected$n_expected_truth, n_expected_receipt = observed,
    status = ifelse(is.na(i), "missing", ifelse(
      is.finite(observed) & observed == expected$n_expected_truth, "match", "mismatch")),
    stringsAsFactors = FALSE)
}

.vr_roster_for_job <- function(job, root = ".", bien_traits_path = "",
                               bien_tree_path = "", n_bins = 10L,
                               regime_contract = NULL, bien_contract = NULL) {
  adapter <- tolower(as.character(.vr_roster_value(job, "adapter", "")))
  scenario <- as.character(.vr_roster_value(job, "scenario", ""))
  method_estimands <- list()
  metrics <- list()
  issue <- NA_character_
  add <- function(method, estimand) {
    method_estimands[[length(method_estimands) + 1L]] <<-
      data.frame(method = as.character(method), estimand = as.character(estimand),
                 stringsAsFactors = FALSE)
  }
  add_many <- function(method, estimands) for (estimand in estimands) add(method, estimand)

  if (adapter %in% c("prediction", "gnn_off", "prediction_reuse")) {
    traits <- .vr_roster_traits(scenario)
    if (is.null(traits)) {
      issue <- paste0("no frozen trait schema for prediction scenario: ", scenario)
    } else {
      arms <- c("gnn_off", "floor")
      if (isTRUE(.vr_roster_value(job, "rphylopars", FALSE)))
        arms <- c(arms, "gnn_off_rphylopars")
      for (arm in arms) {
        for (trait in traits$numeric) {
          add_many(arm, paste(trait,
            c("zRMSE", "width", "interval_score", "conformal_coverage"), sep = "_"))
        }
        for (trait in traits$discrete) {
          add_many(arm, paste(trait,
            c("accuracy", "macroF1", "n_classes_scored", "brier", "ece"), sep = "_"))
        }
      }
      # Bins are a receipt payload requirement, not separate roster metrics.
      # The audit checks the full 1:n_bins set for each method × trait.
      metrics$calibration_bins <- expand.grid(
        method = arms, trait = traits$discrete,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    }
  } else if (adapter %in% c("familyzi", "family_zi", "zi_count", "zi")) {
    traits <- paste0("zi", 1:3)
    for (method in c("mean", "baseline", "pigauto_gnn_off"))
      add_many(method, unlist(lapply(traits, function(trait)
        paste(trait, c("brier", "zero_accuracy", "rmse"), sep = "_"))))
    add_many("pigauto_gnn_off", paste(traits, "positive_magnitude_latent_rmse", sep = "_"))
  } else if (adapter %in% c("composition", "multiproportion", "multi_proportion")) {
    for (method in c("mean", "baseline", "pigauto_gnn_off"))
      add_many(method, paste("comp", c("aitchison", "rmse_clr", "simplex_mae"), sep = "_"))
  } else if (adapter %in% c("workflowmultiobs", "multiobs", "workflow_multiobs")) {
    for (method in c("no_covariates", "covariates")) add(method, "masked_obs_rmse")
  } else if (adapter == "multiobs_mixed") {
    for (method in c("pigauto_hard", "pigauto_soft"))
      add_many(method, c("x1_rmse", "x2_rmse", "y_accuracy", "z_accuracy"))
  } else if (adapter == "covariate") {
    for (variant in c("linear", "nonlinear")) for (method in
      c("pigauto_default", "gnn_on", "gnn_on_covariates"))
      add_many(paste(variant, method, sep = "__"), paste0("trait", 1:3, "_rmse"))
  } else if (adapter %in% c("tree", "trees", "tree_uncertainty")) {
    add_many("multi_tree", paste0(c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length"), "_rmse"))
  } else if (adapter %in% c("posterior", "mi_gls", "posterior_mi")) {
    rid <- suppressWarnings(as.integer(.vr_roster_value(job, "regime_id", sub("^regime_", "", scenario))))
    if (length(rid) != 1L || is.na(rid)) {
      issue <- "posterior roster needs a valid regime_id"
    } else {
      ix <- if (is.null(regime_contract)) NA_integer_ else
        match(rid, regime_contract$regime_id)
      if (is.na(ix)) issue <- if (is.null(regime_contract))
        "posterior regime contract file is absent" else paste("unknown posterior regime_id", rid) else {
        missing <- regime_contract$missing[[ix]]
        add("complete_data", "beta_x")
        add("posterior_full", "beta_x")
        if (identical(missing, "both")) add("posterior_none", "beta_x")
        metrics$posterior_cell_interval <- data.frame(
          trait = if (identical(missing, "both")) c("x", "y") else "x",
          stringsAsFactors = FALSE)
      }
    }
  } else if (adapter %in% c("analysis", "analysislm", "analysis_lm", "analysisglm", "analysis_glm", "analysislmer", "analysis_lmer")) {
    model <- tolower(as.character(.vr_roster_value(job, "model", "")))
    if (!model %in% c("lm", "glm", "lmer")) {
      found <- regmatches(scenario, regexpr("^(lm|glm|lmer)(?=_)" , scenario, perl = TRUE))
      model <- if (length(found)) found else ""
    }
    if (!model %in% c("lm", "glm", "lmer")) issue <- "analysis roster needs model=lm, glm, or lmer" else {
      add("complete_data", "beta_x")
      add(paste0("analysis_", model), "beta_x")
    }
  } else if (adapter == "real") {
    dataset <- tolower(scenario)
    if (dataset == "avonet") {
      add_many("pigauto_gnn_off", paste0(
        c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length"), "_rmse"))
    } else if (dataset == "bien") {
      bien <- bien_contract %||% .vr_roster_bien_traits(bien_traits_path, bien_tree_path)
      if (!isTRUE(bien$available)) {
        # An absent local cache does not imply that a remote worker had no
        # cache. The identity records only cache hashes, not the trait names;
        # refuse to invent an expected roster from successful result rows.
        issue <- bien$reason
      } else add_many("pigauto_gnn_off", paste0(bien$traits, "_rmse"))
    } else if (dataset %in% c("pantheria", "fishbase")) {
      add("unavailable", "attempt")
    } else issue <- paste0("no real-data roster contract for dataset: ", dataset)
  } else {
    issue <- paste0("no frozen expected-key contract for adapter: ", adapter)
  }

  keys <- if (length(method_estimands)) do.call(rbind, method_estimands) else
    data.frame(method = character(), estimand = character(), stringsAsFactors = FALSE)
  keys <- unique(keys)
  if (nrow(keys)) {
    keys$job_id <- as.character(.vr_roster_value(job, "job_id", NA_character_))
    keys$scenario_id <- as.character(.vr_roster_value(job, "scenario_id", NA_character_))
    keys$slice <- as.character(.vr_roster_value(job, "slice", NA_character_))
    keys$adapter <- as.character(.vr_roster_value(job, "adapter", NA_character_))
    keys$rep <- as.integer(.vr_roster_value(job, "rep", NA_integer_))
    keys <- keys[, c("job_id", "scenario_id", "slice", "adapter", "rep", "method", "estimand")]
  }
  if (length(metrics)) for (nm in names(metrics)) {
    metrics[[nm]]$job_id <- as.character(.vr_roster_value(job, "job_id", NA_character_))
    metrics[[nm]]$scenario_id <- as.character(.vr_roster_value(job, "scenario_id", NA_character_))
    metrics[[nm]]$rep <- as.integer(.vr_roster_value(job, "rep", NA_integer_))
    metrics[[nm]]$adapter <- as.character(.vr_roster_value(job, "adapter", NA_character_))
    metrics[[nm]]$metric_family <- nm
  }
  list(keys = keys, metrics = metrics, issue = issue)
}

validation_expected_method_roster <- function(jobs, root = ".",
    bien_traits_path = Sys.getenv("VALIDATION_BIEN_TRAITS", unset = ""),
    bien_tree_path = Sys.getenv("VALIDATION_BIEN_TREE", unset = ""), n_bins = 10L,
    expand_jobs = TRUE, bien_trait_names = NULL) {
  required <- c("job_id", "scenario_id", "adapter", "scenario", "rep")
  if (!is.data.frame(jobs) || !all(required %in% names(jobs)) || anyDuplicated(jobs$job_id))
    stop("jobs must have unique job_id and planned scenario/adapter fields")
  if (n_bins < 1L || n_bins != as.integer(n_bins)) stop("n_bins must be positive integer")
  adapter <- tolower(as.character(jobs$adapter))
  regime_contract <- NULL
  if (any(adapter %in% c("posterior", "mi_gls", "posterior_mi"))) {
    regime_file <- file.path(root, "script", "mi_gls", "regimes.R")
    if (file.exists(regime_file)) {
      regime_env <- new.env(parent = baseenv())
      sys.source(regime_file, envir = regime_env)
      regime_contract <- regime_env$regimes
    }
  }
  bien_contract <- NULL
  if (any(adapter == "real" & tolower(as.character(jobs$scenario)) == "bien")) {
    if (!is.null(bien_trait_names)) {
      if (!is.character(bien_trait_names) || !length(bien_trait_names) ||
          anyNA(bien_trait_names) || any(!nzchar(bien_trait_names)) ||
          anyDuplicated(bien_trait_names)) stop("bien_trait_names must be unique non-empty source-derived names")
      bien_contract <- list(available = TRUE, traits = bien_trait_names)
    } else bien_contract <- .vr_roster_bien_traits(bien_traits_path, bien_tree_path)
  }
  group_fields <- intersect(c("scenario_id", "slice", "adapter", "scenario", "model",
                              "regime_id", "rphylopars"), names(jobs))
  group_key <- do.call(paste, c(lapply(jobs[group_fields], function(x)
    ifelse(is.na(x), "<NA>", as.character(x))), sep = "\034"))
  job_groups <- split(seq_len(nrow(jobs)), group_key, drop = TRUE)
  out <- lapply(job_groups, function(ix) {
    prototype <- .vr_roster_for_job(jobs[ix[[1L]], , drop = FALSE], root,
      bien_traits_path, bien_tree_path, n_bins, regime_contract, bien_contract)
    planned <- jobs[ix, , drop = FALSE]
    if (nrow(prototype$keys)) {
      k <- prototype$keys
      if (expand_jobs) {
        repeat_idx <- rep(seq_len(nrow(k)), times = nrow(planned))
        keys <- k[repeat_idx, , drop = FALSE]
        keys$job_id <- rep(as.character(planned$job_id), each = nrow(k))
        keys$scenario_id <- rep(as.character(planned$scenario_id), each = nrow(k))
        keys$slice <- rep(if ("slice" %in% names(planned)) as.character(planned$slice) else NA_character_, each = nrow(k))
        keys$adapter <- rep(as.character(planned$adapter), each = nrow(k))
        keys$rep <- rep(as.integer(planned$rep), each = nrow(k))
        prototype$keys <- keys[, names(k), drop = FALSE]
      } else {
        compact_fields <- intersect(c("scenario_id", "slice", "adapter", "method", "estimand"), names(k))
        prototype$keys <- unique(k[, compact_fields, drop = FALSE])
        prototype$keys$n_planned <- nrow(planned)
      }
    }
    if (length(prototype$metrics)) for (nm in names(prototype$metrics)) {
      x <- prototype$metrics[[nm]]
      if (expand_jobs) {
        repeat_idx <- rep(seq_len(nrow(x)), times = nrow(planned))
        x <- x[repeat_idx, , drop = FALSE]
        x$job_id <- rep(as.character(planned$job_id), each = nrow(prototype$metrics[[nm]]))
        x$scenario_id <- rep(as.character(planned$scenario_id), each = nrow(prototype$metrics[[nm]]))
        x$rep <- rep(as.integer(planned$rep), each = nrow(prototype$metrics[[nm]]))
        x$adapter <- rep(as.character(planned$adapter), each = nrow(prototype$metrics[[nm]]))
      } else {
        x <- unique(x[, intersect(c("method", "trait"), names(x)), drop = FALSE])
        x$scenario_id <- as.character(planned$scenario_id[[1L]])
        x$slice <- if ("slice" %in% names(planned)) as.character(planned$slice[[1L]]) else NA_character_
        x$adapter <- as.character(planned$adapter[[1L]])
        x$metric_family <- nm
        if (identical(nm, "calibration_bins")) x$n_bins <- as.integer(n_bins)
        x$n_planned <- nrow(planned)
      }
      prototype$metrics[[nm]] <- x
    }
    if (!is.na(prototype$issue)) {
      prototype$issue_row <- data.frame(job_id = planned$job_id[[1L]],
        scenario_id = planned$scenario_id[[1L]], adapter = planned$adapter[[1L]],
        scenario = planned$scenario[[1L]], reason = prototype$issue,
        stringsAsFactors = FALSE)
    } else prototype$issue_row <- NULL
    prototype
  })
  keys <- lapply(out, `[[`, "keys"); keys <- keys[vapply(keys, nrow, integer(1)) > 0L]
  metrics <- unlist(lapply(out, `[[`, "metrics"), recursive = FALSE)
  issues <- lapply(seq_along(out), function(i) {
    out[[i]]$issue_row
  })
  keys <- if (length(keys)) do.call(rbind, keys) else data.frame(
    job_id = character(), scenario_id = character(), slice = character(), adapter = character(),
    rep = integer(), method = character(), estimand = character(), stringsAsFactors = FALSE)
  metrics <- .vr_roster_bind_fill(metrics)
  issues <- Filter(Negate(is.null), issues)
  issues <- if (length(issues)) do.call(rbind, issues) else data.frame(
    job_id = character(), scenario_id = character(), adapter = character(), scenario = character(),
    reason = character(), stringsAsFactors = FALSE)
  list(keys = keys, metrics = metrics, issues = issues)
}

validation_roster_audit <- function(jobs, records, root = ".",
    bien_traits_path = Sys.getenv("VALIDATION_BIEN_TRAITS", unset = ""),
    bien_tree_path = Sys.getenv("VALIDATION_BIEN_TREE", unset = ""), n_bins = 10L,
    expected = NULL) {
  if (is.null(expected)) expected <- validation_expected_method_roster(jobs, root, bien_traits_path,
                                                                         bien_tree_path, n_bins)
  if (!is.list(records)) stop("records must be a list of envelopes")
  result_rows <- lapply(records, function(e) {
    if (!is.list(e) || is.null(e$job_id)) return(data.frame())
    r <- e$results
    if (!is.data.frame(r) || !all(c("method", "estimand", "status") %in% names(r)))
      return(data.frame())
    if (!nrow(r)) return(data.frame())
    r$job_id <- as.character(e$job_id)
    r
  })
  result_rows <- result_rows[vapply(result_rows, nrow, integer(1)) > 0L]
  observed <- if (length(result_rows)) do.call(rbind, result_rows) else data.frame(
    job_id = character(), method = character(), estimand = character(), status = character())
  if (nrow(observed) && (anyNA(observed$status) ||
      any(!observed$status %in% c("ok", "failed", "unavailable"))))
    stop("result rows contain an unrecognized status")
  expected_keys <- expected$keys
  if (nrow(expected_keys)) {
    expected_keys$.key <- paste(expected_keys$job_id, expected_keys$method,
                                expected_keys$estimand, sep = "\034")
    observed$.key <- paste(observed$job_id, observed$method, observed$estimand, sep = "\034")
    observed_expected <- observed[observed$.key %in% expected_keys$.key, , drop = FALSE]
    if (anyDuplicated(observed_expected$.key)) stop("duplicate result keys in roster audit input")
    expected_keys$observed_status <- observed$status[match(expected_keys$.key, observed$.key)]
    expected_keys$present <- !is.na(expected_keys$observed_status)
    group_id <- interaction(expected_keys$scenario_id, expected_keys$slice,
                            expected_keys$adapter, expected_keys$method,
                            expected_keys$estimand, drop = TRUE, lex.order = TRUE)
    groups <- split(seq_len(nrow(expected_keys)), group_id)
    audit <- do.call(rbind, lapply(groups, function(ix) {
      x <- expected_keys[ix, , drop = FALSE]
      present <- x$present
      statuses <- x$observed_status[present]
      data.frame(scenario_id = x$scenario_id[[1L]], slice = x$slice[[1L]],
        adapter = x$adapter[[1L]], method = x$method[[1L]], estimand = x$estimand[[1L]],
        n_planned = nrow(x), n_present = sum(present), n_absent = sum(!present),
        n_ok = sum(statuses == "ok"), n_failed = sum(statuses == "failed"),
        n_unavailable = sum(statuses == "unavailable"),
        absent_all_replicates = !any(present),
        audit_status = if (!any(present)) "absent_all_replicates" else
          if (any(!present)) "partially_absent" else "present_all_replicates",
        stringsAsFactors = FALSE)
    }))
    rownames(audit) <- NULL
  } else audit <- data.frame(scenario_id = character(), slice = character(), adapter = character(),
    method = character(), estimand = character(), n_planned = integer(), n_present = integer(),
    n_absent = integer(), n_ok = integer(), n_failed = integer(), n_unavailable = integer(),
    absent_all_replicates = logical(), audit_status = character())
  expected_ids <- if (nrow(expected_keys)) paste(expected_keys$job_id, expected_keys$method,
    expected_keys$estimand, sep = "\034") else character()
  observed_ids <- if (nrow(observed)) paste(observed$job_id, observed$method,
    observed$estimand, sep = "\034") else character()
  unexpected <- observed[!(observed_ids %in% expected_ids), , drop = FALSE]
  rostered_jobs <- unique(expected$keys$job_id)
  unrostered_jobs <- setdiff(as.character(jobs$job_id), rostered_jobs)
  unrostered_rows <- observed[observed$job_id %in% unrostered_jobs, , drop = FALSE]
  unexpected <- unexpected[!(unexpected$job_id %in% unrostered_jobs), , drop = FALSE]
  attr(audit, "expected_roster") <- expected$keys
  attr(audit, "unexpected_rows") <- unexpected
  attr(audit, "unrostered_rows") <- unrostered_rows
  attr(audit, "roster_issues") <- expected$issues
  audit
}

validation_metric_roster_audit <- function(jobs, records, root = ".",
    bien_traits_path = Sys.getenv("VALIDATION_BIEN_TRAITS", unset = ""),
    bien_tree_path = Sys.getenv("VALIDATION_BIEN_TREE", unset = ""), n_bins = 10L,
    expected = NULL) {
  if (is.null(expected)) expected <- validation_expected_method_roster(jobs, root, bien_traits_path,
                                                                         bien_tree_path, n_bins)
  record_ids <- vapply(records, function(e) as.character(e$job_id), character(1))
  rows <- list()
  metric_keys <- expected$metrics
  if (nrow(metric_keys)) {
    metric_keys <- unique(metric_keys[, c("job_id", "scenario_id", "rep", "adapter",
      "metric_family", intersect(c("method", "trait"), names(metric_keys))), drop = FALSE])
  }
  record_match <- match(metric_keys$job_id, record_ids)
  for (i in seq_len(nrow(metric_keys))) {
    key <- metric_keys[i, , drop = FALSE]
    ri <- record_match[[i]]
    env <- if (!is.na(ri)) records[[ri]] else NULL
    payload <- if (is.null(env) || !is.list(env$metrics)) NULL else env$metrics[[key$metric_family]]
    if (identical(key$metric_family, "calibration_bins")) {
      hit <- if (is.data.frame(payload)) payload[
        payload$method == key$method & payload$trait == key$trait, , drop = FALSE] else data.frame()
      ece <- if (is.data.frame(env$results)) env$results[
        env$results$method == key$method & env$results$estimand == paste(key$trait, "ece", sep = "_") &
          env$results$status %in% c("ok", "failed", "unavailable"), , drop = FALSE] else data.frame()
      payload_status <- .vr_calibration_payload_status(hit, n_bins, ece)
      status <- if (identical(payload_status, "present")) "present" else payload_status
      rows[[length(rows) + 1L]] <- data.frame(job_id = key$job_id,
        scenario_id = key$scenario_id, rep = key$rep, adapter = key$adapter,
        metric_family = key$metric_family, method = key$method, trait = key$trait,
        bin = NA_integer_, n_expected = n_bins, n_present = nrow(hit),
        n_missing_or_nonfinite = NA_integer_, status = status,
        stringsAsFactors = FALSE)
    } else if (identical(key$metric_family, "posterior_cell_interval")) {
      hit <- if (is.data.frame(payload)) payload[payload$trait == key$trait, , drop = FALSE] else data.frame()
      status <- .vr_posterior_interval_status(hit)
      rows[[length(rows) + 1L]] <- data.frame(job_id = key$job_id,
        scenario_id = key$scenario_id, rep = key$rep, adapter = key$adapter,
        metric_family = key$metric_family, method = NA_character_, trait = key$trait,
        bin = NA_integer_, n_expected = if (nrow(hit)) hit$n_expected[[1L]] else NA_integer_,
        n_present = if (nrow(hit)) hit$n_finite_scored[[1L]] else 0L,
        n_missing_or_nonfinite = if (nrow(hit)) hit$n_missing_or_nonfinite[[1L]] else NA_integer_,
        status = status, stringsAsFactors = FALSE)
    }
  }
  if (length(rows)) do.call(rbind, rows) else data.frame(
    job_id = character(), scenario_id = character(), rep = integer(), adapter = character(),
    metric_family = character(), method = character(), trait = character(), bin = integer(),
    n_expected = integer(), n_present = integer(), n_missing_or_nonfinite = integer(), status = character())
}

# Adapter for evaluator-retained CSV receipts. The receipts contain enough
# provenance to reconstruct the per-job envelope shape, but not source identity;
# callers should first validate those tables against their saved identity.
.vr_roster_receipt_table <- function(x, label, optional = FALSE) {
  if (is.null(x) && optional) return(data.frame())
  if (is.data.frame(x)) return(x)
  if (is.character(x) && length(x) == 1L && file.exists(x))
    return(utils::read.csv(x, stringsAsFactors = FALSE, check.names = FALSE))
  stop(label, " must be a data frame or an existing CSV/CSV.GZ path")
}

validation_roster_audit_tables <- function(jobs, result_rows,
    calibration_bins = NULL, posterior_intervals = NULL, attempts = NULL,
    root = ".", bien_traits_path = Sys.getenv("VALIDATION_BIEN_TRAITS", unset = ""),
    bien_tree_path = Sys.getenv("VALIDATION_BIEN_TREE", unset = ""), n_bins = 10L,
    bien_trait_names = NULL) {
  jobs <- .vr_roster_receipt_table(jobs, "jobs")
  result_rows <- .vr_roster_receipt_table(result_rows, "result_rows")
  calibration_bins <- .vr_roster_receipt_table(calibration_bins, "calibration_bins", optional = TRUE)
  posterior_intervals <- .vr_roster_receipt_table(posterior_intervals,
                                                   "posterior_intervals", optional = TRUE)
  attempts <- .vr_roster_receipt_table(attempts, "attempts", optional = TRUE)
  if (!is.data.frame(jobs) || !all(c("job_id", "scenario_id", "rep", "seed") %in% names(jobs)) ||
      anyDuplicated(jobs$job_id)) stop("jobs receipt has invalid job identity")
  if (!is.data.frame(result_rows) || !all(c("job_id", "method", "estimand", "status") %in% names(result_rows)))
    stop("result_rows receipt has an invalid schema")
  if (nrow(result_rows) && (anyNA(result_rows$status) ||
      any(!result_rows$status %in% c("ok", "failed", "unavailable"))))
    stop("result_rows contain an unrecognized status")
  if (any(!result_rows$job_id %in% jobs$job_id)) stop("result_rows contain unplanned job IDs")
  if (anyDuplicated(paste(result_rows$job_id, result_rows$method,
                          result_rows$estimand, sep = "\034")))
    stop("result_rows contain duplicate job/method/estimand keys")
  if (is.null(calibration_bins)) calibration_bins <- data.frame()
  if (!is.data.frame(calibration_bins)) stop("calibration_bins must be a data frame")
  if (nrow(calibration_bins) &&
      !all(c("job_id", "method", "trait", "bin", "n", "status") %in% names(calibration_bins)))
    stop("calibration_bins receipt has an invalid schema")
  if (nrow(calibration_bins) && any(!calibration_bins$job_id %in% jobs$job_id))
    stop("calibration_bins contain unplanned job IDs")
  if (is.null(posterior_intervals)) posterior_intervals <- data.frame()
  if (!is.data.frame(posterior_intervals)) stop("posterior_intervals must be a data frame")
  if (nrow(posterior_intervals) &&
      !all(c("job_id", "trait", "n_expected", "n_finite_scored",
             "n_missing_or_nonfinite", "coverage", "mean_width", "interval_status") %in%
           names(posterior_intervals)))
    stop("posterior_intervals receipt has an invalid schema")
  if (nrow(posterior_intervals) && any(!posterior_intervals$job_id %in% jobs$job_id))
    stop("posterior_intervals contain unplanned job IDs")
  if (is.null(attempts)) attempts <- data.frame()
  if (nrow(attempts) && !all(c("job_id", "status") %in% names(attempts)))
    stop("attempts receipt has an invalid schema")
  if (nrow(attempts) && (anyDuplicated(attempts$job_id) || any(!attempts$job_id %in% jobs$job_id)))
    stop("attempts receipt has duplicate or unplanned job IDs")

  # Partition each large receipt once. Re-filtering million-row result tables
  # once per planned job turns this audit into O(jobs * receipt_rows).
  split_rows <- function(x) {
    if (!nrow(x)) return(list())
    split(seq_len(nrow(x)), as.character(x$job_id), drop = TRUE)
  }
  result_index <- split_rows(result_rows)
  calibration_index <- split_rows(calibration_bins)
  posterior_index <- split_rows(posterior_intervals)
  attempt_match <- if (nrow(attempts)) match(jobs$job_id, attempts$job_id) else rep(NA_integer_, nrow(jobs))
  take_rows <- function(x, index, id) {
    ix <- index[[id]]
    if (is.null(ix)) x[FALSE, , drop = FALSE] else x[ix, , drop = FALSE]
  }

  expected <- validation_expected_method_roster(jobs, root, bien_traits_path,
    bien_tree_path, n_bins, expand_jobs = FALSE, bien_trait_names = bien_trait_names)

  # Summarise expected result keys by scenario contract. This keeps memory and
  # runtime proportional to the number of distinct roster keys, not jobs × keys.
  job_ix <- match(result_rows$job_id, jobs$job_id)
  if ("scenario_id" %in% names(result_rows) &&
      any(as.character(result_rows$scenario_id) != as.character(jobs$scenario_id[job_ix])))
    stop("result_rows scenario_id does not match the planned job table")
  result_rows$scenario_id <- as.character(jobs$scenario_id[job_ix])
  result_rows$slice <- if ("slice" %in% names(jobs)) as.character(jobs$slice[job_ix]) else NA_character_
  result_rows$adapter <- as.character(jobs$adapter[job_ix])
  result_key <- function(x, fields) do.call(paste, c(lapply(x[fields], function(v)
    ifelse(is.na(v), "<NA>", as.character(v))), sep = "\034"))
  key_fields <- c("scenario_id", "slice", "adapter", "method", "estimand")
  expected_keys <- expected$keys
  if (nrow(expected_keys)) {
    expected_keys$.key <- result_key(expected_keys, key_fields)
    observed_group <- split(seq_len(nrow(result_rows)), result_key(result_rows, key_fields), drop = TRUE)
    audit_rows <- lapply(seq_len(nrow(expected_keys)), function(i) {
      x <- expected_keys[i, , drop = FALSE]
      ix <- observed_group[[x$.key]]
      hit <- if (is.null(ix)) result_rows[FALSE, , drop = FALSE] else result_rows[ix, , drop = FALSE]
      n_present <- length(unique(hit$job_id))
      n_planned <- as.integer(x$n_planned)
      data.frame(scenario_id = x$scenario_id, slice = x$slice, adapter = x$adapter,
        method = x$method, estimand = x$estimand, n_planned = n_planned,
        n_present = n_present, n_absent = max(0L, n_planned - n_present),
        n_ok = sum(hit$status == "ok"), n_failed = sum(hit$status == "failed"),
        n_unavailable = sum(hit$status == "unavailable"),
        absent_all_replicates = n_present == 0L,
        audit_status = if (n_present == 0L) "absent_all_replicates" else
          if (n_present < n_planned) "partially_absent" else "present_all_replicates",
        stringsAsFactors = FALSE)
    })
    audit <- do.call(rbind, audit_rows)
  } else audit <- data.frame()
  rownames(audit) <- NULL
  known_keys <- if (nrow(expected_keys)) expected_keys$.key else character()
  observed_keys <- result_key(result_rows, key_fields)
  issue_scenarios <- unique(as.character(expected$issues$scenario_id))
  is_unrostered <- result_rows$scenario_id %in% issue_scenarios
  unrostered_rows <- result_rows[is_unrostered, , drop = FALSE]
  unexpected_rows <- result_rows[!is_unrostered & !observed_keys %in% known_keys, , drop = FALSE]

  # Metric payloads are checked once per job × method × trait (or job × trait
  # for posterior intervals); calibration bins are validated as a set.
  metric_expected <- expected$metrics
  metric_rows <- list()
  if (nrow(metric_expected)) for (m in seq_len(nrow(metric_expected))) {
    spec <- metric_expected[m, , drop = FALSE]
    matching_jobs <- jobs[jobs$scenario_id == spec$scenario_id & jobs$adapter == spec$adapter, , drop = FALSE]
    if ("slice" %in% names(matching_jobs) && !is.na(spec$slice))
      matching_jobs <- matching_jobs[matching_jobs$slice == spec$slice, , drop = FALSE]
    metric_payload <- if (spec$metric_family == "calibration_bins") calibration_bins else posterior_intervals
    for (ji in seq_len(nrow(matching_jobs))) {
      id <- as.character(matching_jobs$job_id[[ji]])
      ix <- if (spec$metric_family == "calibration_bins") calibration_index[[id]] else posterior_index[[id]]
      payload <- if (is.null(ix)) metric_payload[FALSE, , drop = FALSE] else metric_payload[ix, , drop = FALSE]
      if (spec$metric_family == "calibration_bins") {
        hit <- payload[payload$method == spec$method & payload$trait == spec$trait, , drop = FALSE]
        result_ix <- result_index[[id]]
        job_rows <- if (is.null(result_ix)) result_rows[FALSE, , drop = FALSE] else
          result_rows[result_ix, , drop = FALSE]
        ece <- job_rows[job_rows$method == spec$method &
          job_rows$estimand == paste(spec$trait, "ece", sep = "_") &
          job_rows$status %in% c("ok", "failed", "unavailable"), , drop = FALSE]
        payload_status <- .vr_calibration_payload_status(hit, n_bins, ece)
        status <- if (identical(payload_status, "present")) "present" else payload_status
        n_present <- nrow(hit); n_expected <- n_bins; n_missing <- NA_integer_
      } else {
        hit <- payload[payload$trait == spec$trait, , drop = FALSE]
        status <- .vr_posterior_interval_status(hit)
        n_present <- if (nrow(hit)) hit$n_finite_scored[[1L]] else 0L
        n_expected <- if (nrow(hit)) hit$n_expected[[1L]] else NA_integer_
        n_missing <- if (nrow(hit)) hit$n_missing_or_nonfinite[[1L]] else NA_integer_
      }
      metric_rows[[length(metric_rows) + 1L]] <- data.frame(job_id = id,
        scenario_id = spec$scenario_id, rep = as.integer(matching_jobs$rep[[ji]]),
        adapter = spec$adapter, metric_family = spec$metric_family,
        method = if ("method" %in% names(spec)) as.character(spec$method) else NA_character_,
        trait = as.character(spec$trait), bin = NA_integer_, n_expected = n_expected,
        n_present = n_present, n_missing_or_nonfinite = n_missing, status = status,
        stringsAsFactors = FALSE)
    }
  }
  metric_audit <- if (length(metric_rows)) do.call(rbind, metric_rows) else data.frame()
  if (nrow(metric_audit)) {
    group_fields <- c("scenario_id", "adapter", "metric_family", "method", "trait")
    group_id <- result_key(metric_audit, group_fields)
    metric_audit <- do.call(rbind, lapply(split(seq_len(nrow(metric_audit)), group_id), function(ix) {
      x <- metric_audit[ix, , drop = FALSE]
      data.frame(scenario_id = x$scenario_id[[1L]], adapter = x$adapter[[1L]],
        metric_family = x$metric_family[[1L]], method = x$method[[1L]], trait = x$trait[[1L]],
        n_planned = nrow(x), n_present = sum(x$status == "present"),
        n_missing = sum(x$status == "missing"), n_failed = sum(x$status == "failed"),
        n_unavailable = sum(x$status == "unavailable"),
        absent_all_replicates = all(x$status == "missing"), stringsAsFactors = FALSE)
    }))
    rownames(metric_audit) <- NULL
  }
  list(audit = audit, metric_audit = metric_audit,
       expected_roster = expected_keys,
       unexpected_rows = unexpected_rows,
       unrostered_rows = unrostered_rows,
       roster_issues = expected$issues)
}
