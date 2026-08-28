#' Check whether traits and a tree are ready for pigauto
#'
#' Inspect declared inputs without constructing a phylogenetic graph or fitting
#' a model. The check is intended to make input problems and likely compute
#' scale visible before calling [impute()]. Its fingerprint identifies the
#' declared input object for comparison within a workflow; it is not a security
#' hash.
#'
#' @param traits A data.frame of traits, with row names or a species column.
#' @param tree A `phylo` tree.
#' @param species_col Optional column naming species identities.
#' @param trait_types Optional named type overrides passed to
#'   [preprocess_traits()].
#' @param multi_proportion_groups Optional named composition groups passed to
#'   [preprocess_traits()].
#' @return A list of class `pigauto_check` containing a schema version, status,
#'   stable input fingerprint, input summaries, and structured messages.
#' @examples
#' \donttest{
#' tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
#' traits <- data.frame(x = c(1, NA, 3), row.names = tree$tip.label)
#' check_pigauto(traits, tree)
#' }
#' @export
check_pigauto <- function(traits, tree, species_col = NULL,
                          trait_types = NULL,
                          multi_proportion_groups = NULL) {
  .check_pigauto_internal(
    traits, tree, species_col, trait_types, multi_proportion_groups,
    log_transform = TRUE, covariates = NULL
  )$check
}

.check_pigauto_messages <- function() {
  data.frame(code = character(), field = character(), text = character(),
             action = character(), severity = character(),
             stringsAsFactors = FALSE)
}

.check_pigauto_add <- function(messages, severity, code, field, text, action) {
  rbind(messages, data.frame(code = code, field = field, text = text,
                             action = action, severity = severity,
                             stringsAsFactors = FALSE))
}

.check_pigauto_fingerprint <- function(traits, tree, species_col, trait_types,
                                        multi_proportion_groups) {
  canonical <- list(
    traits = traits,
    tree = list(tip.label = tree$tip.label, edge = tree$edge,
                edge.length = tree$edge.length, Nnode = tree$Nnode),
    species_col = species_col, trait_types = trait_types,
    multi_proportion_groups = multi_proportion_groups
  )
  path <- tempfile("pigauto-check-", fileext = ".xdr")
  on.exit(unlink(path), add = TRUE)
  writeBin(serialize(canonical, NULL, version = 2L, xdr = TRUE), path)
  list(algorithm = "md5-xdr-v1", value = unname(tools::md5sum(path)),
       scope = "declared traits, tree topology/lengths, species and type declarations")
}

.check_pigauto_cuda_available <- function() {
  tryCatch(torch::cuda_is_available(), error = function(e) FALSE)
}

.check_pigauto_mps_available <- function() {
  if (!"backends_mps_is_available" %in% getNamespaceExports("torch")) return(FALSE)
  tryCatch(torch::backends_mps_is_available(), error = function(e) FALSE)
}

.check_pigauto_select_device <- function(cuda_available, mps_available) {
  if (isTRUE(cuda_available)) "cuda" else if (isTRUE(mps_available)) "mps" else "cpu"
}

.check_pigauto_probe_device <- function(device, tensor_fun = torch::torch_tensor) {
  probe <- tryCatch(tensor_fun(1, device = device), error = function(e) e)
  if (inherits(probe, "error")) {
    return(list(ok = FALSE, detail = paste0(device, " scalar tensor probe failed: ", conditionMessage(probe))))
  }
  list(ok = TRUE, detail = paste0(device, " scalar tensor probe succeeded"))
}

.check_pigauto_runtime <- function() {
  out <- list(torch_package = requireNamespace("torch", quietly = TRUE),
              torch_is_installed = FALSE, device = "unavailable",
              unavailable = TRUE, detail = "torch package is not installed")
  if (!out$torch_package) return(out)
  installed <- tryCatch(torch::torch_is_installed(), error = function(e) e)
  if (inherits(installed, "error") || !isTRUE(installed)) {
    out$detail <- if (inherits(installed, "error")) conditionMessage(installed) else "torch runtime is not installed"
    return(out)
  }
  cuda <- .check_pigauto_cuda_available()
  mps <- .check_pigauto_mps_available()
  out$device <- .check_pigauto_select_device(cuda, mps)
  probe <- .check_pigauto_probe_device(out$device)
  out$torch_is_installed <- TRUE
  out$unavailable <- !isTRUE(probe$ok)
  out$detail <- probe$detail
  out
}

.check_pigauto_internal <- function(traits, tree, species_col = NULL,
                                    trait_types = NULL,
                                    multi_proportion_groups = NULL,
                                    log_transform = TRUE, covariates = NULL) {
  messages <- .check_pigauto_messages()
  empty <- list(names = character(), examples = character(), count = 0L,
                consequence = "")
  if (!is.data.frame(traits)) {
    messages <- .check_pigauto_add(messages, "error", "traits_not_data_frame", "traits",
                                    "`traits` must be a data.frame.", "Supply a data.frame.")
  }
  if (!inherits(tree, "phylo")) {
    messages <- .check_pigauto_add(messages, "error", "tree_not_phylo", "tree",
                                    "`tree` must inherit from `phylo`.", "Supply an ape phylo tree.")
  }
  if (nrow(messages) > 0L) return(.check_pigauto_finish(messages, empty, empty, empty, empty, empty, NULL))

  fingerprint <- .check_pigauto_fingerprint(traits, tree, species_col, trait_types,
                                             multi_proportion_groups)
  input_n <- nrow(traits)
  if (!is.null(species_col)) {
    if (!species_col %in% names(traits)) {
      messages <- .check_pigauto_add(messages, "error", "species_column_missing", "species_col",
                                      "The requested species column is absent.", "Choose an existing column.")
      ids <- character()
    } else ids <- as.character(traits[[species_col]])
    mode <- "species_col"
  } else {
    ids <- rownames(traits)
    mode <- "rownames"
    if (is.null(ids) || anyDuplicated(ids)) {
      messages <- .check_pigauto_add(messages, "error", "duplicate_or_missing_rownames", "rownames",
                                      "Row names must be present and unique when `species_col` is NULL.",
                                      "Supply unique row names or use `species_col`.")
    }
  }
  if (length(ids) && (anyNA(ids) || any(!nzchar(ids)))) {
    messages <- .check_pigauto_add(messages, "error", "missing_species_id", "species",
                                    "Species identifiers cannot be missing or empty.", "Repair species identifiers.")
  }
  tips <- tree$tip.label %||% character()
  edge <- tree$edge
  edge_length <- tree$edge.length
  if (is.null(edge_length) || length(edge_length) != nrow(edge %||% matrix(, 0, 2))) {
    messages <- .check_pigauto_add(messages, "error", "branch_length_missing_or_wrong", "tree$edge.length",
                                    "Tree branch lengths are missing or have the wrong length.", "Supply one finite length per edge.")
  } else if (any(!is.finite(edge_length)) || any(edge_length < 0)) {
    messages <- .check_pigauto_add(messages, "error", "branch_length_invalid", "tree$edge.length",
                                    "Tree branch lengths must be finite and non-negative.", "Repair branch lengths.")
  } else if (all(edge_length == 0)) {
    messages <- .check_pigauto_add(messages, "error", "branch_length_all_zero", "tree$edge.length",
                                    "All branch lengths are zero.", "Supply positive branch lengths.")
  } else if (any(edge_length == 0)) {
    messages <- .check_pigauto_add(messages, "warning", "branch_length_zero", "tree$edge.length",
                                    "Some branch lengths are zero.", "Inspect whether zero-length branches are intended.")
  }
  if (anyDuplicated(tips)) messages <- .check_pigauto_add(messages, "error", "duplicate_tree_tip", "tree$tip.label",
                                                             "Tree tip labels must be unique.", "Repair duplicate tree labels.")
  repeated <- if (length(ids)) unique(ids[duplicated(ids)]) else character()
  matched <- intersect(unique(ids), tips)
  data_only <- setdiff(unique(ids), tips)
  tree_only <- setdiff(tips, unique(ids))
  species <- list(mode = mode, species_col = species_col, input_rows = input_n,
                  input_species = length(unique(ids)), tree_tips = length(tips),
                  matched = list(names = matched, examples = utils::head(matched, 8L), count = length(matched), consequence = "modeled rows"),
                  data_only = list(names = data_only, examples = utils::head(data_only, 8L), count = length(data_only), consequence = "kept in completed data but not modeled or filled"),
                  tree_only = list(names = tree_only, examples = utils::head(tree_only, 8L), count = length(tree_only), consequence = "internal all-missing rows only"),
                  repeated = list(state = length(repeated) > 0L, count = length(repeated), names = repeated,
                                  examples = utils::head(repeated, 8L), max = if (length(ids)) max(tabulate(match(ids, unique(ids)))) else 0L))
  if (!is.null(species_col) && length(repeated)) messages <- .check_pigauto_add(messages, "info", "repeated_species_rows", species_col,
                                                                                   "Repeated species rows are supported with `species_col`.", "Check that repeated observations are intended.")
  branch_status <- if (is.null(edge_length) || length(edge_length) != nrow(edge %||% matrix(, 0, 2)) ||
                       any(!is.finite(edge_length)) || any(edge_length < 0) || all(edge_length == 0)) "error" else if (any(edge_length == 0)) "warning" else "ok"
  tree_summary <- list(tips = length(tips), edges = if (is.null(edge)) 0L else nrow(edge),
                       rooted = isTRUE(ape::is.rooted(tree)), binary = isTRUE(ape::is.binary(tree)),
                       duplicate_labels = anyDuplicated(tips) > 0L,
                       branch = list(status = branch_status,
                                     count = length(edge_length), min = if (length(edge_length)) min(edge_length) else NA_real_,
                                     max = if (length(edge_length)) max(edge_length) else NA_real_,
                                     total = if (length(edge_length)) sum(edge_length) else NA_real_))
  if (!isTRUE(ape::is.rooted(tree))) messages <- .check_pigauto_add(messages, "info", "tree_unrooted", "tree", "Tree is unrooted.", "Rooting is not required by the current preprocessing path.")
  if (!isTRUE(ape::is.binary(tree))) messages <- .check_pigauto_add(messages, "info", "tree_polytomy", "tree", "Tree has a polytomy.", "Resolve only if your analysis requires a binary tree.")

  prepared <- NULL
  trait_table <- data.frame()
  if (!any(messages$severity == "error")) {
    prep_warnings <- character(); prep_messages <- character()
    prepared <- tryCatch(withCallingHandlers(
      preprocess_traits(traits, tree, species_col = species_col, trait_types = trait_types,
                        multi_proportion_groups = multi_proportion_groups,
                        log_transform = log_transform, covariates = covariates),
      warning = function(w) { prep_warnings <<- c(prep_warnings, conditionMessage(w)); invokeRestart("muffleWarning") },
      message = function(m) { prep_messages <<- c(prep_messages, conditionMessage(m)); invokeRestart("muffleMessage") }),
      error = function(e) e)
    if (inherits(prepared, "error")) {
      messages <- .check_pigauto_add(messages, "error", "preprocess_error", "traits", conditionMessage(prepared), "Repair the declared input and run check_pigauto() again.")
      prepared <- NULL
    } else {
      for (w in prep_warnings) messages <- .check_pigauto_add(messages, "warning", "preprocess_warning", "traits", w, "Review the affected input rows.")
      for (m in prep_messages) messages <- .check_pigauto_add(messages, "info", "preprocess_message", "traits", m, "Review the declared input.")
      trait_table <- .check_pigauto_traits_table(traits, prepared, species_col, trait_types, multi_proportion_groups, matched)
      actionable <- which(nzchar(trait_table$correction))
      for (i in actionable) messages <- .check_pigauto_add(
        messages,
        if (identical(trait_table$type[[i]], "count") && identical(trait_table$type_source[[i]], "auto")) "info" else "info",
        if (identical(trait_table$type[[i]], "count") && identical(trait_table$type_source[[i]], "auto")) "integer_auto_ambiguous" else "trait_declaration_available",
        trait_table$trait[[i]],
        if (identical(trait_table$type[[i]], "count") && identical(trait_table$type_source[[i]], "auto")) paste0("Integer trait `", trait_table$trait[[i]], "` was auto-detected as count and may be a continuous measurement stored as whole numbers.") else paste0("Trait `", trait_table$trait[[i]], "` has an explicit declaration option."),
        trait_table$correction[[i]]
      )
      for (i in seq_len(nrow(trait_table))) {
        nm <- trait_table$trait[[i]]
        if (nm %in% names(traits) && identical(trait_table$type_source[[i]], "auto")) {
          target <- if (is.null(species_col)) rownames(traits) %in% matched else as.character(traits[[species_col]]) %in% matched
          values <- traits[[nm]][target & !is.na(traits[[nm]])]
          if (identical(trait_table$type[[i]], "continuous") && is.numeric(values) && length(values) && all(values >= 0 & values <= 1)) messages <- .check_pigauto_add(messages, "info", "bounded_numeric_may_be_proportion", nm, paste0("Bounded numeric trait `", nm, "` may represent a proportion if its 0-1 scale is semantic."), .check_pigauto_trait_action(nm, "proportion"))
          if (identical(trait_table$type[[i]], "count") && is.integer(values) && any(values == 0L)) messages <- .check_pigauto_add(messages, "info", "zero_integer_may_be_zi_count", nm, paste0("Integer trait `", nm, "` contains zeros and may be zero-inflated if zeros are a separate process."), .check_pigauto_trait_action(nm, "zi_count"))
        }
      }
      if (any(trait_table$all_missing)) messages <- .check_pigauto_add(messages, "error", "trait_all_missing", "traits", "At least one matched trait is all missing.", "Supply observed values for every trait or remove it.")
      if (any(trait_table$constant & !trait_table$all_missing)) messages <- .check_pigauto_add(messages, "warning", "trait_constant", "traits", "At least one matched trait is constant.", "Review whether a constant trait is informative.")
    }
  }
  runtime <- .check_pigauto_runtime()
  if (!isTRUE(runtime$torch_is_installed) || isTRUE(runtime$unavailable)) messages <- .check_pigauto_add(messages, "error", "torch_unavailable", "runtime", runtime$detail, "Install or repair torch before fitting.")
  n_species <- if (!is.null(prepared)) prepared$n_species else length(matched)
  n_obs <- if (!is.null(prepared)) prepared$n_obs else sum(ids %in% matched)
  p_latent <- if (!is.null(prepared)) prepared$p_latent else NA_integer_
  cells <- as.double(n_species) * as.double(n_species)
  bytes <- 4 * 8 * cells
  gib <- bytes / 1024^3
  size <- list(n_species = n_species, n_obs = n_obs, p_latent = p_latent,
               latent_cells = if (is.na(p_latent)) NA_real_ else as.double(n_obs) * p_latent,
               dense_lower_bound_bytes = bytes, dense_lower_bound_gib = gib,
               note = "Lower bound is four dense double n_species by n_species matrices.",
               memory_class = if (gib < .25) "low" else if (gib < 1) "moderate" else if (gib < 4) "high" else "extreme",
               fit_class = if (n_species < 500) "small" else if (n_species < 2000) "moderate" else if (n_species < 5000) "large" else if (n_species <= 10000) "very_large" else "extreme")
  .check_pigauto_finish(messages, species, tree_summary, trait_table, runtime, size,
                         fingerprint, prepared)
}

.check_pigauto_traits_table <- function(traits, pd, species_col, trait_types, groups, matched) {
  source <- function(nm, type) if (type == "multi_proportion") "group" else if (!is.null(trait_types) && nm %in% names(trait_types)) "override" else "auto"
  target <- if (is.null(species_col)) rownames(traits) %in% matched else as.character(traits[[species_col]]) %in% matched
  do.call(rbind, lapply(pd$trait_map, function(tm) {
    cols <- tm$input_cols %||% tm$name
    observed <- traits[target, cols, drop = FALSE]
    usable <- if (length(cols) > 1L) stats::complete.cases(observed) else !is.na(observed[[1L]])
    vals <- unlist(observed, use.names = FALSE)
    correction <- if (length(cols) == 1L && is.integer(traits[[cols]]) && source(tm$name, tm$type) == "auto") {
      .check_pigauto_trait_action(tm$name, "continuous")
    } else if (tm$type %in% c("proportion", "zi_count")) {
      .check_pigauto_trait_action(tm$name, tm$type)
    } else if (tm$type == "multi_proportion") {
      .check_pigauto_group_action(tm$name, cols)
    } else ""
    row_signatures <- if (length(cols) > 1L) apply(observed[usable, , drop = FALSE], 1L, paste, collapse = "\r") else vals[!is.na(vals)]
    transform <- switch(tm$type, continuous = if (isTRUE(tm$log_transform)) "log" else "identity", count = "log1p", ordinal = "z_score", proportion = "logit", zi_count = "binary gate + log1p magnitude", multi_proportion = "CLR", binary = "binary", categorical = "one_hot", tm$type)
    data.frame(trait = tm$name, input_columns = paste(cols, collapse = ","),
               input_class = paste(vapply(observed, function(x) class(x)[1L], character(1)), collapse = ","),
               type = tm$type, type_source = source(tm$name, tm$type),
               transform = transform,
               n_matched_rows = nrow(observed), n_usable = sum(usable),
               n_missing = sum(!usable), missing_fraction = mean(!usable),
               n_species_usable = length(unique((if (is.null(species_col)) rownames(traits) else traits[[species_col]])[target][usable])),
               n_unique_observed = length(unique(vals[!is.na(vals)])),
               all_missing = !any(usable), constant = sum(usable) > 0L && length(unique(row_signatures)) <= 1L,
               correction = correction, stringsAsFactors = FALSE)
  }))
}

.check_pigauto_syntactic_name <- function(x) {
  make.names(x) == x && grepl("^[A-Za-z.]", x)
}

.check_pigauto_string <- function(x) {
  encodeString(x, quote = "\"")
}

.check_pigauto_trait_action <- function(name, type) {
  if (.check_pigauto_syntactic_name(name)) {
    paste0("trait_types = c(", name, " = ", .check_pigauto_string(type), ")")
  } else {
    paste0("trait_types = stats::setNames(", .check_pigauto_string(type), ", ",
           .check_pigauto_string(name), ")")
  }
}

.check_pigauto_group_action <- function(name, cols) {
  values <- paste(vapply(cols, .check_pigauto_string, character(1)), collapse = ", ")
  if (.check_pigauto_syntactic_name(name)) {
    paste0("multi_proportion_groups = list(", name, " = c(", values, "))")
  } else {
    paste0("multi_proportion_groups = stats::setNames(list(c(", values, ")), ",
           .check_pigauto_string(name), ")")
  }
}

.check_pigauto_finish <- function(messages, species, tree, traits, runtime, size,
                                  fingerprint, prepared = NULL) {
  error <- any(messages$severity == "error")
  not_needed <- !error && is.data.frame(traits) && nrow(traits) > 0L && all(traits$n_missing == 0L)
  if (not_needed) messages <- .check_pigauto_add(messages, "info", "fully_observed", "traits",
                                                   "All matched target cells are observed; fitting is not needed.",
                                                   "Use cross_validate() to evaluate predictive performance.")
  status <- if (error) "error" else if (not_needed) "not_needed" else if (any(messages$severity == "warning")) "ready_with_warnings" else "ready"
  check <- structure(list(schema_version = 1L, status = status, ready = status %in% c("ready", "ready_with_warnings"),
                          fingerprint = fingerprint, species = species, tree = tree, traits = traits,
                          runtime = runtime, size = size, messages = messages[, c("code", "field", "text", "action", "severity"), drop = FALSE]),
                     class = c("pigauto_check", "list"))
  list(check = check, prepared = prepared)
}

#' Print a pigauto input check
#'
#' @param x A `pigauto_check` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.pigauto_check <- function(x, ...) {
  cat("pigauto input check:", x$status, "\n")
  if (nrow(x$messages)) {
    cat("  messages:", nrow(x$messages), "\n")
    for (i in seq_len(min(3L, nrow(x$messages)))) cat("  -", x$messages$text[[i]], " Action:", x$messages$action[[i]], "\n")
  }
  invisible(x)
}
