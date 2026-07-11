#!/usr/bin/env Rscript

# Consolidate disjoint, same-SHA calibration result directories into one full
# product-grid campaign without rerunning completed DGP replicates.

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else ""
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path)) else getwd()
source(file.path(script_dir, "0_prepare.R"))

inputs_arg <- .arg_value(args, "inputs", NULL)
output_arg <- .arg_value(args, "output", NULL)
if (is.null(inputs_arg) || is.null(output_arg)) {
  stop("Specify --inputs=dir1,dir2,... and --output=dir.", call. = FALSE)
}
inputs <- trimws(strsplit(inputs_arg, ",", fixed = TRUE)[[1]])
inputs <- normalizePath(inputs, mustWork = TRUE)
output <- normalizePath(output_arg, mustWork = FALSE)
if (output %in% inputs) stop("Output must not be one of the input directories.")
if (dir.exists(output) && length(list.files(output, all.files = TRUE,
                                             no.. = TRUE)) > 0L) {
  stop("Output directory is not empty: ", output, call. = FALSE)
}

manifest_paths <- file.path(inputs, "manifest.rds")
if (any(!file.exists(manifest_paths))) {
  stop("Every input must contain manifest.rds.", call. = FALSE)
}
manifests <- lapply(manifest_paths, readRDS)
manifest <- do.call(rbind, manifests)
rownames(manifest) <- NULL

key <- function(dgp, regime, replicate) {
  paste(dgp, regime, as.integer(replicate), sep = "::")
}
manifest_key <- key(manifest$dgp, manifest$regime, manifest$replicate)
if (anyDuplicated(manifest_key)) {
  stop("Input manifests contain duplicate DGP/regime/replicate cells.",
       call. = FALSE)
}

replicates <- sort(unique(manifest$replicate))
expected <- expand.grid(
  dgp = c("lm", "glm", "lmer"),
  regime = c("phylogeny", "auxiliary"),
  replicate = replicates,
  stringsAsFactors = FALSE
)
expected_key <- key(expected$dgp, expected$regime, expected$replicate)
if (!setequal(manifest_key, expected_key)) {
  missing <- setdiff(expected_key, manifest_key)
  extra <- setdiff(manifest_key, expected_key)
  stop(
    "Inputs do not form the complete 3 x 2 product grid. Missing: ",
    length(missing), "; extra: ", length(extra), ".", call. = FALSE
  )
}

order_index <- order(
  match(manifest$dgp, c("lm", "glm", "lmer")),
  match(manifest$regime, c("phylogeny", "auxiliary")),
  manifest$replicate
)
manifest <- manifest[order_index, , drop = FALSE]
rownames(manifest) <- NULL
manifest$task_id <- seq_len(nrow(manifest))
manifest_key <- key(manifest$dgp, manifest$regime, manifest$replicate)

setting_names <- c("m", "smcfcs_numit", "jomo_nburn", "jomo_nbetween")
if (!all(setting_names %in% names(manifest)) ||
    nrow(unique(manifest[, setting_names, drop = FALSE])) != 1L) {
  stop("Input manifests do not share one M/SMC setting tuple.", call. = FALSE)
}

files <- unlist(lapply(inputs, function(input) {
  list.files(input, pattern = "^task-[0-9]+\\.rds$", full.names = TRUE)
}), use.names = FALSE)
if (length(files) != nrow(manifest)) {
  stop("Found ", length(files), " task files for ", nrow(manifest),
       " manifest rows.", call. = FALSE)
}
results <- lapply(files, readRDS)
meta_value <- function(result, name) unlist(result$meta[[name]], use.names = FALSE)[1]
result_key <- vapply(results, function(result) {
  key(meta_value(result, "dgp"), meta_value(result, "regime"),
      meta_value(result, "replicate"))
}, character(1))
if (anyDuplicated(result_key) || !setequal(result_key, manifest_key)) {
  stop("Task results do not uniquely match the consolidated manifest.",
       call. = FALSE)
}

git_sha <- vapply(results, function(x) x$provenance$git_sha %||% NA_character_,
                   character(1))
git_dirty <- vapply(results, function(x) x$provenance$git_dirty %||% NA,
                    logical(1))
package_version <- vapply(
  results, function(x) x$provenance$package_version %||% NA_character_,
  character(1)
)
if (anyNA(git_sha) || length(unique(git_sha)) != 1L ||
    anyNA(git_dirty) || any(git_dirty) ||
    anyNA(package_version) || length(unique(package_version)) != 1L) {
  stop("Results do not share one clean SHA and package version.", call. = FALSE)
}

dir.create(output, recursive = TRUE, showWarnings = FALSE)
atomic_save_rds(manifest, file.path(output, "manifest.rds"))
for (i in seq_along(results)) {
  task_id <- match(result_key[[i]], manifest_key)
  result <- results[[i]]
  result$meta$task_id <- as.integer(task_id)
  atomic_save_rds(result, file.path(output, sprintf("task-%06d.rds", task_id)))
}

receipt <- c(
  "pigauto MI validation consolidation",
  paste0("created_at: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste0("source_sha: ", unique(git_sha)),
  paste0("package_version: ", unique(package_version)),
  paste0("input_dirs: ", paste(inputs, collapse = ", ")),
  paste0("replicates: ", min(replicates), "-", max(replicates)),
  paste0("tasks: ", nrow(manifest))
)
writeLines(receipt, file.path(output, "CONSOLIDATION.txt"))
cat("Consolidated", nrow(manifest), "tasks at", unique(git_sha), "into",
    output, "\n")
