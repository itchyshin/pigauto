out <- "/private/tmp/pigauto-final-closure-20261017/successor-b2e6054/installed-novice/attempt3"
lib <- "/private/tmp/pigauto-final-closure-20261017/successor-b2e6054/mi-matrix/neither"
.libPaths(c(lib, .libPaths()))
stopifnot(normalizePath(find.package("pigauto")) == normalizePath(file.path(lib, "pigauto")))
stopifnot(as.character(packageVersion("pigauto")) == "0.11.0")

plain_help <- list()
for (topic in c("impute", "pigauto_report")) {
  rd <- utils:::.getHelpFile(utils::help(topic, package = "pigauto", lib.loc = lib))
  txt <- paste(capture.output(tools::Rd2txt(rd)), collapse = "\n")
  txt <- gsub(paste0("_", "\b"), "", txt, fixed = TRUE)
  txt <- gsub("[[:space:]]+", " ", txt)
  plain_help[[topic]] <- txt
  writeLines(txt, file.path(out, paste0("installed-", topic, "-help.txt")))
}
stopifnot(identical(formals(pigauto::impute)$gnn, FALSE))
impute_help <- plain_help$impute
stopifnot(grepl("gnn = FALSE", impute_help, fixed = TRUE))
stopifnot(grepl("baseline alone", impute_help, fixed = TRUE))
stopifnot(grepl("gnn = TRUE", impute_help, fixed = TRUE))
stopifnot(grepl("also train a GNN and calibrate its contribution", impute_help, fixed = TRUE))
report_help <- plain_help$pigauto_report
stopifnot(grepl("For fits created with the default", report_help, fixed = TRUE))
stopifnot(grepl("gnn = FALSE", report_help, fixed = TRUE))
stopifnot(grepl("baseline-only imputation", report_help, fixed = TRUE))
stopifnot(grepl("gnn = TRUE", report_help, fixed = TRUE))
stopifnot(grepl("compares the GNN against the phylogenetic baseline and includes training history", report_help, fixed = TRUE))

run <- file.path(out, "novice-run")
dir.create(run, showWarnings = FALSE)
setwd(run)
browser_requests <- character()
options(browser = function(url, ...) {
  browser_requests <<- c(browser_requests, url)
  invisible(NULL)
})
set.seed(20261018L)
source(system.file("examples", "novice-workflow.R", package = "pigauto"))
stopifnot(isFALSE(result$fit$model_config$gnn))
stopifnot(identical(dim(traits), dim(completed)), identical(names(traits), names(completed)))
for (name in names(traits)) {
  observed <- !is.na(traits[[name]])
  stopifnot(identical(class(traits[[name]]), class(completed[[name]])))
  stopifnot(identical(traits[[name]][observed], completed[[name]][observed]))
}
stopifnot(!anyNA(completed), sum(is.na(traits)) == 2L, sum(is.na(completed)) == 0L)
stopifnot(length(browser_requests) == 1L)
receipt <- list(
  scope = "Installed novice workflow from exact frozen successor archive a93fece7; installed in the successor four-cell matrix's neither library",
  source_commit = "b2e6054685e6bdedb073ba46e058af2c0a619a18",
  package_path = find.package("pigauto"),
  version = as.character(packageVersion("pigauto")),
  impute_gnn_default = formals(pigauto::impute)$gnn,
  impute_help_default_and_opt_in_checked = TRUE,
  report_help_default_and_opt_in_checked = TRUE,
  configured_gnn = result$fit$model_config$gnn,
  torch_namespace_loaded = "torch" %in% loadedNamespaces(),
  gnn_training_disabled_by_config = isFALSE(result$fit$model_config$gnn),
  library_paths = .libPaths(),
  input_missing = sum(is.na(traits)),
  output_missing = sum(is.na(completed)),
  observed_values_and_classes_preserved = TRUE,
  browser_requests = browser_requests,
  browser_scope = "Suppressed only in this R process through options(browser=...); no browser, extension, or system settings changed",
  seed = 20261018L
)
jsonlite::write_json(receipt, file.path(out, "installed-workflow-receipt.json"),
                     auto_unbox = TRUE, pretty = TRUE, null = "null")
capture.output(sessionInfo(), file = file.path(out, "installed-session.txt"))
cat("EXACT_ARCHIVE_INSTALLED_HELP_AND_NOVICE_WORKFLOW_PASS_2_MISSING_FILLED\n")
