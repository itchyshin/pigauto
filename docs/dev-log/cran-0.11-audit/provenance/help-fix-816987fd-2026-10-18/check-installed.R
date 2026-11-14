out <- "/private/tmp/pigauto-final-closure-20261017/help-fix"
lib <- file.path(out, "installed-library")
.libPaths(c(lib, .libPaths()))
stopifnot(normalizePath(find.package("pigauto")) == normalizePath(file.path(lib, "pigauto")))
stopifnot(as.character(packageVersion("pigauto")) == "0.11.0")
for (topic in c("impute", "pigauto_report")) {
  rd <- utils:::.getHelpFile(utils::help(topic, package = "pigauto", lib.loc = lib))
  writeLines(capture.output(tools::Rd2txt(rd)),
             file.path(out, paste0("installed-", topic, "-help.txt")))
}
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
stopifnot(identical(formals(pigauto::impute)$gnn, FALSE))
stopifnot(isFALSE(result$fit$model_config$gnn))
stopifnot(identical(dim(traits), dim(completed)), identical(names(traits), names(completed)))
for (name in names(traits)) {
  observed <- !is.na(traits[[name]])
  stopifnot(identical(class(traits[[name]]), class(completed[[name]])))
  stopifnot(identical(traits[[name]][observed], completed[[name]][observed]))
}
stopifnot(!anyNA(completed), sum(is.na(traits)) == 2L)
receipt <- list(scope = "Provisional source install, not a frozen release archive",
                package_path = find.package("pigauto"), version = as.character(packageVersion("pigauto")),
                gnn_default = formals(pigauto::impute)$gnn, configured_gnn = result$fit$model_config$gnn,
                input_missing = sum(is.na(traits)), output_missing = sum(is.na(completed)),
                observed_values_and_classes_preserved = TRUE,
                browser_requests = browser_requests,
                browser_scope = "Suppressed only in this headless R process; no browser or extension settings changed",
                seed = 20261018L)
jsonlite::write_json(receipt, file.path(out, "installed-workflow-receipt.json"),
                     auto_unbox = TRUE, pretty = TRUE, null = "null")
capture.output(sessionInfo(), file = file.path(out, "installed-session.txt"))
cat("PROVISIONAL_INSTALLED_NOVICE_WORKFLOW_PASS_2_MISSING_FILLED\n")
