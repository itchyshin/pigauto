#!/usr/bin/env Rscript

fail <- function(...) stop(sprintf(...), call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) fail("Usage: verify_pigauto_011_candidate.R [FROZEN_DIR]")
if (!length(args)) {
  pointer <- file.path("checkpoints", "pigauto-011", "CURRENT")
  if (!file.exists(pointer)) fail("No frozen CURRENT candidate exists")
  args <- readLines(pointer, n = 1L, warn = FALSE)
}
frozen <- normalizePath(args[[1L]], mustWork = TRUE)
manifest_path <- file.path(frozen, "manifest.rds")
if (!file.exists(manifest_path)) fail("Frozen candidate manifest is missing")
manifest <- readRDS(manifest_path)

required <- c("package", "version", "source_sha", "tarball", "sha256", "bytes",
              "entry_count", "inventory", "source_status", "logs", "site_inventory",
              "component_ledger", "session_info", "commands", "source_root",
              "staging_root", "library_digest", "library_inventory", "site_digest", "created_at")
missing <- setdiff(required, names(manifest))
if (length(missing)) fail("Manifest field(s) missing: %s", paste(missing, collapse = ", "))
if (!identical(manifest$package, "pigauto") || !identical(manifest$version, "0.11.0")) {
  fail("Manifest package/version is not pigauto 0.11.0")
}
if (nzchar(manifest$source_status)) fail("Frozen source was not clean")

tarball <- normalizePath(file.path(frozen, manifest$tarball), mustWork = TRUE)
sha_line <- system2("shasum", c("-a", "256", shQuote(tarball)), stdout = TRUE)
sha <- strsplit(sha_line[[1L]], "[[:space:]]+")[[1L]][[1L]]
if (!identical(sha, manifest$sha256)) fail("Tarball SHA-256 does not match manifest")
if (!identical(unname(file.info(tarball)$size), manifest$bytes)) fail("Tarball byte size does not match manifest")
entries <- utils::untar(tarball, list = TRUE)
if (!identical(length(entries), manifest$entry_count)) fail("Tarball entry count does not match manifest")
inventory <- readLines(file.path(frozen, manifest$inventory), warn = FALSE)
if (!identical(entries, inventory)) fail("Tarball inventory does not match manifest")
for (field in c("site_inventory", "library_inventory", "component_ledger", "session_info")) {
  path <- file.path(frozen, manifest[[field]])
  if (!file.exists(path) || !file.info(path)$size) fail("Frozen %s is missing or empty", field)
}

sha256_file <- function(path) {
  line <- system2("shasum", c("-a", "256", shQuote(normalizePath(path, mustWork = TRUE))), stdout = TRUE)
  strsplit(line[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
tree_state <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  paths <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  paths <- sort(paths[!file.info(paths)$isdir %in% TRUE])
  relative <- substring(paths, nchar(root) + 2L)
  hashes <- vapply(paths, sha256_file, character(1))
  sizes <- unname(file.info(paths)$size)
  inventory <- paste(relative, hashes, format(sizes, scientific = FALSE, trim = TRUE), sep = "\t")
  digest_line <- system2("shasum", c("-a", "256"), input = inventory, stdout = TRUE)
  digest <- strsplit(digest_line[[1L]], "[[:space:]]+")[[1L]][[1L]]
  list(digest = digest, inventory = inventory)
}

forbidden <- paste0(
  "(^|/)([.]git|[.]github|BACE|script|dev|checkpoints|docs|",
  "avonet|data-raw|useful|pkgdown|misc)(/|$)"
)
if (any(grepl(forbidden, entries, perl = TRUE))) fail("Forbidden source path is present in tarball")

required_ship <- file.path("pigauto", c(
  "data/avonet300.rda", "data/avonet_full.rda", "data/ctmax_sim.rda",
  "data/tree300.rda", "data/tree_full.rda", "data/trees300.rda",
  "inst/NOTICE", "inst/extdata/legacy_fit_v091.rds",
  "inst/extdata/novice_traits.csv", "inst/extdata/novice_tree.nwk",
  "inst/examples/novice-workflow.R", "man/figures/logo.png"
))
missing_ship <- setdiff(required_ship, entries)
if (length(missing_ship)) fail("Required SHIP component(s) absent: %s", paste(missing_ship, collapse = ", "))

expected_markers <- c(
  document = "R1 documentation verified",
  g5 = "G5 MI provenance verified",
  g6 = "G6 preflight verified",
  g7 = "G7 result invariants verified",
  community = "R1 community surface verified",
  full_tests = "R1 full tests verified",
  pkgdown_check = "R1 pkgdown check verified",
  pkgdown_build = "R1 pkgdown build verified",
  build = "* building",
  check = "Status: OK",
  install = "* DONE (pigauto)",
  g4 = "G4 BACE surface verified",
  g8 = "G8 installed novice workflow verified",
  g9 = "G9 claim boundary verified"
)
if (!identical(sort(names(manifest$logs)), sort(names(expected_markers)))) {
  fail("Manifest log set does not match the required candidate commands")
}
if (!identical(sort(names(manifest$commands)), sort(names(expected_markers)))) {
  fail("Manifest command set does not match the required candidate commands")
}
for (name in names(expected_markers)) {
  command <- manifest$commands[[name]]
  if (!identical(command$status, 0L)) fail("Candidate command did not exit zero: %s", name)
  if (!identical(command$source_sha, manifest$source_sha)) fail("Command/source SHA mismatch: %s", name)
  if (!identical(command$head_before, manifest$source_sha) ||
      !identical(command$head_after, manifest$source_sha)) fail("Command HEAD measurement mismatch: %s", name)
  if (!identical(command$source_status_before, "") || !identical(command$source_status_after, "")) {
    fail("Candidate command did not run against a clean source: %s", name)
  }
  path <- file.path(frozen, manifest$logs[[name]])
  if (!file.exists(path) || !file.info(path)$size) fail("Required log is missing or empty: %s", name)
  log <- paste(readLines(path, warn = FALSE), collapse = "\n")
  if (!grepl(expected_markers[[name]], log, fixed = TRUE)) fail("Required marker missing from %s log", name)
  if (!all(c("OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1") %in% command$effective_env)) {
    fail("Thread caps absent from command environment: %s", name)
  }
}

library_path <- file.path(frozen, "candidate-library")
installed_desc <- read.dcf(file.path(library_path, "pigauto", "DESCRIPTION"))
if (!identical(unname(installed_desc[1, "Version"]), manifest$version)) {
  fail("Frozen installed candidate version does not match manifest")
}
library_state <- tree_state(library_path)
if (!identical(library_state$digest, manifest$library_digest)) fail("Frozen candidate-library content digest mismatch")
library_inventory <- readLines(file.path(frozen, manifest$library_inventory), warn = FALSE)
if (!identical(library_state$inventory, library_inventory)) fail("Frozen candidate-library inventory mismatch")
site_state <- tree_state(file.path(frozen, "site"))
if (!identical(site_state$digest, manifest$site_digest)) fail("Frozen rendered-site content digest mismatch")
site_inventory <- readLines(file.path(frozen, manifest$site_inventory), warn = FALSE)
if (!identical(site_state$inventory, site_inventory)) fail("Rendered site inventory mismatch")

for (name in c("check", "install", "g4", "g8", "g9")) {
  binding <- manifest$commands[[name]]$bindings
  if (!identical(binding$tarball_sha256, manifest$sha256)) fail("Tarball binding mismatch: %s", name)
  command_line <- paste(manifest$commands[[name]]$args, collapse = " ")
  if (!grepl(shQuote(binding$tarball_path), command_line, fixed = TRUE) && name %in% c("check", "install")) {
    fail("Recorded command does not name its bound tarball: %s", name)
  }
}
for (name in c("g4", "g8", "g9")) {
  binding <- manifest$commands[[name]]$bindings
  if (!identical(binding$source_root, manifest$source_root)) fail("Source-root binding mismatch: %s", name)
  if (!identical(binding$library_digest, manifest$library_digest)) fail("Library binding mismatch: %s", name)
  if (!identical(binding$library_path, manifest$commands$install$bindings$library_path)) {
    fail("Installed-library path binding mismatch: %s", name)
  }
  expected_env <- c(
    paste0("PIGAUTO_SOURCE=", shQuote(binding$source_root)),
    paste0("PIGAUTO_TARBALL=", shQuote(binding$tarball_path)),
    paste0("PIGAUTO_LIBRARY=", shQuote(binding$library_path))
  )
  if (!all(expected_env %in% manifest$commands[[name]]$effective_env)) {
    fail("Artifact environment binding mismatch: %s", name)
  }
}
if (!identical(manifest$commands$g9$bindings$site_digest, manifest$site_digest)) {
  fail("G9 rendered-site binding mismatch")
}
if (!(paste0("PIGAUTO_SITE=", shQuote(manifest$commands$g9$bindings$site_path)) %in%
      manifest$commands$g9$effective_env)) fail("G9 site environment binding mismatch")
if (!identical(manifest$commands$build$bindings$output_tarball_sha256, manifest$sha256)) {
  fail("Build output tarball binding mismatch")
}
if (!identical(manifest$commands$build$bindings$output_tarball_path,
               manifest$commands$check$bindings$tarball_path) ||
    !identical(manifest$commands$check$bindings$tarball_path,
               manifest$commands$install$bindings$tarball_path)) {
  fail("Build, check and install did not bind the same staging tarball")
}
if (!identical(manifest$commands$install$bindings$library_digest, manifest$library_digest)) {
  fail("Install output library binding mismatch")
}
if (!identical(manifest$commands$pkgdown_build$bindings$site_digest, manifest$site_digest)) {
  fail("Pkgdown output site binding mismatch")
}
if (!identical(manifest$commands$pkgdown_build$bindings$site_path,
               manifest$commands$g9$bindings$site_path)) {
  fail("Pkgdown build and G9 did not bind the same rendered site")
}

expr_arg <- function(x) c("--vanilla", "-e", shQuote(x))
rscript <- file.path(R.home("bin"), "Rscript")
rbin <- file.path(R.home("bin"), "R")
staging_tarball <- manifest$commands$build$bindings$output_tarball_path
staging_library <- manifest$commands$install$bindings$library_path
staging_site <- manifest$commands$pkgdown_build$bindings$site_path
site_expr <- sprintf(
  "pkgdown::build_site(dest_dir=%s, new_process=FALSE, install=FALSE); cat('R1 pkgdown build verified\\n')",
  encodeString(staging_site, quote = "\"")
)
expected_args <- list(
  document = expr_arg("devtools::document(); cat('R1 documentation verified\\n')"),
  g5 = expr_arg("testthat::test_dir('tests/testthat', filter='mi-provenance', stop_on_failure=TRUE, load_package='source', reporter='summary'); cat('G5 MI provenance verified\\n')"),
  g6 = expr_arg("testthat::test_dir('tests/testthat', filter='check-pigauto', stop_on_failure=TRUE, load_package='source', reporter='summary'); cat('G6 preflight verified\\n')"),
  g7 = expr_arg("testthat::test_dir('tests/testthat', filter='result-contract', stop_on_failure=TRUE, load_package='source', reporter='summary'); cat('G7 result invariants verified\\n')"),
  community = expr_arg("devtools::test(filter='community-surface', reporter='summary', stop_on_failure=TRUE); cat('R1 community surface verified\\n')"),
  full_tests = expr_arg("devtools::test(reporter='summary', stop_on_failure=TRUE); cat('R1 full tests verified\\n')"),
  pkgdown_check = expr_arg("pkgdown::check_pkgdown(); cat('R1 pkgdown check verified\\n')"),
  pkgdown_build = expr_arg(site_expr),
  build = c("CMD", "build", shQuote(manifest$source_root)),
  check = c("CMD", "check", "--as-cran", "--run-donttest", shQuote(staging_tarball)),
  install = c("CMD", "INSTALL", paste0("--library=", shQuote(staging_library)), shQuote(staging_tarball)),
  g4 = c("--vanilla", shQuote(file.path(manifest$source_root, "script", "verify_pigauto_011_surface.R"))),
  g8 = c("--vanilla", shQuote(file.path(manifest$source_root, "script", "verify_pigauto_011_workflow.R"))),
  g9 = c("--vanilla", shQuote(file.path(manifest$source_root, "script", "verify_pigauto_011_claims.R")))
)
for (name in names(expected_args)) {
  expected_executable <- if (name %in% c("build", "check", "install")) rbin else rscript
  if (!identical(manifest$commands[[name]]$command, expected_executable) ||
      !identical(manifest$commands[[name]]$args, expected_args[[name]])) {
    fail("Recorded candidate command is not the exact approved operation: %s", name)
  }
}

source_wd_commands <- setdiff(names(expected_args), c("build", "check"))
if (any(vapply(manifest$commands[source_wd_commands], function(x) {
  !identical(x$working_directory, manifest$source_root)
}, logical(1)))) fail("A candidate command used an unexpected source working directory")
if (!identical(manifest$commands$build$working_directory, file.path(manifest$staging_root, "build")) ||
    !identical(manifest$commands$check$working_directory, file.path(manifest$staging_root, "check")) ||
    !identical(dirname(staging_tarball), file.path(manifest$staging_root, "build")) ||
    !identical(staging_library, file.path(manifest$staging_root, "candidate-library")) ||
    !identical(staging_site, file.path(manifest$staging_root, "site"))) {
  fail("Build or check command used an unexpected staging directory")
}
install_library_arg <- paste0("--library=", shQuote(manifest$commands$install$bindings$library_path))
if (!(install_library_arg %in% manifest$commands$install$args)) {
  fail("Install command did not name its bound candidate library")
}

cat("G10 frozen candidate verified\n")
