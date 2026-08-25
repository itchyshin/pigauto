#!/usr/bin/env Rscript

fail <- function(...) stop(sprintf(...), call. = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) fail("Usage: run_pigauto_011_candidate.R [FROZEN_ROOT]")

source_root <- normalizePath(".", mustWork = TRUE)
frozen_root <- if (length(args)) args[[1L]] else file.path(source_root, "checkpoints", "pigauto-011", "frozen")
dir.create(frozen_root, recursive = TRUE, showWarnings = FALSE)
frozen_root <- normalizePath(frozen_root, mustWork = TRUE)

git_out <- function(args) system2("git", args, stdout = TRUE, stderr = TRUE)
head_sha <- function() trimws(git_out(c("rev-parse", "HEAD"))[[1L]])
source_sha <- head_sha()
source_status <- function() paste(git_out(c("status", "--porcelain")), collapse = "\n")
if (nzchar(source_status())) fail("Candidate runner requires a clean source worktree")

desc <- read.dcf(file.path(source_root, "DESCRIPTION"))
if (!identical(unname(desc[1, "Package"]), "pigauto") ||
    !identical(unname(desc[1, "Version"]), "0.11.0")) {
  fail("Source DESCRIPTION is not pigauto 0.11.0")
}

stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
staging <- file.path(source_root, "checkpoints", "pigauto-011", "staging",
                     paste0(stamp, "-", substr(source_sha, 1L, 12L)))
if (file.exists(staging)) fail("Fresh staging directory already exists: %s", staging)
dir.create(staging, recursive = TRUE)
logs_dir <- file.path(staging, "logs")
dir.create(logs_dir)
site_dir <- file.path(staging, "site")
build_dir <- file.path(staging, "build")
check_dir <- file.path(staging, "check")
library_dir <- file.path(staging, "candidate-library")
dir.create(build_dir)
dir.create(check_dir)
dir.create(library_dir)

rscript <- file.path(R.home("bin"), "Rscript")
rbin <- file.path(R.home("bin"), "R")
thread_env <- c("OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1")
commands <- list()

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

run_command <- function(name, command, args, marker, wd = source_root, env = character(), bindings = list()) {
  if (!identical(source_status(), "")) fail("Source became dirty before %s", name)
  head_before <- head_sha()
  if (!identical(head_before, source_sha)) fail("HEAD changed before %s", name)
  log_path <- file.path(logs_dir, paste0(gsub("_", "-", name), ".log"))
  started <- Sys.time()
  cat(sprintf("START %s at %s\n", name, format(started, tz = "UTC", usetz = TRUE)))
  old <- setwd(wd)
  on.exit(setwd(old), add = TRUE)
  status <- system2(command, args, stdout = log_path, stderr = log_path,
                    env = c(thread_env, env), wait = TRUE)
  ended <- Sys.time()
  after <- source_status()
  head_after <- head_sha()
  record <- list(
    status = as.integer(status),
    source_sha = source_sha,
    head_before = head_before,
    head_after = head_after,
    source_status_before = "",
    source_status_after = after,
    started_at = format(started, tz = "UTC", usetz = TRUE),
    ended_at = format(ended, tz = "UTC", usetz = TRUE),
    elapsed_seconds = unname(as.numeric(difftime(ended, started, units = "secs"))),
    command = command,
    args = args,
    working_directory = normalizePath(wd, mustWork = TRUE),
    effective_env = sort(c(thread_env, env)),
    bindings = bindings
  )
  commands[[name]] <<- record
  if (!identical(status, 0L)) fail("Candidate command failed: %s (log: %s)", name, log_path)
  if (nzchar(after)) fail("Candidate command dirtied source: %s", name)
  if (!identical(head_after, source_sha)) fail("HEAD changed during %s", name)
  log_text <- paste(readLines(log_path, warn = FALSE), collapse = "\n")
  if (!grepl(marker, log_text, fixed = TRUE)) fail("Success marker missing after %s", name)
  cat(sprintf("PASS %s in %.1f seconds\n", name, record$elapsed_seconds))
  invisible(log_path)
}

expr_arg <- function(x) c("--vanilla", "-e", shQuote(x))

run_command("document", rscript,
  expr_arg("devtools::document(); cat('R1 documentation verified\\n')"),
  "R1 documentation verified")
run_command("g5", rscript,
  expr_arg("testthat::test_dir('tests/testthat', filter='mi-provenance', stop_on_failure=TRUE, load_package='source', reporter='summary'); cat('G5 MI provenance verified\\n')"),
  "G5 MI provenance verified")
run_command("g6", rscript,
  expr_arg("testthat::test_dir('tests/testthat', filter='check-pigauto', stop_on_failure=TRUE, load_package='source', reporter='summary'); cat('G6 preflight verified\\n')"),
  "G6 preflight verified")
run_command("g7", rscript,
  expr_arg("testthat::test_dir('tests/testthat', filter='result-contract', stop_on_failure=TRUE, load_package='source', reporter='summary'); cat('G7 result invariants verified\\n')"),
  "G7 result invariants verified")
run_command("community", rscript,
  expr_arg("devtools::test(filter='community-surface', reporter='summary', stop_on_failure=TRUE); cat('R1 community surface verified\\n')"),
  "R1 community surface verified")
run_command("full_tests", rscript,
  expr_arg("devtools::test(reporter='summary', stop_on_failure=TRUE); cat('R1 full tests verified\\n')"),
  "R1 full tests verified")
run_command("pkgdown_check", rscript,
  expr_arg("pkgdown::check_pkgdown(); cat('R1 pkgdown check verified\\n')"),
  "R1 pkgdown check verified")

site_expr <- sprintf(
  "pkgdown::build_site(dest_dir=%s, new_process=FALSE, install=FALSE); cat('R1 pkgdown build verified\\n')",
  encodeString(site_dir, quote = "\"")
)
run_command("pkgdown_build", rscript, expr_arg(site_expr), "R1 pkgdown build verified")

run_command("build", rbin, c("CMD", "build", shQuote(source_root)), "* building", wd = build_dir)
tarballs <- list.files(build_dir, pattern = "^pigauto_0[.]11[.]0[.]tar[.]gz$", full.names = TRUE)
if (length(tarballs) != 1L) fail("Build did not produce exactly one pigauto_0.11.0.tar.gz")
tarball <- normalizePath(tarballs[[1L]], mustWork = TRUE)
if (file.info(tarball)$mtime < as.POSIXct(commands$build$started_at)) {
  fail("Built tarball predates the recorded build command")
}

tarball_sha256 <- sha256_file(tarball)
commands$build$bindings <- list(output_tarball_path = tarball,
                                output_tarball_sha256 = tarball_sha256)
run_command("check", rbin,
  c("CMD", "check", "--as-cran", "--run-donttest", shQuote(tarball)),
  "Status: OK", wd = check_dir,
  bindings = list(tarball_path = tarball, tarball_sha256 = tarball_sha256))
run_command("install", rbin,
  c("CMD", "INSTALL", paste0("--library=", shQuote(library_dir)), shQuote(tarball)),
  "* DONE (pigauto)",
  bindings = list(tarball_path = tarball, tarball_sha256 = tarball_sha256,
                  library_path = library_dir))

library_before_gates <- tree_state(library_dir)
site_before_gates <- tree_state(site_dir)
commands$install$bindings$library_digest <- library_before_gates$digest
commands$pkgdown_build$bindings <- list(site_path = site_dir,
                                        site_digest = site_before_gates$digest)

shared_env <- c(
  paste0("PIGAUTO_SOURCE=", shQuote(source_root)),
  paste0("PIGAUTO_TARBALL=", shQuote(tarball)),
  paste0("PIGAUTO_LIBRARY=", shQuote(library_dir))
)
artifact_bindings <- list(
  source_root = source_root,
  tarball_path = tarball,
  tarball_sha256 = tarball_sha256,
  library_path = library_dir,
  library_digest = library_before_gates$digest
)
run_command("g4", rscript, c("--vanilla", shQuote(file.path(source_root, "script", "verify_pigauto_011_surface.R"))),
            "G4 BACE surface verified", env = shared_env, bindings = artifact_bindings)
run_command("g8", rscript, c("--vanilla", shQuote(file.path(source_root, "script", "verify_pigauto_011_workflow.R"))),
            "G8 installed novice workflow verified", env = shared_env, bindings = artifact_bindings)
run_command("g9", rscript, c("--vanilla", shQuote(file.path(source_root, "script", "verify_pigauto_011_claims.R"))),
            "G9 claim boundary verified", env = c(shared_env, paste0("PIGAUTO_SITE=", shQuote(site_dir))),
            bindings = c(artifact_bindings, list(site_path = site_dir,
                                                 site_digest = site_before_gates$digest)))

if (!identical(source_status(), "")) fail("Source was not clean at candidate freeze")
if (!identical(head_sha(), source_sha)) fail("HEAD changed before candidate freeze")
library_final <- tree_state(library_dir)
site_final <- tree_state(site_dir)
if (!identical(library_final$digest, library_before_gates$digest)) fail("Candidate library changed during artifact gates")
if (!identical(site_final$digest, site_before_gates$digest)) fail("Rendered site changed during artifact gates")
sha256 <- sha256_file(tarball)
if (!identical(sha256, tarball_sha256)) fail("Candidate tarball changed during artifact gates")
entries <- utils::untar(tarball, list = TRUE)
target <- file.path(frozen_root, substr(source_sha, 1L, 12L), sha256)
if (file.exists(target)) fail("Frozen candidate directory already exists: %s", target)
partial <- paste0(target, ".partial-", Sys.getpid())
if (file.exists(partial)) fail("Partial freeze path already exists: %s", partial)
dir.create(partial, recursive = TRUE)

copy_tree <- function(from, to) {
  dir.create(to, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(from, recursive = TRUE, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  dirs <- files[file.info(files)$isdir %in% TRUE]
  for (dir in dirs) dir.create(file.path(to, substring(dir, nchar(from) + 2L)), recursive = TRUE, showWarnings = FALSE)
  regular <- files[!file.info(files)$isdir %in% TRUE]
  for (file in regular) {
    dest <- file.path(to, substring(file, nchar(from) + 2L))
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(file, dest, copy.mode = TRUE, copy.date = TRUE)) fail("Could not copy %s", file)
  }
}

dir.create(file.path(partial, "logs"))
log_names <- c(
  document = "document.log", g5 = "g5.log", g6 = "g6.log", g7 = "g7.log",
  community = "community.log", full_tests = "full-tests.log",
  pkgdown_check = "pkgdown-check.log", pkgdown_build = "pkgdown-build.log",
  build = "build.log", check = "check.log", install = "install.log",
  g4 = "g4.log", g8 = "g8.log", g9 = "g9.log"
)
for (name in names(log_names)) {
  from <- file.path(logs_dir, paste0(gsub("_", "-", name), ".log"))
  if (!file.copy(from, file.path(partial, "logs", log_names[[name]]))) fail("Could not freeze %s log", name)
}
if (!file.copy(tarball, file.path(partial, basename(tarball)))) fail("Could not freeze tarball")
copy_tree(library_dir, file.path(partial, "candidate-library"))
copy_tree(site_dir, file.path(partial, "site"))
writeLines(entries, file.path(partial, "inventory.txt"))
writeLines(library_final$inventory, file.path(partial, "candidate-library-inventory.txt"))
writeLines(site_final$inventory, file.path(partial, "pkgdown-site-inventory.txt"))
if (!file.copy(file.path(source_root, "docs", "dev-log", "2026-08-25-pigauto-011-component-ledger.md"),
               file.path(partial, "component-ledger.md"))) fail("Could not freeze component ledger")
writeLines(capture.output(sessionInfo()), file.path(partial, "sessionInfo.txt"))

manifest <- list(
  package = "pigauto", version = "0.11.0", source_sha = source_sha,
  tarball = basename(tarball), sha256 = sha256,
  bytes = unname(file.info(tarball)$size), entry_count = length(entries),
  inventory = "inventory.txt", source_status = "",
  logs = file.path("logs", log_names), commands = commands,
  source_root = source_root, staging_root = staging,
  library_digest = library_final$digest,
  library_inventory = "candidate-library-inventory.txt",
  site_digest = site_final$digest,
  site_inventory = "pkgdown-site-inventory.txt",
  component_ledger = "component-ledger.md", session_info = "sessionInfo.txt",
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
saveRDS(manifest, file.path(partial, "manifest.rds"), version = 3)
receipt <- c(
  "# pigauto 0.11.0 frozen local candidate", "",
  sprintf("- Source commit: `%s`", source_sha),
  sprintf("- Tarball: `%s`", basename(tarball)),
  sprintf("- SHA-256: `%s`", sha256),
  sprintf("- Bytes: %s", format(manifest$bytes, scientific = FALSE)),
  sprintf("- Archive entries: %d", manifest$entry_count),
  sprintf("- Frozen at: %s", manifest$created_at),
  "- Source status at freeze: clean",
  "- Highest proven CRAN rung: tarball-clean on this local macOS host",
  "- Next unproven rung: platform-clean",
  "- Not authorized or claimed: submission-ready, submitted, accepted, or released", ""
)
writeLines(receipt, file.path(partial, "RECEIPT.md"))

dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
if (!file.rename(partial, target)) fail("Atomic candidate freeze rename failed")
pointer <- file.path(source_root, "checkpoints", "pigauto-011", "CURRENT")
writeLines(normalizePath(target), pointer)
cat("FROZEN ", normalizePath(target), "\n", sep = "")
