coverage_workflow_safe <- function(producer_text, uploader_text) {
  tryCatch({
    p <- yaml::yaml.load(producer_text)
    u <- yaml::yaml.load(uploader_text)
    event_key <- if ("on" %in% names(u)) "on" else "TRUE"
    event <- u[[event_key]]
    job <- u$jobs$upload
    norm <- function(x) gsub("[[:space:]]+", " ", trimws(x))
    exact_names <- function(x, expected) length(names(x)) == length(expected) && setequal(names(x), expected)
    step <- function(steps, name) {
      found <- Filter(function(x) identical(x$name, name), steps)
      if (length(found) != 1L) stop("Expected one named step")
      found[[1L]]
    }
    contains_key <- function(x, key) {
      is.list(x) && (key %in% names(x) || any(vapply(x, contains_key, logical(1), key = key)))
    }
    steps <- p$jobs[["test-coverage"]]$steps
    gate <- step(steps, "Check public source surfaces")
    artifact <- step(steps, "Save coverage report for trusted upload")
    checkout <- step(job$steps, "Check out the tested commit")
    download <- step(job$steps, "Download the trusted coverage artifact")
    upload <- step(job$steps, "Upload coverage to Codecov")
    expected_if <- paste(
      "github.event.workflow_run.event == 'push' &&",
      "github.event.workflow_run.head_branch == 'main' &&",
      "github.event.workflow_run.head_repository.full_name == github.repository &&",
      "github.event.workflow_run.conclusion == 'success'"
    )
    identical(names(event), "workflow_run") &&
      identical(event$workflow_run$workflows, "test-coverage") &&
      identical(event$workflow_run$types, "completed") &&
      identical(event$workflow_run$branches, "main") &&
      exact_names(u, c("name", event_key, "permissions", "jobs")) &&
      length(u$permissions) == 0L && identical(names(u$jobs), "upload") &&
      exact_names(job, c("if", "runs-on", "permissions", "steps")) &&
      identical(job$`runs-on`, "ubuntu-latest") &&
      identical(norm(job[["if"]]), expected_if) &&
      identical(job$permissions, list(actions = "read", contents = "read", `id-token` = "write")) &&
      !contains_key(p, "id-token") &&
      !any(grepl("secrets.CODECOV_TOKEN", unlist(p, use.names = FALSE), fixed = TRUE)) &&
      !any(vapply(steps, function(x) !is.null(x$uses) && grepl("^codecov/", x$uses), logical(1))) &&
      identical(gate$env$PIGAUTO_REQUIRE_SOURCE, "true") &&
      grepl('community-surface|coverage-workflow', gate$run, fixed = TRUE) &&
      grepl('stop_on_failure = TRUE', gate$run, fixed = TRUE) &&
      identical(artifact[["if"]], "github.event_name == 'push' && github.ref == 'refs/heads/main'") &&
      exact_names(artifact, c("name", "if", "uses", "with")) &&
      exact_names(artifact$with, c("name", "path", "if-no-files-found", "retention-days")) &&
      identical(artifact$with$name, "main-coverage") &&
      identical(artifact$with$path, "cobertura.xml") &&
      identical(artifact$with$`if-no-files-found`, "error") &&
      all(vapply(job$steps, exact_names, logical(1), expected = c("name", "uses", "with"))) &&
      exact_names(checkout$with, c("ref", "persist-credentials")) &&
      exact_names(download$with, c("name", "path", "run-id", "github-token")) &&
      exact_names(upload$with, c("use_oidc", "fail_ci_if_error", "files", "disable_search", "plugins", "version", "override_branch", "override_commit", "override_build_url")) &&
      identical(checkout$with$ref, "${{ github.event.workflow_run.head_sha }}") &&
      identical(checkout$with$`persist-credentials`, FALSE) &&
      identical(download$with$name, "main-coverage") &&
      identical(download$with$path, "${{ runner.temp }}/coverage") &&
      identical(download$with$`run-id`, "${{ github.event.workflow_run.id }}") &&
      identical(download$with$`github-token`, "${{ github.token }}") &&
      identical(upload$with$override_commit, "${{ github.event.workflow_run.head_sha }}") &&
      identical(upload$with$override_branch, "main") &&
      identical(upload$with$override_build_url, "${{ github.event.workflow_run.html_url }}") &&
      identical(upload$with$use_oidc, TRUE) && identical(upload$with$fail_ci_if_error, TRUE) &&
      identical(upload$with$disable_search, TRUE) && identical(upload$with$plugins, "noop") &&
      identical(upload$with$version, "v11.3.2") &&
      identical(upload$with$files, "${{ runner.temp }}/coverage/cobertura.xml") &&
      identical(vapply(job$steps, `[[`, character(1), "name"),
                c("Check out the tested commit", "Download the trusted coverage artifact", "Upload coverage to Codecov")) &&
      identical(vapply(job$steps, `[[`, character(1), "uses"), c(
        "actions/checkout@11d5960a326750d5838078e36cf38b85af677262",
        "actions/download-artifact@d3f86a106a0bac45b974a628896c90dbdf5c8093",
        "codecov/codecov-action@0fb7174895f61a3b6b78fc075e0cd60383518dac")) &&
      identical(artifact$uses, "actions/upload-artifact@ea165f8d65b6e75b540449e92b4886f43607fa02")
  }, error = function(e) FALSE)
}

test_that("coverage uploads require a successful trusted main push", {
  root <- testthat::test_path("..", "..")
  producer <- file.path(root, ".github", "workflows", "test-coverage.yaml")
  uploader <- file.path(root, ".github", "workflows", "codecov-upload.yaml")
  if (!file.exists(producer) || !file.exists(uploader)) {
    if (file.exists(file.path(root, "R", "report.R"))) fail("coverage workflow file missing from checkout")
    skip("workflow sources are unavailable in the installed package")
  }
  skip_if_not_installed("yaml")
  p <- paste(readLines(producer, warn = FALSE), collapse = "\n")
  u <- paste(readLines(uploader, warn = FALSE), collapse = "\n")
  expect_true(coverage_workflow_safe(p, u))
  # A comment containing the removed guard must never satisfy the check.
  missing_push <- sub("github.event.workflow_run.event == 'push' &&", "", u, fixed = TRUE)
  missing_push <- paste(missing_push, "# github.event.workflow_run.event == 'push' &&", sep = "\n")
  expect_false(coverage_workflow_safe(p, missing_push))
  missing_success <- sub("github.event.workflow_run.conclusion == 'success'", "true", u, fixed = TRUE)
  missing_success <- paste(missing_success, "# github.event.workflow_run.conclusion == 'success'", sep = "\n")
  expect_false(coverage_workflow_safe(p, missing_success))
  extra_trigger <- sub("  workflow_run:", "  pull_request:\n  workflow_run:", u, fixed = TRUE)
  expect_false(coverage_workflow_safe(p, extra_trigger))
  privileged_producer <- sub("permissions:\n  contents: read", "permissions:\n  contents: read\n  id-token: write", p, fixed = TRUE)
  expect_false(coverage_workflow_safe(privileged_producer, u))
  secret_producer <- sub("GITHUB_PAT:", "CODECOV: ${{ secrets.CODECOV_TOKEN }}\n      GITHUB_PAT:", p, fixed = TRUE)
  expect_false(coverage_workflow_safe(secret_producer, u))
  evil_uploader <- sub("codecov/codecov-action@", "evil/uploader@", u, fixed = TRUE)
  expect_false(coverage_workflow_safe(p, evil_uploader))
  evil_downloader <- sub("actions/download-artifact@", "evil/downloader@", u, fixed = TRUE)
  expect_false(coverage_workflow_safe(p, evil_downloader))
  evil_checkout <- sub("actions/checkout@", "evil/checkout@", u, fixed = TRUE)
  expect_false(coverage_workflow_safe(p, evil_checkout))
  evil_artifact <- sub("actions/upload-artifact@ea165f8d65b6e75b540449e92b4886f43607fa02", "evil/artifact@ea165f8d65b6e75b540449e92b4886f43607fa02", p, fixed = TRUE)
  expect_false(coverage_workflow_safe(evil_artifact, u))
  extra_step <- sub("      - name: Upload coverage to Codecov", "      - name: Untrusted step\n        run: echo bad\n      - name: Upload coverage to Codecov", u, fixed = TRUE)
  expect_false(coverage_workflow_safe(p, extra_step))
  searches_files <- sub("disable_search: true", "disable_search: false", u, fixed = TRUE)
  expect_false(coverage_workflow_safe(p, searches_files))
  altered_plugins <- sub("plugins: noop", "plugins: gcov", u, fixed = TRUE)
  expect_false(coverage_workflow_safe(p, altered_plugins))
  wrong_producer_path <- sub("path: cobertura.xml", "path: other.xml", p, fixed = TRUE)
  expect_false(coverage_workflow_safe(wrong_producer_path, u))
  wrong_download_path <- sub("path: ${{ runner.temp }}/coverage", "path: unexpected-dir", u, fixed = TRUE)
  expect_false(coverage_workflow_safe(p, wrong_download_path))
  malicious_url <- sub("use_oidc: true", "url: https://example.invalid\n          use_oidc: true", u, fixed = TRUE)
  expect_false(coverage_workflow_safe(p, malicious_url))
  malicious_env <- sub("    runs-on: ubuntu-latest", "    env:\n      UNTRUSTED: yes\n    runs-on: ubuntu-latest", u, fixed = TRUE)
  expect_false(coverage_workflow_safe(p, malicious_env))
  other_repo <- sub("          ref: ${{ github.event.workflow_run.head_sha }}", "          repository: example/other\n          ref: ${{ github.event.workflow_run.head_sha }}", u, fixed = TRUE)
  expect_false(coverage_workflow_safe(p, other_repo))
})
