coverage_workflow_safe <- function(producer_text, uploader_text) {
  tryCatch({
    p <- yaml::yaml.load(producer_text)
    u <- yaml::yaml.load(uploader_text)
    event <- u[["on"]] %||% u[["TRUE"]]
    job <- u$jobs$upload
    norm <- function(x) gsub("[[:space:]]+", " ", trimws(x))
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
      length(u$permissions) == 0L && identical(names(u$jobs), "upload") &&
      identical(norm(job[["if"]]), expected_if) &&
      identical(job$permissions, list(actions = "read", contents = "read", `id-token` = "write")) &&
      !contains_key(p, "id-token") &&
      !any(grepl("secrets.CODECOV_TOKEN", unlist(p, use.names = FALSE), fixed = TRUE)) &&
      !any(vapply(steps, function(x) !is.null(x$uses) && grepl("^codecov/", x$uses), logical(1))) &&
      identical(gate$env$PIGAUTO_REQUIRE_SOURCE, "true") &&
      grepl('community-surface|coverage-workflow', gate$run, fixed = TRUE) &&
      grepl('stop_on_failure = TRUE', gate$run, fixed = TRUE) &&
      identical(artifact[["if"]], "github.event_name == 'push' && github.ref == 'refs/heads/main'") &&
      identical(artifact$with$name, "main-coverage") &&
      identical(artifact$with$`if-no-files-found`, "error") &&
      identical(checkout$with$ref, "${{ github.event.workflow_run.head_sha }}") &&
      identical(checkout$with$`persist-credentials`, FALSE) &&
      identical(download$with$name, "main-coverage") &&
      identical(download$with$`run-id`, "${{ github.event.workflow_run.id }}") &&
      identical(upload$with$override_commit, "${{ github.event.workflow_run.head_sha }}") &&
      identical(upload$with$override_branch, "main") &&
      identical(upload$with$use_oidc, TRUE) && identical(upload$with$fail_ci_if_error, TRUE) &&
      identical(upload$with$files, "${{ runner.temp }}/coverage/cobertura.xml") &&
      all(vapply(job$steps, function(x) grepl("@[0-9a-f]{40}$", x$uses), logical(1))) &&
      grepl("@[0-9a-f]{40}$", artifact$uses)
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
})
