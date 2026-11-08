test_that("coverage uploads require a successful trusted main push", {
  root <- testthat::test_path("..", "..")
  producer <- file.path(root, ".github", "workflows", "test-coverage.yaml")
  uploader <- file.path(root, ".github", "workflows", "codecov-upload.yaml")
  if (!file.exists(producer) || !file.exists(uploader)) {
    if (file.exists(file.path(root, "R", "report.R"))) fail("coverage workflow file missing from checkout")
    skip("workflow sources are unavailable in the installed package")
  }
  p <- paste(readLines(producer, warn = FALSE), collapse = "\n")
  u <- paste(readLines(uploader, warn = FALSE), collapse = "\n")
  safe_producer <- function(x) {
    grepl("github.event_name == 'push' && github.ref == 'refs/heads/main'", x, fixed = TRUE) &&
      grepl("if-no-files-found: error", x, fixed = TRUE) &&
      !grepl("id-token:", x, fixed = TRUE) && !grepl("codecov/codecov-action@", x, fixed = TRUE)
  }
  safe_uploader <- function(x) {
    !grepl("pull_request:", x, fixed = TRUE) &&
      !grepl("pull_request_target:", x, fixed = TRUE) &&
      !grepl("workflow_dispatch:", x, fixed = TRUE) &&
      all(vapply(c("workflow_run:", "workflows: [test-coverage]", "types: [completed]", "branches: [main]",
      "github.event.workflow_run.event == 'push'", "github.event.workflow_run.head_branch == 'main'",
      "github.event.workflow_run.head_repository.full_name == github.repository",
      "github.event.workflow_run.conclusion == 'success'", "permissions: {}", "id-token: write",
      "run-id: ${{ github.event.workflow_run.id }}", "override_commit: ${{ github.event.workflow_run.head_sha }}",
      "use_oidc: true", "fail_ci_if_error: true", "actions: read", "contents: read",
      "persist-credentials: false", "ref: ${{ github.event.workflow_run.head_sha }}"), grepl, logical(1), x = x, fixed = TRUE))
  }
  expect_true(safe_producer(p))
  expect_true(safe_uploader(u))
  expect_false(safe_producer(paste(p, "id-token: write")))
  expect_false(safe_uploader(paste(u, "pull_request:")))
  expect_false(safe_uploader(sub("github.event.workflow_run.event == 'push'", "true", u, fixed = TRUE)))
  expect_false(safe_uploader(sub("github.event.workflow_run.conclusion == 'success'", "true", u, fixed = TRUE)))
})
