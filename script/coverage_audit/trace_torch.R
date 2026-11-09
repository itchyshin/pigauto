# covr does not discover torch's R6 generator hidden on a factory attribute.
# Expose the actual generator to covr's supported R6 environment traversal.
coverage_trace_graph <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  env <- new.env(parent = asNamespace("testthat"))
  sys.source(file.path(root, "R", "graph_transformer_block.R"),
             envir = env, keep.source = TRUE)
  env$.coverage_graph_generator <- attr(env$GraphTransformerBlock, "module")
  stopifnot(inherits(env$.coverage_graph_generator, "R6ClassGenerator"))
  tests <- file.path(root, "tests", "testthat",
                    c("helper-libtorch.R", "test-graph-transformer-block.R"))
  extra <- file.path(root, "tests", "testthat", "test-coverage-100-transformer.R")
  if (file.exists(extra)) tests <- c(tests, extra)
  stopifnot(all(file.exists(tests)))
  result <- covr::environment_coverage(env, test_files = tests)
  attr(result, "root") <- root
  tally <- covr::tally_coverage(result, by = "line")
  stopifnot(nrow(tally) > 0L, any(tally$value > 0))
  files <- unique(tally$filename)
  stopifnot(identical(files, "R/graph_transformer_block.R"))
  result
}

coverage_combine_graph <- function(primary, supplemental, root) {
  a <- covr::tally_coverage(primary, by = "line")
  b <- covr::tally_coverage(supplemental, by = "line")
  # An upstream covr fix may discover the class itself. Refuse overlapping
  # instrumentation rather than count expressions twice or discard a counter.
  stopifnot(length(intersect(unique(a$filename), unique(b$filename))) == 0L)
  result <- c(primary, supplemental)
  class(result) <- "coverage"
  attr(result, "root") <- normalizePath(root)
  attr(result, "package") <- attr(primary, "package")
  result
}
