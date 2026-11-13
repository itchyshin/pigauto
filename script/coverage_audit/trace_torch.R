# The shipped helper registers torch's hidden generators during the full suite.
# Independently require counters for each model method in the returned coverage.
coverage_torch_inventory <- function(cov, root) {
  root <- normalizePath(root, mustWork = TRUE)
  pkgload::load_all(root, quiet = TRUE, helpers = FALSE, export_all = FALSE)
  ns <- asNamespace("pigauto")
  tally <- covr::tally_coverage(cov, by = "line")
  methods <- list()
  for (factory in c("GraphTransformerBlock", "ResidualPhyloDAE")) {
    generator <- attr(get(factory, envir = ns), "module")
    stopifnot(inherits(generator, "R6ClassGenerator"))
    for (method in c("initialize", "forward")) {
      ref <- attr(generator$public_methods[[method]], "srcref")
      stopifnot(inherits(ref, "srcref"))
      file <- paste0("R/", basename(attr(ref, "srcfile")$filename))
      first <- as.integer(ref[1]); last <- as.integer(ref[3])
      idx <- tally$filename == file & tally$line >= first & tally$line <= last
      stopifnot(any(idx))
      methods[[length(methods) + 1L]] <- list(factory = factory, method = method,
        file = file, first_line = first, last_line = last,
        tracked_lines = sort(unique(tally$line[idx])),
        covered_lines = sort(unique(tally$line[idx & tally$value > 0])))
    }
  }
  methods
}
