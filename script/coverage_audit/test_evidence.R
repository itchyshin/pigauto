# Extract named test witnesses from covr's retained per-test counter records.
coverage_test_evidence <- function(cov, output) {
  tests <- attr(cov, "tests")
  if (is.null(tests)) stop("Coverage has no per-test records")
  describe <- function(trace) {
    for (call in trace) {
      if (!is.call(call) || length(call) < 2L) next
      head <- paste(deparse(call[[1]]), collapse = "")
      if (head %in% c("test_that", "testthat::test_that") &&
          is.character(call[[2]]) && length(call[[2]]) == 1L) {
        source <- getSrcFilename(call, full.names = FALSE)
        if (!length(source)) source <- "unknown-test-file"
        return(paste0(basename(source), "::", call[[2]]))
      }
    }
    "<setup-or-unattributed-call>"
  }
  labels <- vapply(tests, describe, character(1))
  raw <- as.data.frame(cov, sort = FALSE)
  full <- getFromNamespace("tally_coverage", "covr")(cov, by = "line")
  full <- aggregate(value ~ filename + line, data = full, FUN = max)
  full <- full[full$value > 0, ]
  keys <- paste(full$filename, full$line, sep = ":")
  witnesses <- setNames(vector("list", length(keys)), keys)
  for (i in seq_along(cov)) {
    if (raw$value[i] <= 0 || is.na(raw$first_line[i])) next
    tally <- cov[[i]]$tests
    if (is.null(tally) || !nrow(tally)) next
    ids <- unique(tally[, "test"])
    ids <- ids[!is.na(ids) & ids >= 1L & ids <= length(labels)]
    if (!length(ids)) next
    covered_keys <- intersect(keys, paste(raw$filename[i],
      seq.int(raw$first_line[i], raw$last_line[i]), sep = ":"))
    for (key in covered_keys) witnesses[[key]] <- union(witnesses[[key]], labels[ids])
  }
  full$test_evidence <- vapply(witnesses, function(x) paste(sort(x), collapse = " | "), character(1))
  full$named_test <- vapply(witnesses, function(x)
    any(x != "<setup-or-unattributed-call>"), logical(1))
  write.csv(full, output, row.names = FALSE)
  invisible(full)
}
