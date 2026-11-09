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
  stopifnot(nrow(raw) == length(cov),
    isTRUE(all.equal(raw$value, vapply(cov, function(z) z$value, numeric(1)),
                     check.attributes = FALSE)),
    identical(as.integer(raw$first_line),
      unname(vapply(cov, function(z) as.integer(z$srcref[1]), integer(1)))))
  full <- getFromNamespace("tally_coverage", "covr")(cov, by = "line")
  full <- full[full$value > 0, ]
  group_key <- function(filename, fun, line) paste(filename, fun, line, sep = "\r")
  keys <- group_key(full$filename, full$functions, full$line)
  groups <- setNames(vector("list", length(keys)), keys)
  counter_tests <- vector("list", length(cov))
  for (i in seq_along(cov)) {
    if (is.na(raw$first_line[i])) next
    tally <- cov[[i]]$tests
    ids <- if (is.null(tally) || !nrow(tally)) integer() else unique(tally[, "test"])
    counter_tests[[i]] <- ids[!is.na(ids) & ids >= 1L & ids <= length(labels)]
    covered_keys <- intersect(keys, group_key(raw$filename[i], raw$functions[i],
      seq.int(raw$first_line[i], raw$last_line[i])))
    for (key in covered_keys) groups[[key]] <- c(groups[[key]], i)
  }
  full$single_test_evidence <- full$collective_test_evidence <- ""
  full$named_test <- FALSE
  for (i in seq_len(nrow(full))) {
    counters <- groups[[keys[i]]]
    sets <- counter_tests[counters]
    if (!length(sets)) next
    # Covr takes the minimum across overlapping counters within a function.
    # A single-test witness must hit every such counter, not just its parent.
    single <- Reduce(intersect, sets)
    collective <- Reduce(union, sets)
    full$single_test_evidence[i] <- paste(sort(unique(labels[single])), collapse = " | ")
    full$collective_test_evidence[i] <- paste(sort(unique(labels[collective])), collapse = " | ")
    full$named_test[i] <- any(labels[single] != "<setup-or-unattributed-call>")
  }
  write.csv(full, output, row.names = FALSE)
  invisible(full)
}
