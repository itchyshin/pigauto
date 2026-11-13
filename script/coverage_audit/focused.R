args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
filters <- c(presentation = "^(plot-contract|coverage-100-presentation)$",
             inputs = "^coverage-100-inputs$",
             prediction = "^coverage-100-prediction$",
             gnn = "^coverage-100-(gnn|transformer)$", mi = "^coverage-100-mi$")
stopifnot(args[1] %in% names(filters))
r <- devtools::test(filter = filters[[args[1]]], stop_on_failure = TRUE,
                    reporter = "summary")
counts <- as.data.frame(r)
stopifnot(nrow(counts) > 0L, sum(counts$nb) > 0L,
          !any(counts$failed > 0L), !any(counts$error),
          !any(counts$warning > 0L), !any(counts$skipped))
cat("FOCUSED_PASS", args[1], "assertions", sum(counts$nb), "\n")
