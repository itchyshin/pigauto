args <- commandArgs(TRUE)
options(useFancyQuotes = FALSE)
objects <- readRDS(args[[1]])
row <- function(x, id) {
  src <- attr(x$srcref, "srcfile")
  calls <- x$trace$call
  list(id = id, message = x$message, test = x$test,
       file = if (is.null(src$filename)) NULL else basename(src$filename),
       line = if (is.null(x$srcref)) NULL else as.integer(x$srcref[[1]]),
       last_call = if (!length(calls)) NULL else
         paste(deparse(tail(calls, 1L)[[1]]), collapse = " "),
       calls = unname(lapply(calls, function(call) paste(deparse(call), collapse = " "))))
}
out <- lapply(c("warnings", "skips", "problems"), function(key) {
  values <- objects[[key]]
  unname(lapply(seq_along(values), function(i) row(values[[i]], i)))
})
names(out) <- c("warnings", "skips", "problems")
cat(jsonlite::toJSON(out, auto_unbox = TRUE, null = "null", na = "null", digits = NA))
