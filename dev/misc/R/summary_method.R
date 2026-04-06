
#' @export
summary.phyloimpute_result <- function(object, ...) {
  cat('phyloimpute result\n')
  cat(sprintf('Rows: %d  Columns: %d\n',
              nrow(object$completed_data), ncol(object$completed_data)))
  if (!is.null(object$note)) cat('Note:', object$note, '\n')
}

