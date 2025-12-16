
#' Create boolean mask of observed entries
#' @noRd
mask_matrix <- function(mat) {
  !is.na(mat)
}

