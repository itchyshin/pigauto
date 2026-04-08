#' Save a fitted pigauto model
#'
#' Serialises a \code{pigauto_fit} object to disk, including the torch
#' model state.  The saved file can be loaded on any machine with the
#' same pigauto version.
#'
#' @details
#' Standard \code{saveRDS()} does not work for pigauto fits because
#' torch tensor states cannot be serialised by R's native mechanism.
#' This function converts the model state to a portable raw format
#' before saving.
#'
#' @param fit pigauto_fit object.
#' @param path character.  File path to save to (recommended extension:
#'   \code{.pigauto}).
#' @param compress logical.  Use gzip compression (default \code{TRUE}).
#' @return Invisible \code{path}.
#' @examples
#' \dontrun{
#' save_pigauto(fit, "my_model.pigauto")
#' fit2 <- load_pigauto("my_model.pigauto")
#' }
#' @export
save_pigauto <- function(fit, path, compress = TRUE) {
  if (!inherits(fit, "pigauto_fit")) {
    stop("'fit' must be a pigauto_fit object.")
  }

  # Convert torch state_dict to a serialisable list of raw vectors
  state <- fit$model_state
  serial_state <- lapply(state, function(t) {
    torch::torch_serialize(t)
  })

  # Build a portable list (no torch tensors)
  portable <- fit
  portable$model_state <- serial_state
  portable$.pigauto_version <- as.character(utils::packageVersion("pigauto"))

  saveRDS(portable, file = path, compress = compress)
  if (compress) {
    message("Model saved to ", path, " (",
            format(file.size(path), big.mark = ","), " bytes)")
  }
  invisible(path)
}


#' Load a saved pigauto model
#'
#' Reads a \code{pigauto_fit} object previously saved with
#' \code{\link{save_pigauto}}.
#'
#' @param path character.  File path to load from.
#' @return An object of class \code{"pigauto_fit"}.
#' @examples
#' \dontrun{
#' fit <- load_pigauto("my_model.pigauto")
#' pred <- predict(fit)
#' }
#' @export
load_pigauto <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)

  portable <- readRDS(path)

  if (!inherits(portable, "pigauto_fit")) {
    stop("File does not contain a pigauto_fit object.")
  }

  # Deserialise torch state_dict (raw vectors back to tensors)
  portable$model_state <- lapply(portable$model_state, function(raw) {
    torch::torch_load(raw)
  })

  # Version check
  saved_ver <- portable$.pigauto_version
  current_ver <- as.character(utils::packageVersion("pigauto"))
  if (!is.null(saved_ver) && saved_ver != current_ver) {
    message("Note: model was saved with pigauto ", saved_ver,
            "; current version is ", current_ver)
  }
  portable$.pigauto_version <- NULL

  portable
}
