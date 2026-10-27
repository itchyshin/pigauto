#' Save a fitted pigauto model
#'
#' Serialises a \code{pigauto_fit} object to disk, including the torch
#' model state, for later loading with \code{load_pigauto()}. Matching
#' pigauto versions alone do not guarantee cross-version or cross-platform
#' compatibility.
#'
#' @details
#' Use \code{save_pigauto()} and \code{load_pigauto()} as the package's
#' save/load interface. For GNN fits containing torch tensors, native
#' \code{saveRDS()} / \code{readRDS()} does not preserve a usable tensor
#' state; this function serialises that state before saving. pigauto checks
#' its package version after loading the tensors, but does not check R,
#' torch/LibTorch, or operating-system compatibility. Retain the software
#' environment used to save the fit and verify loading and predictions in
#' the destination environment. GNN-off fits have an empty model state and
#' do not use this tensor serialization path.
#'
#' @param fit pigauto_fit object, or a \code{pigauto_result} from
#'   \code{\link{impute}} (its \code{$fit} is saved).
#' @param path character.  File path to save to (recommended extension:
#'   \code{.pigauto}).
#' @param compress logical.  Use gzip compression (default \code{TRUE}).
#' @return Invisible \code{path}.
#' @examples
#' \donttest{
#' data(avonet300, tree300)
#' tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
#' traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
#'                      c("Mass", "Wing.Length"), drop = FALSE]
#' rownames(traits) <- tree$tip.label
#' traits$Mass[seq_len(3L)] <- NA_real_
#' fit <- impute(traits, tree, epochs = 5L, verbose = FALSE)$fit
#' path <- tempfile(fileext = ".pigauto")
#' suppressMessages(save_pigauto(fit, path))
#' fit2 <- load_pigauto(path)
#' }
#' @export
save_pigauto <- function(fit, path, compress = TRUE) {
  if (inherits(fit, "pigauto_result")) {
    fit <- fit$fit
  }
  if (!inherits(fit, "pigauto_fit")) {
    stop("'fit' must be a pigauto_fit or pigauto_result object.")
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
#' \donttest{
#' data(avonet300, tree300)
#' tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
#' traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
#'                      c("Mass", "Wing.Length"), drop = FALSE]
#' rownames(traits) <- tree$tip.label
#' traits$Mass[seq_len(3L)] <- NA_real_
#' fit <- impute(traits, tree, epochs = 5L, verbose = FALSE)$fit
#' path <- tempfile(fileext = ".pigauto")
#' suppressMessages(save_pigauto(fit, path))
#' pred <- predict(load_pigauto(path))
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
