#' Read trait data from a CSV file or data frame
#'
#' Loads trait data and sets species names as row names.  If a CSV path is
#' supplied, it is read with \code{read.csv}.  By default, all non-species
#' columns are returned (numeric, factor, integer, etc.).  Use
#' \code{trait_cols} to select a subset.
#'
#' @param x character path to a CSV file, or a \code{data.frame}.
#' @param species_col character. Name of the column containing species names
#'   (default \code{"species"}).
#' @param trait_cols character vector of column names to include. If
#'   \code{NULL} (default), all columns except \code{species_col} are used.
#' @param factor_cols character vector of columns to coerce to \code{factor}.
#' @param ordered_cols character vector of columns to coerce to \code{ordered}.
#' @return A \code{data.frame} with species names as row names.
#' @examples
#' df <- data.frame(species = c("Sp_a", "Sp_b"), mass = c(10, 20))
#' traits <- read_traits(df, species_col = "species")
#' @importFrom utils read.csv
#' @export
read_traits <- function(x, species_col = "species", trait_cols = NULL,
                        factor_cols = NULL, ordered_cols = NULL) {
  if (is.character(x)) {
    if (!file.exists(x)) stop("Trait file not found: ", x)
    x <- read.csv(x, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(x)) {
    stop("'x' must be a data.frame or a path to a CSV file.")
  }
  if (!(species_col %in% names(x))) {
    stop("Column '", species_col, "' not found. Available columns: ",
         paste(names(x), collapse = ", "))
  }

  species <- as.character(x[[species_col]])
  if (anyDuplicated(species)) {
    stop("Duplicate species names found in column '", species_col, "'.")
  }

  if (is.null(trait_cols)) {
    trait_cols <- setdiff(names(x), species_col)
    if (length(trait_cols) == 0) {
      stop("No trait columns found (excluding '", species_col, "').")
    }
  } else {
    missing_cols <- setdiff(trait_cols, names(x))
    if (length(missing_cols) > 0) {
      stop("Trait columns not found: ", paste(missing_cols, collapse = ", "))
    }
  }

  out <- x[, trait_cols, drop = FALSE]
  rownames(out) <- species

  # Coerce to factor / ordered as requested
  for (nm in intersect(factor_cols %||% character(0), names(out))) {
    out[[nm]] <- factor(out[[nm]])
  }
  for (nm in intersect(ordered_cols %||% character(0), names(out))) {
    out[[nm]] <- factor(out[[nm]], ordered = TRUE)
  }

  out
}
