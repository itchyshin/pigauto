#' Read trait data from a CSV file or data frame
#'
#' Loads trait data and sets species names as row names. If a CSV path is
#' supplied, it is read with \code{read.csv}. All species identifiers must be
#' present and unique because row names identify one trait row per species. For
#' repeated observations, use \code{read.csv()} and pass the resulting data
#' frame with \code{species_col} to \code{check_pigauto()} and \code{impute()}.
#' By default, all non-species columns are returned (numeric, factor, integer,
#' etc.). Use \code{trait_cols} to select a subset.
#'
#' @param x character path to a CSV file, or a \code{data.frame}.
#' @param species_col character. Name of the column containing species names
#'   (default \code{"species"}).
#' @param trait_cols character vector of column names to include. If
#'   \code{NULL} (default), all columns except \code{species_col} are used.
#' @param factor_cols character vector of columns to coerce to \code{factor}.
#' @param ordered_cols character vector of columns to coerce to \code{ordered}.
#'   Without \code{ordinal_levels}, levels follow R's default sorting order.
#' @param ordinal_levels named list specifying the scientific order of ordinal
#'   levels, for example \code{list(status = c("LC", "NT", "VU", "EN", "CR"))}.
#'   Each named column is converted to an ordered factor. Unlisted observed
#'   values are an error, so a misspelling cannot silently become missing.
#' @return A \code{data.frame} with unique species names as row names.
#' @examples
#' df <- data.frame(species = c("Sp_a", "Sp_b"), mass = c(10, 20))
#' traits <- read_traits(df, species_col = "species")
#' @importFrom utils read.csv
#' @export
read_traits <- function(x, species_col = "species", trait_cols = NULL,
                        factor_cols = NULL, ordered_cols = NULL,
                        ordinal_levels = NULL) {
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

  out <- .apply_ordinal_levels(out, ordinal_levels)

  out
}


#' Internal: coerce declared ordinal columns to ordered factors
#'
#' Shared by read_traits(), preprocess_traits(), impute(), multi_impute(),
#' multi_impute_trees() and check_pigauto() so a declared level order means
#' the same thing on every entry point. Observed values absent from the
#' declared levels are an error, so a misspelling cannot silently become
#' missing.
#'
#' @keywords internal
#' @noRd
.apply_ordinal_levels <- function(traits, ordinal_levels) {
  if (is.null(ordinal_levels)) return(traits)
  if (!is.list(ordinal_levels) || is.null(names(ordinal_levels)) ||
      anyNA(names(ordinal_levels)) || any(!nzchar(names(ordinal_levels))) ||
      anyDuplicated(names(ordinal_levels))) {
    stop("'ordinal_levels' must be a uniquely named list.", call. = FALSE)
  }
  absent <- setdiff(names(ordinal_levels), names(traits))
  if (length(absent)) stop("Ordinal columns not found: ",
                           paste(absent, collapse = ", "), call. = FALSE)
  for (nm in names(ordinal_levels)) {
    lev <- ordinal_levels[[nm]]
    if (!is.character(lev) || length(lev) < 2L || anyNA(lev) ||
        any(!nzchar(lev)) || anyDuplicated(lev)) {
      stop("Levels for '", nm, "' must be unique, nonmissing strings (at least two).",
           call. = FALSE)
    }
    unknown <- setdiff(as.character(traits[[nm]][!is.na(traits[[nm]])]), lev)
    if (length(unknown)) stop("Values in '", nm, "' not included in ordinal_levels: ",
                              paste(unknown, collapse = ", "), call. = FALSE)
    traits[[nm]] <- ordered(as.character(traits[[nm]]), levels = lev)
  }
  traits
}
