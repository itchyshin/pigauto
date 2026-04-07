#' Preprocess trait data: align to tree, encode into latent space
#'
#' Aligns species in the trait data frame to the tree, detects or accepts
#' trait types (continuous, binary, categorical, ordinal, count), and
#' encodes each trait into a continuous latent matrix.  The output rows
#' are reordered to match \code{tree$tip.label} order, which is required
#' for downstream graph operations.
#'
#' @details
#' **Type detection** (when \code{trait_types = NULL}) follows the column
#' class of each trait:
#' \describe{
#'   \item{\code{numeric} (non-integer)}{continuous}
#'   \item{\code{integer}}{count}
#'   \item{\code{factor} with 2 levels}{binary}
#'   \item{\code{factor} (unordered) with >2 levels}{categorical}
#'   \item{\code{ordered} factor}{ordinal}
#'   \item{\code{character}}{converted to \code{factor}, then classified}
#'   \item{\code{logical}}{converted to binary (FALSE = 0, TRUE = 1)}
#' }
#'
#' **Latent encoding** per type:
#' \describe{
#'   \item{continuous}{optional \code{log()}, then z-score (1 latent column)}
#'   \item{binary}{0/1 encoding (1 latent column)}
#'   \item{count}{\code{log1p()}, then z-score (1 latent column)}
#'   \item{ordinal}{integer coding (0 to K-1), then z-score (1 latent column)}
#'   \item{categorical}{one-hot encoding (K latent columns)}
#' }
#'
#' @param traits \code{data.frame} with species as row names.  Columns may
#'   be numeric, integer, factor, ordered, character, or logical.
#' @param tree object of class \code{"phylo"}.
#' @param trait_types named character vector overriding auto-detection,
#'   e.g. \code{c(Mass = "continuous", Diet = "categorical")}.
#'   Valid types: \code{"continuous"}, \code{"binary"}, \code{"categorical"},
#'   \code{"ordinal"}, \code{"count"}.  Unspecified traits are auto-detected.
#' @param log_cols character vector of continuous trait names to
#'   log-transform.  Default \code{NULL} means auto-detect (log if all
#'   observed values are positive).  Set to \code{character(0)} to disable.
#' @param log_transform logical.  Legacy parameter: if \code{TRUE} and
#'   \code{log_cols} is \code{NULL}, log-transform all continuous traits
#'   with all-positive values.  Overridden by \code{log_cols} when both
#'   are supplied.
#' @param center logical. Subtract column means for continuous/count/ordinal?
#'   Default \code{TRUE}.
#' @param scale logical. Divide by column SDs for continuous/count/ordinal?
#'   Default \code{TRUE}.
#' @return A list of class \code{"pigauto_data"} with components:
#'   \describe{
#'     \item{X_scaled}{Numeric matrix (n_species x p_latent), latent encoding.}
#'     \item{X_raw}{Numeric matrix of continuous traits after optional log
#'       but before z-scoring (for backward compatibility).}
#'     \item{X_original}{Original data.frame (aligned to tree, before encoding).}
#'     \item{means}{Named numeric vector of column means used for z-scoring
#'       (continuous/count/ordinal traits only).}
#'     \item{sds}{Named numeric vector of column SDs.}
#'     \item{species_names}{Character vector matching \code{tree$tip.label}
#'       order.}
#'     \item{trait_names}{Character vector of original trait names.}
#'     \item{latent_names}{Character vector of latent column names.}
#'     \item{trait_map}{List of trait descriptors (see Details).}
#'     \item{p_latent}{Integer, total number of latent columns.}
#'     \item{log_transform}{Logical, legacy field (TRUE if any continuous
#'       trait was log-transformed).}
#'   }
#' @examples
#' # All-numeric (backward compatible)
#' data(avonet300, tree300, package = "pigauto")
#' traits <- avonet300
#' rownames(traits) <- traits$Species_Key
#' traits$Species_Key <- NULL
#' pd <- preprocess_traits(traits, tree300)
#' dim(pd$X_scaled)   # 300 x 4
#'
#' @importFrom ape keep.tip
#' @export
preprocess_traits <- function(traits, tree, trait_types = NULL,
                              log_cols = NULL, log_transform = TRUE,
                              center = TRUE, scale = TRUE) {
  if (!is.data.frame(traits)) stop("'traits' must be a data.frame.")
  if (!inherits(tree, "phylo")) stop("'tree' must be a phylo object.")
  if (is.null(rownames(traits))) stop("'traits' must have species as row names.")

  # ---- Align species to tree ------------------------------------------------
  in_tree  <- rownames(traits) %in% tree$tip.label
  in_data  <- tree$tip.label %in% rownames(traits)

  n_dropped <- sum(!in_tree)
  if (n_dropped > 0) {
    warning(n_dropped, " species in traits not found in tree -- dropped.")
    traits <- traits[in_tree, , drop = FALSE]
  }
  n_missing_from_data <- sum(!in_data)
  if (n_missing_from_data > 0) {
    message(n_missing_from_data,
            " tree tip(s) have no trait data and will have all-NA rows.")
  }

  idx <- match(tree$tip.label, rownames(traits))
  traits <- traits[idx, , drop = FALSE]
  rownames(traits) <- tree$tip.label

  # ---- Convert convenience types -------------------------------------------
  for (nm in names(traits)) {
    if (is.character(traits[[nm]])) traits[[nm]] <- factor(traits[[nm]])
    if (is.logical(traits[[nm]]))   traits[[nm]] <- factor(traits[[nm]],
                                                           levels = c(FALSE, TRUE))
  }

  X_original <- traits  # save before encoding

  # ---- Detect trait types ---------------------------------------------------
  types <- detect_trait_types(traits, trait_types)

  # ---- Determine which continuous traits to log-transform -------------------
  cont_names <- names(types)[types == "continuous"]
  if (!is.null(log_cols)) {
    # Explicit list overrides everything
    log_set <- intersect(log_cols, cont_names)
  } else if (log_transform) {
    # Legacy: log all continuous traits with all-positive observed values
    log_set <- vapply(cont_names, function(nm) {
      vals <- traits[[nm]][!is.na(traits[[nm]])]
      length(vals) > 0 && all(vals > 0)
    }, logical(1))
    log_set <- cont_names[log_set]
  } else {
    log_set <- character(0)
  }

  # ---- Build trait map and encode to latent ---------------------------------
  trait_map <- build_trait_map(traits, types, log_set, center, scale)
  latent    <- encode_to_latent(traits, trait_map)

  # ---- Backward-compat X_raw (continuous traits only) -----------------------
  cont_cols <- which(types[colnames(traits)] == "continuous")
  if (length(cont_cols) > 0) {
    X_raw <- as.matrix(traits[, cont_cols, drop = FALSE])
    storage.mode(X_raw) <- "double"
    for (nm in intersect(log_set, colnames(X_raw))) {
      X_raw[, nm] <- log(X_raw[, nm])
    }
  } else {
    X_raw <- matrix(nrow = nrow(traits), ncol = 0)
  }

  # ---- Backward-compat means/sds vectors ------------------------------------
  all_means <- vapply(trait_map, function(tm) {
    if (!is.null(tm$mean) && !is.na(tm$mean)) tm$mean else NA_real_
  }, numeric(1))
  names(all_means) <- vapply(trait_map, "[[", character(1), "name")

  all_sds <- vapply(trait_map, function(tm) {
    if (!is.null(tm$sd) && !is.na(tm$sd)) tm$sd else NA_real_
  }, numeric(1))
  names(all_sds) <- names(all_means)

  structure(
    list(
      X_scaled      = latent$X,
      X_raw         = X_raw,
      X_original    = X_original,
      means         = all_means,
      sds           = all_sds,
      species_names = tree$tip.label,
      trait_names   = colnames(traits),
      latent_names  = colnames(latent$X),
      trait_map     = trait_map,
      p_latent      = ncol(latent$X),
      log_transform = length(log_set) > 0
    ),
    class = "pigauto_data"
  )
}


# ---- Internal: detect trait types -------------------------------------------

#' @keywords internal
detect_trait_types <- function(traits, overrides = NULL) {
  types <- character(ncol(traits))
  names(types) <- colnames(traits)

  for (nm in colnames(traits)) {
    x <- traits[[nm]]
    if (!is.null(overrides) && nm %in% names(overrides)) {
      types[nm] <- match.arg(overrides[nm],
                             c("continuous", "binary", "categorical",
                               "ordinal", "count"))
      next
    }
    if (is.ordered(x)) {
      nlev <- nlevels(x)
      types[nm] <- if (nlev == 2) "binary" else "ordinal"
    } else if (is.factor(x)) {
      nlev <- nlevels(x)
      types[nm] <- if (nlev == 2) "binary" else "categorical"
    } else if (is.integer(x)) {
      types[nm] <- "count"
    } else if (is.numeric(x)) {
      types[nm] <- "continuous"
    } else {
      stop("Cannot determine type for column '", nm, "' (class: ",
           paste(class(x), collapse = "/"), ").")
    }
  }
  types
}


# ---- Internal: build trait map ----------------------------------------------

#' @keywords internal
build_trait_map <- function(traits, types, log_set, center, scale) {
  tmap <- list()
  col_offset <- 0L

  for (nm in colnames(traits)) {
    tp <- unname(types[nm])
    x  <- traits[[nm]][!is.na(traits[[nm]])]

    entry <- list(name = nm, type = tp)

    if (tp == "continuous") {
      do_log <- nm %in% log_set
      vals <- if (do_log) log(x) else as.numeric(x)
      m <- if (center && length(vals) > 0) mean(vals) else 0
      s <- if (scale && length(vals) > 1) stats::sd(vals) else 1
      if (s == 0) s <- 1
      entry$n_latent      <- 1L
      entry$latent_cols   <- col_offset + 1L
      entry$levels        <- NULL
      entry$log_transform <- do_log
      entry$mean          <- m
      entry$sd            <- s

    } else if (tp == "count") {
      vals <- log1p(as.numeric(x))
      m <- if (center && length(vals) > 0) mean(vals) else 0
      s <- if (scale && length(vals) > 1) stats::sd(vals) else 1
      if (s == 0) s <- 1
      entry$n_latent      <- 1L
      entry$latent_cols   <- col_offset + 1L
      entry$levels        <- NULL
      entry$log_transform <- FALSE
      entry$mean          <- m
      entry$sd            <- s

    } else if (tp == "ordinal") {
      int_vals <- as.integer(x) - 1L  # 0-based
      vals <- as.numeric(int_vals)
      m <- if (center && length(vals) > 0) mean(vals) else 0
      s <- if (scale && length(vals) > 1) stats::sd(vals) else 1
      if (s == 0) s <- 1
      entry$n_latent      <- 1L
      entry$latent_cols   <- col_offset + 1L
      entry$levels        <- levels(traits[[nm]])
      entry$log_transform <- FALSE
      entry$mean          <- m
      entry$sd            <- s

    } else if (tp == "binary") {
      levs <- levels(traits[[nm]])
      if (is.null(levs)) levs <- sort(unique(x))
      entry$n_latent      <- 1L
      entry$latent_cols   <- col_offset + 1L
      entry$levels        <- levs
      entry$log_transform <- FALSE
      entry$mean          <- NA_real_
      entry$sd            <- NA_real_

    } else if (tp == "categorical") {
      levs <- levels(traits[[nm]])
      if (is.null(levs)) levs <- sort(unique(x))
      K <- length(levs)
      entry$n_latent      <- K
      entry$latent_cols   <- col_offset + seq_len(K)
      entry$levels        <- levs
      entry$log_transform <- FALSE
      entry$mean          <- NA_real_
      entry$sd            <- NA_real_
    }

    col_offset <- col_offset + entry$n_latent
    tmap[[nm]] <- entry
  }
  tmap
}


# ---- Internal: encode to latent matrix --------------------------------------

#' @keywords internal
encode_to_latent <- function(traits, trait_map) {
  n       <- nrow(traits)
  p_latent <- sum(vapply(trait_map, "[[", integer(1), "n_latent"))
  X       <- matrix(NA_real_, nrow = n, ncol = p_latent)

  lat_names <- character(p_latent)

  for (tm in trait_map) {
    nm  <- tm$name
    col <- traits[[nm]]

    if (tm$type == "continuous") {
      vals <- as.numeric(col)
      if (tm$log_transform) vals <- log(vals)
      vals <- (vals - tm$mean) / tm$sd
      X[, tm$latent_cols] <- vals
      lat_names[tm$latent_cols] <- nm

    } else if (tm$type == "count") {
      vals <- log1p(as.numeric(col))
      vals <- (vals - tm$mean) / tm$sd
      X[, tm$latent_cols] <- vals
      lat_names[tm$latent_cols] <- nm

    } else if (tm$type == "ordinal") {
      vals <- as.numeric(as.integer(col) - 1L)
      vals <- (vals - tm$mean) / tm$sd
      X[, tm$latent_cols] <- vals
      lat_names[tm$latent_cols] <- nm

    } else if (tm$type == "binary") {
      # Second level = 1, first level = 0
      vals <- as.numeric(col == tm$levels[2])
      vals[is.na(col)] <- NA_real_
      X[, tm$latent_cols] <- vals
      lat_names[tm$latent_cols] <- nm

    } else if (tm$type == "categorical") {
      K    <- tm$n_latent
      levs <- tm$levels
      for (k in seq_len(K)) {
        vals <- as.numeric(col == levs[k])
        vals[is.na(col)] <- NA_real_
        X[, tm$latent_cols[k]] <- vals
        lat_names[tm$latent_cols[k]] <- paste0(nm, "=", levs[k])
      }
    }
  }

  rownames(X) <- rownames(traits)
  colnames(X) <- lat_names
  list(X = X)
}


#' @export
print.pigauto_data <- function(x, ...) {
  cat("pigauto_data\n")
  cat("  Species:", length(x$species_names), "\n")
  cat("  Traits: ", length(x$trait_names), "\n")

  if (!is.null(x$trait_map)) {
    types <- vapply(x$trait_map, "[[", character(1), "type")
    type_tab <- table(types)
    cat("  Types:  ",
        paste(names(type_tab), type_tab, sep = "=", collapse = ", "), "\n")
    cat("  Latent columns:", x$p_latent, "\n")
  } else {
    cat("  Columns:", paste(x$trait_names, collapse = ", "), "\n")
    cat("  Log-transformed:", x$log_transform, "\n")
  }

  missing_frac <- mean(is.na(x$X_scaled))
  if (missing_frac > 0) {
    cat("  Missing values:", round(100 * missing_frac, 1), "%\n")
  }
  invisible(x)
}
