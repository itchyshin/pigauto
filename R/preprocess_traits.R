#' Preprocess trait data: align to tree, encode into latent space
#'
#' Aligns species in the trait data frame to the tree, detects or accepts
#' trait types (continuous, binary, categorical, ordinal, count, proportion,
#' zi_count), and encodes each trait into a continuous latent matrix.
#'
#' When each species has one observation (the default), output rows match
#' \code{tree$tip.label} order.  When \code{species_col} is supplied,
#' multiple observations per species are supported: the output matrix has
#' one row per observation, plus an \code{obs_to_species} mapping for the
#' GNN (which operates at species level).
#'
#' @details
#' **Automatic type detection** (when \code{trait_types = NULL}) follows the
#' R class of each column — no user input is required for most data:
#' \tabular{ll}{
#'   \strong{R class} \tab \strong{pigauto type} \cr
#'   \code{numeric} (non-integer) \tab continuous \cr
#'   \code{integer} \tab count \cr
#'   \code{factor} with 2 levels \tab binary \cr
#'   \code{factor} (unordered) with >2 levels \tab categorical \cr
#'   \code{ordered} / \code{factor(..., ordered = TRUE)} \tab ordinal \cr
#'   \code{character} \tab converted to \code{factor}, then binary or categorical \cr
#'   \code{logical} \tab binary (FALSE = 0, TRUE = 1) \cr
#' }
#'
#' **Two types require an explicit override** because they cannot be
#' distinguished from their R class alone:
#' \describe{
#'   \item{\code{"proportion"}}{A \code{numeric} column bounded 0–1 looks
#'     identical to continuous.  Declare it explicitly:
#'     \code{trait_types = c(SurvivalRate = "proportion")}.
#'     Encoded via \code{qlogis(clamp(x, 0.001, 0.999))}.}
#'   \item{\code{"zi_count"}}{An \code{integer} column with excess zeros
#'     looks identical to count.  Declare it explicitly:
#'     \code{trait_types = c(Parasites = "zi_count")}.
#'     Encoded as a binary zero/non-zero gate plus \code{log1p}-z magnitude.}
#' }
#'
#' **Practical examples of type assignment:**
#' \itemize{
#'   \item Body mass (\code{numeric}, all positive) → \code{continuous},
#'     auto-log-transformed.
#'   \item Clutch size (\code{integer}) → \code{count}.
#'   \item Migratory (\code{factor} with levels "Yes"/"No") → \code{binary}.
#'   \item Diet (\code{factor} with >2 levels) → \code{categorical}.
#'   \item IUCN status (\code{ordered} factor, LC < NT < VU < EN < CR) →
#'     \code{ordinal}.  If left as an unordered factor, it becomes
#'     \code{categorical} — both are valid depending on the question.
#'   \item Parasite load (\code{integer} with many zeros) → needs
#'     \code{trait_types = c(Parasites = "zi_count")}.
#'   \item Survival rate (\code{numeric}, values in 0 to 1) → needs
#'     \code{trait_types = c(Survival = "proportion")}.
#' }
#'
#' **Latent encoding** per type:
#' \describe{
#'   \item{continuous}{optional \code{log()}, then z-score (1 latent column)}
#'   \item{binary}{0/1 encoding (1 latent column)}
#'   \item{count}{\code{log1p()}, then z-score (1 latent column)}
#'   \item{ordinal}{integer coding (0 to K-1), then z-score (1 latent column)}
#'   \item{categorical}{one-hot encoding (K latent columns)}
#'   \item{proportion}{\code{qlogis(clamp(x, 0.001, 0.999))}, then z-score (1 latent column)}
#'   \item{zi_count}{gate (0/1) + \code{log1p}-z of non-zeros (2 latent columns)}
#'   \item{multi_proportion}{centred log-ratio (CLR) + per-component z-score
#'     (K latent columns per group).  Rows must sum to 1.  Declared via the
#'     \code{multi_proportion_groups} argument, not \code{trait_types}.}
#' }
#'
#' @param traits \code{data.frame} with species as row names (one row per
#'   species), or with a species column identified by \code{species_col}
#'   (potentially multiple rows per species).  Columns may be numeric,
#'   integer, factor, ordered, character, or logical.
#' @param tree object of class \code{"phylo"}.
#' @param species_col character.  Name of the column in \code{traits} that
#'   identifies species.  When supplied, \code{traits} may have multiple
#'   rows per species.  The column is removed from trait columns before
#'   encoding.  Default \code{NULL} uses row names (one row per species).
#' @param trait_types named character vector overriding auto-detection,
#'   e.g. \code{c(Mass = "continuous", Diet = "categorical")}.
#'   Valid types: \code{"continuous"}, \code{"binary"}, \code{"categorical"},
#'   \code{"ordinal"}, \code{"count"}, \code{"proportion"}, \code{"zi_count"}.
#'   Proportion and zi_count are override-only (not auto-detected).
#'   Unspecified traits are auto-detected.  Note that
#'   \code{"multi_proportion"} is NOT set here — use
#'   \code{multi_proportion_groups} instead.
#' @param multi_proportion_groups named list declaring compositional
#'   (multi-proportion) trait groups.  Each element is a character vector
#'   of column names whose row-wise values sum to 1 (e.g.
#'   \code{list(diet = c("plants", "insects", "fish"))}).  The group
#'   name becomes a single trait in the output, encoded via centred
#'   log-ratio (CLR) + per-component z-score.  Group names must NOT
#'   match any column in \code{traits}.  Default \code{NULL}.
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
#' @param covariates data.frame or numeric matrix of environmental covariates.
#'   Covariates are conditioners: they inform imputation but are not themselves
#'   imputed, so they must be \strong{fully observed} (no NAs — if a variable
#'   has missing values, put it in \code{traits} instead).  Must have the same
#'   number of rows as \code{traits} after alignment to the tree.
#'   \describe{
#'     \item{Numeric / integer columns}{z-scored automatically.}
#'     \item{Factor / ordered columns}{one-hot encoded (K binary columns per
#'       factor with K levels).  Column names become \code{"var.level"}.}
#'     \item{Character / logical columns}{coerced to factor, then one-hot.}
#'   }
#'   Default \code{NULL} (no covariates).
#' @return A list of class \code{"pigauto_data"} with components:
#'   \describe{
#'     \item{X_scaled}{Numeric matrix (n_obs x p_latent), latent encoding.
#'       When \code{species_col} is \code{NULL}, n_obs = n_species.}
#'     \item{X_raw}{Numeric matrix of continuous traits after optional log
#'       but before z-scoring (for backward compatibility).}
#'     \item{X_original}{Original data.frame (aligned to tree, before encoding).}
#'     \item{means}{Named numeric vector of column means used for z-scoring
#'       (continuous/count/ordinal traits only).}
#'     \item{sds}{Named numeric vector of column SDs.}
#'     \item{species_names}{Character vector of unique species matching
#'       \code{tree$tip.label} order (length = n_species).}
#'     \item{obs_species}{Character vector of species labels per observation
#'       (length = n_obs).  When multi-obs, can have duplicates.
#'       \code{NULL} when \code{species_col} is \code{NULL}.}
#'     \item{obs_to_species}{Integer vector (length = n_obs) mapping each
#'       observation to its species index in \code{species_names}.
#'       \code{NULL} when \code{species_col} is \code{NULL}.}
#'     \item{n_species}{Integer, number of unique species.}
#'     \item{n_obs}{Integer, number of observations (= n_species when
#'       single-obs).}
#'     \item{multi_obs}{Logical, \code{TRUE} when multiple observations
#'       per species are present.}
#'     \item{trait_names}{Character vector of original trait names.}
#'     \item{latent_names}{Character vector of latent column names.}
#'     \item{trait_map}{List of trait descriptors (see Details).}
#'     \item{p_latent}{Integer, total number of latent columns.}
#'     \item{log_transform}{Logical, legacy field (TRUE if any continuous
#'       trait was log-transformed).}
#'   }
#' @examples
#' # Single-obs per species (backward compatible)
#' data(avonet300, tree300, package = "pigauto")
#' traits <- avonet300
#' rownames(traits) <- traits$Species_Key
#' traits$Species_Key <- NULL
#' pd <- preprocess_traits(traits, tree300)
#' dim(pd$X_scaled)   # 300 x p_latent
#'
#' # Multi-obs per species (via species_col)
#' pd2 <- preprocess_traits(avonet300, tree300, species_col = "Species_Key")
#' pd2$n_obs      # number of observations
#' pd2$n_species  # number of unique species
#'
#' @importFrom ape keep.tip
#' @export
preprocess_traits <- function(traits, tree, species_col = NULL,
                              trait_types = NULL,
                              multi_proportion_groups = NULL,
                              log_cols = NULL, log_transform = TRUE,
                              center = TRUE, scale = TRUE,
                              covariates = NULL) {
  if (!is.data.frame(traits)) stop("'traits' must be a data.frame.")
  if (!inherits(tree, "phylo")) stop("'tree' must be a phylo object.")

  # ---- Validate multi_proportion_groups -------------------------------------
  if (!is.null(multi_proportion_groups)) {
    if (!is.list(multi_proportion_groups) ||
        is.null(names(multi_proportion_groups)) ||
        any(names(multi_proportion_groups) == "")) {
      stop("'multi_proportion_groups' must be a named list, e.g. ",
           "list(colour = c('black','blue','red')).", call. = FALSE)
    }
    all_group_cols <- unlist(multi_proportion_groups, use.names = FALSE)
    if (anyDuplicated(all_group_cols)) {
      stop("Columns appear in more than one multi_proportion group: ",
           paste(all_group_cols[duplicated(all_group_cols)], collapse = ", "),
           call. = FALSE)
    }
    if (any(names(multi_proportion_groups) %in% colnames(traits))) {
      bad <- intersect(names(multi_proportion_groups), colnames(traits))
      stop("multi_proportion group name(s) collide with existing column(s): ",
           paste(bad, collapse = ", "),
           ". Pick group names that are NOT column names in `traits`.",
           call. = FALSE)
    }
    missing_cols <- setdiff(all_group_cols, colnames(traits))
    if (length(missing_cols) > 0) {
      stop("multi_proportion_groups references columns not in `traits`: ",
           paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    for (gnm in names(multi_proportion_groups)) {
      gcols <- multi_proportion_groups[[gnm]]
      if (length(gcols) < 2L) {
        stop("multi_proportion group '", gnm, "' needs >= 2 columns.",
             call. = FALSE)
      }
      sub <- traits[, gcols, drop = FALSE]
      if (!all(vapply(sub, is.numeric, logical(1)))) {
        stop("multi_proportion group '", gnm, "' has non-numeric columns.",
             call. = FALSE)
      }
    }
  }

  # ---- Determine species identity per row -----------------------------------
  multi_obs <- !is.null(species_col)
  if (multi_obs) {
    if (!species_col %in% names(traits)) {
      stop("species_col '", species_col, "' not found in traits.")
    }
    obs_species_raw <- as.character(traits[[species_col]])
    traits[[species_col]] <- NULL  # remove from trait columns
  } else {
    if (is.null(rownames(traits))) {
      stop("'traits' must have species as row names (or supply species_col).")
    }
    obs_species_raw <- rownames(traits)
  }

  # ---- Align species to tree ------------------------------------------------
  unique_species <- unique(obs_species_raw)
  in_tree  <- unique_species %in% tree$tip.label
  in_data  <- tree$tip.label %in% unique_species

  n_dropped_sp <- sum(!in_tree)
  if (n_dropped_sp > 0) {
    warning(n_dropped_sp, " species in traits not found in tree -- dropped.")
    keep_rows <- obs_species_raw %in% tree$tip.label
    traits <- traits[keep_rows, , drop = FALSE]
    obs_species_raw <- obs_species_raw[keep_rows]
    unique_species <- unique(obs_species_raw)
  }
  n_missing_from_data <- sum(!in_data)
  if (n_missing_from_data > 0) {
    message(n_missing_from_data,
            " tree tip(s) have no trait data and will have all-NA rows.")
  }

  if (multi_obs) {
    # Order by tree tip order; within species, preserve original row order
    sp_order <- tree$tip.label[tree$tip.label %in% unique_species]
    sp_rank  <- setNames(seq_along(sp_order), sp_order)
    row_order <- order(sp_rank[obs_species_raw])
    traits <- traits[row_order, , drop = FALSE]
    obs_species_raw <- obs_species_raw[row_order]
    rownames(traits) <- NULL

    # `input_row_order[k]` = original input-row index that was placed at
    # internal position k.  Used by `build_completed()` (and any downstream
    # caller that needs to re-align internal-order results back to the
    # user's input row order).
    input_row_order <- as.integer(row_order)

    # Add all-NA rows for tree tips missing from data
    missing_tips <- tree$tip.label[!in_data]
    if (length(missing_tips) > 0) {
      na_rows <- traits[rep(NA, length(missing_tips)), , drop = FALSE]
      rownames(na_rows) <- NULL
      traits <- rbind(traits, na_rows)
      obs_species_raw <- c(obs_species_raw, missing_tips)
      # Pad input_row_order with NAs for these synthetic rows (no original
      # input row corresponds to them).
      input_row_order <- c(input_row_order,
                            rep(NA_integer_, length(missing_tips)))
    }

    # Build species-level outputs
    species_names <- tree$tip.label[tree$tip.label %in%
                                      c(unique_species, missing_tips)]
    obs_to_species <- match(obs_species_raw, species_names)
    n_obs <- nrow(traits)
    n_species <- length(species_names)
    obs_species <- obs_species_raw
  } else {
    # Single-obs path: reorder rows to match tree$tip.label exactly
    idx <- match(tree$tip.label, rownames(traits))
    # `input_row_order[k]` = original input-row index that was placed at
    # internal position k.  `idx[k]` is the row of `traits` (= input row
    # index, since rownames are preserved at this point) that maps to
    # tree$tip.label[k].  When a tip has no input row, idx[k] = NA.
    input_row_order <- as.integer(idx)
    traits <- traits[idx, , drop = FALSE]
    rownames(traits) <- tree$tip.label

    species_names <- tree$tip.label
    obs_species <- NULL
    obs_to_species <- NULL
    n_obs <- nrow(traits)
    n_species <- n_obs
  }

  # ---- Convert convenience types -------------------------------------------
  for (nm in names(traits)) {
    if (is.character(traits[[nm]])) traits[[nm]] <- factor(traits[[nm]])
    if (is.logical(traits[[nm]]))   traits[[nm]] <- factor(traits[[nm]],
                                                           levels = c(FALSE, TRUE))
  }

  X_original <- traits  # save before encoding

  # ---- Detect trait types ---------------------------------------------------
  # Columns that belong to a multi_proportion group are excluded from
  # the per-column type detection — they're handled as a single group entry
  # in build_trait_map below.
  group_member_cols <- unlist(multi_proportion_groups, use.names = FALSE)
  standalone_cols   <- setdiff(colnames(traits), group_member_cols)
  types <- detect_trait_types(traits[, standalone_cols, drop = FALSE],
                              trait_types)

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
  # Pass multi_proportion_groups so build_trait_map can add one entry per group
  trait_map <- build_trait_map(traits, types, log_set, center, scale,
                               multi_proportion_groups = multi_proportion_groups)
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
  # Note: multi_proportion traits have K-vector mean/sd (one per component),
  # not a scalar — we store NA_real_ in these scalar-per-trait back-compat
  # vectors since no downstream code legitimately consumes a scalar mean/sd
  # for multi_proportion. The full mean/sd lives in trait_map[[gnm]]$mean.
  all_means <- vapply(trait_map, function(tm) {
    if (!is.null(tm$mean) && length(tm$mean) == 1L && !is.na(tm$mean)) {
      tm$mean
    } else {
      NA_real_
    }
  }, numeric(1))
  names(all_means) <- vapply(trait_map, "[[", character(1), "name")

  all_sds <- vapply(trait_map, function(tm) {
    if (!is.null(tm$sd) && length(tm$sd) == 1L && !is.na(tm$sd)) {
      tm$sd
    } else {
      NA_real_
    }
  }, numeric(1))
  names(all_sds) <- names(all_means)

  # ---- Environmental covariates (conditioners, not imputed) ------------------
  # Supports: numeric columns (z-scored), factor/ordered columns (one-hot).
  # No NAs allowed — covariates are conditioning variables, not imputation
  # targets.  If a variable has missing values, put it in `traits` instead.
  cov_scaled <- NULL
  cov_means  <- NULL
  cov_sds    <- NULL
  cov_names  <- NULL

  if (!is.null(covariates)) {
    if (!is.data.frame(covariates) && !is.matrix(covariates)) {
      stop("`covariates` must be a data.frame or matrix.", call. = FALSE)
    }
    # Coerce matrix to data.frame so column classes are preserved
    if (is.matrix(covariates)) {
      covariates <- as.data.frame(covariates)
    }

    if (nrow(covariates) != n_obs) {
      stop("`covariates` has ", nrow(covariates), " rows but traits has ",
           n_obs, " rows after alignment. They must match.", call. = FALSE)
    }
    if (any(vapply(covariates, anyNA, logical(1)))) {
      stop("`covariates` must be fully observed (no NAs). ",
           "If a covariate has missing values, include it in `traits` ",
           "instead so pigauto can impute it jointly.",
           call. = FALSE)
    }

    # Encode each covariate column -----------------------------------------
    # numeric/integer  → z-score (1 column)
    # factor/ordered   → one-hot (K columns, column names = "var.levelK")
    # character/logical → coerced to factor then one-hot
    cov_parts      <- list()
    cov_means_list <- list()   # only for numeric cols
    cov_sds_list   <- list()

    raw_names <- colnames(covariates)
    if (is.null(raw_names))
      raw_names <- paste0("cov", seq_len(ncol(covariates)))

    for (j in seq_along(raw_names)) {
      nm  <- raw_names[j]
      col <- covariates[[j]]

      if (is.character(col) || is.logical(col)) col <- factor(col)

      if (is.factor(col) || is.ordered(col)) {
        # One-hot encoding (drop-none; model learns to handle collinearity)
        levs   <- levels(col)
        K      <- length(levs)
        oh_mat <- matrix(0.0, nrow = n_obs, ncol = K)
        colnames(oh_mat) <- paste0(nm, ".", levs)
        for (k in seq_len(K)) oh_mat[, k] <- as.numeric(col == levs[k])
        cov_parts[[nm]] <- oh_mat
        # No z-score means/sds for one-hot columns
        cov_means_list[[nm]] <- rep(NA_real_, K)
        cov_sds_list[[nm]]   <- rep(NA_real_, K)

      } else if (is.numeric(col) || is.integer(col)) {
        col   <- as.double(col)
        cmean <- mean(col)
        csd   <- stats::sd(col)
        if (!is.finite(csd) || csd < 1e-12) csd <- 1.0  # constant column
        scaled_col <- (col - cmean) / csd
        mat <- matrix(scaled_col, ncol = 1, dimnames = list(NULL, nm))
        cov_parts[[nm]]      <- mat
        cov_means_list[[nm]] <- cmean
        cov_sds_list[[nm]]   <- csd

      } else {
        stop("Column '", nm, "' in `covariates` has unsupported class '",
             paste(class(col), collapse = "/"), "'. ",
             "Use numeric, integer, factor, ordered, character, or logical.",
             call. = FALSE)
      }
    }

    cov_scaled <- do.call(cbind, cov_parts)
    cov_names  <- colnames(cov_scaled)
    cov_means  <- unlist(cov_means_list)
    cov_sds    <- unlist(cov_sds_list)
    names(cov_means) <- cov_names
    names(cov_sds)   <- cov_names
  }

  structure(
    list(
      X_scaled       = latent$X,
      X_raw          = X_raw,
      X_original     = X_original,
      means          = all_means,
      sds            = all_sds,
      species_names  = species_names,
      obs_species    = obs_species,
      obs_to_species = obs_to_species,
      input_row_order = input_row_order,
      n_species      = n_species,
      n_obs          = n_obs,
      multi_obs      = multi_obs,
      trait_names    = colnames(traits),
      latent_names   = colnames(latent$X),
      trait_map      = trait_map,
      p_latent       = ncol(latent$X),
      log_transform  = length(log_set) > 0,
      covariates     = cov_scaled,
      cov_means      = cov_means,
      cov_sds        = cov_sds,
      cov_names      = cov_names
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
                               "ordinal", "count", "proportion",
                               "zi_count"))
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

#' Internal: build trait map from type labels
#' @keywords internal
#' @noRd
build_trait_map <- function(traits, types, log_set, center, scale,
                            multi_proportion_groups = NULL) {
  tmap <- list()
  col_offset <- 0L
  group_member_cols <- unlist(multi_proportion_groups, use.names = FALSE)

  # Iterate over named columns of `types` (which excludes group members)
  for (nm in names(types)) {
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
      # Phase G (2026-05-01): observed range in original units, used by
      # the optional `clamp_outliers` argument to predict.pigauto_fit() to
      # cap tail-extrapolation amplified by the back-transform.  Only
      # populated for log-transformed traits where expm1/exp() can blow up
      # latent overshoots; left NA for untransformed continuous (which
      # decode linearly and don't need clamping).
      if (do_log && length(stats::na.omit(as.numeric(x))) > 0L) {
        entry$obs_max <- max(as.numeric(x), na.rm = TRUE)
        entry$obs_min <- min(as.numeric(x), na.rm = TRUE)
      } else {
        entry$obs_max <- NA_real_
        entry$obs_min <- NA_real_
      }

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
      # Phase G: count is always log1p-z; expm1 back-transform amplifies
      # tail overshoots.  Clamp at obs_max * clamp_factor when the user
      # opts in via predict(clamp_outliers = TRUE).
      if (length(stats::na.omit(as.numeric(x))) > 0L) {
        entry$obs_max <- max(as.numeric(x), na.rm = TRUE)
        entry$obs_min <- min(as.numeric(x), na.rm = TRUE)
      } else {
        entry$obs_max <- NA_real_
        entry$obs_min <- NA_real_
      }

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

    } else if (tp == "proportion") {
      # Logit transform: qlogis(clamp(x, 0.001, 0.999)), then z-score
      vals <- stats::qlogis(pmin(pmax(as.numeric(x), 0.001), 0.999))
      m <- if (center && length(vals) > 0) mean(vals) else 0
      s <- if (scale && length(vals) > 1) stats::sd(vals) else 1
      if (s == 0) s <- 1
      entry$n_latent      <- 1L
      entry$latent_cols   <- col_offset + 1L
      entry$levels        <- NULL
      entry$log_transform <- FALSE
      entry$mean          <- m
      entry$sd            <- s

    } else if (tp == "zi_count") {
      # Two latent columns: gate (binary 0/1) + magnitude (log1p-z of non-zeros)
      nz <- x[x > 0]
      nz_log <- if (length(nz) > 0) log1p(as.numeric(nz)) else numeric(0)
      m_nz <- if (center && length(nz_log) > 0) mean(nz_log) else 0
      s_nz <- if (scale && length(nz_log) > 1) stats::sd(nz_log) else 1
      if (s_nz == 0) s_nz <- 1
      entry$n_latent      <- 2L
      entry$latent_cols   <- col_offset + 1:2
      entry$levels        <- NULL
      entry$log_transform <- FALSE
      entry$mean          <- m_nz      # for magnitude column z-scoring
      entry$sd            <- s_nz
      entry$zero_frac     <- sum(x == 0) / length(x)
      # Phase G: same expm1 amplification hazard as count.  Use the
      # non-zero observations to derive the clamp range (zero values
      # would push obs_min to 0 and make the lower clamp meaningless).
      if (length(nz) > 0L) {
        entry$obs_max <- max(as.numeric(nz), na.rm = TRUE)
        entry$obs_min <- min(as.numeric(nz), na.rm = TRUE)
      } else {
        entry$obs_max <- NA_real_
        entry$obs_min <- NA_real_
      }
    }

    col_offset <- col_offset + entry$n_latent
    tmap[[nm]] <- entry
  }

  # ---- Multi-proportion group entries (CLR + per-component z-scoring) ----
  # Each group becomes ONE trait_map entry with K latent columns.
  if (!is.null(multi_proportion_groups)) {
    eps <- 1e-6   # for log-safety on observed zeros
    for (gnm in names(multi_proportion_groups)) {
      gcols <- multi_proportion_groups[[gnm]]
      K <- length(gcols)
      mat_raw <- as.matrix(traits[, gcols, drop = FALSE])
      storage.mode(mat_raw) <- "double"

      # Only use rows with NO NAs across the group for fitting mean/sd
      complete_rows <- stats::complete.cases(mat_raw)
      mat_fit <- mat_raw[complete_rows, , drop = FALSE]

      # Handle zeros: replace with eps, re-normalise so rows still sum to ~1
      if (nrow(mat_fit) > 0L) {
        mat_fit <- pmax(mat_fit, eps)
        rs <- rowSums(mat_fit)
        mat_fit <- mat_fit / rs
        # CLR: log(x) - mean(log(x)) row-wise
        log_mat <- log(mat_fit)
        gmean_row <- rowMeans(log_mat)
        clr_mat <- log_mat - gmean_row
      } else {
        clr_mat <- matrix(NA_real_, 0L, K)
      }

      mu <- if (center && nrow(clr_mat) > 0L) colMeans(clr_mat) else rep(0, K)
      sd_vec <- if (scale && nrow(clr_mat) > 1L) {
        apply(clr_mat, 2L, stats::sd)
      } else {
        rep(1, K)
      }
      sd_vec[sd_vec == 0] <- 1

      entry <- list(
        name          = gnm,
        type          = "multi_proportion",
        n_latent      = K,
        latent_cols   = col_offset + seq_len(K),
        levels        = gcols,      # component names
        input_cols    = gcols,      # original column names in traits
        log_transform = FALSE,
        mean          = mu,         # K-vector
        sd            = sd_vec,     # K-vector
        epsilon       = eps
      )
      col_offset <- col_offset + K
      tmap[[gnm]] <- entry
    }
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

    } else if (tm$type == "proportion") {
      vals <- as.numeric(col)
      vals <- stats::qlogis(pmin(pmax(vals, 0.001), 0.999))
      vals <- (vals - tm$mean) / tm$sd
      X[, tm$latent_cols] <- vals
      lat_names[tm$latent_cols] <- nm

    } else if (tm$type == "multi_proportion") {
      # K CLR columns, per-component z-score, small epsilon for zeros
      mat_raw <- as.matrix(traits[, tm$input_cols, drop = FALSE])
      storage.mode(mat_raw) <- "double"
      K  <- tm$n_latent
      lc <- tm$latent_cols
      eps <- tm$epsilon

      complete_rows <- stats::complete.cases(mat_raw)

      # Replace zeros/small values with epsilon and re-normalise
      mat <- mat_raw
      mat[complete_rows, ] <- pmax(mat_raw[complete_rows, , drop = FALSE], eps)
      rs <- rowSums(mat[complete_rows, , drop = FALSE])
      mat[complete_rows, ] <- mat[complete_rows, , drop = FALSE] / rs

      # CLR
      clr <- matrix(NA_real_, nrow(mat_raw), K)
      if (any(complete_rows)) {
        log_mat <- log(mat[complete_rows, , drop = FALSE])
        clr[complete_rows, ] <- log_mat - rowMeans(log_mat)
      }
      # Per-component z-score
      for (k in seq_len(K)) {
        clr[, k] <- (clr[, k] - tm$mean[k]) / tm$sd[k]
      }
      X[, lc] <- clr
      lat_names[lc] <- paste0(tm$name, "=", tm$levels)
      next

    } else if (tm$type == "zi_count") {
      lc <- tm$latent_cols
      raw <- as.numeric(col)
      # Column 1 (gate): 0 = zero, 1 = non-zero, NA = missing
      gate <- rep(NA_real_, length(raw))
      gate[!is.na(raw) & raw == 0] <- 0
      gate[!is.na(raw) & raw > 0]  <- 1
      X[, lc[1]] <- gate
      lat_names[lc[1]] <- paste0(nm, "_gate")
      # Column 2 (magnitude): log1p-z of non-zero values, NA otherwise
      mag <- rep(NA_real_, length(raw))
      nz_idx <- !is.na(raw) & raw > 0
      if (any(nz_idx)) {
        mag[nz_idx] <- (log1p(raw[nz_idx]) - tm$mean) / tm$sd
      }
      X[, lc[2]] <- mag
      lat_names[lc[2]] <- paste0(nm, "_mag")
    }
  }

  rownames(X) <- rownames(traits)
  colnames(X) <- lat_names
  list(X = X)
}


#' @export
print.pigauto_data <- function(x, ...) {
  cat("pigauto_data\n")
  cat("  Species:", x$n_species %||% length(x$species_names), "\n")
  if (isTRUE(x$multi_obs)) {
    cat("  Observations:", x$n_obs, "\n")
    obs_per_sp <- table(x$obs_to_species)
    cat("  Obs/species: min=", min(obs_per_sp), " median=",
        stats::median(obs_per_sp), " max=", max(obs_per_sp), "\n", sep = "")
  }
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
