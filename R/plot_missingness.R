#' Plot observed missingness patterns
#'
#' Display missingness by trait, observation and phylogenetic tip order, or
#' display the fraction of observations missing both traits. Patterns can
#' suggest structure to investigate; they cannot establish whether data are
#' missing completely at random. Trait groups expanded into multiple latent
#' columns count as missing when all their columns are missing.
#'
#' @param x A preprocessed `pigauto_data` object.
#' @param type One of `"summary"`, `"matrix"`, `"phylo"`, or `"co_missing"`.
#' @param tree A `phylo` tree, required for `type = "phylo"`. Rows are
#'   ordered by its tip order; repeated observations appear together.
#' @return Invisibly, a logical observation-by-trait missingness matrix.
#' @examples
#' data(avonet300, tree300)
#' tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30)])
#' d <- avonet300[match(tree$tip.label, avonet300$Species_Key),
#'                c("Mass", "Wing.Length"), drop = FALSE]
#' rownames(d) <- tree$tip.label
#' d$Mass[1:3] <- NA
#' x <- preprocess_traits(d, tree)
#' plot_missingness(x, "summary")
#' plot_missingness(x, "phylo", tree)
#' @export
plot_missingness <- function(x, type = c("summary", "matrix", "phylo", "co_missing"), tree = NULL) {
  if (!inherits(x, "pigauto_data")) stop("'x' must be a pigauto_data object.")
  type <- match.arg(type)
  miss <- vapply(x$trait_map, function(tm) {
    rowSums(!is.na(x$X_scaled[, tm$latent_cols, drop = FALSE])) == 0L
  }, logical(nrow(x$X_scaled)))
  miss <- matrix(miss, nrow = nrow(x$X_scaled),
                 dimnames = list(rownames(x$X_scaled), vapply(x$trait_map, `[[`, character(1), "name")))
  types <- vapply(x$trait_map, `[[`, character(1), "type")
  labels <- paste0(colnames(miss), " (", types, ")")
  # Reserve room for rotated trait labels, and restore the caller's graphics state.
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  label_lines <- max(5, min(12, ceiling(max(graphics::strwidth(labels, units = "inches")) /
                                              graphics::par("cin")[2]) + 2))
  margins <- graphics::par("mar")
  margins[1] <- max(margins[1], label_lines)
  if (type == "co_missing") margins[2] <- max(margins[2], label_lines)
  graphics::par(mar = margins)
  if (type == "summary") {
    graphics::barplot(colMeans(miss), names.arg = labels, las = 2,
                      ylim = c(0, 1), ylab = "Fraction missing", col = "#397a9b")
  } else if (type == "co_missing") {
    joint <- crossprod(miss) / nrow(miss)
    graphics::image(seq_len(ncol(miss)), seq_len(ncol(miss)), joint,
                    zlim = c(0, 1), col = grDevices::colorRampPalette(c("white", "#275d78"))(50),
                    axes = FALSE, xlab = "", ylab = "")
    graphics::axis(1, seq_along(labels), labels, las = 2)
    graphics::axis(2, seq_along(labels), labels, las = 2)
    for (i in seq_len(ncol(miss))) {
      for (j in seq_len(ncol(miss))) {
        graphics::text(i, j, sprintf("%.2f", joint[i, j]),
                       col = if (joint[i, j] > 0.5) "white" else "black")
      }
    }
    graphics::title("Fraction jointly missing")
  } else {
    if (type == "phylo") {
      if (!inherits(tree, "phylo")) stop("'tree' must be a phylo object for phylo ordering.")
      species <- if (is.null(x$obs_to_species)) x$species_names else x$species_names[x$obs_to_species]
      if (anyNA(match(species, tree$tip.label))) stop("All data species must occur in 'tree'.")
      # tip.label order need not be traversal order after tree manipulation
      ordered_tree <- ape::reorder.phylo(tree, "cladewise")
      tips <- ordered_tree$edge[ordered_tree$edge[, 2] <= length(tree$tip.label), 2]
      miss <- miss[order(match(species, tree$tip.label[tips])), , drop = FALSE]
    }
    graphics::image(seq_len(ncol(miss)), seq_len(nrow(miss)), t(miss * 1),
                    zlim = c(0, 1), col = c("#eeeeee", "#397a9b"), axes = FALSE,
                    xlab = "", ylab = if (type == "phylo") "Observations in tree order" else "Observation")
    graphics::axis(1, seq_along(labels), labels, las = 2)
    graphics::legend("topright", legend = c("Observed", "Missing"), fill = c("#eeeeee", "#397a9b"), bty = "n")
  }
  invisible(miss)
}
