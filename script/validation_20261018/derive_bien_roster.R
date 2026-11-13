#!/usr/bin/env Rscript

traits_path <- Sys.getenv("VALIDATION_BIEN_TRAITS", unset = "")
tree_path <- Sys.getenv("VALIDATION_BIEN_TREE", unset = "")
if (!nzchar(traits_path) || !nzchar(tree_path) ||
    !file.exists(traits_path) || !file.exists(tree_path))
  stop("both registered BIEN cache paths must exist")

trait_means <- readRDS(traits_path)
tree_pkg <- readRDS(tree_path)
if (!is.list(trait_means) || is.null(names(trait_means)) ||
    !is.list(tree_pkg) || is.null(tree_pkg[["scenario.3"]]) ||
    !inherits(tree_pkg[["scenario.3"]], "phylo"))
  stop("BIEN cache schema does not match the frozen runner contract")

normalize_species <- function(x) {
  gsub("_", " ", gsub("[[:space:]]+", " ", trimws(as.character(x))))
}
species_lists <- lapply(trait_means, function(d) {
  if (is.data.frame(d) && "species" %in% names(d)) normalize_species(d[["species"]]) else character()
})
all_species <- Reduce(union, species_lists, init = character())
observed <- matrix(FALSE, nrow = length(all_species), ncol = length(trait_means),
                   dimnames = list(NULL, names(trait_means)))
for (trait in names(trait_means)) {
  d <- trait_means[[trait]]
  if (!is.data.frame(d) || !("species" %in% names(d))) next
  value_col <- intersect(c("mean_value", "trait_value", "value"), names(d))[1L]
  if (is.na(value_col)) next
  idx <- match(all_species, normalize_species(d[["species"]]))
  observed[, trait] <- is.finite(suppressWarnings(as.numeric(d[[value_col]][idx])))
}
eligible <- all_species[rowSums(observed) >= 2L]
tree <- tree_pkg[["scenario.3"]]
matched <- intersect(normalize_species(tree[["tip.label"]]), eligible)
if (length(matched) < 300L) stop("fewer than 300 eligible species match the BIEN tree")

out <- data.frame(trait = names(trait_means), matched_species = length(matched),
  subset_n = 300L, subset_seed = 20261018L,
  traits_cache_md5 = unname(tools::md5sum(traits_path)),
  tree_cache_md5 = unname(tools::md5sum(tree_path)), stringsAsFactors = FALSE)
utils::write.csv(out, stdout(), row.names = FALSE, quote = FALSE)
