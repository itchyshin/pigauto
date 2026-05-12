#' Safe Impute Wrapper
#' Automatically detects K >= 3 ordinal traits and adjusts defaults 
#' to prevent class collapse (known Phase B3).
#'
#' @param traits A dataframe containing the species traits.
#' @param tree A phylogenetic tree of class phylo.
#' @param n_imputations Number of imputations (defaults to 1L).
#' @param pool_method Pooling method (defaults to "median").
#' @param ... Additional arguments passed to pigauto::impute()
#' @export
safe_impute <- function(traits, tree, n_imputations = 1L, pool_method = "median", ...) {
  
  # Scan the dataset for any column that is a factor with 3 or more levels
  has_k3_ordinal <- any(sapply(traits, function(col) is.factor(col) && length(levels(col)) >= 3))
  
  # If a K>=3 trait is found AND the user is using the default settings:
  if (has_k3_ordinal && n_imputations == 1L && pool_method == "median") {
    
    # print a message so the user knows what is happening
    message("🛡️  [pigauto: safe_impute] K >= 3 ordinal trait detected.")
    message("    -> Automatically overriding default settings.")
    message("    -> Setting n_imputations = 20L and pool_method = 'mode'.")
    
    # Override the variables
    n_imputations <- 20L
    pool_method <- "mode"
  }
  
  # Pass the (now safe) variables down to the original impute function
  pigauto::impute(
    traits = traits,
    tree = tree, 
    n_imputations = n_imputations, 
    pool_method = pool_method, 
    ...
  )
}