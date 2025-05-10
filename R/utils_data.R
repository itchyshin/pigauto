
#' Detect trait column types
#' @param df data.frame of traits
#' @return list of column names by type
#' @noRd
detect_trait_types <- function(df) {
  cont    <- names(df)[sapply(df, is.numeric)]
  ordered <- names(df)[sapply(df, is.ordered)]
  factor_cols <- names(df)[sapply(df, is.factor)]
  binary  <- factor_cols[sapply(df[factor_cols], function(x) nlevels(x)==2)]
  nominal <- setdiff(factor_cols, c(ordered, binary))
  list(continuous = cont,
       ordinal    = ordered,
       binary     = binary,
       nominal    = nominal)
}

