
#' Plot uncertainty heatmap
#' @param result phyloimpute_result
#' @export
plot_uncertainty <- function(result) {
  if (is.null(result$uncertainty))
    stop('Uncertainty not stored (placeholder).')
  df <- as.data.frame(result$uncertainty)
  df$obs <- seq_len(nrow(df))
  df_long <- tidyr::pivot_longer(df, -obs, names_to='trait', values_to='sd')
  ggplot(df_long, aes(x=obs, y=trait, fill=sd)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_minimal()
}

