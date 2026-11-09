# Stored-draw sensitivity only. No model is fitted and no private inputs ship.
# Input is a named list of trees. Each element has completed_draws (a named
# list of equal-length numeric rho_E draw vectors) and reference_draws.
# Correlations must be on their raw [-1,1] scale, not Fisher-z transformed.
subset_intervals <- function(trees, m_grid=c(5L,10L,20L), repetitions=200L, seed=2511109L) {
  valid <- function(x) is.numeric(x) && length(x)>1L && all(is.finite(x)) && all(abs(x)<=1)
  if (!is.list(trees) || !length(trees) || is.null(names(trees)) ||
      any(!nzchar(names(trees))) || anyDuplicated(names(trees))) stop("Name each tree uniquely")
  if (!is.numeric(m_grid) || !length(m_grid) || anyNA(m_grid) || any(!is.finite(m_grid)) ||
      any(m_grid<2 | m_grid!=as.integer(m_grid))) stop("m_grid must contain integers >=2")
  if (length(repetitions)!=1 || !is.finite(repetitions) || repetitions<1 ||
      repetitions!=as.integer(repetitions)) stop("Invalid repetitions")
  set.seed(seed)
  output <- list()
  for (tree in names(trees)) {
    x <- trees[[tree]]
    draws <- x$completed_draws
    if (!valid(x$reference_draws) || !is.list(draws) || length(draws)<2 ||
        is.null(names(draws)) || any(!nzchar(names(draws))) || anyDuplicated(names(draws)) ||
        !all(vapply(draws, valid, logical(1))) || length(unique(lengths(draws)))!=1L)
      stop("Invalid or unequal-length draws for tree ", tree)
    if (any(m_grid>length(draws))) stop("Requested m exceeds stored fit count for tree ", tree,
      "; choose a fully available grid. No levels were silently dropped.")
    ref <- unname(quantile(x$reference_draws, c(.025,.975), type=7))
    ref_width <- diff(ref)
    if (ref_width<=0) stop("Reference interval has zero width")
    for (m in unique(m_grid)) {
      for (rep in seq_len(if(m==length(draws)) 1L else repetitions)) {
        selected <- if(m==length(draws)) seq_along(draws) else sample.int(length(draws), m)
        pooled <- unlist(draws[selected], use.names=FALSE)
        interval <- unname(quantile(pooled, c(.025,.975), type=7))
        output[[length(output)+1L]] <- data.frame(tree=tree,m=m,rep=rep,
          completed_ids=paste(names(draws)[selected],collapse=";"),
          mean=mean(pooled),lower=interval[1],upper=interval[2],
          reference_width=ref_width,width_ratio=diff(interval)/ref_width)
      }
    }
  }
  if(!length(output)) stop("No requested m is available; cannot extrapolate beyond stored fits")
  do.call(rbind,output)
}

if (sys.nframe()==0L) {
  args <- commandArgs(TRUE)
  if(length(args)!=2L) stop("Usage: Rscript subset_intervals.R input.rds output.csv")
  write.csv(subset_intervals(readRDS(args[1])),args[2],row.names=FALSE)
}
