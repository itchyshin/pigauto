# No refitting. Input: list(traits=<original single-observation frame>,
# tree=<matching phylo>, result=<saved pigauto_result>, trait="log_ls").
diagnose_ties <- function(bundle, digits=8L, max_pairs=20L) {
  d <- bundle$traits
  tree <- bundle$tree
  result <- bundle$result
  trait <- bundle$trait
  if(!is.data.frame(d)||!inherits(tree,"phylo")||!is.character(trait)||length(trait)!=1L||
     !trait %in% names(d)||!is.numeric(d[[trait]])||anyDuplicated(rownames(d)))
    stop("Supply a numeric single-observation trait and uniquely named species rows")
  if(!is.data.frame(result$completed)||!trait %in% names(result$completed)||
     !identical(rownames(d),rownames(result$completed)))
    stop("Completed rows must match original input rows exactly")
  if(isTRUE(result$data$multi_obs)) stop("This diagnostic is single-observation only")
  if(!identical(tree,result$tree)) stop("Use the exact tree saved in result$tree")
  observed <- !is.na(d[[trait]])
  if(any(!is.finite(d[[trait]][observed])) ||
     !identical(as.numeric(d[[trait]][observed]),
                as.numeric(result$completed[[trait]][observed])))
    stop("Observed values must agree with the saved completion")
  if(length(digits)!=1L||!is.finite(digits)||digits<0||digits!=as.integer(digits)||
     length(max_pairs)!=1L||!is.finite(max_pairs)||max_pairs<1||max_pairs!=as.integer(max_pairs))
    stop("Invalid digits or max_pairs")
  modeled <- rownames(d) %in% tree$tip.label
  missing <- which(is.na(d[[trait]]) & modeled)
  values <- result$completed[[trait]][missing]
  if(!is.numeric(values)||any(!is.finite(values))) stop("Missing-cell completions must be finite")
  tied <- function(x) duplicated(x)|duplicated(x,fromLast=TRUE)
  exact <- tied(values)
  rounded <- tied(round(values,digits))
  summary <- data.frame(n_modeled_missing=length(missing),
    n_unmodeled_missing=sum(is.na(d[[trait]]) & !modeled),
    n_unique_exact=length(unique(values)),n_exact_tied=sum(exact),
    n_unique_rounded=length(unique(round(values,digits))),n_rounded_tied=sum(rounded))
  # A bounded list of exact-tied pairs, with row indices rather than names.
  # Correlation-row identity is only a possible univariate-BM mechanism.
  groups <- split(seq_along(values),match(values,unique(values)))
  pairs <- do.call(rbind,lapply(groups,function(g) if(length(g)>1L)
    cbind(g[1],g[-1]) else NULL))
  checks <- data.frame()
  if(!is.null(pairs)) {
    pairs <- pairs[seq_len(min(nrow(pairs),max_pairs)),,drop=FALSE]
    observed <- rownames(d)[!is.na(d[[trait]]) & modeled]
    if(!length(observed)) stop("No observed trait values for covariance comparison")
    R <- stats::cov2cor(ape::vcv.phylo(tree))
    checks <- do.call(rbind,lapply(seq_len(nrow(pairs)),function(k) {
      rows <- missing[pairs[k,]]
      difference <- max(abs(R[rownames(d)[rows[1]],observed]-R[rownames(d)[rows[2]],observed]))
      data.frame(input_row_1=rows[1],input_row_2=rows[2],
                 max_observed_correlation_difference=difference)
    }))
  }
  list(summary=summary,pairs=checks,
    route=result$fit$baseline_full$predict_method_by_trait,
    lambda=result$fit$baseline_full$lambda_per_trait,
    interpretation="Covariance-row identity can explain equal univariate BM means; it does not diagnose joint-model, decoding or other causes.")
}

if(sys.nframe()==0L) {
  args <- commandArgs(TRUE)
  if(length(args)!=2L) stop("Usage: Rscript diagnose_ties.R bundle.rds report.rds")
  saveRDS(diagnose_ties(readRDS(args[1])),args[2])
}
