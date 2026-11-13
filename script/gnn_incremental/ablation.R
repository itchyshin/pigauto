# Same-fit counterfactual: remove only the calibrated neural contribution.

ablate_gnn_in_fit <- function(fit) {
  if (!inherits(fit, "pigauto_fit")) {
    stop("fit must be a pigauto_fit object.", call. = FALSE)
  }
  if (!isTRUE(fit$model_config$gnn)) {
    stop("The supplied fit must have been trained with gnn = TRUE.",
         call. = FALSE)
  }
  slots <- c("r_cal_bm", "r_cal_gnn", "r_cal_mean")
  if (any(vapply(slots, function(nm) is.null(fit[[nm]]), logical(1)))) {
    stop("The fit is missing three-way calibrated gate state.", call. = FALSE)
  }
  if (length(fit$r_cal_bm) != length(fit$r_cal_gnn) ||
      length(fit$r_cal_bm) != length(fit$r_cal_mean) ||
      any(!is.finite(fit$r_cal_bm)) || any(!is.finite(fit$r_cal_gnn)) ||
      any(!is.finite(fit$r_cal_mean))) {
    stop("The fit's calibrated gate vectors are malformed.", call. = FALSE)
  }
  if (any(fit$r_cal_mean != 0) && is.null(fit$mean_baseline_per_col)) {
    stop("The fit has nonzero mean weights but no stored mean baseline.",
         call. = FALSE)
  }

  ablated <- fit
  ablated$r_cal_bm <- fit$r_cal_bm + fit$r_cal_gnn
  ablated$r_cal_gnn <- rep(0, length(fit$r_cal_gnn))
  ablated$calibrated_gates <- ablated$r_cal_gnn
  ablated$r_cal <- ablated$r_cal_gnn
  ablated
}
