#' Generate an HTML imputation report from a pigauto fit
#'
#' Produces an HTML file with charts. Baseline-only figures are embedded as
#' SVG and work offline; GNN charts load Chart.js from its public CDN. The
#' report compares the fitted baseline with the calibrated blend, which may
#' include the neural component and other baseline contributions. This is not
#' a comparison isolating the neural component. With \code{gnn = FALSE},
#' no neural network is fitted. Per-trait metrics, gate values, and held-out
#' diagnostic conformal coverage are included when available.
#'
#' @param fit A \code{pigauto_fit} object (or a \code{pigauto_result} from
#'   \code{\link{impute}}).
#' @param data Optional \code{pigauto_data} object.  Extracted automatically
#'   when \code{fit} is a \code{pigauto_result}, or reconstructed from a fit
#'   that retains its scaled observations. When supplied it must match
#'   the fit (species and traits).
#' @param splits Optional splits object from \code{\link{make_missing_splits}}.
#'   Extracted automatically from the result or fit when available.
#'   When supplied it must match the fit.
#' @param output_path Character of length 1.  File path for the HTML report
#'   (default \code{"pigauto_report.html"} in the working directory).
#' @param title Character of length 1.  Report title (HTML-escaped on write).
#' @param open Logical of length 1.  Open the report in a browser when done
#'   (default \code{TRUE}).
#' @return The output path (invisibly).
#' @examples
#' \donttest{
#' data(avonet300, tree300)
#' tree <- ape::keep.tip(tree300, tree300$tip.label[seq_len(30L)])
#' traits <- avonet300[match(tree$tip.label, avonet300$Species_Key),
#'                      c("Mass", "Wing.Length"), drop = FALSE]
#' rownames(traits) <- tree$tip.label
#' traits$Mass[seq_len(3L)] <- NA_real_
#' result <- impute(traits, tree, epochs = 5L, verbose = FALSE)
#' pigauto_report(result, output_path = tempfile(fileext = ".html"),
#'                open = FALSE)
#' }
#' @export
pigauto_report <- function(fit, data = NULL, splits = NULL,
                           output_path = "pigauto_report.html",
                           title = "pigauto Imputation Report",
                           open = TRUE) {

  .check_report_title(title)
  .check_report_output_path(output_path)
  .check_report_open(open)
  .check_report_data(data)
  .check_report_splits(splits)

  # Unwrap pigauto_result

  if (inherits(fit, "pigauto_result")) {
    if (is.null(data))   data   <- fit$data
    if (is.null(splits)) splits <- fit$splits
    pred   <- fit$prediction
    fit_obj <- fit$fit
    saved_check <- fit$check
  } else if (inherits(fit, "pigauto_fit")) {
    fit_obj <- fit
    pred <- predict(fit_obj, return_se = TRUE)
    saved_check <- NULL
  } else {
    stop("'fit' must be a pigauto_fit or pigauto_result object.")
  }

  # Current fits retain their input matrix and holdout split. Older saved fits
  # may not; in that case the report remains available without test metrics.
  if (is.null(data) && is.matrix(fit_obj$X_scaled)) {
    data <- structure(list(
      X_scaled = fit_obj$X_scaled, species_names = fit_obj$species_names,
      trait_names = fit_obj$trait_names, trait_map = fit_obj$trait_map,
      obs_species = fit_obj$obs_species, obs_to_species = fit_obj$obs_to_species
    ), class = "pigauto_data")
  }
  if (is.null(splits)) splits <- fit_obj$splits

  .check_report_data(data)
  .check_report_splits(splits)
  .check_report_data_matches_fit(data, fit_obj)
  .check_report_splits_match_fit(splits, fit_obj, data)

  trait_map <- fit_obj$trait_map
  if (is.null(trait_map)) stop("Report requires a trait_map (mixed-type fit).")

  # ---- Collect metrics -------------------------------------------------------
  # Metrics must score held-out cells without feeding their truth back into
  # the DAE context. Keep `pred` unchanged for the report's user-facing
  # prediction payload, but use a masked prediction for split metrics.
  metric_pred <- pred
  if (!is.null(data) && !is.null(splits) && length(splits$test_idx) > 0L) {
    metric_pred <- predict(
      fit_obj, return_se = TRUE,
      .mask_observed_idx = c(splits$val_idx, splits$test_idx)
    )
  }
  metrics <- collect_report_metrics(fit_obj, metric_pred, data, splits)

  # ---- Build HTML ------------------------------------------------------------
  html <- build_report_html(
    title     = .escape_report_html(title),
    fit_obj   = fit_obj,
    pred      = pred,
    metrics   = metrics,
    trait_map = trait_map
  )
  html <- sub("</body>", paste0(.report_contract_html(saved_check, fit), "</body>"),
              html, fixed = TRUE)

  writeLines(html, output_path)
  if (open) utils::browseURL(output_path)
  invisible(output_path)
}

.escape_report_html <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  gsub("'", "&#39;", x, fixed = TRUE)
}

# Encode labels as JavaScript strings without allowing a closing script tag.
.report_js_string <- function(x) {
  chars <- utf8ToInt(enc2utf8(x))
  escaped <- vapply(chars, function(cp) {
    if (cp <= 0xffffL) return(sprintf("\\u%04x", cp))
    cp <- cp - 0x10000L
    paste0(sprintf("\\u%04x", 0xd800L + cp %/% 1024L),
           sprintf("\\u%04x", 0xdc00L + cp %% 1024L))
  }, character(1))
  paste0('"', paste0(escaped, collapse = ""), '"')
}

.check_report_title <- function(title) {
  if (!is.character(title) || length(title) != 1L || is.na(title)) {
    stop("'title' must be a single character string.")
  }
  invisible(title)
}

.check_report_output_path <- function(output_path) {
  if (!is.character(output_path) || length(output_path) != 1L ||
      is.na(output_path) || !nzchar(output_path)) {
    stop("'output_path' must be a single non-empty character string.")
  }
  invisible(output_path)
}

.check_report_open <- function(open) {
  if (!is.logical(open) || length(open) != 1L || is.na(open)) {
    stop("'open' must be TRUE or FALSE.")
  }
  invisible(open)
}

.check_report_data <- function(data) {
  if (is.null(data)) return(invisible(NULL))
  if (!inherits(data, "pigauto_data")) {
    stop("'data' must be a pigauto_data object.")
  }
  invisible(data)
}

.check_report_splits <- function(splits) {
  if (is.null(splits)) return(invisible(NULL))
  ok <- is.list(splits) && !is.data.frame(splits) &&
    !is.null(splits$val_idx) && !is.null(splits$test_idx)
  if (!ok) {
    stop("'splits' must be a list from make_missing_splits().")
  }
  invisible(splits)
}

.report_latent_dims <- function(fit_obj) {
  n <- if (isTRUE(fit_obj$multi_obs)) {
    length(fit_obj$obs_to_species)
  } else {
    length(fit_obj$species_names)
  }
  p <- if (!is.null(fit_obj$baseline$mu)) {
    ncol(as.matrix(fit_obj$baseline$mu))
  } else {
    sum(vapply(fit_obj$trait_map, function(tm) length(tm$latent_cols), integer(1)))
  }
  c(n = as.integer(n), p = as.integer(p))
}

.check_report_data_matches_fit <- function(data, fit_obj) {
  if (is.null(data)) return(invisible(NULL))
  same_sp <- identical(as.character(data$species_names),
                       as.character(fit_obj$species_names))
  same_tr <- identical(as.character(data$trait_names),
                       as.character(fit_obj$trait_names))
  dims <- .report_latent_dims(fit_obj)
  same_dim <- is.matrix(data$X_scaled) &&
    nrow(data$X_scaled) == dims[["n"]] &&
    ncol(data$X_scaled) == dims[["p"]]
  if (!isTRUE(same_sp) || !isTRUE(same_tr) || !isTRUE(same_dim)) {
    stop("'data' does not match fit (species/traits).")
  }
  invisible(data)
}

.check_report_splits_match_fit <- function(splits, fit_obj, data = NULL) {
  if (is.null(splits)) return(invisible(NULL))
  dims <- .report_latent_dims(fit_obj)
  if (!is.null(data) && is.matrix(data$X_scaled)) {
    dims <- c(n = as.integer(nrow(data$X_scaled)),
              p = as.integer(ncol(data$X_scaled)))
  }
  n_cells <- dims[["n"]] * dims[["p"]]
  idx <- c(splits$val_idx, splits$test_idx)
  if (!is.numeric(idx) || anyNA(idx) ||
      length(idx) && (any(idx != as.integer(idx)) ||
                      min(idx) < 1L || max(idx) > n_cells)) {
    stop("'splits' do not match fit.")
  }
  if (!is.null(splits$n) && as.integer(splits$n) != dims[["n"]]) {
    stop("'splits' do not match fit.")
  }
  if (!is.null(splits$p) && as.integer(splits$p) != dims[["p"]]) {
    stop("'splits' do not match fit.")
  }
  invisible(splits)
}

.report_contract_html <- function(check, result = NULL) {
  input <- if (is.null(check)) {
    "<p>Preflight check unavailable.</p>"
  } else {
    msgs <- check$messages
    rows <- if (is.null(msgs) || !nrow(msgs)) "<tr><td colspan=\"3\">No messages</td></tr>" else paste(vapply(seq_len(nrow(msgs)), function(i) {
      sprintf("<tr><td>%s</td><td>%s</td><td>%s</td></tr>",
              .escape_report_html(msgs$code[[i]]), .escape_report_html(msgs$text[[i]]),
              .escape_report_html(msgs$action[[i]]))
    }, character(1)), collapse = "")
    paste0("<p>Status: <strong>", .escape_report_html(check$status),
           "</strong></p><table><thead><tr><th>Code</th><th>Message</th><th>Action</th></tr></thead><tbody>", rows, "</tbody></table>")
  }
  roles <- if (inherits(result, "pigauto_result")) .describe_pigauto_result(result) else data.frame(role = "legacy", text = "Preflight check unavailable for a bare fit.")
  role_rows <- paste(vapply(seq_len(nrow(roles)), function(i) sprintf("<li><strong>%s:</strong> %s</li>", .escape_report_html(roles$role[[i]]), .escape_report_html(roles$text[[i]])), character(1)), collapse = "")
  paste0("<h2>Input checks</h2>", input, "<h2>How to use outputs</h2><ul>", role_rows, "</ul>")
}


# ---- Internal: collect per-trait metrics ------------------------------------

collect_report_metrics <- function(fit_obj, pred, data, splits) {
  trait_map <- fit_obj$trait_map
  baseline <- fit_obj$baseline
  results <- list()
  for (tm in trait_map) {
    nm <- tm$name
    tp <- tm$type
    entry <- list(name = nm, type = tp)
    if (!is.null(fit_obj$calibrated_gates)) {
      entry$gate <- mean(fit_obj$calibrated_gates[tm$latent_cols])
    }
    if (!is.null(fit_obj$conformal_scores) && nm %in% names(fit_obj$conformal_scores)) {
      entry$conformal_score <- fit_obj$conformal_scores[nm]
    }
    if (!is.null(data) && !is.null(splits) && length(splits$test_idx) > 0) {
      X_truth <- data$X_scaled
      n <- nrow(X_truth)
      lc <- tm$latent_cols
      test_mat <- matrix(FALSE, n, ncol(X_truth))
      test_mat[splits$test_idx] <- TRUE
      test_rows <- which(test_mat[, lc[1]])
      if (length(test_rows) > 0) {
        entry$n_test <- length(test_rows)
        bm <- baseline$mu
        if (!is.null(bm) && isTRUE(fit_obj$multi_obs)) {
          bm <- bm[fit_obj$obs_to_species, , drop = FALSE]
        }
        if (tp %in% c("continuous", "count", "ordinal", "proportion")) {
          truth <- X_truth[test_rows, lc[1]]
          gnn <- if (is.null(pred$imputed_latent)) rep(NA_real_, length(test_rows)) else pred$imputed_latent[test_rows, lc[1]]
          base <- if (is.null(bm)) rep(NA_real_, length(test_rows)) else bm[test_rows, lc[1]]
          valid_gnn <- is.finite(truth) & is.finite(gnn)
          valid_bm <- is.finite(truth) & is.finite(base)
          paired <- valid_gnn & valid_bm
          entry$n_gnn <- sum(valid_gnn); entry$n_bm <- sum(valid_bm); entry$n_paired <- sum(paired)
          if (entry$n_gnn) {
            entry$rmse_gnn <- sqrt(mean((truth[valid_gnn] - gnn[valid_gnn])^2))
            if (entry$n_gnn >= 2L) entry$r_gnn <- suppressWarnings(stats::cor(truth[valid_gnn], gnn[valid_gnn]))
          }
          if (entry$n_bm) {
            entry$rmse_bm <- sqrt(mean((truth[valid_bm] - base[valid_bm])^2))
            if (entry$n_bm >= 2L) entry$r_bm <- suppressWarnings(stats::cor(truth[valid_bm], base[valid_bm]))
          }
          if (entry$n_paired) {
            entry$rmse_paired_bm <- sqrt(mean((truth[paired] - base[paired])^2))
            entry$rmse_paired_gnn <- sqrt(mean((truth[paired] - gnn[paired])^2))
          }
          entry$n_coverage <- 0L
          if (!is.null(pred$conformal_lower) && !is.null(pred$conformal_upper) &&
              nm %in% colnames(pred$conformal_lower) && nm %in% colnames(pred$conformal_upper)) {
            low <- pred$conformal_lower[test_rows, nm]
            high <- pred$conformal_upper[test_rows, nm]
            truth_orig <- truth * tm$sd + tm$mean
            if (tp == "continuous" && isTRUE(tm$log_transform)) truth_orig <- exp(truth_orig)
            if (tp == "count") truth_orig <- round(expm1(truth_orig))
            if (tp == "ordinal") truth_orig <- round(truth_orig)
            if (tp == "proportion") truth_orig <- stats::plogis(truth_orig)
            valid <- is.finite(low) & is.finite(high) & is.finite(truth_orig) & low <= high
            entry$n_coverage <- sum(valid)
            if (entry$n_coverage) entry$conformal_coverage <- mean(truth_orig[valid] >= low[valid] & truth_orig[valid] <= high[valid])
          }
        } else if (tp %in% c("binary", "categorical")) {
          if (tp == "binary") {
            truth <- X_truth[test_rows, lc[1]]
            prob <- pred$probabilities[[nm]]
            gnn <- if (is.null(prob)) rep(NA_real_, length(test_rows)) else ifelse(is.finite(prob[test_rows]), as.numeric(prob[test_rows] >= 0.5), NA_real_)
            base <- if (is.null(bm)) rep(NA_real_, length(test_rows)) else ifelse(is.finite(bm[test_rows, lc[1]]), as.numeric(bm[test_rows, lc[1]] >= 0), NA_real_)
          } else {
            # A row with an unavailable class score cannot be assigned a class.
            classes <- function(x) vapply(seq_len(nrow(x)), function(i) {
              if (any(!is.finite(x[i, ]))) NA_integer_ else which.max(x[i, ])
            }, integer(1))
            truth <- classes(X_truth[test_rows, lc, drop = FALSE])
            probs <- pred$probabilities[[nm]]
            gnn <- if (is.null(probs) || !is.matrix(probs)) rep(NA_integer_, length(test_rows)) else classes(probs[test_rows, , drop = FALSE])
            base <- if (is.null(bm)) rep(NA_integer_, length(test_rows)) else classes(bm[test_rows, lc, drop = FALSE])
          }
          valid_gnn <- is.finite(truth) & is.finite(gnn)
          valid_bm <- is.finite(truth) & is.finite(base)
          paired <- valid_gnn & valid_bm
          entry$n_gnn <- sum(valid_gnn); entry$n_bm <- sum(valid_bm); entry$n_paired <- sum(paired)
          if (entry$n_gnn) entry$accuracy_gnn <- mean(gnn[valid_gnn] == truth[valid_gnn])
          if (entry$n_bm) entry$accuracy_bm <- mean(base[valid_bm] == truth[valid_bm])
          if (entry$n_paired) {
            entry$accuracy_paired_bm <- mean(base[paired] == truth[paired])
            entry$accuracy_paired_gnn <- mean(gnn[paired] == truth[paired])
          }
          if (tp == "binary" && !is.null(prob)) {
            valid_prob <- is.finite(truth) & is.finite(prob[test_rows])
            if (any(valid_prob)) entry$brier_gnn <- mean((prob[test_rows][valid_prob] - truth[valid_prob])^2)
          }
        }
      }
    }
    results[[nm]] <- entry
  }
  results
}

# ---- Internal: build HTML string -------------------------------------------

build_report_html <- function(title, fit_obj, pred, metrics, trait_map) {
  n_species <- length(fit_obj$species_names)
  n_traits  <- length(trait_map)
  types     <- vapply(trait_map, "[[", character(1), "type")
  type_tab  <- table(types)
  type_str  <- paste(names(type_tab), type_tab, sep = "=", collapse = ", ")

  has_attention  <- isTRUE(fit_obj$model_config$use_attention)
  has_calibration <- !is.null(fit_obj$calibrated_gates)
  has_conformal   <- !is.null(fit_obj$conformal_scores)
  k_eigen <- fit_obj$model_config$k_eigen

  # The final prediction is a calibrated blend, not an isolated GNN prediction.
  # With gnn = FALSE the blend may still differ from the baseline.
  gnn_off <- isFALSE(fit_obj$model_config$gnn)
  # Build metrics table rows
  trait_rows <- ""
  for (m in metrics) {
    nm <- m$name
    tp <- m$type

    gate_val <- if (!is.null(m$gate)) sprintf("%.3f", m$gate) else "&mdash;"

    if (tp %in% c("continuous", "count", "ordinal", "proportion")) {
      bm_rmse  <- if (isTRUE(is.finite(m$rmse_bm))) sprintf("%.4f (n=%d)", m$rmse_bm, m$n_bm) else "&mdash;"
      gnn_rmse <- if (isTRUE(is.finite(m$rmse_gnn))) sprintf("%.4f (n=%d)", m$rmse_gnn, m$n_gnn) else "&mdash;"
      bm_r     <- if (isTRUE(is.finite(m$r_bm))) sprintf("%.3f", m$r_bm) else "&mdash;"
      gnn_r    <- if (isTRUE(is.finite(m$r_gnn))) sprintf("%.3f", m$r_gnn) else "&mdash;"
      n_test   <- if (!is.null(m$n_test))   m$n_test                     else "&mdash;"

      pct <- ""
      verdict <- '<span class="pill pill-neutral">no data</span>'
      if (isTRUE(is.finite(m$rmse_paired_bm)) && isTRUE(is.finite(m$rmse_paired_gnn)) && m$rmse_paired_bm > 0) {
        pct_val <- (m$rmse_paired_bm - m$rmse_paired_gnn) / m$rmse_paired_bm * 100
        if (pct_val > 0.5) {
          pct <- sprintf('<span class="pos">+%.1f%%</span>', pct_val)
          verdict <- '<span class="pill pill-green">Blend better</span>'
        } else if (pct_val < -0.5) {
          pct <- sprintf('<span class="neg">%.1f%%</span>', pct_val)
          verdict <- '<span class="pill pill-red">Baseline better</span>'
        } else {
          pct <- sprintf('<span class="zero">%.1f%%</span>', pct_val)
          verdict <- '<span class="pill pill-neutral">tied</span>'
        }
        pct <- paste0(pct, sprintf(" (paired n=%d)", m$n_paired))
        if (gnn_off) verdict <- '<span class="pill pill-neutral">GNN off</span>'
      }

      cov <- "&mdash;"
      if (isTRUE(is.finite(m$conformal_coverage))) {
        cov_val <- m$conformal_coverage * 100
        cov <- sprintf('%.1f%% (coverage n=%d)', cov_val, m$n_coverage)
      }

      trait_rows <- paste0(trait_rows, sprintf(
        '<tr><td>%s</td><td>%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td>%s</td></tr>\n',
        .escape_report_html(nm), .escape_report_html(tp), n_test, bm_rmse, gnn_rmse, pct, bm_r, gnn_r, cov, verdict
      ))

    } else {
      bm_acc  <- if (isTRUE(is.finite(m$accuracy_bm))) sprintf("%.1f%% (n=%d)", m$accuracy_bm * 100, m$n_bm) else "&mdash;"
      gnn_acc <- if (isTRUE(is.finite(m$accuracy_gnn))) sprintf("%.1f%% (n=%d)", m$accuracy_gnn * 100, m$n_gnn) else "&mdash;"
      n_test  <- if (!is.null(m$n_test))       m$n_test                                 else "&mdash;"

      pct <- ""
      verdict <- '<span class="pill pill-neutral">no data</span>'
      if (isTRUE(is.finite(m$accuracy_paired_bm)) && isTRUE(is.finite(m$accuracy_paired_gnn))) {
        diff <- (m$accuracy_paired_gnn - m$accuracy_paired_bm) * 100
        if (diff > 0.5) {
          pct <- sprintf('<span class="pos">+%.1f pp</span>', diff)
          verdict <- '<span class="pill pill-green">Blend better</span>'
        } else if (diff < -0.5) {
          pct <- sprintf('<span class="neg">%.1f pp</span>', diff)
          verdict <- '<span class="pill pill-red">Baseline better</span>'
        } else {
          pct <- sprintf('<span class="zero">%.1f pp</span>', diff)
          verdict <- '<span class="pill pill-neutral">tied</span>'
        }
        pct <- paste0(pct, sprintf(" (paired n=%d)", m$n_paired))
        if (gnn_off) verdict <- '<span class="pill pill-neutral">GNN off</span>'
      }

      trait_rows <- paste0(trait_rows, sprintf(
        '<tr><td>%s</td><td>%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num" colspan="2">&mdash;</td><td class="num">&mdash;</td><td>%s</td></tr>\n',
        .escape_report_html(nm), .escape_report_html(tp), n_test, bm_acc, gnn_acc, pct, verdict
      ))
    }
  }

  # Build gate chart data
  gate_names  <- vapply(trait_map, "[[", character(1), "name")
  gate_values <- if (has_calibration) {
    vapply(metrics, function(m) if (!is.null(m$gate)) m$gate else 0, numeric(1))
  } else {
    rep(0, n_traits)
  }
  gate_types <- vapply(trait_map, "[[", character(1), "type")
  gate_colors <- ifelse(gate_types == "continuous", "'rgba(25,135,84,0.7)'",
                 ifelse(gate_types == "count",      "'rgba(13,110,253,0.7)'",
                 ifelse(gate_types == "ordinal",    "'rgba(255,193,7,0.7)'",
                 ifelse(gate_types == "binary",     "'rgba(220,53,69,0.7)'",
                                                    "'rgba(111,66,193,0.7)'"))))

  gate_labels_js <- paste0("[", paste(vapply(gate_names, .report_js_string, character(1)), collapse = ","), "]")
  gate_values_js <- paste0("[", paste(sprintf("%.4f", gate_values), collapse = ","), "]")
  gate_colors_js <- paste0("[", paste(gate_colors, collapse = ","), "]")

  # Training history chart data
  hist <- fit_obj$history
  if (!is.null(hist) && nrow(hist) > 0) {
    hist_epochs <- paste(hist$epoch, collapse = ",")
    hist_rec    <- paste(sprintf("%.4f", hist$loss_rec), collapse = ",")
    hist_val    <- paste(ifelse(is.na(hist$val_loss), "null",
                                sprintf("%.4f", hist$val_loss)), collapse = ",")
  } else {
    hist_epochs <- ""
    hist_rec    <- ""
    hist_val    <- ""
  }

  # Static figures for the default baseline-only fit. Each measure has its
  # own scale: continuous latent-scale RMSE and discrete accuracy are not
  # directly comparable, so they are displayed in separate panels.
  baseline_figures <- if (gnn_off) {
    paste0('<div class="chart-grid">',
      .report_svg_metrics(metrics, "rmse_gnn", "Held-out RMSE (latent scale)"),
      .report_svg_metrics(metrics, "accuracy_gnn", "Held-out accuracy", percent = TRUE),
      .report_svg_metrics(metrics, "conformal_score", "Conformal scores (latent scale)"),
      '</div>')
  } else NULL

  # Features list
  features <- c(
    if (gnn_off) "GNN off: phylogenetic baseline only (no torch)" else
      if (has_attention) "Attention message passing" else "Standard message passing",
    if (has_calibration) "Validation-calibrated gates",
    if (has_conformal)   "Conformal prediction intervals",
    sprintf("Adaptive spectral encoding (k=%d)", k_eigen),
    "Phylogenetic label propagation (discrete baselines)"
  )
  features_html <- paste0("<li>", features, "</li>", collapse = "\n")

  timestamp <- format(Sys.time(), "%d %B %Y, %H:%M")

  # Assemble HTML
  html <- sprintf('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>%s</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.4/dist/chart.umd.min.js"></script>
<style>
:root{--bg:#f8f9fa;--card:#fff;--border:#dee2e6;--text:#212529;--muted:#6c757d;--accent:#0d6efd;--green:#198754;--red:#dc3545;}
*{box-sizing:border-box;margin:0;padding:0;}
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;background:var(--bg);color:var(--text);line-height:1.6;padding:2rem;max-width:1300px;margin:0 auto;}
h1{font-size:1.8rem;margin-bottom:.3rem;}
h2{font-size:1.3rem;margin:2rem 0 1rem;border-bottom:2px solid var(--accent);padding-bottom:.3rem;}
h3{font-size:1.05rem;margin:1rem 0 .5rem;color:var(--muted);}
.subtitle{color:var(--muted);margin-bottom:1.5rem;font-size:.95rem;}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:1rem;margin-bottom:2rem;}
.card{background:var(--card);border:1px solid var(--border);border-radius:8px;padding:1.2rem;text-align:center;}
.card .label{font-size:.78rem;color:var(--muted);text-transform:uppercase;letter-spacing:.05em;}
.card .value{font-size:1.8rem;font-weight:700;margin:.2rem 0;}
.card .detail{font-size:.82rem;color:var(--muted);}
.card.green .value{color:var(--green);}
.card.red .value{color:var(--red);}
.card.blue .value{color:var(--accent);}
table{width:100%%;border-collapse:collapse;background:var(--card);border-radius:8px;overflow:hidden;margin-bottom:1.5rem;font-size:.85rem;}
thead{background:#f1f3f5;}
th{padding:.55rem .7rem;text-align:left;font-weight:600;font-size:.78rem;text-transform:uppercase;letter-spacing:.03em;color:var(--muted);border-bottom:2px solid var(--border);}
td{padding:.5rem .7rem;border-bottom:1px solid #f1f3f5;}
tr:hover{background:#f8f9fa;}
.num{text-align:right;font-variant-numeric:tabular-nums;}
.pos{color:var(--green);font-weight:600;}
.neg{color:var(--red);font-weight:600;}
.zero{color:var(--muted);}
.pill{display:inline-block;padding:.12rem .45rem;border-radius:10px;font-size:.72rem;font-weight:600;}
.pill-green{background:#d1e7dd;color:var(--green);}
.pill-red{background:#f8d7da;color:var(--red);}
.pill-neutral{background:#e9ecef;color:var(--muted);}
.chart-grid{display:grid;grid-template-columns:1fr 1fr;gap:1.5rem;margin-bottom:2rem;}
.chart-box{background:var(--card);border:1px solid var(--border);border-radius:8px;padding:1rem;}
.features{background:var(--card);border:1px solid var(--border);border-radius:8px;padding:1.2rem;margin-bottom:1.5rem;}
.features ul{margin-left:1.2rem;}
.features li{margin-bottom:.3rem;}
.footer{text-align:center;color:var(--muted);font-size:.8rem;margin-top:3rem;padding-top:1rem;border-top:1px solid var(--border);}
@media(max-width:800px){.chart-grid{grid-template-columns:1fr;}}
</style>
</head>
<body>
<h1>%s</h1>
<p class="subtitle">%d species &bull; %d traits (%s) &bull; generated %s</p>

<div class="cards">
<div class="card blue"><div class="label">Species</div><div class="value">%d</div><div class="detail">in phylogeny</div></div>
<div class="card blue"><div class="label">Traits</div><div class="value">%d</div><div class="detail">%s</div></div>
<div class="card green"><div class="label">Val loss</div><div class="value">%.4f</div><div class="detail">best validation</div></div>
<div class="card %s"><div class="label">Test loss</div><div class="value">%s</div><div class="detail">held-out</div></div>
<div class="card blue"><div class="label">Attention</div><div class="value">%s</div><div class="detail">message passing</div></div>
<div class="card blue"><div class="label">k_eigen</div><div class="value">%d</div><div class="detail">spectral features</div></div>
</div>

<div class="features">
<h3>Enabled features</h3>
<ul>%s</ul>
</div>

<h2>Per-Trait Performance</h2>
<p>The baseline column shows the fitted baseline; the calibrated blend may include neural and other baseline contributions. This comparison does not isolate the neural component.</p>
<p>Held-out cells are known values hidden during fitting. The root mean squared error (RMSE) measures prediction error on the transformed, standardised scale; lower is better. Accuracy is the percentage of held-out binary or categorical values predicted correctly. Each metric shows the number of finite values scored. Differences compare only cells scored by both methods. Coverage is the percentage of held-out values within the prediction bounds; its count may differ from the prediction counts.</p>
<table>
<thead><tr><th>Trait</th><th>Type</th><th class="num">n held out</th><th class="num">Baseline</th><th class="num">ARM_LABEL</th><th class="num">%%&Delta;</th><th class="num">r (Baseline)</th><th class="num">r (ARM_LABEL)</th><th class="num">Held-out coverage (diagnostic)</th><th>Verdict</th></tr></thead>
<tbody>%s</tbody>
</table>

<div class="chart-grid">
<div class="chart-box"><h3>Calibrated gate values per trait</h3><canvas id="gateChart"></canvas></div>
<div class="chart-box"><h3>Training history</h3><canvas id="histChart"></canvas></div>
</div>

<div class="footer">
pigauto &bull; Report generated %s &bull; <a href="https://github.com/itchyshin/pigauto">github.com/itchyshin/pigauto</a>
</div>

<script>
Chart.defaults.font.family="-apple-system,BlinkMacSystemFont,sans-serif";
Chart.defaults.font.size=12;
new Chart(document.getElementById("gateChart"),{type:"bar",data:{labels:%s,datasets:[{data:%s,backgroundColor:%s,borderWidth:0}]},options:{responsive:true,plugins:{legend:{display:false}},scales:{y:{title:{display:true,text:"Gate value (0=baseline, cap=full GNN)"},beginAtZero:true,max:0.85},x:{grid:{display:false}}}}});
var hEp=[%s],hRec=[%s],hVal=[%s];
if(hEp.length>0){new Chart(document.getElementById("histChart"),{type:"line",data:{labels:hEp,datasets:[{label:"Rec loss",data:hRec,borderColor:"rgba(13,110,253,.8)",backgroundColor:"transparent",pointRadius:2,tension:.3},{label:"Val loss",data:hVal,borderColor:"rgba(220,53,69,.8)",backgroundColor:"transparent",pointRadius:2,tension:.3}]},options:{responsive:true,plugins:{legend:{position:"top"}},scales:{y:{title:{display:true,text:"Loss"},grid:{color:"#eee"}},x:{title:{display:true,text:"Epoch"},grid:{display:false}}}}});}
</script>
</body></html>',
    title,                                               # page title
    title,                                               # h1
    n_species, n_traits, type_str, timestamp,            # subtitle
    n_species,                                           # card 1
    n_traits, type_str,                                  # card 2
    fit_obj$val_rmse,                                    # card 3
    if (is.na(fit_obj$test_rmse)) "blue" else "green",  # card 4 color
    if (is.na(fit_obj$test_rmse)) "&mdash;" else sprintf("%.4f", fit_obj$test_rmse),
    if (has_attention) "ON" else "OFF",                  # card 5
    k_eigen,                                             # card 6
    features_html,                                       # features
    trait_rows,                                           # table
    timestamp,                                            # footer
    gate_labels_js, gate_values_js, gate_colors_js,      # gate chart
    hist_epochs, hist_rec, hist_val                       # history chart
  )
  if (gnn_off) {
    start <- regexpr('<div class="chart-grid">', html, fixed = TRUE)[1]
    end <- regexpr('<div class="footer">', html, fixed = TRUE)[1]
    html <- paste0(substr(html, 1, start - 1L), baseline_figures, substr(html, end, nchar(html)))
    html <- sub('<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.4/dist/chart.umd.min.js"></script>', '', html, fixed = TRUE)
    script <- regexpr('<script>\nChart.defaults', html)[1]
    if (script > 0L) html <- paste0(substr(html, 1, script - 1L), '</body></html>')
  }
  gsub("ARM_LABEL", if (gnn_off) "Calibrated blend (GNN off)" else "Calibrated blend",
       html, fixed = TRUE)
}

.report_svg_metrics <- function(metrics, field, heading, percent = FALSE) {
  vals <- vapply(metrics, function(m) {
    v <- m[[field]]
    if (is.null(v) || length(v) != 1L || !is.finite(v)) NA_real_ else as.numeric(v)
  }, numeric(1))
  keep <- which(is.finite(vals) & vals >= 0)
  if (!length(keep)) return(paste0('<div class="chart-box"><h3>', heading,
                                    '</h3><p>No data available.</p></div>'))
  vals <- vals[keep]
  labs <- vapply(metrics[keep], function(m) .escape_report_html(m$name), character(1))
  scale <- if (percent) 1 else max(vals, 1e-12)
  height <- 35 + 30 * length(vals)
  bars <- vapply(seq_along(vals), function(i) {
    y <- 15 + (i - 1L) * 30
    sprintf('<text x="5" y="%d">%s</text><rect x="160" y="%d" width="%.2f" height="16" fill="#397a9b"/><text x="%.2f" y="%d">%s</text>',
            y + 13, labs[i], y, 240 * vals[i] / scale, 165 + 240 * vals[i] / scale,
            y + 13, if (percent) sprintf("%.1f%%", 100 * vals[i]) else sprintf("%.3f", vals[i]))
  }, character(1))
  paste0('<div class="chart-box"><h3>', heading, '</h3><svg role="img" aria-label="', heading,
         '" viewBox="0 0 480 ', height, '" xmlns="http://www.w3.org/2000/svg">',
         paste(bars, collapse = ""), '</svg></div>')
}
