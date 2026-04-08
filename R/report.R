#' Generate an HTML benchmark report from a pigauto fit
#'
#' Produces a self-contained HTML file with interactive charts comparing the
#' GNN against the phylogenetic baseline.  The report includes per-trait
#' metrics, gate values, conformal coverage, and training history.
#'
#' @param fit A \code{pigauto_fit} object (or a \code{pigauto_result} from
#'   \code{\link{impute}}).
#' @param data Optional \code{pigauto_data} object.  Extracted automatically
#'   when \code{fit} is a \code{pigauto_result}.
#' @param splits Optional splits object.  Extracted automatically when
#'   \code{fit} is a \code{pigauto_result}.
#' @param output_path Character.  File path for the HTML report (default
#'   \code{"pigauto_report.html"} in the working directory).
#' @param title Character.  Report title.
#' @param open Logical.  Open the report in a browser when done (default
#'   \code{TRUE}).
#' @return The output path (invisibly).
#' @export
pigauto_report <- function(fit, data = NULL, splits = NULL,
                           output_path = "pigauto_report.html",
                           title = "pigauto Imputation Report",
                           open = TRUE) {

  # Unwrap pigauto_result

  if (inherits(fit, "pigauto_result")) {
    if (is.null(data))   data   <- fit$data
    if (is.null(splits)) splits <- fit$splits
    pred   <- fit$prediction
    fit_obj <- fit$fit
  } else if (inherits(fit, "pigauto_fit")) {
    fit_obj <- fit
    pred <- predict(fit_obj, return_se = TRUE)
  } else {
    stop("'fit' must be a pigauto_fit or pigauto_result object.")
  }

  trait_map <- fit_obj$trait_map
  if (is.null(trait_map)) stop("Report requires a trait_map (mixed-type fit).")

  # ---- Collect metrics -------------------------------------------------------
  metrics <- collect_report_metrics(fit_obj, pred, data, splits)

  # ---- Build HTML ------------------------------------------------------------
  html <- build_report_html(
    title     = title,
    fit_obj   = fit_obj,
    pred      = pred,
    metrics   = metrics,
    trait_map = trait_map
  )

  writeLines(html, output_path)
  if (open) utils::browseURL(output_path)
  invisible(output_path)
}


# ---- Internal: collect per-trait metrics ------------------------------------

collect_report_metrics <- function(fit_obj, pred, data, splits) {
  trait_map <- fit_obj$trait_map
  baseline  <- fit_obj$baseline

  results <- list()
  for (tm in trait_map) {
    nm <- tm$name
    tp <- tm$type
    entry <- list(name = nm, type = tp)

    # Gate value
    if (!is.null(fit_obj$calibrated_gates)) {
      entry$gate <- mean(fit_obj$calibrated_gates[tm$latent_cols])
    }

    # Conformal score
    if (!is.null(fit_obj$conformal_scores) && nm %in% names(fit_obj$conformal_scores)) {
      entry$conformal_score <- fit_obj$conformal_scores[nm]
    }

    # Test-set evaluation (when data + splits available)
    if (!is.null(data) && !is.null(splits) && length(splits$test_idx) > 0) {
      X_truth <- data$X_scaled
      n <- nrow(X_truth)
      p <- ncol(X_truth)
      lc <- tm$latent_cols

      # Extract test cells for this trait
      test_mat <- matrix(FALSE, n, p)
      test_mat[splits$test_idx] <- TRUE
      test_rows <- which(test_mat[, lc[1]])

      if (length(test_rows) > 0) {
        entry$n_test <- length(test_rows)

        if (tp %in% c("continuous", "count", "ordinal")) {
          truth_latent <- X_truth[test_rows, lc[1]]
          pred_latent  <- pred$imputed_latent[test_rows, lc[1]]

          # BM baseline latent
          bm_latent <- baseline$mu
          if (isTRUE(fit_obj$multi_obs)) {
            bm_latent <- bm_latent[fit_obj$obs_to_species, , drop = FALSE]
          }
          bm_pred <- bm_latent[test_rows, lc[1]]

          entry$rmse_gnn <- rmse_vec(truth_latent, pred_latent)
          entry$rmse_bm  <- rmse_vec(truth_latent, bm_pred)
          entry$r_gnn    <- tryCatch(stats::cor(truth_latent, pred_latent,
                                                 use = "complete.obs"),
                                      error = function(e) NA_real_)
          entry$r_bm     <- tryCatch(stats::cor(truth_latent, bm_pred,
                                                 use = "complete.obs"),
                                      error = function(e) NA_real_)

          # Conformal coverage
          if (!is.null(pred$conformal_lower) && nm %in% colnames(pred$conformal_lower)) {
            low  <- pred$conformal_lower[test_rows, nm]
            high <- pred$conformal_upper[test_rows, nm]
            # Need truth in original scale
            if (tp == "continuous") {
              truth_orig <- truth_latent * tm$sd + tm$mean
              if (isTRUE(tm$log_transform)) truth_orig <- exp(truth_orig)
            } else {
              truth_orig <- truth_latent * tm$sd + tm$mean
            }
            valid <- !is.na(low) & !is.na(high) & !is.na(truth_orig)
            if (sum(valid) > 0) {
              entry$conformal_coverage <- mean(
                truth_orig[valid] >= low[valid] & truth_orig[valid] <= high[valid]
              )
            }
          }

        } else if (tp == "binary") {
          truth_val <- X_truth[test_rows, lc[1]]
          pred_prob <- pred$probabilities[[nm]][test_rows]
          if (!is.null(pred_prob)) {
            pred_class <- ifelse(pred_prob >= 0.5, 1, 0)
            entry$accuracy_gnn <- mean(pred_class == truth_val, na.rm = TRUE)
            entry$brier_gnn    <- mean((pred_prob - truth_val)^2, na.rm = TRUE)
          }
          # BM baseline accuracy
          bm_logit <- baseline$mu[, lc[1]]
          if (isTRUE(fit_obj$multi_obs)) {
            bm_logit <- bm_logit[fit_obj$obs_to_species]
          }
          bm_prob  <- 1 / (1 + exp(-bm_logit[test_rows]))
          bm_class <- ifelse(bm_prob >= 0.5, 1, 0)
          entry$accuracy_bm <- mean(bm_class == truth_val, na.rm = TRUE)

        } else if (tp == "categorical") {
          truth_oh <- X_truth[test_rows, lc, drop = FALSE]
          truth_class <- apply(truth_oh, 1, which.max)
          pred_probs  <- pred$probabilities[[nm]]
          if (!is.null(pred_probs) && is.matrix(pred_probs)) {
            pred_class <- apply(pred_probs[test_rows, , drop = FALSE], 1, which.max)
            entry$accuracy_gnn <- mean(pred_class == truth_class, na.rm = TRUE)
          }
          # BM baseline
          bm_logprobs <- baseline$mu[, lc, drop = FALSE]
          if (isTRUE(fit_obj$multi_obs)) {
            bm_logprobs <- bm_logprobs[fit_obj$obs_to_species, , drop = FALSE]
          }
          bm_class <- apply(bm_logprobs[test_rows, , drop = FALSE], 1, which.max)
          entry$accuracy_bm <- mean(bm_class == truth_class, na.rm = TRUE)
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

  # Build metrics table rows
  trait_rows <- ""
  for (m in metrics) {
    nm <- m$name
    tp <- m$type

    gate_val <- if (!is.null(m$gate)) sprintf("%.3f", m$gate) else "&mdash;"

    if (tp %in% c("continuous", "count", "ordinal")) {
      bm_rmse  <- if (!is.null(m$rmse_bm))  sprintf("%.4f", m$rmse_bm)  else "&mdash;"
      gnn_rmse <- if (!is.null(m$rmse_gnn)) sprintf("%.4f", m$rmse_gnn) else "&mdash;"
      bm_r     <- if (!is.null(m$r_bm))     sprintf("%.3f", m$r_bm)     else "&mdash;"
      gnn_r    <- if (!is.null(m$r_gnn))    sprintf("%.3f", m$r_gnn)    else "&mdash;"
      n_test   <- if (!is.null(m$n_test))   m$n_test                     else "&mdash;"

      pct <- ""
      verdict <- '<span class="pill pill-neutral">no data</span>'
      if (!is.null(m$rmse_bm) && !is.null(m$rmse_gnn) && m$rmse_bm > 0) {
        pct_val <- (m$rmse_bm - m$rmse_gnn) / m$rmse_bm * 100
        if (pct_val > 0.5) {
          pct <- sprintf('<span class="pos">+%.1f%%</span>', pct_val)
          verdict <- '<span class="pill pill-green">GNN better</span>'
        } else if (pct_val < -0.5) {
          pct <- sprintf('<span class="neg">%.1f%%</span>', pct_val)
          verdict <- '<span class="pill pill-red">BM better</span>'
        } else {
          pct <- sprintf('<span class="zero">%.1f%%</span>', pct_val)
          verdict <- '<span class="pill pill-neutral">tied</span>'
        }
      }

      cov <- "&mdash;"
      if (!is.null(m$conformal_coverage)) {
        cov_val <- m$conformal_coverage * 100
        cov_class <- if (cov_val >= 93) "pos" else if (cov_val >= 85) "zero" else "neg"
        cov <- sprintf('<span class="%s">%.1f%%</span>', cov_class, cov_val)
      }

      trait_rows <- paste0(trait_rows, sprintf(
        '<tr><td>%s</td><td>%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td>%s</td></tr>\n',
        nm, tp, n_test, bm_rmse, gnn_rmse, pct, bm_r, gnn_r, cov, verdict
      ))

    } else {
      bm_acc  <- if (!is.null(m$accuracy_bm))  sprintf("%.1f%%", m$accuracy_bm * 100)  else "&mdash;"
      gnn_acc <- if (!is.null(m$accuracy_gnn)) sprintf("%.1f%%", m$accuracy_gnn * 100) else "&mdash;"
      n_test  <- if (!is.null(m$n_test))       m$n_test                                 else "&mdash;"

      pct <- ""
      verdict <- '<span class="pill pill-neutral">no data</span>'
      if (!is.null(m$accuracy_bm) && !is.null(m$accuracy_gnn)) {
        diff <- (m$accuracy_gnn - m$accuracy_bm) * 100
        if (diff > 0.5) {
          pct <- sprintf('<span class="pos">+%.1f pp</span>', diff)
          verdict <- '<span class="pill pill-green">GNN better</span>'
        } else if (diff < -0.5) {
          pct <- sprintf('<span class="neg">%.1f pp</span>', diff)
          verdict <- '<span class="pill pill-red">BM better</span>'
        } else {
          pct <- '<span class="zero">0.0 pp</span>'
          verdict <- '<span class="pill pill-neutral">tied</span>'
        }
      }

      trait_rows <- paste0(trait_rows, sprintf(
        '<tr><td>%s</td><td>%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td class="num" colspan="2">&mdash;</td><td class="num">&mdash;</td><td>%s</td></tr>\n',
        nm, tp, n_test, bm_acc, gnn_acc, pct, verdict
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

  gate_labels_js <- paste0("[", paste0("'", gate_names, "'", collapse = ","), "]")
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

  # Features list
  features <- c(
    if (has_attention)   "Attention message passing" else "Standard message passing",
    if (has_calibration) "Validation-calibrated gates",
    if (has_conformal)   "Conformal prediction intervals",
    sprintf("Adaptive spectral encoding (k=%d)", k_eigen),
    "Phylogenetic label propagation (discrete baselines)"
  )
  features_html <- paste0("<li>", features, "</li>", collapse = "\n")

  timestamp <- format(Sys.time(), "%d %B %Y, %H:%M")

  # Assemble HTML
  sprintf('<!DOCTYPE html>
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
<table>
<thead><tr><th>Trait</th><th>Type</th><th class="num">n test</th><th class="num">BM</th><th class="num">GNN</th><th class="num">%%&Delta;</th><th class="num">r (BM)</th><th class="num">r (GNN)</th><th class="num">Coverage</th><th>Verdict</th></tr></thead>
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
}
