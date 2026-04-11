#!/usr/bin/env Rscript
#
# script/make_walkthrough_multi_obs_html.R
#
# Build a self-contained HTML walkthrough:
#   "Multi-observation per species with observation-level covariates"
#
# This walkthrough demonstrates pigauto's multi-obs + covariate support
# using the bundled ctmax_sim dataset. It runs a small imputation inline
# (no precompute needed) on 300 species.
#
# Writes two identical copies:
#   inst/doc/pigauto_walkthrough_multi_obs.html
#   pkgdown/assets/pigauto_walkthrough_multi_obs.html

suppressPackageStartupMessages({
  library(ape)
  devtools::load_all(
    "/Users/z3437171/Dropbox/Github Local/pigauto",
    quiet = TRUE
  )
})

timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M")

# ---- Run a quick imputation for display ------------------------------------
cat("Running imputation on ctmax_sim...\n")

data(ctmax_sim, tree300, package = "pigauto")

# Truth (stored as attribute by make_ctmax_sim.R)
ctmax_truth <- attr(ctmax_sim, "ctmax_truth")

# Prepare data
traits <- ctmax_sim[, c("species", "CTmax")]
covs   <- data.frame(acclim_temp = ctmax_sim$acclim_temp)

# Impute without covariates
t0 <- proc.time()["elapsed"]
res_no_cov <- tryCatch(
  impute(traits, tree300, species_col = "species",
         epochs = 200L, verbose = FALSE, seed = 42L,
         eval_every = 25L, patience = 15L),
  error = function(e) NULL
)
t1 <- proc.time()["elapsed"]
time_no_cov <- round(t1 - t0, 1)

# Impute with covariates
t0 <- proc.time()["elapsed"]
res_cov <- tryCatch(
  impute(traits, tree300, species_col = "species",
         covariates = covs,
         epochs = 200L, verbose = FALSE, seed = 42L,
         eval_every = 25L, patience = 15L),
  error = function(e) NULL
)
t1 <- proc.time()["elapsed"]
time_cov <- round(t1 - t0, 1)

# ---- Compute metrics -------------------------------------------------------
missing_mask <- is.na(ctmax_sim$CTmax)

compute_metrics <- function(result, label) {
  if (is.null(result)) return(data.frame(method = label, obs_rmse = NA,
                                          pearson_r = NA, n_imputed = NA))
  preds <- result$prediction$imputed$CTmax
  ok <- missing_mask & is.finite(preds) & is.finite(ctmax_truth)
  truth_held <- ctmax_truth[ok]
  preds_held <- preds[ok]
  data.frame(
    method    = label,
    obs_rmse  = round(sqrt(mean((truth_held - preds_held)^2)), 3),
    pearson_r = round(cor(truth_held, preds_held), 3),
    n_imputed = sum(ok)
  )
}

metrics <- rbind(
  compute_metrics(res_no_cov, "pigauto (no covariates)"),
  compute_metrics(res_cov,    "pigauto (with covariates)")
)

# Check within-species variation (only for the covariates model)
within_sp_varies <- FALSE
if (!is.null(res_cov)) {
  preds_cov <- res_cov$prediction$imputed$CTmax
  sp_list <- unique(ctmax_sim$species)
  n_varies <- 0L
  for (sp in sp_list[1:50]) {
    rows <- which(ctmax_sim$species == sp)
    if (length(rows) > 1) {
      rng <- max(preds_cov[rows]) - min(preds_cov[rows])
      if (rng > 0.001) n_varies <- n_varies + 1L
    }
  }
  within_sp_varies <- n_varies > 0L
}

cat(sprintf("Metrics computed. Within-species variation: %s\n",
            if (within_sp_varies) "YES" else "NO"))
print(metrics)

# ---- Data summary ----------------------------------------------------------
n_total   <- nrow(ctmax_sim)
n_species <- length(unique(ctmax_sim$species))
n_missing <- sum(missing_mask)
n_unobs   <- length(setdiff(tree300$tip.label,
                             unique(ctmax_sim$species[!missing_mask])))
obs_per_sp <- table(ctmax_sim$species)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "&ndash;", formatC(x, digits = digits, format = "f"))
}

# ---------------------------------------------------------------------------
# Build HTML
# ---------------------------------------------------------------------------

html <- character()
h <- function(...) html <<- c(html, paste0(...))

h('<!DOCTYPE html>')
h('<html lang="en">')
h('<head>')
h('<meta charset="utf-8">')
h('<meta name="viewport" content="width=device-width, initial-scale=1">')
h('<title>Multi-observation per species &mdash; pigauto walkthrough</title>')
h('<style>')
h('
  :root { --primary: #059669; --bg: #f8fafb; --text: #1e293b; }
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body { font-family: "Inter", system-ui, sans-serif; color: var(--text);
         background: var(--bg); line-height: 1.65; padding: 2rem 1rem; }
  .container { max-width: 56rem; margin: 0 auto; }
  h1 { font-size: 1.85rem; color: var(--primary); margin-bottom: .3rem; }
  h2 { font-size: 1.35rem; color: var(--primary); margin: 2.2rem 0 .8rem;
       border-bottom: 2px solid #d1fae5; padding-bottom: .25rem; }
  h3 { font-size: 1.1rem; color: #374151; margin: 1.2rem 0 .4rem; }
  p, li { margin-bottom: .75rem; }
  ul, ol { padding-left: 1.4rem; }
  code { font-family: "JetBrains Mono", monospace; font-size: 0.88em;
         background: #e2e8f0; padding: 0.15em 0.4em; border-radius: 3px; }
  pre { background: #1e293b; color: #e2e8f0; padding: 1rem 1.2rem;
        border-radius: 6px; overflow-x: auto; margin: .8rem 0 1.2rem;
        font-size: 0.85rem; line-height: 1.5; }
  pre code { background: none; padding: 0; color: inherit; }
  table { border-collapse: collapse; width: 100%; margin: 1rem 0; font-size: 0.9rem; }
  th { background: var(--primary); color: #fff; padding: .55rem .7rem; text-align: left; }
  td { padding: .45rem .7rem; border-bottom: 1px solid #e2e8f0; }
  tr:nth-child(even) { background: #f1f5f9; }
  .subtitle { color: #64748b; font-size: 0.95rem; margin-bottom: 1.5rem; }
  .note { background: #ecfdf5; border-left: 4px solid var(--primary);
          padding: .8rem 1rem; margin: 1rem 0; border-radius: 0 4px 4px 0;
          font-size: 0.92rem; }
  .warning { background: #fef9c3; border-left: 4px solid #eab308;
             padding: .8rem 1rem; margin: 1rem 0; border-radius: 0 4px 4px 0;
             font-size: 0.92rem; }
  .footer { margin-top: 3rem; padding-top: 1rem; border-top: 1px solid #e2e8f0;
            font-size: 0.82rem; color: #94a3b8; }
  .highlight { color: var(--primary); font-weight: 600; }
  .diagram { background: #f1f5f9; padding: 1rem; border-radius: 6px;
             margin: 1rem 0; font-family: "JetBrains Mono", monospace;
             font-size: 0.82rem; line-height: 1.7; white-space: pre; }
')
h('</style>')
h('</head>')
h('<body>')
h('<div class="container">')

# ---- Title -----------------------------------------------------------------
h('<h1>Multi-observation per species</h1>')
h('<p class="subtitle">pigauto v0.6.0 &mdash; observation-level covariate',
  ' imputation &mdash; generated ', timestamp, '</p>')

# ---- Section 1: The problem ------------------------------------------------
h('<h2>1. The problem</h2>')
h('<p>Comparative datasets increasingly have <strong>multiple data points per',
  ' species</strong> measured under different conditions. Examples include:',
  ' critical thermal maximum (CTmax) at different acclimation temperatures',
  ' (Pottier et al. 2025, <em>Nature</em>), metabolic rate at different',
  ' body temperatures, or growth rate across populations.</p>')
h('<p>Missing data is ubiquitous in these datasets &mdash; many species are',
  ' entirely unobserved, and even measured species have incomplete coverage.',
  ' The imputation challenge is two-fold:</p>')
h('<ol>')
h('<li>Species-level: use the phylogenetic tree to predict trait values',
  ' for unobserved species (the classical comparative imputation problem).</li>')
h('<li>Observation-level: use <em>observation-level covariates</em>',
  ' (e.g., acclimation temperature) to produce covariate-conditional',
  ' predictions within species. A species measured at 20&deg;C and 30&deg;C',
  ' acclimation should get different imputed CTmax values.</li>')
h('</ol>')
h('<p>pigauto (v0.6.0+) handles both levels through its',
  ' <strong>observation-level refinement</strong> architecture.</p>')

# ---- Section 2: The data ---------------------------------------------------
h('<h2>2. The bundled dataset: <code>ctmax_sim</code></h2>')
h('<p>pigauto ships a simulated multi-observation CTmax dataset for',
  ' tutorials and testing. The data-generating process is:</p>')
h('<pre><code>CTmax_ij = 38 + phylo_i + 0.50 * acclim_temp_ij + epsilon_ij</code></pre>')
h('<p>where <code>phylo_i</code> follows Brownian motion on',
  ' <code>tree300</code> and <code>epsilon ~ N(0, 1.5)</code>.</p>')

h('<table>')
h('<tr><th>Property</th><th>Value</th></tr>')
h(sprintf('<tr><td>Species</td><td>%d</td></tr>', n_species))
h(sprintf('<tr><td>Total observations</td><td>%s</td></tr>',
          format(n_total, big.mark = ",")))
h(sprintf('<tr><td>Obs per species</td><td>%.1f (range %d&ndash;%d)</td></tr>',
          mean(obs_per_sp), min(obs_per_sp), max(obs_per_sp)))
h(sprintf('<tr><td>Missing CTmax cells</td><td>%d (%.1f%%)</td></tr>',
          n_missing, 100 * n_missing / n_total))
h(sprintf('<tr><td>Entirely unobserved species</td><td>%d (%.0f%%)</td></tr>',
          n_unobs, 100 * n_unobs / n_species))
h('<tr><td>Observation-level covariate</td><td><code>acclim_temp</code></td></tr>')
h('</table>')

h('<pre><code>data(ctmax_sim, tree300)
head(ctmax_sim)
#>   species acclim_temp CTmax
#> 1     t28        22.1 39.82
#> 2     t28        17.5    NA
#> 3     t28        30.8 41.03
#> ...</code></pre>')

# ---- Section 3: Imputation -------------------------------------------------
h('<h2>3. Imputation: with and without covariates</h2>')

h('<h3>Without covariates</h3>')
h('<pre><code>traits &lt;- ctmax_sim[, c("species", "CTmax")]

result_no_cov &lt;- impute(traits, tree300, species_col = "species")</code></pre>')
h(sprintf('<p>Runtime: %.1f seconds (200 epochs, 300 species, %s observations).</p>',
          time_no_cov, format(n_total, big.mark = ",")))

h('<h3>With observation-level covariates</h3>')
h('<pre><code>covs &lt;- data.frame(acclim_temp = ctmax_sim$acclim_temp)

result_cov &lt;- impute(traits, tree300, species_col = "species",
                     covariates = covs)</code></pre>')
h(sprintf('<p>Runtime: %.1f seconds.</p>', time_cov))

h('<h3>RMSE comparison on held-out cells</h3>')
h('<table>')
h('<tr><th>Method</th><th>Obs-level RMSE</th><th>Pearson r</th><th>N imputed</th></tr>')
for (i in seq_len(nrow(metrics))) {
  bold <- if (i == which.min(metrics$obs_rmse)) ' style="font-weight:600;"' else ''
  h(sprintf('<tr%s><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>',
            bold,
            metrics$method[i],
            fmt(metrics$obs_rmse[i]),
            fmt(metrics$pearson_r[i]),
            format(metrics$n_imputed[i], big.mark = ",")))
}
h('</table>')

# ---- Section 4: How it works -----------------------------------------------
h('<h2>4. Architecture: observation-level refinement</h2>')
h('<p>The key challenge is that the GNN operates at <strong>species level</strong>:',
  ' observations are aggregated via <code>scatter_mean</code> before phylogenetic',
  ' message passing, then broadcast back. Without the refinement layer,',
  ' all observations of the same species receive the identical prediction.</p>')

h('<div class="diagram">')
h('Observations (n_obs x p)')
h('  |')
h('  v')
h('Encoder: concat(X, spectral_coords, covariates) -> h    [obs-level]')
h('  |')
h('  v')
h('scatter_mean(h, species_index) -> h_species              [AGGREGATION]')
h('  |')
h('  v')
h('GNN message passing (2 layers) -> h_species              [species-level]')
h('  |')
h('  v')
h('broadcast(h_species, species_index) -> h_obs             [species -> obs]')
h('  |')
h('  v')
h('<span class="highlight">h_obs + obs_refine(cat(h_obs, user_covs)) -> h_refined</span>    [NEW]')
h('  |')
h('  v')
h('Decoder -> delta (obs-specific predictions)')
h('</div>')

h('<p>The <code>obs_refine</code> MLP (hidden_dim + n_cov &rarr; hidden_dim',
  ' &rarr; hidden_dim) uses a <strong>residual connection</strong>: it',
  ' adjusts the species-level representation rather than replacing it.',
  ' The phylogenetic signal from message passing is preserved; the',
  ' covariate adjustment is additive.</p>')

h('<div class="note">')
h('<strong>Backward compatible.</strong> When no covariates are supplied or',
  ' in single-obs mode, no <code>obs_refine</code> MLP is created. The model',
  ' architecture is identical to v0.5.0 in those cases.')
h('</div>')

if (within_sp_varies) {
  h('<div class="note">')
  h('<strong>Verification.</strong> In this run, observations of the same',
    ' species with different acclimation temperatures received',
    ' <em>different</em> predicted CTmax values &mdash; confirming that',
    ' the refinement layer produces covariate-conditional predictions.')
  h('</div>')
}

# ---- Section 5: Usage pattern -----------------------------------------------
h('<h2>5. Typical usage pattern</h2>')
h('<pre><code>library(pigauto)
library(glmmTMB)
library(ape)

data(ctmax_sim, tree300)

# 1. Separate traits from covariates
traits &lt;- ctmax_sim[, c("species", "CTmax")]
covs   &lt;- data.frame(acclim_temp = ctmax_sim$acclim_temp)

# 2. Multiple imputation (M = 50)
mi &lt;- multi_impute(traits, tree300, m = 50L,
                   species_col = "species",
                   covariates = covs)

# 3. Downstream regression on each imputed dataset
Vphy &lt;- cov2cor(ape::vcv(tree300))
fits &lt;- with_imputations(mi, function(d) {
  d$acclim   &lt;- covs$acclim_temp
  d$sp       &lt;- factor(d$species, levels = rownames(Vphy))
  d$dummy    &lt;- factor(1)
  glmmTMB(
    CTmax ~ acclim +
      propto(0 + sp | dummy, Vphy),
    data = d
  )
})

# 4. Pool with Rubin\'s rules
pool_mi(fits,
        coef_fun = function(f) fixef(f)$cond,
        vcov_fun = function(f) vcov(f)$cond)</code></pre>')

# ---- Section 6: Comparison with BACE ----------------------------------------
h('<h2>6. When to use pigauto vs BACE</h2>')
h('<table>')
h('<tr><th>Feature</th><th>pigauto</th><th>BACE</th></tr>')
h('<tr><td>Paradigm</td><td>Frequentist (GNN + Rubin\'s rules)</td>')
h('    <td>Bayesian (MCMCglmm)</td></tr>')
h('<tr><td>Speed (1,000 species)</td><td>Minutes</td><td>Hours</td></tr>')
h('<tr><td>Speed (5,000+ species)</td><td>5&ndash;15 min</td><td>Days</td></tr>')
h('<tr><td>Trait types</td><td>All 7 types in one call</td>')
h('    <td>Continuous + binary (Gaussian/ordinal family)</td></tr>')
h('<tr><td>Obs-level covariates</td><td>Yes (v0.6.0+, refinement MLP)</td>')
h('    <td>Yes (MCMCglmm fixed effects)</td></tr>')
h('<tr><td>Phylo uncertainty</td><td><code>multi_impute_trees()</code></td>')
h('    <td>Via posterior tree sampling</td></tr>')
h('<tr><td>Downstream inference</td><td>Rubin\'s rules via <code>pool_mi()</code></td>')
h('    <td>Posterior concatenation</td></tr>')
h('</table>')

h('<p><strong>Recommendation:</strong> Use pigauto when you need speed,',
  ' mixed trait types, or datasets with &gt;1,000 species. Use BACE when',
  ' you want a fully Bayesian workflow or need integrated uncertainty',
  ' quantification from the imputation model itself.</p>')

# ---- Section 7: Summary -----------------------------------------------------
h('<h2>7. Summary</h2>')
h('<ul>')
h('<li>Comparative datasets with multiple observations per species are',
  ' increasingly common (AmphiTherm, AnAge, metabolic databases).</li>')
h('<li>pigauto\'s <code>species_col</code> parameter handles multi-obs data;',
  ' adding <code>covariates</code> enables observation-level predictions.</li>')
h('<li>The <code>obs_refine</code> MLP re-injects covariates after phylogenetic',
  ' message passing via a residual connection, preserving between-species',
  ' phylogenetic signal while enabling within-species covariate-conditional',
  ' predictions.</li>')
h('<li>The full MI workflow (<code>multi_impute()</code> &rarr;',
  ' <code>with_imputations()</code> &rarr; <code>pool_mi()</code>) works',
  ' seamlessly with multi-obs data.</li>')
h('</ul>')

# ---- Footer -----------------------------------------------------------------
h('<div class="footer">')
h(sprintf('<p>Generated by <code>script/make_walkthrough_multi_obs_html.R</code> on %s.<br>',
          timestamp))
h('pigauto v0.6.0 &middot; Nakagawa S (2026).<br>')
h('Data: bundled <code>ctmax_sim</code> + <code>tree300</code>.</p>')
h('</div>')

h('</div>')  # container
h('</body>')
h('</html>')

# ---- Write HTML --------------------------------------------------------------
here <- "/Users/z3437171/Dropbox/Github Local/pigauto"

out1 <- file.path(here, "inst", "doc", "pigauto_walkthrough_multi_obs.html")
out2 <- file.path(here, "pkgdown", "assets", "pigauto_walkthrough_multi_obs.html")

dir.create(dirname(out1), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out2), recursive = TRUE, showWarnings = FALSE)

writeLines(html, out1)
writeLines(html, out2)

cat(sprintf("Wrote %s (%s bytes)\n", out1, file.size(out1)))
cat(sprintf("Wrote %s (%s bytes)\n", out2, file.size(out2)))
cat("Done.\n")
