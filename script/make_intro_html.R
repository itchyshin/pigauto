#!/usr/bin/env Rscript
# Build a self-contained HTML introduction to pigauto, drawn from README.md
# and the architecture summary in CLAUDE.md. Writes to inst/doc/ so it ships
# with the installed package and can be found via:
#   system.file("doc", "pigauto_intro.html", package = "pigauto")

out <- "inst/doc/pigauto_intro.html"
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

# Pull version from DESCRIPTION so the badge stays in sync with the source of
# truth without re-parsing at load time.
desc    <- readLines("DESCRIPTION")
ver_line <- grep("^Version:", desc, value = TRUE)
version <- if (length(ver_line)) sub("^Version:\\s*", "", ver_line) else "?"

timestamp <- format(Sys.time(), "%Y-%m-%d")

html <- paste0(
'<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>pigauto &mdash; introduction</title>
<style>
  :root { --fg:#111827; --muted:#6b7280; --line:#e5e7eb; --soft:#f3f4f6;
          --accent:#059669; --code-bg:#f3f4f6; --hi:#eef2ff; --hi-fg:#3730a3; }
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         max-width: 880px; margin: 2em auto; padding: 0 1.5em; color: var(--fg);
         line-height: 1.6; }
  h1 { border-bottom: 2px solid var(--fg); padding-bottom: .25em; margin-bottom: .2em; }
  h2 { margin-top: 2em; color: #1f2937; border-bottom: 1px solid var(--line);
       padding-bottom: .15em; }
  h3 { color: #111827; }
  .tagline { color: var(--muted); font-size: 16px; margin-top: 0; }
  .badge { display: inline-block; background: var(--hi); color: var(--hi-fg);
           padding: 2px 10px; border-radius: 999px; font-size: 13px;
           font-weight: 600; margin-left: 8px; vertical-align: middle; }
  pre { background: var(--code-bg); padding: 1em 1.2em; border-radius: 6px;
        overflow-x: auto; font-size: 13px; line-height: 1.5;
        border-left: 3px solid var(--accent); }
  code { background: var(--code-bg); padding: 1px 6px; border-radius: 3px;
         font-size: 13px; }
  pre code { background: transparent; padding: 0; }
  table { border-collapse: collapse; width: 100%; margin: 1em 0; font-size: 14px; }
  th, td { padding: 8px 12px; border-bottom: 1px solid var(--line);
           text-align: left; vertical-align: top; }
  th { background: var(--soft); font-weight: 600; }
  .note { background: #fef3c7; border-left: 4px solid #d97706; padding: 1em 1.2em;
          border-radius: 4px; margin: 1.5em 0; }
  .note b { color: #92400e; }
  ul li { margin: 6px 0; }
  .meta { color: var(--muted); font-size: 13px; }
  footer { color: var(--muted); font-size: 12px; margin-top: 3em;
           border-top: 1px solid var(--line); padding-top: 1em; }
</style>
</head>
<body>

<h1>pigauto <span class="badge">v', version, '</span></h1>
<p class="tagline">Phylogenetic trait imputation via graph neural network.</p>
<p class="meta">Generated ', timestamp, ' from <code>README.md</code> and <code>CLAUDE.md</code>.</p>

<h2>What it does</h2>
<p>
<code>pigauto</code> fits a Residual Phylogenetic Denoising Autoencoder
(ResidualPhyloDAE) that combines a Brownian-motion phylogenetic baseline
(via <code>Rphylopars</code>) with a residual graph autoencoder. The GNN
learns corrections from tree topology and inter-trait correlations.
Continuous, binary, categorical, ordinal, and count traits all live in a
single unified latent space.
</p>
<ul>
<li><b>Mixed-type traits</b> &mdash; continuous, binary, categorical, ordinal, and count in one model.</li>
<li><b>Attention-based message passing</b> with a learnable phylogenetic adjacency bias &mdash; the model attends to informative neighbours rather than treating all species equally.</li>
<li><b>Validation-calibrated gates</b> &mdash; post-training gate optimisation prevents the GNN from adding noise when the baseline is already strong.</li>
<li><b>Conformal prediction intervals</b> &mdash; distribution-free 95% coverage guarantees for continuous/count/ordinal traits.</li>
<li><b>Adaptive spectral encoding</b> &mdash; <code>k_eigen</code> scales automatically with tree size (4 for ~30 tips up to 32 for 640+).</li>
<li><b>Multiple observations per species</b> &mdash; aggregates individual-level data before phylogenetic message passing.</li>
</ul>

<h2>Installation</h2>
<pre><code># Install from local source
devtools::install()

# First-time torch setup (required)
torch::install_torch()</code></pre>

<h2>Quick start</h2>
<p>One call does preprocessing, baseline fitting, GNN training, calibration, and prediction:</p>
<pre><code>library(pigauto)
data(avonet300, tree300)

# Set species as rownames
df &lt;- avonet300
rownames(df) &lt;- df$Species_Key
df$Species_Key &lt;- NULL

# One-call imputation (attention + calibration + conformal by default)
result &lt;- impute(df, tree300)

# Access results
result$completed                             # user df with NAs filled in
result$prediction$imputed$Mass               # imputed values (one trait)
result$prediction$se[, "Mass"]               # standard errors
result$prediction$conformal_lower[, "Mass"]  # 95% CI lower
result$prediction$conformal_upper[, "Mass"]  # 95% CI upper
result$prediction$probabilities$Trophic.Level  # categorical probabilities

# Auto-generate an HTML report
pigauto_report(result)</code></pre>

<h2>Trait types supported</h2>
<p>All five types coexist in one model. pigauto auto-detects from the R column class; you can override with <code>trait_types</code>.</p>
<table>
<thead><tr>
<th>R class</th>
<th>pigauto type</th>
<th>Encoding</th>
<th>Loss</th>
<th>Baseline</th>
</tr></thead>
<tbody>
<tr><td><code>numeric</code></td>   <td>continuous</td><td>optional log + z-score</td><td>MSE</td>          <td>Rphylopars BM</td></tr>
<tr><td><code>integer</code></td>   <td>count</td>     <td><code>log1p</code> + z-score</td><td>MSE</td> <td>Rphylopars BM</td></tr>
<tr><td><code>factor(2)</code></td> <td>binary</td>    <td>0/1</td>              <td>BCE</td>          <td>Phylo label propagation</td></tr>
<tr><td><code>factor(&gt;2)</code></td><td>categorical</td><td>one-hot (K cols)</td><td>cross-entropy</td><td>Phylo label propagation</td></tr>
<tr><td><code>ordered</code></td>   <td>ordinal</td>   <td>integer + z-score</td><td>MSE</td>          <td>Rphylopars BM</td></tr>
</tbody>
</table>

<h2>Architecture in one paragraph</h2>
<p>
Every prediction is
</p>
<pre><code>pred = (1 - r_cal) * BM_baseline + r_cal * GNN_delta</code></pre>
<p>
where <code>r_cal</code> is a per-trait gate learned during training and
<b>re-calibrated on the validation set</b> afterward. When the Brownian-motion
baseline is already optimal, the calibrated gate closes to zero and the GNN
becomes a no-op &mdash; so pigauto is guaranteed never to underperform the
baseline on held-out data. The GNN uses two message-passing layers with
scaled dot-product attention and a learnable log-adjacency bias, so the
model starts close to the phylogenetic prior and learns deviations only
where the data supports them.
</p>

<h2>Pipeline functions for fine-grained control</h2>
<p>If you want to swap a component (custom baseline, different split, cross-validation), call the pipeline directly:</p>
<pre><code>pd       &lt;- preprocess_traits(df, tree300, log_transform = TRUE)
splits   &lt;- make_missing_splits(pd$X_scaled, missing_frac = 0.25, seed = 42)
graph    &lt;- build_phylo_graph(tree300)  # adaptive k_eigen
baseline &lt;- fit_baseline(pd, tree300, splits = splits)
fit      &lt;- fit_pigauto(pd, tree300, splits = splits, graph = graph,
                        baseline = baseline, epochs = 2000L)
pred     &lt;- predict(fit, return_se = TRUE, n_imputations = 10L)

# Evaluate on held-out test set
eval_df  &lt;- evaluate(fit, data = pd, splits = splits)

# k-fold cross-validation
cv       &lt;- cross_validate(pd, tree300, k = 5, seeds = 1:3)</code></pre>

<h2>Multiple observations per species</h2>
<p>When your data has multiple rows per species (individual-level measurements, multi-specimen museum records), use <code>species_col</code>:</p>
<pre><code># traits_df has columns: species, mass, wing_length, ...
# Multiple rows per species are allowed
result &lt;- impute(traits_df, tree, species_col = "species")

# Or via the lower-level API
pd &lt;- preprocess_traits(traits_df, tree, species_col = "species")
pd$n_obs      # number of observations
pd$n_species  # number of unique species</code></pre>
<p>
The GNN aggregates observations to species level before phylogenetic
message passing, then broadcasts back to observation level. Rphylopars
handles within-species replication natively in the baseline.
</p>

<h2>What&rsquo;s in the box</h2>
<ul>
<li><code>avonet300</code> &mdash; 300 bird species with 4 continuous morphometric traits
  (Mass, Beak.Length_Culmen, Tarsus.Length, Wing.Length) plus 3 ecological traits
  (Trophic.Level [categorical], Primary.Lifestyle [categorical], Migration [ordinal]).</li>
<li><code>tree300</code> &mdash; pruned BirdTree phylogeny matching <code>avonet300</code>.</li>
</ul>

<div class="note">
<b>When should I use pigauto vs alternatives?</b>
Prefer pigauto for mid-sized phylogenetic datasets (50&ndash;5000 species) where
you want a single calibrated model covering all trait types and do not
want to wait hours for MCMC. For formal Bayesian multiple imputation with
posterior pooling, consider the BACE package. For small continuous-only
datasets where interpretability of a specific evolutionary model matters,
Rphylopars is the right tool.
</div>

<h2>Further reading</h2>
<ul>
<li><code>README.md</code> &mdash; full user guide, benchmarks, and plotting examples.</li>
<li><code>CLAUDE.md</code> &mdash; architecture reference for contributors.</li>
<li><code>vignettes/</code> &mdash; getting-started and mixed-types vignettes.</li>
<li>Citation: Nakagawa S (2026). <i>pigauto: Phylogenetic Imputation via Graph Autoencoder.</i> R package version ', version, '.</li>
</ul>

<footer>
Generator: <code>script/make_intro_html.R</code> &middot;
Installed path: <code>system.file("doc", "pigauto_intro.html", package = "pigauto")</code>
</footer>

</body>
</html>
')

writeLines(html, out)
cat("Wrote ", out, "\n", sep = "")
