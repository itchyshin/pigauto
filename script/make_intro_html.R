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
<p class="tagline">Fill in missing species traits using a phylogenetic tree.</p>
<p class="meta">Generated ', timestamp, ' from <code>README.md</code> and <code>CLAUDE.md</code>.</p>

<h2>What it does</h2>
<p>
Comparative analyses often fail because trait databases are incomplete.
<code>pigauto</code> fills in the gaps by combining three sources of
information: <b>(1)</b> the phylogenetic tree, which tells us that
closely related species tend to share similar traits, <b>(2)</b>
correlations among the traits themselves, which let observed traits
inform predictions of missing ones, and <b>(3)</b> optional
environmental covariates (climate, habitat, geography), which capture
trait variation that phylogeny alone cannot explain. The package
handles continuous measurements, counts, binary variables, ordered
categories, and unordered categories &mdash; all in a single call.
</p>
<h3>How it works</h3>
<p>
Under the hood, pigauto blends a <b>phylogenetic baseline</b>
(Brownian-motion imputation for continuous traits, phylogenetic label
propagation for discrete traits) with a <b>graph neural network</b>
that learns additional patterns from tree structure and cross-trait
correlations. A per-trait gate calibrated on held-out data controls how
much the GNN contributes: when the baseline is already good enough, the
gate closes and the network stays out of the way. The final prediction
is the blend
<code>(1 &minus; r) &times; baseline + r &times; GNN</code>.
</p>
<div class="note">
<b>Technical note.</b> The internal torch module is named
<code>ResidualPhyloDAE</code>; the &ldquo;Residual&rdquo; refers to
ResNet-style skip connections inside the GNN layers, not to a
statistical residual. The GNN&rsquo;s output is a full prediction
trained end-to-end, not <code>y &minus; baseline</code>.
</div>
<ul>
<li><b>Multiple trait types</b> &mdash; continuous, binary, categorical, ordinal, count, and proportion in one model.</li>
<li><b>Uses the phylogeny</b> &mdash; closely related species inform predictions, as you would expect from a comparative method.</li>
<li><b>Cross-trait patterns</b> &mdash; if body mass predicts beak length, observed masses help impute missing beak lengths.</li>
<li><b>Safe by default</b> &mdash; the per-trait gate prevents the neural network from degrading traits the baseline already handles well.</li>
<li><b>Uncertainty quantification</b> &mdash; conformal prediction intervals give 95% coverage for continuous/count/ordinal traits.</li>
<li><b>Multiple imputation</b> &mdash; <code>multi_impute()</code> &rarr; <code>with_imputations()</code> &rarr; <code>pool_mi()</code> for downstream inference via Rubin&rsquo;s rules.</li>
<li><b>Scales to large trees</b> &mdash; tested up to 10,000 species on a laptop.</li>
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
<tr><td><code>numeric</code></td>   <td>continuous</td><td>optional log + z-score</td><td>MSE</td>          <td>Phylogenetic BM</td></tr>
<tr><td><code>integer</code></td>   <td>count</td>     <td><code>log1p</code> + z-score</td><td>MSE</td> <td>Phylogenetic BM</td></tr>
<tr><td><code>factor(2)</code></td> <td>binary</td>    <td>0/1</td>              <td>BCE</td>          <td>Phylo label propagation</td></tr>
<tr><td><code>factor(&gt;2)</code></td><td>categorical</td><td>one-hot (K cols)</td><td>cross-entropy</td><td>Phylo label propagation</td></tr>
<tr><td><code>ordered</code></td>   <td>ordinal</td>   <td>integer + z-score</td><td>MSE</td>          <td>Phylogenetic BM</td></tr>
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
message passing, then broadcasts back to observation level.
</p>

<h2>Why impute? Downstream analysis with multiple imputation</h2>
<p>
Imputation is a <i>means</i>, not an end. The reason to fill in missing
cells is almost always to run a downstream analysis &mdash; a phylogenetic
regression, a comparative study of trait evolution, a trait-environment
test. Plugging point-estimate imputations directly into a regression
underestimates standard errors, because it treats imputed cells as if
they were observed. The consequences are well documented for ecological
and phylogenetic datasets in Nakagawa &amp; Freckleton
(<a href="https://doi.org/10.1016/j.tree.2008.06.014">2008, <i>TREE</i></a>)
and Nakagawa &amp; Freckleton
(<a href="https://doi.org/10.1007/s00265-010-1044-7">2011, <i>Behav Ecol Sociobiol</i></a>).
</p>
<p>
The remedy, due to Rubin (1987), is <b>multiple imputation</b>: generate
<code>M</code> stochastic completions of the trait matrix, fit the
downstream model on each, and pool the <code>M</code> model outputs via
Rubin&rsquo;s rules. pigauto exposes this workflow as three functions:
</p>
<pre><code>library(pigauto)
library(glmmTMB)   # attach so propto() is visible to the formula parser
library(ape)

data(avonet300, tree300)
df &lt;- avonet300; rownames(df) &lt;- df$Species_Key; df$Species_Key &lt;- NULL

# 1. Generate 100 complete datasets (MC-dropout sampling from the GNN)
mi &lt;- multi_impute(df, tree300, m = 100)

# 2. Fit a phylogenetic mixed model on each imputation.
#    Vphy must be a CORRELATION matrix (diag = 1), not a covariance:
#    glmmTMB&rsquo;s propto() estimates sigma^2 freely, so passing
#    vcv(tree) (diagonal = tree height) would rescale the variance
#    component.
Vphy &lt;- cov2cor(ape::vcv(tree300))
fits &lt;- with_imputations(mi, function(d) {
  d$species &lt;- factor(rownames(d), levels = rownames(Vphy))
  d$dummy   &lt;- factor(1)
  glmmTMB(
    log(Mass) ~ log(Wing.Length) + Trophic.Level +
      propto(0 + species | dummy, Vphy),
    data = d
  )
})

# 3. Pool with Rubin&rsquo;s rules (Barnard-Rubin df optional)
pool_mi(
  fits,
  coef_fun = function(f) fixef(f)$cond,
  vcov_fun = function(f) vcov(f)$cond
)
#&gt;               term estimate std.error    df statistic p.value  conf.low  conf.high    fmi   riv
#&gt;        (Intercept)   ...
#&gt;   log(Wing.Length)   ...
#&gt;     Trophic.LevelC   ...</code></pre>
<p>
The pooled table carries a few diagnostics worth reading. <code>fmi</code>
is the <i>fraction of missing information</i> for each coefficient &mdash;
the share of total variance that comes from imputation noise rather than
sampling. <code>riv</code> is the relative increase in variance due to
non-response. A good rule of thumb from Rubin is
<code>M &ge; 100 &times; fmi</code>; when <code>fmi</code> is small
(say &lt; 0.1), even <code>M = 20</code> is plenty, but if some
coefficients show <code>fmi &gt; 0.3</code> the pooled SEs continue to
drift until <code>M</code> is in the hundreds.
</p>
<div class="note">
<b>Bayesian alternative.</b> <code>pool_mi()</code> uses Rubin&rsquo;s rules,
which are the right tool for frequentist fits (<code>lm</code>,
<code>nlme::gls</code>, <code>glmmTMB</code>, <code>phylolm</code>, etc.).
For a fully Bayesian workflow you want to concatenate posterior samples
across imputations rather than decompose variance. The companion
<b>BACE</b> package (bundled in this repo) implements that via
<code>pool_posteriors()</code> on top of <code>MCMCglmm</code>; see
<code>BACE/R/pool_posteriors.R</code>. <code>pool_mi()</code> deliberately
rejects <code>MCMCglmm</code> fits so you do not accidentally mix the two
paradigms.
</div>

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
a dedicated phylogenetic comparative-methods package may be more appropriate.
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

# Dual-write: one copy ships with the installed package (inst/doc/),
# one copy lives under pkgdown/assets/ so pkgdown::build_site() picks
# it up and exposes it on the web site.
targets <- c(out, "pkgdown/assets/pigauto_intro.html")
for (t in targets) {
  dir.create(dirname(t), showWarnings = FALSE, recursive = TRUE)
  writeLines(html, t)
  cat("Wrote ", t, "\n", sep = "")
}
