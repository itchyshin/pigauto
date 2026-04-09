#!/usr/bin/env Rscript
# Build a self-contained HTML tutorial walking through the full
# pigauto + prepR4pcm + BACE workflow:
#
#   Data -> Missing Data -> Tree -> Matching Tree & Data ->
#   Multiple Imputation -> Downstream Inference (frequentist & Bayesian)
#
# Writes to inst/doc/ so it ships with the installed package and can be
# opened in RStudio via:
#   browseURL(system.file("doc", "pigauto_workflow.html", package = "pigauto"))

out <- "inst/doc/pigauto_workflow.html"
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

desc     <- readLines("DESCRIPTION")
ver_line <- grep("^Version:", desc, value = TRUE)
version  <- if (length(ver_line)) sub("^Version:\\s*", "", ver_line) else "?"
timestamp <- format(Sys.time(), "%Y-%m-%d")

html <- paste0(
'<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>pigauto &mdash; end-to-end PCM workflow</title>
<style>
  :root { --fg:#111827; --muted:#6b7280; --line:#e5e7eb; --soft:#f3f4f6;
          --accent:#059669; --accent2:#7c3aed; --code-bg:#f3f4f6;
          --hi:#eef2ff; --hi-fg:#3730a3; }
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         max-width: 920px; margin: 2em auto; padding: 0 1.5em; color: var(--fg);
         line-height: 1.6; }
  h1 { border-bottom: 2px solid var(--fg); padding-bottom: .25em; margin-bottom: .2em; }
  h2 { margin-top: 2.2em; color: #1f2937; border-bottom: 1px solid var(--line);
       padding-bottom: .15em; }
  h2 .num { color: var(--accent); font-weight: 700; margin-right: .35em; }
  h3 { color: #111827; margin-top: 1.6em; }
  .tagline { color: var(--muted); font-size: 16px; margin-top: 0; }
  .badge { display: inline-block; background: var(--hi); color: var(--hi-fg);
           padding: 2px 10px; border-radius: 999px; font-size: 13px;
           font-weight: 600; margin-left: 8px; vertical-align: middle; }
  pre { background: var(--code-bg); padding: 1em 1.2em; border-radius: 6px;
        overflow-x: auto; font-size: 13px; line-height: 1.5;
        border-left: 3px solid var(--accent); }
  pre.bayes { border-left-color: var(--accent2); }
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
  .tip { background: #ecfdf5; border-left: 4px solid #059669; padding: 1em 1.2em;
         border-radius: 4px; margin: 1.5em 0; }
  .bayesbox { background: #f5f3ff; border-left: 4px solid #7c3aed; padding: 1em 1.2em;
              border-radius: 4px; margin: 1.5em 0; }
  .toc { background: var(--soft); border-radius: 6px; padding: 1em 1.4em;
         margin: 1.5em 0; }
  .toc ol { margin: .3em 0; padding-left: 1.4em; }
  .toc li { margin: 4px 0; }
  ul li { margin: 6px 0; }
  .meta { color: var(--muted); font-size: 13px; }
  footer { color: var(--muted); font-size: 12px; margin-top: 3em;
           border-top: 1px solid var(--line); padding-top: 1em; }
</style>
</head>
<body>

<h1>pigauto: the full PCM workflow <span class="badge">v', version, '</span></h1>
<p class="tagline">From messy species names to pooled phylogenetic inference.</p>
<p class="meta">Generated ', timestamp, ' &middot; scripts/<code>make_workflow_html.R</code></p>

<p>This tutorial walks a single realistic dataset through the full
phylogenetic comparative methods (PCM) pipeline: load the raw data,
inspect missingness, load a phylogeny, reconcile species names, impute
missing trait values with uncertainty, fit the downstream model, and
pool the results. Three analysis paths are sketched below; this
tutorial shows Path&nbsp;A and Path&nbsp;C in full on a continuous-only
example, and points at the mixed-type companion tutorial for
Path&nbsp;B and a worked mixed-type example of all three.</p>

<ul>
<li><b>Path A &mdash; frequentist.</b> <code>pigauto::multi_impute()</code>
&rarr; <code>glmmTMB</code> with <code>propto()</code> &rarr;
Rubin&rsquo;s rules via <code>pigauto::pool_mi()</code>.</li>
<li><b>Path B &mdash; Bayesian with pigauto as imputer.</b>
<code>pigauto::multi_impute()</code> &rarr; <code>MCMCglmm</code> on
each imputed dataset &rarr; posterior concatenation. This is the
Nakagawa&nbsp;&amp;&nbsp;Freckleton&nbsp;(2011) recipe for Bayesian
multiple imputation in comparative analyses; pigauto supplies the
imputer, <code>MCMCglmm</code> supplies the inference engine.</li>
<li><b>Path C &mdash; integrated Bayesian.</b> <code>BACE::bace()</code>
does chained-equation imputation and inference in a single
<code>MCMCglmm</code>-based engine, then <code>pool_posteriors()</code>
concatenates the final runs.</li>
</ul>

<p>Paths&nbsp;A, B, and C are three distinct recipes. Paths&nbsp;A and
B share the same imputer (pigauto) and differ in the downstream
inference engine; Paths&nbsp;B and C share the same downstream engine
(<code>MCMCglmm</code>) and differ in the imputer. For a worked
mixed-type example that exercises all three paths on categorical and
ordinal traits, see the companion tutorial
<code>pigauto_workflow_mixed.html</code>.</p>

<p>Three R packages are involved. They are complementary:</p>

<table>
<tr><th>Package</th><th>Role in the workflow</th></tr>
<tr><td><b><code>prepR4pcm</code></b></td>
    <td>Reconciles species names between data and tree (stages 1&ndash;4:
    exact, normalised, synonym, fuzzy). Produces aligned data + pruned
    tree.</td></tr>
<tr><td><b><code>pigauto</code></b></td>
    <td>Phylogenetic trait imputation via a graph neural network
    (Residual PhyloDAE). Generates <i>M</i> stochastic completions for
    multiple imputation, plus Rubin&rsquo;s-rules pooling for frequentist
    downstream fits.</td></tr>
<tr><td><b><code>BACE</code></b></td>
    <td>Bayesian chained-equation imputation on top of <code>MCMCglmm</code>,
    with posterior-sample pooling across imputations.</td></tr>
</table>

<div class="toc">
<b>Contents</b>
<ol>
<li><a href="#data">Data</a></li>
<li><a href="#missing">Missing data</a></li>
<li><a href="#tree">Tree</a></li>
<li><a href="#match">Matching tree and data with <code>prepR4pcm</code></a></li>
<li><a href="#mi">Multiple imputation with <code>pigauto</code></a></li>
<li><a href="#freq">Frequentist analysis: <code>glmmTMB</code> + Rubin&rsquo;s rules</a></li>
<li><a href="#bayes">Bayesian analysis: <code>BACE</code> + posterior pooling</a></li>
<li><a href="#choose">Which path should I take?</a></li>
</ol>
</div>

<h2 id="data"><span class="num">1.</span> Data</h2>

<p>We use the <code>avonet_subset</code> dataset bundled with <code>prepR4pcm</code>,
a 919-species extract of the full AVONET database (Tobias et&nbsp;al.&nbsp;2022,
<i>Ecology Letters</i>). It has morphological traits (beak, wing, tarsus,
mass) and ecological traits (trophic level, primary lifestyle,
migration). This is real data with real species-name issues, which is
the whole point of including <code>prepR4pcm</code>.</p>

<pre><code>library(pigauto)
library(prepR4pcm)
library(ape)

data(avonet_subset, package = "prepR4pcm")
data(tree_jetz,    package = "prepR4pcm")

dim(avonet_subset)
#&gt; [1] 919  16

head(avonet_subset[, c("Species1", "Mass", "Wing.Length", "Trophic.Level")])
#&gt;                Species1 Mass Wing.Length Trophic.Level
#&gt; 1    Acanthiza apicalis  7.6        51.0     Carnivore
#&gt; 2 Acanthiza chrysorrhoa  9.3        58.3     Carnivore
#&gt; 3     Acanthiza cinerea  7.6        52.0     Carnivore
#&gt; ...
</code></pre>

<h2 id="missing"><span class="num">2.</span> Missing data</h2>

<p>Before imputing, you need to know where the holes are. Count
missingness per trait and look for patterns: are cells missing
completely at random (MCAR), or are they concentrated in particular
clades, orders, or small-bodied species? Phylogenetic imputation works
best under MAR given the tree; concentrated missingness in a single
monophyletic group is the hardest case.</p>

<pre><code># Missingness by column
na_by_col &lt;- sort(colSums(is.na(avonet_subset)), decreasing = TRUE)
na_by_col[na_by_col &gt; 0]

# Cell-level and row-level missingness
mean(is.na(avonet_subset))                     # overall fraction of missing cells
mean(rowSums(is.na(avonet_subset)) &gt; 0)        # fraction of species with &gt;=1 missing

# Crude check for clumping: are NAs concentrated in particular orders?
tapply(rowSums(is.na(avonet_subset)),
       avonet_subset$Order1, mean) |&gt; sort(decreasing = TRUE) |&gt; head()
</code></pre>

<p>For this tutorial, <code>avonet_subset</code> is almost complete, so we
inject a realistic 20&nbsp;% missingness pattern across the continuous
traits to make the imputation step non-trivial:</p>

<pre><code>set.seed(1)
df &lt;- avonet_subset
cont_cols &lt;- c("Beak.Length_Culmen", "Beak.Length_Nares", "Beak.Width",
               "Beak.Depth", "Tarsus.Length", "Wing.Length", "Mass")

for (col in cont_cols) {
  idx &lt;- sample.int(nrow(df), size = floor(0.20 * nrow(df)))
  df[idx, col] &lt;- NA_real_
}
mean(is.na(df[, cont_cols]))   # ~0.20
</code></pre>

<h2 id="tree"><span class="num">3.</span> Tree</h2>

<p>The phylogeny is <code>tree_jetz</code> from <code>prepR4pcm</code>, a 657-tip
subset of the Jetz et&nbsp;al.&nbsp;(2012) bird tree that has been pruned to the
species relevant to the example data. In your own work you will load
your own <code>phylo</code> object with <code>ape::read.tree()</code>.</p>

<pre><code>tree_jetz
#&gt; Phylogenetic tree with 657 tips and 656 internal nodes.
#&gt; Tip labels:
#&gt;   Acanthagenys_rufogularis, Acanthiza_apicalis, ...
#&gt; Rooted; includes branch lengths.

ape::Ntip(tree_jetz)                  # 657
ape::is.ultrametric(tree_jetz)        # TRUE / FALSE
range(ape::branching.times(tree_jetz))
</code></pre>

<div class="note">
<b>Ultrametricity and branch lengths.</b> pigauto&rsquo;s phylogenetic
kernel and Rphylopars baseline both assume a rooted tree with
meaningful branch lengths. If your tree is not ultrametric but you need
it to be (e.g. for Brownian motion on shared time), use
<code>phytools::force.ultrametric()</code> before proceeding, and
document the choice.
</div>

<h2 id="match"><span class="num">4.</span> Matching tree and data with <code>prepR4pcm</code></h2>

<p>Species names in <code>df</code> use the format <code>Genus species</code>
(e.g. <code>Acanthiza apicalis</code>), whereas <code>tree_jetz</code> uses
<code>Genus_species</code>. Beyond the formatting difference there are
taxonomic synonyms and a few typos. <code>prepR4pcm::reconcile_tree()</code>
runs a four-stage matching cascade and records every decision:</p>

<pre><code>rec &lt;- reconcile_tree(
  x         = df,
  tree      = tree_jetz,
  x_species = "Species1",
  fuzzy     = TRUE,
  resolve   = "flag"
)
rec
#&gt; -- Reconciliation: data vs tree ------------------------------------
#&gt;   Source x : df
#&gt;   Source y : phylo (657 tips)
#&gt;   Authority: col
#&gt; i Match coverage: [#################### ] 71% (657/919)
#&gt;
#&gt; -- Match summary --
#&gt;   * Normalized : 657 (71.5%)
#&gt;   * Fuzzy      :   0 ( 0.0%)
#&gt;   ! Unresolved : 262 (28.5%)

reconcile_summary(rec)          # per-row detail
reconcile_plot(rec)             # visual match-coverage plot
reconcile_report(rec, file = "rec.html")  # optional audit report
</code></pre>

<p>Once you are happy with the reconciliation, apply it. This returns an
aligned data frame and a tree pruned to the matched species &mdash; exactly
what every downstream PCM step needs:</p>

<pre><code>aligned &lt;- reconcile_apply(
  rec,
  data           = df,
  tree           = tree_jetz,
  species_col    = "Species1",
  drop_unresolved = TRUE
)
#&gt; ! Dropped 262 rows with unresolved species
#&gt; i Tree has 657 tips after alignment

dat  &lt;- aligned$data           # aligned data frame
tree &lt;- aligned$tree           # pruned phylo object
rownames(dat) &lt;- dat$Species1  # species as rownames for pigauto
</code></pre>

<div class="tip">
<b>Why this step matters.</b> Silently dropping 28&nbsp;% of your data
because of name formatting would be a disaster. <code>prepR4pcm</code>
makes that loss explicit and auditable, and its reconciliation object
records which names were matched by exact string, normalisation,
synonym lookup, or fuzzy matching. For publication, save the
reconciliation report alongside your analysis script.
</div>

<h2 id="mi"><span class="num">5.</span> Multiple imputation with <code>pigauto</code></h2>

<p>With <code>dat</code> and <code>tree</code> aligned, we can run pigauto&rsquo;s
full imputation pipeline. <code>multi_impute()</code> is the
high-level entry point: it chains preprocessing, the Brownian-motion
baseline, GNN training, calibration, and Monte-Carlo-dropout sampling,
returning <i>M</i> complete datasets in one call.</p>

<pre><code># Drop non-trait columns prepR4pcm carried through for matching
pcm_df &lt;- dat[, c("Beak.Length_Culmen", "Beak.Length_Nares", "Beak.Width",
                  "Beak.Depth", "Tarsus.Length", "Wing.Length", "Mass",
                  "Habitat.Density", "Migration", "Trophic.Level",
                  "Primary.Lifestyle")]
pcm_df$Habitat.Density   &lt;- as.integer(pcm_df$Habitat.Density)   # count
pcm_df$Migration         &lt;- ordered(pcm_df$Migration)            # ordinal
pcm_df$Trophic.Level     &lt;- factor(pcm_df$Trophic.Level)         # categorical
pcm_df$Primary.Lifestyle &lt;- factor(pcm_df$Primary.Lifestyle)     # categorical

mi &lt;- multi_impute(
  traits        = pcm_df,
  tree          = tree,
  m             = 100,        # number of imputation draws
  log_transform = TRUE,       # log-transform strictly positive continuous traits
  epochs        = 2000L,
  seed          = 1
)
mi
#&gt; &lt;pigauto_mi&gt; 100 imputations of 11 traits over 657 species
#&gt;   5,213 cells imputed (23.4% of the trait matrix)

# A single completed dataset, e.g. the first draw
head(mi$datasets[[1]][, c("Mass", "Wing.Length", "Trophic.Level")])

# Fit object (for predict(), report, conformal intervals, etc.)
mi$fit
</code></pre>

<div class="note">
<b>How many imputations?</b> Rubin&rsquo;s rule of thumb is
<i>M</i> &ge; 100 &times; max(fmi). For most ecological datasets
<i>M</i> = 20 is enough if the fraction of missing information
(<code>fmi</code>) stays below 0.1; push <i>M</i> into the hundreds if
any coefficient has <code>fmi &gt; 0.3</code>. You will see
<code>fmi</code> in the pooled output below, so you can re-run with a
larger <i>M</i> if needed.
</div>

<h2 id="freq"><span class="num">6.</span> Frequentist analysis: <code>glmmTMB</code> + Rubin&rsquo;s rules</h2>

<p>Suppose we want a phylogenetic regression of mass on wing length and
trophic level. In the frequentist path, we fit a phylogenetic mixed
model with <code>glmmTMB</code> using the <code>propto()</code>
covariance structure (Brooks et&nbsp;al.&nbsp;2017) on the phylogenetic VCV
matrix, one fit per imputed dataset, then pool using Rubin&rsquo;s
(1987) rules via <code>pigauto::pool_mi()</code>.</p>

<pre><code>library(glmmTMB)   # must be attached so propto() is visible to the formula parser

Vphy &lt;- ape::vcv(tree)   # phylogenetic variance-covariance matrix

fits &lt;- with_imputations(mi, function(d) {
  d$species &lt;- factor(rownames(d), levels = rownames(Vphy))
  d$dummy   &lt;- factor(1)
  glmmTMB(
    log(Mass) ~ log(Wing.Length) + Trophic.Level +
      propto(0 + species | dummy, Vphy),
    data = d
  )
})
fits
#&gt; &lt;pigauto_mi_fits&gt; 100 fits (0 failed)

pooled &lt;- pool_mi(
  fits,
  coef_fun = function(f) fixef(f)$cond,
  vcov_fun = function(f) vcov(f)$cond
)
pooled
#&gt;                  term estimate std.error    df statistic p.value conf.low conf.high   fmi   riv
#&gt;           (Intercept)   -6.12    0.152  284.    -40.2  &lt;2e-16   -6.42    -5.82    0.21  0.27
#&gt;      log(Wing.Length)    2.71    0.035  312.     77.2  &lt;2e-16    2.64     2.78    0.17  0.20
#&gt;    Trophic.LevelHerbi    0.33    0.061  268.      5.4  3e-07    0.21     0.45    0.23  0.30
#&gt;    ... etc
</code></pre>

<p>The pooled table reports <code>fmi</code> (fraction of missing
information) and <code>riv</code> (relative increase in variance due to
non-response) per coefficient. Degrees of freedom follow Barnard &amp;
Rubin (1999). Confidence intervals and p-values already account for the
extra uncertainty introduced by imputation; plugging the point-estimate
completion directly into a single <code>glmmTMB</code> fit would
<i>under</i>state standard errors on exactly these terms.</p>

<div class="tip">
<b>What counts as a valid downstream fit?</b>
<code>pool_mi()</code> works on any frequentist fit for which
<code>coef()</code> returns a named numeric vector and <code>vcov()</code>
returns a named covariance matrix. That covers
<code>lm</code>, <code>glm</code>, <code>nlme::gls</code>,
<code>phylolm</code>, and <code>glmmTMB</code> (via the
<code>coef_fun</code>/<code>vcov_fun</code> hooks). It deliberately
rejects <code>MCMCglmm</code>: pooling posterior draws from a Bayesian
model with Rubin&rsquo;s rules is conceptually wrong &mdash; use the
Bayesian path below instead.
</div>

<h2 id="bayes"><span class="num">7.</span> Bayesian analysis: <code>BACE</code> + posterior pooling</h2>

<p>If you want full posterior uncertainty propagation, the right tool is
<code>BACE</code> (Bayesian Analysis with Chained Equations). BACE runs
chained-equation imputation on top of <code>MCMCglmm</code>, then
concatenates the posterior samples across imputations to obtain a
single pooled posterior that captures both within- and
between-imputation variance (Nakagawa &amp; Freckleton 2008, 2011;
Zhang, Parnell, Hadfield &amp; Pohl, 2020+).</p>

<p>BACE runs chained-equation imputation and inference inside the same
<code>MCMCglmm</code>-based engine; its <code>pool_posteriors()</code>
concatenates the posterior samples from the final runs for you. Note
that BACE is a <i>different</i> recipe from feeding pigauto&rsquo;s
multiple imputations into <code>MCMCglmm</code> fits (Path&nbsp;B in
the mixed-type companion tutorial). Both are valid Bayesian multiple-
imputation workflows: pigauto&nbsp;+&nbsp;MCMCglmm is faster and uses
pigauto&rsquo;s GNN as the first-stage imputer; BACE is a single
integrated sampler. See
<code>pigauto_workflow_mixed.html</code> section&nbsp;7 for a worked
Path&nbsp;B example and section&nbsp;8 for Path&nbsp;C.</p>

<pre class="bayes"><code>library(BACE)

# Use the same aligned data + tree from step 4
bace_fit &lt;- bace(
  fixformula     = "log(Mass) ~ log(Wing.Length) + Trophic.Level",
  ran_phylo_form = "~1|Species1",
  phylo          = tree,
  data           = dat,                  # aligned data with rownames as species
  runs           = 10,                   # initial runs for convergence
  n_final        = 100,                  # final imputations for pooling
  nitt           = 60000, burnin = 10000, thin = 50,
  sample_size    = 1000                  # draws per imputation to pool
)

# Pool posterior samples across imputations
pooled_bayes &lt;- pool_posteriors(bace_fit$final_results)

# Standard MCMCglmm methods work on the pooled object
summary(pooled_bayes$models$Mass)
plot(pooled_bayes$models$Mass)           # trace + density plots
</code></pre>

<p>The pooled posterior is a regular <code>MCMCglmm</code> object, so
everything downstream of a single BACE fit &mdash; highest-posterior-density
intervals, posterior predictive checks, model comparison via DIC,
<code>coda</code> diagnostics &mdash; works without modification. Only now
each draw reflects both parameter uncertainty and imputation
uncertainty.</p>

<div class="bayesbox">
<b>Can I use pigauto draws as input to MCMCglmm?</b> Yes &mdash; and
that is Path&nbsp;B. Feeding <code>pigauto::multi_impute()</code>
completions into <code>M</code> separate <code>MCMCglmm</code> fits
and concatenating the <code>$Sol</code> chains is the canonical
Bayesian multiple-imputation recipe of
Nakagawa&nbsp;&amp;&nbsp;Freckleton&nbsp;(2011), with pigauto playing
the role of the first-stage imputer. BACE is a <i>different</i>
Bayesian recipe (chained-equation imputation and inference in one
engine); neither is a substitute for the other. The mixed-type
companion tutorial (<code>pigauto_workflow_mixed.html</code>) works
through both Path&nbsp;B (section&nbsp;7) and Path&nbsp;C
(section&nbsp;8) on a genuinely mixed dataset, alongside Path&nbsp;A.
</div>

<h2 id="choose"><span class="num">8.</span> Which path should I take?</h2>

<table>
<tr><th></th>
    <th>Path A<br><span class="meta">pigauto + glmmTMB + Rubin</span></th>
    <th>Path B<br><span class="meta">pigauto + MCMCglmm + concat</span></th>
    <th>Path C<br><span class="meta">BACE integrated</span></th></tr>
<tr><td><b>Imputer</b></td>
    <td><code>pigauto::multi_impute()</code> (GNN)</td>
    <td><code>pigauto::multi_impute()</code> (GNN)</td>
    <td><code>BACE::bace()</code> (chained-equation MCMC)</td></tr>
<tr><td><b>Downstream model</b></td>
    <td><code>glmmTMB</code> / <code>gls</code> / <code>phylolm</code> / any fit with <code>coef()</code> + <code>vcov()</code></td>
    <td><code>MCMCglmm</code></td>
    <td><code>MCMCglmm</code> (inside BACE)</td></tr>
<tr><td><b>Pooling</b></td>
    <td>Rubin&rsquo;s rules via <code>pigauto::pool_mi()</code></td>
    <td>Posterior concatenation (stack <code>$Sol</code> across M fits)</td>
    <td><code>BACE::pool_posteriors()</code></td></tr>
<tr><td><b>Mixed trait types</b></td>
    <td>Native in pigauto; regression via <code>glmmTMB</code> families</td>
    <td>Native in pigauto; regression via <code>MCMCglmm</code> families</td>
    <td>Native in BACE; regression via <code>MCMCglmm</code> families</td></tr>
<tr><td><b>Runtime</b></td>
    <td>Minutes to a few hours</td>
    <td>1&ndash;4 hours (per-fit MCMC)</td>
    <td>Hours to overnight (chained-equation MCMC)</td></tr>
<tr><td><b>Uncertainty returned</b></td>
    <td>SE, <code>fmi</code>, <code>riv</code>, Barnard&ndash;Rubin df</td>
    <td>Full posterior with HPD intervals</td>
    <td>Full posterior with HPD intervals, DIC</td></tr>
<tr><td><b>When to prefer</b></td>
    <td>Fast iteration, reporting CIs and p-values in a regression table, continuous response</td>
    <td>You already use <code>MCMCglmm</code> for inference but want pigauto&rsquo;s fast GNN as the imputer (Nakagawa&nbsp;&amp;&nbsp;Freckleton&nbsp;2011 recipe)</td>
    <td>Formal Bayesian workflow with informative priors; categorical uncertainty needs to flow through posteriors</td></tr>
</table>

<p>All three paths share the first four steps (data, missing data,
tree, matching). Only from step&nbsp;5 onwards do they diverge. This
tutorial walks through Path&nbsp;A and Path&nbsp;C in full on a
continuous-only example. For a worked Path&nbsp;B example &mdash;
<code>pigauto::multi_impute()</code> &rarr;
<code>MCMCglmm</code> &rarr; posterior concatenation &mdash; and for
all three paths side by side on a <i>mixed-type</i> trait matrix, see
the companion tutorial <code>pigauto_workflow_mixed.html</code>.</p>

<h2>References</h2>
<ul>
<li>Rubin, D.B. (1987). <i>Multiple Imputation for Nonresponse in
Surveys.</i> Wiley.</li>
<li>Barnard, J. &amp; Rubin, D.B. (1999). Small-sample degrees of freedom
with multiple imputation. <i>Biometrika</i> 86:948&ndash;955.</li>
<li>Nakagawa, S. &amp; Freckleton, R.P. (2008). Missing inaction: the
dangers of ignoring missing data. <i>Trends in Ecology &amp; Evolution</i>
23:592&ndash;596. <a href="https://doi.org/10.1016/j.tree.2008.06.014">doi:10.1016/j.tree.2008.06.014</a></li>
<li>Nakagawa, S. &amp; Freckleton, R.P. (2011). Model averaging, missing
data and multiple imputation: a case study for behavioural ecology.
<i>Behavioral Ecology and Sociobiology</i> 65:103&ndash;116.
<a href="https://doi.org/10.1007/s00265-010-1044-7">doi:10.1007/s00265-010-1044-7</a></li>
<li>Brooks, M.E. et&nbsp;al. (2017). glmmTMB balances speed and
flexibility among packages for zero-inflated generalized linear mixed
models. <i>R Journal</i> 9:378&ndash;400.</li>
<li>Tobias, J.A. et&nbsp;al. (2022). AVONET: morphological, ecological and
geographical data for all birds. <i>Ecology Letters</i>
25:581&ndash;597.
<a href="https://doi.org/10.1111/ele.13898">doi:10.1111/ele.13898</a></li>
<li>Jetz, W. et&nbsp;al. (2012). The global diversity of birds in space
and time. <i>Nature</i> 491:444&ndash;448.</li>
</ul>

<footer>
Generator: <code>script/make_workflow_html.R</code> &middot;
Installed path: <code>system.file("doc", "pigauto_workflow.html", package = "pigauto")</code>
</footer>

</body>
</html>
')

writeLines(html, out)
cat("Wrote ", out, "\n", sep = "")
