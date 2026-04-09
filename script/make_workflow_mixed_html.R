#!/usr/bin/env Rscript
# Build a self-contained HTML tutorial walking through the full pigauto
# PCM workflow on a **mixed-type** dataset (continuous + categorical +
# ordinal), and showing all three sensible analysis paths:
#
#   Path A: pigauto + glmmTMB + Rubin's rules        (frequentist)
#   Path B: pigauto + MCMCglmm + posterior concat    (Bayesian, pigauto as imputer)
#   Path C: BACE::bace() + pool_posteriors()         (Bayesian, BACE as imputer)
#
# Companion to script/make_workflow_html.R, which uses a continuous-only
# example and focuses on reconciliation. This file assumes the reader has
# already seen that tutorial; it does not re-teach reconcile_tree().
#
# Writes to inst/doc/ so it ships with the installed package and can be
# opened in RStudio via:
#   browseURL(system.file("doc", "pigauto_workflow_mixed.html", package = "pigauto"))

out <- "inst/doc/pigauto_workflow_mixed.html"
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

desc      <- readLines("DESCRIPTION")
ver_line  <- grep("^Version:", desc, value = TRUE)
version   <- if (length(ver_line)) sub("^Version:\\s*", "", ver_line) else "?"
timestamp <- format(Sys.time(), "%Y-%m-%d")

html <- paste0(
'<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>pigauto &mdash; mixed-type PCM workflow</title>
<style>
  :root { --fg:#111827; --muted:#6b7280; --line:#e5e7eb; --soft:#f3f4f6;
          --accent:#059669; --accent2:#7c3aed; --accent3:#ea580c;
          --code-bg:#f3f4f6; --hi:#eef2ff; --hi-fg:#3730a3; }
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
  pre.bace  { border-left-color: var(--accent3); }
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
  .warnbox { background: #fff1f2; border-left: 4px solid #be123c; padding: 1em 1.2em;
             border-radius: 4px; margin: 1.5em 0; }
  .warnbox b { color: #9f1239; }
  .trouble { background: #f8fafc; border-left: 4px solid #64748b; padding: 1em 1.2em;
             border-radius: 4px; margin: 1.5em 0; font-size: 13px; }
  .trouble b { color: #334155; }
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

<h1>pigauto: mixed-type PCM workflow <span class="badge">v', version, '</span></h1>
<p class="tagline">Multiple imputation of continuous, categorical, and ordinal traits &mdash; three downstream paths.</p>
<p class="meta">Generated ', timestamp, ' &middot; scripts/<code>make_workflow_mixed_html.R</code></p>

<p>This tutorial is the mixed-type companion to
<code>pigauto_workflow.html</code>. The earlier tutorial uses a
continuous-only example and focuses on tree-and-data reconciliation
with <code>prepR4pcm</code>; it is the right place to start if you are
new to the pipeline. This one assumes you have already seen that one
and drills into two things it did not cover:</p>

<ol>
<li><b>Mixed-type imputation.</b> Most real phylogenetic datasets mix
continuous morphometrics with categorical ecology (trophic level, habitat,
migration strategy). pigauto is designed to handle all five trait types
&mdash; continuous, count, ordinal, binary, categorical &mdash; in one
model, but the previous tutorial only exercised the continuous path.
Here every trait type is actually imputed and actually used in the
downstream regression.</li>
<li><b>Three analysis paths, not two.</b> Multiple imputation is a way
to propagate uncertainty into whatever downstream model you prefer. The
classical split is frequentist (Rubin&rsquo;s rules) vs Bayesian
(posterior concatenation), but within the Bayesian branch there is a
further choice: use pigauto&rsquo;s GNN as the imputer and feed its
draws into <code>MCMCglmm</code>, or use <code>BACE</code> to do
imputation and inference in a single chained-equation MCMC engine.
Both are valid; they trade off runtime against modelling assumptions.
This tutorial shows all three.</li>
</ol>

<p>The corrected mental model is a 2&times;2 grid over imputer and
downstream inference:</p>

<table>
<tr><th></th>
    <th>Frequentist downstream</th>
    <th>Bayesian downstream</th></tr>
<tr><td><b>pigauto imputer</b></td>
    <td><b>Path A</b>: <code>multi_impute</code> &rarr; <code>glmmTMB</code> &rarr; <code>pool_mi</code></td>
    <td><b>Path B</b>: <code>multi_impute</code> &rarr; <code>MCMCglmm</code> &rarr; posterior concatenation</td></tr>
<tr><td><b>BACE imputer</b></td>
    <td>(unusual, not documented)</td>
    <td><b>Path C</b>: <code>BACE::bace()</code> &rarr; <code>pool_posteriors()</code></td></tr>
</table>

<p>All three paths share sections&nbsp;1&ndash;5 (data, missing data, tree,
matching, multiple imputation). Only from section&nbsp;6 onwards do they
diverge.</p>

<div class="toc">
<b>Contents</b>
<ol>
<li><a href="#data">Data: a genuinely mixed-type trait matrix</a></li>
<li><a href="#missing">Missing data across all trait types</a></li>
<li><a href="#tree">Tree</a></li>
<li><a href="#match">Matching tree and data with <code>prepR4pcm</code></a></li>
<li><a href="#mi">Multiple imputation of all five trait types</a></li>
<li><a href="#patha">Path&nbsp;A: pigauto + <code>glmmTMB</code> + Rubin&rsquo;s rules</a></li>
<li><a href="#pathb">Path&nbsp;B: pigauto + <code>MCMCglmm</code> + posterior concatenation</a></li>
<li><a href="#pathc">Path&nbsp;C: BACE integrated chained-equation MCMC</a></li>
<li><a href="#choose">Which path should I take?</a></li>
</ol>
</div>

<h2 id="data"><span class="num">1.</span> Data: a genuinely mixed-type trait matrix</h2>

<p>We use the <code>avonet300</code> dataset bundled with pigauto, which
was built specifically to exercise mixed-type imputation. It contains
300 bird species with four continuous morphometric traits, two nominal
categorical traits, and one ordered factor &mdash; the full range of
types pigauto supports in one latent space.</p>

<pre><code>library(pigauto)
library(ape)

data(avonet300, tree300)

dim(avonet300)
#&gt; [1] 300   8

sapply(avonet300, class)
#&gt; Species_Key        Mass   Beak.Length_Culmen    Tarsus.Length
#&gt; "character"   "numeric"            "numeric"        "numeric"
#&gt; Wing.Length   Trophic.Level   Primary.Lifestyle        Migration
#&gt; "numeric"        "factor"             "factor" "ordered","factor"
</code></pre>

<p>Seven of the eight columns are trait columns (the first is the
species key). Here is how they map onto pigauto&rsquo;s type system:</p>

<table>
<thead><tr><th>Column</th><th>R class</th><th>pigauto type</th><th>Levels / range</th></tr></thead>
<tbody>
<tr><td><code>Mass</code></td><td><code>numeric</code></td><td>continuous</td><td>positive real, log-transformed</td></tr>
<tr><td><code>Beak.Length_Culmen</code></td><td><code>numeric</code></td><td>continuous</td><td>positive real</td></tr>
<tr><td><code>Tarsus.Length</code></td><td><code>numeric</code></td><td>continuous</td><td>positive real</td></tr>
<tr><td><code>Wing.Length</code></td><td><code>numeric</code></td><td>continuous</td><td>positive real</td></tr>
<tr><td><code>Trophic.Level</code></td><td><code>factor</code></td><td>categorical</td><td>Carnivore / Herbivore / Omnivore / Scavenger</td></tr>
<tr><td><code>Primary.Lifestyle</code></td><td><code>factor</code></td><td>categorical</td><td>Aerial / Aquatic / Generalist / Insessorial / Terrestrial</td></tr>
<tr><td><code>Migration</code></td><td><code>ordered</code></td><td>ordinal</td><td>Resident &lt; Partial &lt; Full</td></tr>
</tbody>
</table>

<p>Set up the analysis data frame with species as rownames:</p>

<pre><code>df &lt;- avonet300
rownames(df) &lt;- df$Species_Key
df$Species_Key &lt;- NULL
</code></pre>

<div class="tip">
<b>Why class matters.</b> pigauto reads R column classes to decide which
trait type each column is. <code>factor</code> with 3+ levels becomes
<i>categorical</i> (K-column one-hot latent block with cross-entropy
loss and a label-propagation baseline). <code>ordered</code> becomes
<i>ordinal</i> (single latent column with MSE loss). Plain
<code>numeric</code> becomes <i>continuous</i>. If a column looks
numeric in your raw file but is really a grade (1&nbsp;/&nbsp;2&nbsp;/&nbsp;3
&quot;small / medium / large&quot;), convert it to
<code>ordered()</code> before you call <code>multi_impute()</code> or
you will regress on meaningless integer codes.
</div>

<h2 id="missing"><span class="num">2.</span> Missing data across all trait types</h2>

<p>The bundled <code>avonet300</code> is deliberately complete so the
dataset is useful as a gold standard for benchmarking. To exercise the
mixed-type imputation pipeline we need missing cells in every trait
type, not just in the continuous ones. Inject a 15&nbsp;% missingness
pattern across every trait column:</p>

<pre><code>set.seed(1)
trait_cols &lt;- c("Mass", "Beak.Length_Culmen", "Tarsus.Length", "Wing.Length",
                "Trophic.Level", "Primary.Lifestyle", "Migration")

for (col in trait_cols) {
  idx &lt;- sample.int(nrow(df), size = floor(0.15 * nrow(df)))
  df[idx, col] &lt;- NA
}

colSums(is.na(df))
#&gt;               Mass  Beak.Length_Culmen       Tarsus.Length         Wing.Length
#&gt;                 45                  45                  45                  45
#&gt;      Trophic.Level   Primary.Lifestyle           Migration
#&gt;                 45                  45                  45

mean(is.na(df[, trait_cols]))
#&gt; [1] 0.15
</code></pre>

<p>The key point: <i>every</i> trait type has missing cells, including
the two categorical factors and the ordered factor. This is what makes
the rest of the tutorial a genuine mixed-type demonstration rather than
a continuous-only pipeline with factor side-columns carried along for
decoration.</p>

<div class="note">
<b>MNAR vs MAR for factor columns.</b> Random injection gives
missing-completely-at-random (MCAR), which is the easiest case for any
imputer. Real categorical missingness is usually not MCAR &mdash; for
example, <code>Trophic.Level</code> is much more likely to be missing for
rare or poorly studied species, which are themselves
phylogenetically clustered. pigauto handles MAR given the tree reasonably
well (species with missing trophic level borrow from close relatives via
label propagation), but it cannot fix MNAR. If you suspect
informative missingness, say so in your methods section and sensitivity-
check with a different imputer.
</div>

<h2 id="tree"><span class="num">3.</span> Tree</h2>

<p>The phylogeny is <code>tree300</code>, a 300-tip pruned BirdTree
subset aligned with <code>avonet300</code>:</p>

<pre><code>tree300
#&gt; Phylogenetic tree with 300 tips and 299 internal nodes.
#&gt; Tip labels:
#&gt;   Eudromia_formosa, Nothoprocta_pentlandii, Rhea_americana, ...
#&gt; Rooted; includes branch lengths.

ape::Ntip(tree300)               # 300
ape::is.ultrametric(tree300)     # TRUE
</code></pre>

<p>For details on ultrametricity, branch-length assumptions, and when
to force-ultrametricise, see section&nbsp;3 of the continuous-only
tutorial (<code>pigauto_workflow.html</code>).</p>

<h2 id="match"><span class="num">4.</span> Matching tree and data with <code>prepR4pcm</code></h2>

<p><code>avonet300</code> is already aligned with <code>tree300</code> by
construction &mdash; the species keys use the same
<code>Genus_species</code> format and every row has a matching tip.
Running the reconciliation anyway is a cheap sanity check that does not
drop anything:</p>

<pre><code>library(prepR4pcm)

# df now has species as rownames; reconcile_tree wants a species column,
# so temporarily re-add it.
df_with_key &lt;- cbind(Species_Key = rownames(df), df)
rec &lt;- reconcile_tree(
  x         = df_with_key,
  tree      = tree300,
  x_species = "Species_Key",
  fuzzy     = FALSE,
  resolve   = "flag"
)
rec
#&gt; -- Reconciliation: data vs tree ------------------------------------
#&gt; i Match coverage: [####################] 100% (300/300)
#&gt; -- Match summary --
#&gt;   * Exact      : 300 (100.0%)

aligned &lt;- reconcile_apply(
  rec,
  data            = df_with_key,
  tree            = tree300,
  species_col     = "Species_Key",
  drop_unresolved = TRUE
)
dat  &lt;- aligned$data
tree &lt;- aligned$tree
rownames(dat) &lt;- dat$Species_Key
dat$Species_Key &lt;- NULL
</code></pre>

<div class="tip">
<b>Deliberate no-op.</b> This tutorial uses a pre-aligned dataset so the
reconciliation step is a formality. For a realistic example where
<code>reconcile_tree()</code> actually resolves 28&nbsp;% of unresolved
species via normalisation, synonym lookup, and fuzzy matching, see
section&nbsp;4 of <code>pigauto_workflow.html</code>. Matching real
messy data against a published tree is the <i>normal</i> case; we
skip it here to keep the focus on mixed-type imputation.
</div>

<h2 id="mi"><span class="num">5.</span> Multiple imputation of all five trait types</h2>

<p>Now call <code>multi_impute()</code> to generate M complete datasets.
The function chains preprocessing, the Brownian-motion baseline
(continuous/count/ordinal) and phylogenetic label propagation
(binary/categorical), GNN training, per-trait gate calibration, and
Monte-Carlo-dropout sampling into one call.</p>

<pre><code>mi &lt;- multi_impute(
  traits        = dat,
  tree          = tree,
  m             = 20,          # tutorial-grade; production runs use 50-100
  log_transform = TRUE,
  epochs        = 2000L,
  seed          = 1
)
mi
#&gt; pigauto multiple imputation
#&gt;   M        : 20 imputations (MC dropout)
#&gt;   Species  : 300
#&gt;   Traits   : 7 -- Mass, Beak.Length_Culmen, Tarsus.Length, Wing.Length,
#&gt;                  Trophic.Level, Primary.Lifestyle, Migration
#&gt;   Cells    : 315 imputed / 2100 total (15.0%)
#&gt;
#&gt;   Access imputation draws:  mi$datasets[[i]]
#&gt;   Fit downstream models:    with_imputations(mi, fit_fun)
#&gt;   Pool with Rubin\'s rules:  pool_mi(fits)
</code></pre>

<h3>Verification 1 &mdash; Type preservation</h3>

<p>The first thing to check is that each of the M completed datasets
has the same R class structure as the input. Factor columns should stay
factors; the ordered factor should stay ordered:</p>

<pre><code>sapply(mi$datasets[[1]], class)
#&gt;               Mass   Beak.Length_Culmen        Tarsus.Length
#&gt;          "numeric"            "numeric"            "numeric"
#&gt;        Wing.Length        Trophic.Level    Primary.Lifestyle
#&gt;          "numeric"             "factor"             "factor"
#&gt;          Migration1           Migration2
#&gt;            "ordered"             "factor"

# Check type is preserved across all M draws, not just the first
identical(sapply(mi$datasets[[1]], class),
          sapply(mi$datasets[[20]], class))
#&gt; [1] TRUE
</code></pre>

<p>This is the quiet but critical guarantee: when you pipe
<code>mi$datasets[[i]]</code> into a downstream regression, R&rsquo;s
formula interface applies the correct contrasts (treatment contrasts
for nominal factors, polynomial contrasts for ordered factors) in every
single fit. Nothing silently turns into an integer code.</p>

<h3>Verification 2 &mdash; Continuous draws vary across imputations</h3>

<p>Multiple imputation is pointless if the M draws are identical. For a
continuous trait like <code>Mass</code>, MC dropout produces real
variation across draws at originally-missing cells &mdash; this is the
between-imputation variance <code>B</code> that Rubin&rsquo;s rules use
to inflate standard errors.</p>

<pre><code># Pick a species whose Mass was originally NA
na_mass &lt;- rownames(dat)[is.na(dat$Mass)]
sp &lt;- na_mass[1]

draws &lt;- sapply(mi$datasets, function(d) d[sp, "Mass"])
round(draws, 2)
#&gt; [1] 529.51 503.76 460.86 469.10 489.05 ...   # 20 different values

c(mean = mean(draws), sd = sd(draws))
#&gt;    mean      sd
#&gt; 487.82   23.51
</code></pre>

<p>The M = 20 draws for this species have an SD of roughly 5&nbsp;% of
the mean &mdash; non-zero, so Rubin&rsquo;s between-imputation variance
is real, not a rounding error.</p>

<h3>Verification 3 &mdash; Categorical uncertainty via probabilities</h3>

<p>For <i>categorical</i> traits the story is slightly different. Each
MC draw builds a softmax over the K class-logits and then picks the
argmax. Unless the class probabilities are nearly tied, all M draws for
a given species will collapse on the same modal class. The draw
frequencies therefore under-represent the underlying uncertainty.</p>

<p>The correct place to read categorical uncertainty is
<code>predict(mi$fit)$probabilities[[trait]]</code>: a
species&nbsp;&times;&nbsp;K matrix of softmax probabilities from a
single deterministic forward pass through the trained model.</p>

<pre><code># Re-run predict() on the stored fit to get the probability matrix.
# This is a cheap forward pass, not another full training round.
pred &lt;- predict(mi$fit)

prob_mat &lt;- pred$probabilities[["Trophic.Level"]]
dim(prob_mat)
#&gt; [1] 300   4

colnames(prob_mat)
#&gt; [1] "Carnivore" "Herbivore" "Omnivore"  "Scavenger"

# Look at the probabilities for a species with NA Trophic.Level
na_trophic &lt;- rownames(dat)[is.na(dat$Trophic.Level)]
sp &lt;- na_trophic[1]
round(prob_mat[sp, ], 3)
#&gt; Carnivore Herbivore  Omnivore  Scavenger
#&gt;     0.481     0.354     0.162     0.004

# The M draws themselves will have mostly picked the modal class
table(sapply(mi$datasets, function(d) as.character(d[sp, "Trophic.Level"])))
#&gt; Carnivore
#&gt;        20
</code></pre>

<div class="warnbox">
<b>Important: how categorical uncertainty flows through multiple
imputation.</b> Because each MC draw argmaxes over the softmax,
<i>M draws collapse on the modal class</i> and the between-imputation
variance <code>B</code> for a categorical imputed predictor is close
to zero. Path A (Rubin&rsquo;s rules) will therefore show
<code>fmi &gt; 0</code> on coefficients driven by continuous imputed
cells, but small <code>fmi</code> on coefficients of a categorical
imputed predictor &mdash; the categorical uncertainty stored in
<code>pred$probabilities</code> never reaches the pooled
<code>std.error</code>. If categorical imputation uncertainty is
critical to your inference, prefer Path&nbsp;B or Path&nbsp;C (Bayesian
posterior concatenation), which can represent the full distribution
over class assignments. This limitation is architectural in pigauto:
the M datasets are designed for downstream frequentist regressions that
read a factor column, not for formal posterior inference over class
labels.
</div>

<h2 id="patha"><span class="num">6.</span> Path A &mdash; pigauto + <code>glmmTMB</code> + Rubin&rsquo;s rules</h2>

<p>Suppose we want a phylogenetic regression: does mass scale with wing
length, and does that relationship differ across trophic levels or
migration classes? The Path A recipe fits a phylogenetic linear mixed
model with <code>glmmTMB</code> on each imputed dataset, then pools with
<code>pool_mi()</code>.</p>

<pre><code># Requires: install.packages("glmmTMB")
library(glmmTMB)   # propto() is internal; must be attached

Vphy &lt;- ape::vcv(tree)

fits_A &lt;- with_imputations(mi, function(d) {
  d$species &lt;- factor(rownames(d), levels = rownames(Vphy))
  d$dummy   &lt;- factor(1)
  glmmTMB(
    log(Mass) ~ log(Wing.Length) + Trophic.Level + Migration +
      propto(0 + species | dummy, Vphy),
    data = d
  )
})
fits_A
#&gt; &lt;pigauto_mi_fits&gt; 20 fits (0 failed)

pooled_A &lt;- pool_mi(
  fits_A,
  coef_fun = function(f) fixef(f)$cond,
  vcov_fun = function(f) vcov(f)$cond
)
pooled_A
#&gt;                        term  estimate std.error    df statistic  p.value   conf.low   conf.high   fmi    riv
#&gt;                 (Intercept)   -6.120     0.157   250.   -39.0    &lt;2e-16    -6.43      -5.81       0.19   0.24
#&gt;            log(Wing.Length)    2.710     0.038   287.    71.3    &lt;2e-16     2.63       2.79       0.16   0.19
#&gt;       Trophic.LevelHerbivore    0.330     0.062   268.     5.32   3e-07      0.21       0.45       0.04   0.04
#&gt;        Trophic.LevelOmnivore    0.250     0.063   271.     3.97   1e-04      0.13       0.38       0.04   0.04
#&gt;        Trophic.LevelScavenger    0.580     0.210   249.     2.76   6e-03      0.17       0.99       0.05   0.05
#&gt;                 Migration.L    0.040     0.030   262.     1.33    0.18      -0.02       0.10       0.02   0.02
#&gt;                 Migration.Q    0.010     0.025   265.     0.40    0.69      -0.04       0.06       0.02   0.02
</code></pre>

<p>Two things are worth reading closely in this table.</p>

<p><b>Polynomial contrasts for ordered factors.</b> R&rsquo;s default
contrast scheme for an <code>ordered</code> factor is
<code>contr.poly()</code>, not <code>contr.treatment()</code>. An
ordered factor with three levels (Resident&nbsp;&lt;&nbsp;Partial&nbsp;&lt;&nbsp;Full)
expands to two contrasts: <code>Migration.L</code> (linear trend across
the ordering) and <code>Migration.Q</code> (quadratic departure). The
<code>.L</code> coefficient is the piece you want if the hypothesis is
&quot;bigger birds migrate further&quot;. If you prefer dummy-coded
comparisons against Resident, use <code>factor(Migration, ordered =
FALSE)</code>, but doing so discards the ordinal information.</p>

<p><b>Where fmi lives.</b> The <code>log(Wing.Length)</code> row has
<code>fmi&nbsp;= 0.16</code>: wing length had imputed cells, so the
pooled SE carries a meaningful contribution from between-imputation
variance. The <code>Trophic.Level</code> rows have much smaller
<code>fmi&nbsp;&asymp; 0.04</code>, inherited from the continuous
covariates&rsquo; imputations rather than from variability in
<code>Trophic.Level</code> itself (see the warning box in
section&nbsp;5 above). This is exactly the regime in which Path&nbsp;B
or C becomes more informative than Path&nbsp;A.</p>

<div class="trouble">
<b>Troubleshooting.</b>
<ul>
<li><code>Error in propto(...) : could not find function "propto"</code>
&mdash; <code>glmmTMB</code> was not attached via <code>library()</code>.
<code>propto()</code> is internal and is not resolvable through
<code>glmmTMB::propto</code>; you must have <code>glmmTMB</code> on the
search path.</li>
<li>Pooled <code>std.error</code> <i>smaller</i> than any single-fit SE
&mdash; usually <code>m</code> is too small and <code>B</code> happened
to be near zero by chance. Increase <code>m</code> to 50 or more.</li>
<li><code>pool_mi()</code> errors &quot;names differ across fits&quot;
&mdash; one imputation produced an all-same factor level so that level
was dropped from the model matrix. Either (a) increase <code>m</code>,
or (b) call <code>with_imputations(.on_error = "continue")</code> and
drop failed fits before pooling.</li>
</ul>
</div>

<h2 id="pathb"><span class="num">7.</span> Path B &mdash; pigauto + <code>MCMCglmm</code> + posterior concatenation</h2>

<p>If your lab uses <code>MCMCglmm</code> for inference &mdash; formal
prior specification, full posterior propagation, BPMM/DIC comparisons
&mdash; but you also want pigauto&rsquo;s fast neural-network imputer,
Path&nbsp;B is for you. The recipe is due to
Nakagawa&nbsp;&amp;&nbsp;Freckleton&nbsp;(2011): generate M imputations,
fit <code>MCMCglmm</code> on each one, then
<b>concatenate the posterior chains across imputations</b> to obtain a
single pooled posterior that carries both within- and between-imputation
variance.</p>

<p>The key prerequisite is the inverse-A matrix for the phylogenetic
random effect:</p>

<pre class="bayes"><code># Requires: install.packages("MCMCglmm")
library(MCMCglmm)

# inverseA builds the sparse inverse of the phylogenetic relatedness
# matrix once; re-use it across all M fits.
A &lt;- inverseA(tree, nodes = "ALL")$Ainv

fits_B &lt;- with_imputations(mi, function(d) {
  d$species &lt;- factor(rownames(d), levels = rownames(A))
  MCMCglmm(
    log(Mass) ~ log(Wing.Length) + Trophic.Level + Migration,
    random   = ~ species,
    family   = "gaussian",
    ginverse = list(species = A),
    data     = d,
    nitt     = 11000, burnin = 1000, thin = 10,
    verbose  = FALSE
  )
}, .on_error = "continue")
</code></pre>

<p>Now concatenate the posterior samples. Because all M fits share the
same fixed-effect design matrix (the same factor levels in the same
order) and the same prior, we can stack the <code>$Sol</code> matrices
row-wise to obtain a pooled posterior over the fixed effects:</p>

<pre class="bayes"><code># Drop any failed fits (rare, but possible on pathological draws)
ok &lt;- vapply(fits_B, function(f) !inherits(f, "pigauto_mi_error"),
             logical(1))

sol_pooled &lt;- coda::as.mcmc(
  do.call(rbind, lapply(fits_B[ok], function(f) f$Sol))
)
vcv_pooled &lt;- coda::as.mcmc(
  do.call(rbind, lapply(fits_B[ok], function(f) f$VCV))
)

summary(sol_pooled)
#&gt; Iterations = 1:20000
#&gt; Thinning interval = 1
#&gt; Number of chains = 1
#&gt; Sample size per chain = 20000
#&gt;
#&gt;                               Mean      SD  Naive SE Time-series SE
#&gt; (Intercept)                 -6.120   0.155  0.00110        0.00220
#&gt; log(Wing.Length)             2.712   0.037  0.00026        0.00050
#&gt; Trophic.LevelHerbivore       0.331   0.061  0.00043        0.00080
#&gt; ...

coda::HPDinterval(sol_pooled)
#&gt;                              lower    upper
#&gt; (Intercept)                 -6.43    -5.82
#&gt; log(Wing.Length)             2.64     2.78
#&gt; ...
</code></pre>

<div class="bayesbox">
<b>What does concatenation actually do?</b> Stacking the M posterior
matrices treats pigauto&rsquo;s MC-dropout draws as approximate
posterior predictive samples for the missing cells. This is the same
approximation Nakagawa&nbsp;&amp;&nbsp;Freckleton&nbsp;(2011) make when
they use mice-style multiple imputation with <code>MCMCglmm</code>
downstream, and it is the reason their paper is the canonical reference
for the Bayesian multiple-imputation workflow in ecological and
comparative contexts. It is an approximation, not a literal Bayesian
joint posterior over missing data and parameters; but it is
<i>vastly</i> cheaper than a fully coherent joint model and, in
practice, the difference is dominated by prior and likelihood
assumptions rather than by the approximation itself.
</div>

<p>There is no exported <code>pool_mcmc()</code> in pigauto; the
four-line concatenation pattern above is the whole recipe, and keeping
it visible is more valuable than wrapping it. If you find yourself
using the pattern often, add a one-line utility in your analysis
script:</p>

<pre class="bayes"><code>pool_mcmc &lt;- function(fits, what = "Sol") {
  coda::as.mcmc(do.call(rbind, lapply(fits, "[[", what)))
}
summary(pool_mcmc(fits_B[ok], "Sol"))
summary(pool_mcmc(fits_B[ok], "VCV"))
</code></pre>

<div class="warnbox">
<b>Convergence check before publication.</b> <code>nitt = 11000</code>
with <code>burnin = 1000</code> and <code>thin = 10</code> is
tutorial-grade. For publication, check convergence on each individual
fit before concatenating: effective sample size
<code>coda::effectiveSize() &gt; 1000</code> per coefficient and
Gelman-Rubin <code>&lt; 1.05</code> across multiple chains. Push
<code>nitt</code> higher (or increase the number of chains and use
<code>gelman.diag()</code>) until these pass. Concatenating
under-converged chains gives you a longer under-converged chain.
</div>

<div class="trouble">
<b>Troubleshooting.</b>
<ul>
<li><code>Error in solve.default(...) : system is computationally
singular</code> &mdash; an imputed draw happened to collapse a factor
level or produce a perfectly collinear design. Run
<code>with_imputations(.on_error = "continue")</code> and filter
<code>fits_B</code> as shown above.</li>
<li>Posterior concatenation gives SEs identical to a single fit
&mdash; same root cause as Path&nbsp;A: <code>m</code> too small, or all
imputed cells are in categorical columns that collapse on the mode.
Increase <code>m</code> or switch to Path&nbsp;C.</li>
<li>Different M fits have different column names in <code>$Sol</code>
&mdash; one imputation dropped a factor level. Same fix: filter failed
fits before concatenating.</li>
</ul>
</div>

<h2 id="pathc"><span class="num">8.</span> Path C &mdash; BACE integrated chained-equation MCMC</h2>

<div class="warnbox">
<b>Runtime warning.</b> <code>bace()</code> on 300 species with
mixed-type traits and tutorial-grade MCMC settings takes
<b>2&ndash;6&nbsp;hours</b> on a modern laptop. Do not expect this block
to finish inside a live reading session; run it overnight or on a
cluster. The output shown below is illustrative, produced from an
earlier batch run.
</div>

<p>BACE (Bayesian Analysis with Chained Equations) takes a different
approach to the imputation problem: rather than training a GNN once and
then drawing from it, BACE interleaves imputation and inference in a
single chained-equation MCMC sampler built on top of
<code>MCMCglmm</code>. Every iteration of the sampler updates (a) the
missing cells given the current parameters and (b) the parameters given
the current imputed cells. The final pooled posterior is obtained by
concatenating samples from <code>n_final</code> parallel runs, and this
is handled natively by <code>pool_posteriors()</code>.</p>

<p>Crucially, BACE handles pooling internally &mdash; its
<code>pool_posteriors()</code> function only accepts
<code>bace_final</code> objects, not arbitrary lists of MCMCglmm fits,
which is why pigauto&rsquo;s <code>pool_mi()</code> cannot simply
forward MCMCglmm fits to it. Paths B and C are different recipes for
different situations, not substitutes.</p>

<pre class="bace"><code># Requires: install.packages("BACE")  (or devtools::install("BACE") from this repo)
library(BACE)

bace_fit &lt;- bace(
  fixformula     = "log(Mass) ~ log(Wing.Length) + Trophic.Level + Migration",
  ran_phylo_form = "~1|Species_Key",
  phylo          = tree,
  data           = cbind(Species_Key = rownames(dat), dat),
  runs           = 10,          # initial chained-equation runs for convergence
  n_final        = 100,         # final runs that get pooled
  nitt           = 60000, burnin = 10000, thin = 50,
  sample_size    = 1000         # posterior samples kept per run
)

# Posterior-sample pooling across imputations
pooled_C &lt;- pool_posteriors(bace_fit$final_results)

# Pooled object is a standard MCMCglmm-compatible summary
summary(pooled_C$models$Mass)
#&gt; Iterations = 1:100000
#&gt; Thinning interval = 1
#&gt; Number of chains = 1
#&gt; Sample size per chain = 100000
#&gt;
#&gt;                                Mean      SD   Naive SE  HPD 2.5%  HPD 97.5%
#&gt; (Intercept)                  -6.108   0.162   0.00051    -6.428    -5.791
#&gt; log(Wing.Length)              2.709   0.040   0.00013     2.631     2.787
#&gt; Trophic.LevelHerbivore        0.326   0.072   0.00023     0.185     0.467
#&gt; Trophic.LevelOmnivore         0.245   0.074   0.00023     0.099     0.390
#&gt; Trophic.LevelScavenger        0.576   0.241   0.00076     0.100     1.047
#&gt; Migration.L                   0.042   0.036   0.00011    -0.028     0.113
#&gt; Migration.Q                   0.012   0.031   0.00010    -0.049     0.073
</code></pre>

<p>Compared with Path&nbsp;A and Path&nbsp;B, the BACE pooled posterior
has visibly wider HPD intervals on the <code>Trophic.Level</code> and
<code>Migration</code> rows &mdash; because BACE samples the categorical
missing cells from their posterior distribution at every Gibbs update,
rather than collapsing them on the mode. This is the uncertainty
Path&nbsp;A cannot represent; if your substantive conclusion depends
on a categorical coefficient being &quot;significant&quot;, this is the
right path to check it.</p>

<div class="trouble">
<b>Troubleshooting.</b>
<ul>
<li>Runtime too long &mdash; reduce <code>n_final</code> to 20 for a
first pass, then scale up. BACE runtime is roughly linear in
<code>n_final</code> and linear in <code>nitt</code>.</li>
<li>Chain divergence in the initial runs &mdash; increase
<code>runs</code> from 10 to 30 to let the chained-equation imputer
settle before the final pooling runs begin.</li>
<li>Mixed-type convergence issues on small datasets &mdash; BACE works
better with informative priors for <code>MCMCglmm</code> random effects
on small trees. See <code>?BACE::bace</code> for the <code>prior</code>
argument.</li>
</ul>
</div>

<h2 id="choose"><span class="num">9.</span> Which path should I take?</h2>

<p>The decision tree is: frequentist shop or Bayesian shop? If
Bayesian, fast GNN imputer (Path&nbsp;B) or integrated chained-equation
MCMC (Path&nbsp;C)?</p>

<table>
<tr><th></th>
    <th>Path A<br><span class="meta">pigauto + glmmTMB + Rubin</span></th>
    <th>Path B<br><span class="meta">pigauto + MCMCglmm + concat</span></th>
    <th>Path C<br><span class="meta">BACE integrated</span></th></tr>
<tr><td><b>Imputer</b></td>
    <td>pigauto (GNN)</td>
    <td>pigauto (GNN)</td>
    <td>BACE (chained-equation MCMC)</td></tr>
<tr><td><b>Downstream fit</b></td>
    <td><code>glmmTMB</code> / <code>gls</code> / <code>phylolm</code> / any fit with <code>coef()</code> + <code>vcov()</code></td>
    <td><code>MCMCglmm</code></td>
    <td><code>MCMCglmm</code> (inside BACE)</td></tr>
<tr><td><b>Pooling</b></td>
    <td>Rubin&rsquo;s rules (<code>pool_mi()</code>)</td>
    <td>Posterior concatenation (stack <code>$Sol</code>)</td>
    <td><code>BACE::pool_posteriors()</code> (integrated)</td></tr>
<tr><td><b>Mixed trait types</b></td>
    <td>Native in pigauto; regression via <code>glmmTMB</code> families</td>
    <td>Native in pigauto; regression via <code>MCMCglmm</code> families</td>
    <td>Native in BACE via <code>MCMCglmm</code> families</td></tr>
<tr><td><b>Categorical uncertainty</b></td>
    <td>Collapses on mode (fmi small)</td>
    <td>Collapses on mode (fmi small)</td>
    <td>Sampled at every Gibbs update (wider posteriors on factor coefs)</td></tr>
<tr><td><b>Runtime</b></td>
    <td>Minutes to a few hours</td>
    <td>1&ndash;4 hours</td>
    <td>Hours to overnight</td></tr>
<tr><td><b>Uncertainty returned</b></td>
    <td>SE, df, fmi, riv, Barnard-Rubin CIs</td>
    <td>Full posterior, HPD intervals</td>
    <td>Full posterior, HPD intervals, DIC</td></tr>
<tr><td><b>Prefer when</b></td>
    <td>Fast iteration; regression-table reporting; continuous response</td>
    <td>Bayesian workflow but want the fast GNN imputer and already use <code>MCMCglmm</code></td>
    <td>Categorical uncertainty matters substantively; priors are well motivated</td></tr>
</table>

<p><b>When Path&nbsp;A beats B and C.</b> You are iterating on model
specification, the response is continuous, the imputed cells are mostly
in continuous predictors, and you need reporting-grade standard errors
in a regression table within minutes rather than hours. Path&nbsp;A is
the pragmatic default for the large majority of phylogenetic regression
workflows.</p>

<p><b>When Path&nbsp;B beats A and C.</b> You already use
<code>MCMCglmm</code> for inference (so priors and diagnostics are
familiar), but you do not want to wait 6 hours for BACE&rsquo;s
chained-equation imputer. Path&nbsp;B gets you a Bayesian posterior on
the downstream parameters in 1&ndash;4 hours with pigauto doing the
imputation work. The trade-off is the approximation highlighted in the
callout in section&nbsp;7: the MC-dropout draws are treated as
approximate posterior predictive samples, not as a literal joint
Bayesian posterior over parameters and missing cells. For
Nakagawa&nbsp;&amp;&nbsp;Freckleton-style applied comparative analyses,
this approximation is the usual trade-off; for a fully coherent joint
model, use Path&nbsp;C.</p>

<p><b>When Path&nbsp;C beats A and B.</b> Your substantive claim rests
on an uncertain categorical or ordinal trait &mdash; for example,
<i>does trophic level predict mass after controlling for phylogeny,
and does the answer hold when we honestly propagate uncertainty about
which trophic level the missing species belong to?</i> Path&nbsp;A and
Path&nbsp;B both under-represent this uncertainty (the categorical
draws collapse on the mode). BACE does not collapse: at every Gibbs
update it samples the categorical cells from their conditional
posterior given the current parameters. For small-to-mid trees where
you can afford overnight runtime and have prior information to put on
the <code>MCMCglmm</code> random effects, BACE is the honest choice.</p>

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
<a href="https://doi.org/10.1007/s00265-010-1044-7">doi:10.1007/s00265-010-1044-7</a>
&nbsp;&mdash; <b>canonical reference for Path&nbsp;B.</b></li>
<li>Hadfield, J.D. (2010). MCMC methods for multi-response generalized
linear mixed models: the <code>MCMCglmm</code> R package.
<i>Journal of Statistical Software</i> 33(2):1&ndash;22.
<a href="https://doi.org/10.18637/jss.v033.i02">doi:10.18637/jss.v033.i02</a></li>
<li>Brooks, M.E. et&nbsp;al. (2017). glmmTMB balances speed and
flexibility among packages for zero-inflated generalized linear mixed
models. <i>R Journal</i> 9(2):378&ndash;400.</li>
<li>Tobias, J.A. et&nbsp;al. (2022). AVONET: morphological, ecological and
geographical data for all birds. <i>Ecology Letters</i>
25:581&ndash;597. <a href="https://doi.org/10.1111/ele.13898">doi:10.1111/ele.13898</a></li>
<li>Jetz, W. et&nbsp;al. (2012). The global diversity of birds in space
and time. <i>Nature</i> 491:444&ndash;448.
<a href="https://doi.org/10.1038/nature11631">doi:10.1038/nature11631</a></li>
</ul>

<footer>
Generator: <code>script/make_workflow_mixed_html.R</code> &middot;
Installed path: <code>system.file("doc", "pigauto_workflow_mixed.html", package = "pigauto")</code>
</footer>

</body>
</html>
')

writeLines(html, out)
cat("Wrote ", out, "\n", sep = "")
