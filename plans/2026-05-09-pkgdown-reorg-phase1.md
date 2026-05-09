# pkgdown reorganisation Phase 1 — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reorganise the pigauto pkgdown site so a Bhavya-class new user (priority audience A from the spec) lands on a clean Articles menu, finds Common Pitfalls in one click, and can reach methodology evidence under its own dropdown — without changing any R code.

**Architecture:** Three independent units in one PR off `main`. Unit 1 rewrites `_pkgdown.yml` (navbar slim + new Methodology dropdown + DOCS.md home-link removed). Unit 2 ships `vignettes/common-pitfalls.Rmd` covering five pitfalls with a fixed Symptom/Why/Diagnose/Fix/See-also template. Unit 3 retires `DOCS.md` and adds one live-site link line to `README.md`. Plus a `NEWS.md` entry.

**Tech Stack:** R 4.5, pkgdown 2.x, knitr/rmarkdown, devtools, no new dependencies. Bundled data (`avonet300`, `tree300`) for all vignette examples.

**Spec:** `specs/2026-05-09-pkgdown-reorg-phase1-design.md`

---

## Task 0: Verify the prerequisite

**Files:** none (read-only check).

The Common Pitfalls vignette quotes language from PR #69's `?impute` "What gets imputed (read this first)" section. PR #69 (`docs/clarify-imputed-vs-observed`) **must** be merged to `main` before this work starts; otherwise the vignette's "See also" links and re-used phrasing will reference docstring text that doesn't exist on `main`.

- [ ] **Step 1: Check PR #69 state**

Run: `gh pr view 69 --json state,mergedAt`
Expected: `"state":"MERGED"` and a non-null `mergedAt`. If `OPEN`, stop and merge PR #69 first.

- [ ] **Step 2: Sync local main**

Run: `cd /Users/z3437171/Dropbox/Github\ Local/pigauto && git fetch origin main && git checkout main && git pull --ff-only`
Expected: HEAD includes the PR #69 merge commit.

- [ ] **Step 3: Confirm `R/impute.R` has the new section**

Run: `grep -n "What gets imputed (read this first)" R/impute.R`
Expected: One match. If zero matches, PR #69 has not landed; stop.

---

## Task 1: Branch + pre-build verify the 23 benchmark HTMLs

**Files:** none (verification + branch creation).

The Methodology dropdown will link to 23 existing HTMLs under `docs/dev/`. If any of them are missing on a fresh clone, `pkgdown::build_site_local()` will succeed but the navbar links will 404. Catch that now.

- [ ] **Step 1: Create the working branch**

Run: `git checkout -b docs/pkgdown-reorg-phase1 origin/main`
Expected: `Switched to a new branch 'docs/pkgdown-reorg-phase1'`.

- [ ] **Step 2: Verify all 23 benchmark HTMLs exist**

Run:

```bash
for f in \
  bench_continuous bench_binary bench_ordinal bench_count bench_categorical \
  bench_proportion bench_zi_count bench_multi_proportion \
  bench_scaling_v090 bench_avonet_missingness bench_missingness_mechanism \
  bench_tree_uncertainty pantheria_summary bench_delhey bench_multi_obs \
  phase8_summary bench_signal_sweep bench_correlation_sweep \
  bench_evo_model_sweep bench_clade_missingness bench_covariate_sim \
  bench_bace_avonet_head_to_head bench_pantheria_bace_head_to_head; do
  if [[ ! -f "docs/dev/$f.html" ]]; then echo "MISSING: docs/dev/$f.html"; fi
done
```

Expected: no output (silent = all 23 exist). If any are missing, run the corresponding `script/make_*_html.R` driver to regenerate before continuing. Do not proceed to Task 2 until all 23 are present.

---

## Task 2: Navbar redesign in `_pkgdown.yml`

**Files:**
- Modify: `_pkgdown.yml` (the entire `navbar:` block lines ~38–105 in current file; plus removing one entry from `home.links` lines ~28–32)

The full navbar block has three changes: (a) add `methodology` to `navbar.structure.left`, (b) collapse Articles to 8 items, (c) add a new `methodology` component with 23 working items + a Phase-2 placeholder comment. Plus delete the `home.links` entry pointing at `DOCS.md`.

- [ ] **Step 1: Replace the navbar block**

Replace the entire `navbar:` block in `_pkgdown.yml` with:

```yaml
navbar:
  structure:
    left:  [intro, articles, methodology, reference, news]
    right: [search, github]
  components:
    intro:
      text: Get started
      href: articles/getting-started.html
    articles:
      text: Articles
      menu:
        - text: "Knitted vignettes"
        - text: "Get started"
          href: articles/getting-started.html
        - text: "Common pitfalls / FAQ"
          href: articles/common-pitfalls.html
        - text: "Mixed-type traits"
          href: articles/mixed-types.html
        - text: "Propagating tree uncertainty"
          href: articles/tree-uncertainty.html
        - text: "-------"
        - text: "Walk-throughs"
        - text: "Architecture overview"
          href: pigauto_intro.html
        - text: "Mixed-type workflow (Paths A / B / C)"
          href: pigauto_workflow_mixed.html
        - text: "Comparative study with covariates"
          href: pigauto_walkthrough_covariates.html
        - text: "Multi-observation per species"
          href: pigauto_walkthrough_multi_obs.html
    methodology:
      text: Methodology
      menu:
        - text: "Per-trait benches"
        - text: "Continuous (BM / OU / regime / nonlinear)"
          href: dev/bench_continuous.html
        - text: "Binary (signal strength sweep)"
          href: dev/bench_binary.html
        - text: "Ordinal (level count sweep)"
          href: dev/bench_ordinal.html
        - text: "Count (Poisson / NegBin)"
          href: dev/bench_count.html
        - text: "Categorical (K-level sweep)"
          href: dev/bench_categorical.html
        - text: "Proportion (signal sweep)"
          href: dev/bench_proportion.html
        - text: "Zero-inflated counts (zero fraction sweep)"
          href: dev/bench_zi_count.html
        - text: "Multi-proportion / compositional"
          href: dev/bench_multi_proportion.html
        - text: "-------"
        - text: "Cross-dataset benches"
        - text: "Scaling up to ~10,000 species (AVONET full)"
          href: dev/bench_scaling_v090.html
        - text: "AVONET missingness sweep (20 / 50 / 80 %)"
          href: dev/bench_avonet_missingness.html
        - text: "Missingness mechanisms (MCAR / MAR / MNAR)"
          href: dev/bench_missingness_mechanism.html
        - text: "Tree uncertainty (Rubin's rules, 10 posterior trees)"
          href: dev/bench_tree_uncertainty.html
        - text: "PanTHERIA mammals"
          href: dev/pantheria_summary.html
        - text: "Delhey 5,809-species birds"
          href: dev/bench_delhey.html
        - text: "Multi-observation imputation (CTmax-like)"
          href: dev/bench_multi_obs.html
        - text: "-------"
        - text: "Phase 8 simulations + head-to-heads"
        - text: "Phase 8 summary"
          href: dev/phase8_summary.html
        - text: "Phylogenetic signal sweep (Pagel's λ)"
          href: dev/bench_signal_sweep.html
        - text: "Cross-trait correlation sweep"
          href: dev/bench_correlation_sweep.html
        - text: "Evolutionary-model sweep (BM / OU / regime / nonlinear)"
          href: dev/bench_evo_model_sweep.html
        - text: "Clade-correlated missingness (realistic MAR)"
          href: dev/bench_clade_missingness.html
        - text: "Covariate effectiveness (simulation sweep)"
          href: dev/bench_covariate_sim.html
        - text: "AVONET 300 head-to-head (pigauto vs BACE)"
          href: dev/bench_bace_avonet_head_to_head.html
        - text: "PanTHERIA head-to-head (pigauto vs BACE)"
          href: dev/bench_pantheria_bace_head_to_head.html
        # -------- Phase 2 placeholder --------
        # When Phase 2 ships, add an entry here for the math vignette:
        #   - text: "Mathematical & algorithmic description"
        #     href: articles/methodology.html
```

- [ ] **Step 2: Remove the DOCS.md entry from `home.links`**

In `_pkgdown.yml`, find the `home.links` block. Delete only the entry titled `Documentation index (DOCS.md)`. Keep the remaining two links (`Source code`, `Report a bug`). The expected result:

```yaml
  links:
    - text: Source code
      href: https://github.com/itchyshin/pigauto
    - text: Report a bug
      href: https://github.com/itchyshin/pigauto/issues
```

- [ ] **Step 3: Verify YAML still parses**

Run: `Rscript -e 'invisible(yaml::read_yaml("_pkgdown.yml"))'`
Expected: no error, no output. If a `Scanner error` or similar appears, indentation drifted; re-paste the block above with two-space indentation throughout.

- [ ] **Step 4: Build site locally**

Run: `Rscript -e 'pkgdown::build_site_local(preview = FALSE)'`
Expected: completes within ~30 s with no errors. Warnings about the missing `articles/common-pitfalls.html` are expected at this stage and resolve in Task 3.

- [ ] **Step 5: Commit**

```bash
git add _pkgdown.yml
git commit -m "docs(pkgdown): slim Articles dropdown, add Methodology dropdown, drop DOCS.md home link

Phase 1 of the pkgdown reorganisation (specs/2026-05-09-pkgdown-reorg-phase1-design.md):

- navbar.left now [intro, articles, methodology, reference, news]
- Articles dropdown trimmed to 8 first-class items (3 vignettes + Common
  Pitfalls + 4 walk-throughs); the existing 'Logo gallery' and 'Test
  catalogue' developer-internal entries are removed from the navbar (they
  remain accessible via build outputs)
- New Methodology dropdown owns the 23 benchmark HTMLs in three sub-sections:
  per-trait (8), cross-dataset (7), Phase 8 sims + head-to-heads (8)
- Phase-2 placeholder comment in YAML reserves a slot for the math vignette
- home.links 'Documentation index (DOCS.md)' entry removed (DOCS.md retired
  in a follow-up commit)
- The Common Pitfalls Articles entry will 404 until Task 3 of this plan
  ships the vignette in the same PR
"
```

---

## Task 3: Common Pitfalls vignette skeleton + YAML frontmatter

**Files:**
- Create: `vignettes/common-pitfalls.Rmd`

Create the vignette file with YAML frontmatter, opening prose, and an empty section per pitfall. Sections are filled in Tasks 4–8 (one per pitfall).

- [ ] **Step 1: Create the file with frontmatter + opening + 5 empty sections**

Write `vignettes/common-pitfalls.Rmd`:

```markdown
---
title: "Common pitfalls / FAQ"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Common pitfalls / FAQ}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6,
  fig.height = 4
)
```

This vignette collects the questions early users most often hit when running
pigauto on their own data. Each section follows the same template:

- **Symptom** — the surprising output, in user voice.
- **Why this happens** — the mechanism in 2–3 sentences.
- **Diagnose** — 1–3 R commands that confirm the pattern.
- **Fix** — concrete code change with a short explanation.
- **See also** — links to the relevant `?function`, design memo, or
  Methodology bench under the Methodology navbar dropdown.

If you hit something that isn't here and feels surprising, please
[open an issue](https://github.com/itchyshin/pigauto/issues/new) — most
of the items below were added because real users tripped on them.

## 1. "I called `impute()` and `result$prediction$imputed` looks like my input"

*(written in Task 4)*

## 2. "My ordinal trait predicted 100 % majority class"

*(written in Task 5)*

## 3. "The gate stays closed and the GNN seems to do nothing"

*(written in Task 6)*

## 4. "How do I know if my dataset has enough phylogenetic signal?"

*(written in Task 7)*

## 5. "Predictions are way bigger than anything I observed"

*(written in Task 8)*
```

- [ ] **Step 2: Verify the skeleton knits**

Run: `Rscript -e 'rmarkdown::render("vignettes/common-pitfalls.Rmd", quiet = TRUE)'`
Expected: completes in <5 s; produces `vignettes/common-pitfalls.html`. The HTML has an "(written in Task N)" placeholder under each section header — that's correct at this stage.

- [ ] **Step 3: Commit the skeleton**

```bash
git add vignettes/common-pitfalls.Rmd
git commit -m "docs(vignette): scaffold Common pitfalls / FAQ vignette

Skeleton with YAML frontmatter, fixed Symptom/Why/Diagnose/Fix/See-also
template description, and 5 empty sections (one per pitfall).  Sections
are filled in subsequent commits.
"
```

---

## Task 4: Pitfall section 1 — `prediction$imputed` semantics

**Files:**
- Modify: `vignettes/common-pitfalls.Rmd:34-36` (replace the section-1 placeholder)

Address Bhavya issue #67. Reuses the language from PR #69's `?impute` "What gets imputed (read this first)" section. The runnable example masks 30 cells of `avonet300$Mass` and shows truth vs imputed.

- [ ] **Step 1: Replace the section-1 placeholder with the worked content**

Replace the line `## 1. "I called \`impute()\` and \`result$prediction$imputed\` looks like my input"` and its placeholder line `*(written in Task 4)*` with:

````markdown
## 1. "I called `impute()` and `result$prediction$imputed` looks like my input"

**Symptom.** You run `result <- impute(df, tree)` on a fully-observed
dataset (e.g. the bundled `avonet300`) and read
`result$prediction$imputed$Mass` expecting "the imputed values" — but
the values look exactly like your input data, including legitimately
huge ones (a 25 kg rhea, a 12 kg vulture).

**Why this happens.** `impute()` only *imputes* cells that are `NA`
in the input. Your input was fully observed, so nothing was imputed:
`result$completed` equals the input, `sum(result$imputed_mask)` is
zero, and `result$prediction$imputed` contains the model's prediction
for **every** cell — observed and missing alike. For observed cells,
the well-calibrated gate keeps the prediction close to the input
value, so what comes back is essentially the original data passed
through. The slot is intended for diagnostics (checking calibration on
training cells), not as the imputed-values output.

**Diagnose.**

```{r section1-diag, eval = FALSE}
library(pigauto)
data(avonet300, tree300)
df <- avonet300
rownames(df) <- df$Species_Key
df$Species_Key <- NULL
sum(is.na(df))                  # if 0, there's nothing for impute() to do
```

**Fix.** Mask some cells before calling `impute()`, then evaluate
predictions only on the held-out cells:

```{r section1-fix, eval = FALSE}
set.seed(1L)
hide   <- sample(which(!is.na(df$Mass)), 30L)
df_obs <- df
df_obs$Mass[hide] <- NA          # hide 30 mass values

result <- impute(df_obs, tree300)

result$completed$Mass[hide]      # pigauto's imputations
df$Mass[hide]                    # held-out truth, for comparison
sum(result$imputed_mask[, "Mass"])  # 30
```

For your own data with real `NA`s, the imputed values you actually
care about are `result$completed[result$imputed_mask]`, not
`result$prediction$imputed`.

**See also.** `?impute` ("What gets imputed (read this first)"),
[issue #67](https://github.com/itchyshin/pigauto/issues/67).
````

- [ ] **Step 2: Re-knit to verify Section 1 renders**

Run: `Rscript -e 'rmarkdown::render("vignettes/common-pitfalls.Rmd", quiet = TRUE)'`
Expected: completes in <10 s; HTML now shows Section 1 fully rendered with code blocks and `eval = FALSE` chunks shown but not executed.

- [ ] **Step 3: Commit**

```bash
git add vignettes/common-pitfalls.Rmd
git commit -m "docs(vignette): pitfall 1 — prediction\$imputed vs completed semantics

Mirrors PR #69's ?impute 'What gets imputed (read this first)' section.
Worked example on avonet300: hide 30 Mass cells, impute, compare to truth.
Closes the documentation gap behind issue #67.
"
```

---

## Task 5: Pitfall section 2 — K=3 ordinal majority-class collapse

**Files:**
- Modify: `vignettes/common-pitfalls.Rmd` (replace the section-2 placeholder)

Address Bhavya issue #68. Worked example on `avonet300$Migration`: default settings vs `n_imputations = 20L, pool_method = "mode"`.

- [ ] **Step 1: Replace the section-2 placeholder**

Replace `## 2. "My ordinal trait predicted 100 % majority class"` and its placeholder `*(written in Task 5)*` with:

````markdown
## 2. "My ordinal trait predicted 100 % majority class"

**Symptom.** You impute an ordinal trait and the prediction is the
majority class for every species. For example, on `avonet300$Migration`
(K = 3 ordinal: Resident / Partial / Full), 300/0/0.

**Why this happens.** Two things compound:

1. If your input has no `NA`s in that column, there's nothing to impute
   externally (see Pitfall 1) — `result$prediction$imputed$Migration`
   reflects the model's calibrated-gate output, not new imputations.
2. At default settings (`n_imputations = 1L`, `pool_method = "median"`),
   a small ordinal trait whose marginal distribution is heavily skewed
   (AVONET `Migration` is ~78 % Resident / 14 % Partial / 8 % Full at
   n = 300) can have its calibrated gate snap to a corner that returns
   the majority class for every species.

**Diagnose.**

```{r section2-diag, eval = FALSE}
library(pigauto)
data(avonet300, tree300)
df <- avonet300; rownames(df) <- df$Species_Key; df$Species_Key <- NULL

table(df$Migration)              # check the marginal distribution
result <- impute(df, tree300, verbose = FALSE)
table(result$prediction$imputed$Migration)
```

**Fix.** For imbalanced K-class ordinal traits, increase
`n_imputations` and switch to `pool_method = "mode"` (Phase H). On
the AVONET multi-seed bench this gave +6.6 percentage-point accuracy
on `Migration` (K = 3) versus the default median pool.

```{r section2-fix, eval = FALSE}
set.seed(1L)
hide   <- sample(which(!is.na(df$Migration)), 30L)
df_obs <- df
df_obs$Migration[hide] <- NA

# Default settings: prone to majority-class collapse on imbalanced K = 3
result <- impute(df_obs, tree300, verbose = FALSE)
table(result$completed$Migration[hide], df$Migration[hide])

# Recommended for K = 3 ordinal: more draws + mode pooling
result_mode <- impute(df_obs, tree300, n_imputations = 20L,
                      pool_method = "mode", verbose = FALSE)
table(result_mode$completed$Migration[hide], df$Migration[hide])
```

**See also.** `?impute` ("Imbalanced K-class traits"),
[Phase H memo](https://github.com/itchyshin/pigauto/blob/main/useful/MEMO_2026-05-01_phase_h_results.md),
[issue #68](https://github.com/itchyshin/pigauto/issues/68).
````

- [ ] **Step 2: Re-knit**

Run: `Rscript -e 'rmarkdown::render("vignettes/common-pitfalls.Rmd", quiet = TRUE)'`
Expected: <10 s; Section 2 renders with both diagnostic and fix code blocks.

- [ ] **Step 3: Commit**

```bash
git add vignettes/common-pitfalls.Rmd
git commit -m "docs(vignette): pitfall 2 — K=3 ordinal majority-class collapse

Documents the recommended n_imputations=20 + pool_method='mode' config
for imbalanced K=3 ordinal traits.  Phase H memo recorded +6.6 pp on
AVONET Migration; this vignette section makes that recipe discoverable.
Closes the documentation gap behind issue #68.
"
```

---

## Task 6: Pitfall section 3 — safety floor / closed gate

**Files:**
- Modify: `vignettes/common-pitfalls.Rmd` (replace the section-3 placeholder)

Explain why `r_cal_gnn` may end up close to zero on low-signal traits, and how to inspect it.

- [ ] **Step 1: Replace the section-3 placeholder**

Replace `## 3. "The gate stays closed and the GNN seems to do nothing"` and its placeholder with:

````markdown
## 3. "The gate stays closed and the GNN seems to do nothing"

**Symptom.** You expected the GNN to dominate, but inspecting the
fitted model shows the calibrated gate is fully or near-fully closed
(`r_cal_gnn` ≈ 0) — predictions equal the BM baseline.

**Why this happens.** This is the safety-floor design behaviour, not a
bug. After training, pigauto picks the per-latent-column gate that
minimises validation loss across the simplex
$r_\text{BM} \cdot \text{BM} + r_\text{GNN} \cdot \text{GNN} + r_\text{MEAN} \cdot \text{MEAN}$.
When the GNN cannot beat BM on the held-out validation set, the
optimum is `r_cal_gnn = 0` — so the calibrated prediction is
guaranteed never to be worse than BM (or grand mean, whichever is
better).  This is what the package was designed to do on
high-phylogenetic-signal traits where BM is already optimal.

**Diagnose.**

```{r section3-diag, eval = FALSE}
library(pigauto)
data(avonet300, tree300)
df <- avonet300; rownames(df) <- df$Species_Key; df$Species_Key <- NULL
fit <- impute(df, tree300, verbose = FALSE)$fit

# Per-latent-column calibrated gates (since v0.9.1.9002):
fit$r_cal_bm        # r assigned to the BM baseline
fit$r_cal_gnn       # r assigned to the GNN delta
fit$r_cal_mean      # r assigned to the grand mean
```

A row where `r_cal_gnn` is small (< 0.1) means the gate has effectively
closed for that latent column.

**Fix.** Often there is nothing to fix — the closed gate is
**evidence of high phylogenetic signal**, not a problem. If you suspect
the GNN should be helping (e.g. you've added covariates, or the trait
has known cross-trait structure) but the gate is closed:

- Check the validation set is not pathologically small (the
  "small validation set" warning during fitting is a red flag — see
  `?fit_pigauto` "Calibration at small n").
- Verify covariates are not all-NA or constant after preprocessing.
- For ordinal / categorical traits, see Pitfall 2 — the gate may be
  closing onto a majority-class corner that mode pooling resolves.

**See also.** `?fit_pigauto` ("Safety floor"), `?phylo_signal`,
[design spec](https://github.com/itchyshin/pigauto/blob/main/specs/2026-04-23-safety-floor-mean-gate-design.md).
````

- [ ] **Step 2: Re-knit and commit**

```bash
Rscript -e 'rmarkdown::render("vignettes/common-pitfalls.Rmd", quiet = TRUE)'
git add vignettes/common-pitfalls.Rmd
git commit -m "docs(vignette): pitfall 3 — closed gate is evidence, not a bug

Explains the safety-floor design: r_cal_gnn = 0 means BM beat the GNN
on validation, by construction.  Documents how to inspect the per-
latent-column calibrated gates (fit\$r_cal_bm / r_cal_gnn / r_cal_mean).
"
```

---

## Task 7: Pitfall section 4 — phylogenetic signal diagnosis

**Files:**
- Modify: `vignettes/common-pitfalls.Rmd` (replace the section-4 placeholder)

Show how to compute phylogenetic signal on a trait before deciding pigauto is the right tool.

- [ ] **Step 1: Replace the section-4 placeholder**

Replace `## 4. "How do I know if my dataset has enough phylogenetic signal?"` and its placeholder with:

````markdown
## 4. "How do I know if my dataset has enough phylogenetic signal?"

**Symptom.** You aren't sure whether pigauto's BM kriging baseline
will outperform a simple mean impute on your dataset.

**Why this matters.** pigauto's BM baseline buys you accuracy in
proportion to phylogenetic signal in the trait. At Pagel's λ ≈ 0
(no signal), BM kriging reduces to the species mean and pigauto won't
beat a simple mean baseline; at λ ≈ 1 (strong signal), BM kriging
materially outperforms the mean. The Phase 8 signal-strength sweep
(see the **Methodology → Phylogenetic signal sweep** dropdown) shows
the crossover empirically.

**Diagnose.** Estimate λ on the observed cells of each trait before
imputing:

```{r section4-diag, eval = FALSE}
library(pigauto)
data(avonet300, tree300)
df <- avonet300; rownames(df) <- df$Species_Key; df$Species_Key <- NULL

# Per-trait Pagel's lambda on observed values
phylo_signal(df, tree300)
```

The output reports lambda per continuous trait. Discrete traits use
the binary / categorical analogue.

**Fix.** Use the lambda estimate to set expectations:

- λ > 0.7: pigauto's BM baseline alone usually beats grand-mean.
  GNN may add little if the trait is well-conserved.
- 0.3 < λ < 0.7: pigauto's GNN typically helps on top of BM,
  especially if you have covariates.
- λ < 0.3: phylogenetic information is weak. Consider whether a
  phylogenetic imputation method is the right tool at all. A simple
  mean impute or a covariate-only regression may do as well.

**See also.** `?phylo_signal`,
[Phase 8 signal sweep](https://itchyshin.github.io/pigauto/dev/bench_signal_sweep.html)
under the Methodology dropdown.
````

- [ ] **Step 2: Re-knit and commit**

```bash
Rscript -e 'rmarkdown::render("vignettes/common-pitfalls.Rmd", quiet = TRUE)'
git add vignettes/common-pitfalls.Rmd
git commit -m "docs(vignette): pitfall 4 — phylo signal diagnosis before imputing

Points users at phylo_signal() to estimate Pagel's lambda on each trait
before assuming pigauto is the right tool; cross-references the Phase 8
signal sweep under the Methodology dropdown.
"
```

---

## Task 8: Pitfall section 5 — clamp_outliers tail safety

**Files:**
- Modify: `vignettes/common-pitfalls.Rmd` (replace the section-5 placeholder)

Show the AVONET Casuarius / FishBase weight outlier symptom and the `clamp_outliers = TRUE` fix.

- [ ] **Step 1: Replace the section-5 placeholder**

Replace `## 5. "Predictions are way bigger than anything I observed"` and its placeholder with:

````markdown
## 5. "Predictions are way bigger than anything I observed"

**Symptom.** A masked log-transformed continuous trait (body mass,
seed mass, fish weight) predicts a value 50–100× larger than anything
observed. On AVONET, the canonical case is the cassowary: truth ≈ 35 kg,
predicted up to ~540 kg.

**Why this happens.** For log-transformed traits, the GNN's
MC-dropout draws are on the log scale. A latent ~+3-4σ above the
training distribution survives as a ~50-100× value error after
`expm1()` back-transformation. With `n_imputations = 1`, a single
unlucky dropout pattern can produce this; with `pool_method = "median"`
(default) the median of M draws is robust to one bad draw, but a small
M (≤ 5) on the long tail of the latent distribution can still mis-pool.

**Diagnose.**

```{r section5-diag, eval = FALSE}
# After running impute(), check whether any imputation exceeds the
# observed maximum by an unrealistic factor:
predicted_mass <- result$completed$Mass[result$imputed_mask[, "Mass"]]
obs_max <- max(df$Mass, na.rm = TRUE)
sum(predicted_mass > 5 * obs_max)
```

A non-zero count is a signal of tail extrapolation.

**Fix.** Phase G `clamp_outliers = TRUE` caps post-back-transform
predictions for log-transformed continuous, count, and `zi_count`
magnitude traits at `obs_max * clamp_factor` (default 5). This is
opt-in because for legitimate growth-curve datasets where 5× the
observed maximum is plausible, you don't want it.

```{r section5-fix, eval = FALSE}
result <- impute(df_obs, tree300,
                 clamp_outliers = TRUE,
                 clamp_factor   = 5,    # Tukey-style outlier cap
                 verbose = FALSE)
```

**See also.** `?impute` (`clamp_outliers`, `clamp_factor` arguments),
[AVONET Mass diagnosis memo](https://github.com/itchyshin/pigauto/blob/main/useful/MEMO_2026-05-01_avonet_mass_diag.md),
[Phase G results](https://github.com/itchyshin/pigauto/blob/main/useful/MEMO_2026-05-01_phase_g_results.md).
````

- [ ] **Step 2: Re-knit and verify total knit time stays inside the budget**

Run: `time Rscript -e 'rmarkdown::render("vignettes/common-pitfalls.Rmd", quiet = TRUE)'`
Expected: total wall time ≤ 30 s. If it exceeds 30 s, the `eval = FALSE` chunks must already be inert; the only realistic culprit would be an unintended `eval = TRUE` chunk — search for `eval = TRUE` in the Rmd and re-mark to `eval = FALSE`. (At time of writing, every chunk in this vignette is `eval = FALSE`, so the budget is comfortable.)

- [ ] **Step 3: Commit**

```bash
git add vignettes/common-pitfalls.Rmd
git commit -m "docs(vignette): pitfall 5 — clamp_outliers tail safety

Documents the Phase G opt-in clamp for tail-extrapolated predictions
on log-transformed continuous, count, and zi_count magnitude traits.
Includes the AVONET Casuarius case as the canonical example.
"
```

---

## Task 9: Index consolidation — retire `DOCS.md`, link the live site from `README.md`

**Files:**
- Delete: `DOCS.md`
- Modify: `README.md` (add one line near the top)

- [ ] **Step 1: Delete `DOCS.md`**

Run: `git rm DOCS.md`
Expected: `rm 'DOCS.md'`

- [ ] **Step 2: Verify nothing in tracked source still references `DOCS.md`**

Run: `git grep -n DOCS.md`
Expected: no hits (we already removed the `home.links` entry in Task 2). If a hit appears, edit the offending file to remove the dead reference.

- [ ] **Step 3: Add the live-site link line to `README.md`**

Open `README.md`. Immediately after line 2 (`**Missing trait data should not stop a comparative analysis.**`), insert a blank line and then:

```markdown
> Live documentation: <https://itchyshin.github.io/pigauto>
```

The first 6 lines of `README.md` should now read:

```markdown
# pigauto: Phylogenetic Imputation via Graph AUTO-encoders <img src="man/figures/logo.png" align="right" height="139" alt="pigauto logo"/>

**Missing trait data should not stop a comparative analysis.**

> Live documentation: <https://itchyshin.github.io/pigauto>

pigauto fills gaps in species trait matrices by combining the phylogenetic
```

- [ ] **Step 4: Commit**

```bash
git add README.md DOCS.md
git commit -m "docs: retire DOCS.md, point README at live pkgdown site

DOCS.md (282 lines) was a parallel index of repo artefacts that
duplicated and drifted from README + pkgdown's auto-generated
navigation.  Single source of truth: README is canonical for repo
visitors; pkgdown navbar is canonical for live-site visitors.
"
```

---

## Task 10: NEWS.md entry

**Files:**
- Modify: `NEWS.md` (insert at the top, under the existing top-most version section header)

- [ ] **Step 1: Identify the current top-of-file version**

Run: `head -3 NEWS.md`
Expected: the first line is a `# pigauto X.Y.Z (dev)` header. Note `X.Y.Z`.

- [ ] **Step 2: Insert the new section**

Immediately after the top-most version's section header line, insert a blank line then:

```markdown
## Documentation: pkgdown reorganisation Phase 1 (2026-05-09)

User-facing docs reorganisation following the issues by @b1805 (#67, #68)
that surfaced confusion the README and `?impute` should have prevented.
PR #69 fixed the function-level docstring; this release reorganises the
pkgdown site itself:

- The Articles dropdown is slimmed from 39 entries to 8 first-class items
  (3 vignettes + Common Pitfalls + 4 walk-throughs).
- A new top-level Methodology dropdown owns the 23 benchmark HTMLs in
  three sub-sections: per-trait benches (8), cross-dataset benches (7),
  and Phase 8 simulations + head-to-heads (8).
- New `vignettes/common-pitfalls.Rmd`: five sections covering
  `prediction$imputed` semantics, K=3 ordinal majority-class collapse,
  the safety-floor closed-gate behaviour, phylogenetic signal diagnosis,
  and `clamp_outliers` tail safety.  Each follows a fixed
  Symptom / Why / Diagnose / Fix / See-also template.
- `DOCS.md` retired; README + live pkgdown site are the single source of
  truth for documentation navigation.

No code changes, no DESCRIPTION bump, no test changes.
```

- [ ] **Step 3: Commit**

```bash
git add NEWS.md
git commit -m "docs(news): pkgdown reorganisation Phase 1 entry"
```

---

## Task 11: Full pkgdown rebuild + R CMD check

**Files:** none (verification only).

- [ ] **Step 1: Rebuild the site from scratch**

Run: `Rscript -e 'pkgdown::clean_site(); pkgdown::build_site_local(preview = FALSE)'`
Expected: completes in ~60 s with no errors and no warnings about missing target HTMLs (the `articles/common-pitfalls.html` warning from Task 2 is gone now).

- [ ] **Step 2: Open the site index and sanity-check the navbar**

Run: `Rscript -e 'utils::browseURL(file.path("docs", "index.html"))'`
(Or open `docs/index.html` in your browser.)

Visual check:

- Navbar shows: `[Get started] [Articles] [Methodology] [Reference] [Changelog]`.
- Articles dropdown has 8 items: Get started, Common pitfalls / FAQ, Mixed-type traits, Propagating tree uncertainty, Architecture overview, Mixed-type workflow, Comparative study with covariates, Multi-observation per species.
- Methodology dropdown has three sub-sections totalling 23 working items + the Phase-2 placeholder is intentionally absent (it's a YAML comment).
- The home page no longer has the "Documentation index (DOCS.md)" link.
- The home page has the new "Live documentation:" line in the README rendering.

If anything is wrong, fix the YAML / vignette / README and re-run Step 1.

- [ ] **Step 3: Click each Methodology link**

Cycle through all 23 Methodology entries. Each should resolve to a non-empty HTML page (the existing benchmark output). If any link 404s, regenerate the corresponding bench HTML using the relevant `script/make_bench_*_html.R` driver and re-run Step 1.

- [ ] **Step 4: Click "Common pitfalls / FAQ"**

Verify all five sections render with code chunks visible (chunks are `eval = FALSE`, so the code is shown as text without an output block). Verify the "See also" links resolve: `?impute`, the issues #67 / #68, the memo links, and the Methodology bench link in Pitfall 4.

- [ ] **Step 5: Run R CMD check**

Run: `Rscript -e 'devtools::check(args = c("--no-manual"), error_on = "warning")'`
Expected: `Status: OK`. The new vignette is part of the vignette-build step. If any warning appears, fix and re-run.

- [ ] **Step 6: Final lint of NEWS.md formatting**

Run: `head -40 NEWS.md`
Expected: clean Markdown — no stray characters, the new section appears under the top-most `# pigauto X.Y.Z (dev)` header.

---

## Task 12: Open the PR

**Files:** none.

- [ ] **Step 1: Push the branch**

Run: `git push -u origin docs/pkgdown-reorg-phase1`

- [ ] **Step 2: Open the PR**

Run:

```bash
gh pr create --title "docs(pkgdown): reorganisation phase 1 (closes confusion behind #67, #68)" --body "$(cat <<'EOF'
## Summary

Phase 1 of the pkgdown reorganisation specced in
\`specs/2026-05-09-pkgdown-reorg-phase1-design.md\`.  Targets the
Bhavya-class new-user audience (priority A) by making the site easier
to navigate; methodology and developer audiences (priorities B and C)
are unchanged in this PR (Phase 2 will lift the math vignette).

## Changes

- \`_pkgdown.yml\`: Articles dropdown slimmed to 8 first-class items;
  new top-level Methodology dropdown owns 23 benchmark HTMLs in three
  sub-sections (per-trait, cross-dataset, Phase 8 sims); home.links
  DOCS.md entry removed.
- \`vignettes/common-pitfalls.Rmd\`: NEW.  Five sections following a
  fixed Symptom/Why/Diagnose/Fix/See-also template, covering the
  documented user confusions surfaced by issues #67 and #68 plus three
  more: closed-gate behaviour, phylo-signal diagnosis, clamp_outliers
  tail safety.
- \`DOCS.md\`: deleted (282 lines of parallel index that drifted from
  README and pkgdown).
- \`README.md\`: one new line near the top — \`> Live documentation:
  <https://itchyshin.github.io/pigauto>\` — making the live site the
  canonical doc destination for repo visitors.
- \`NEWS.md\`: new section under the top-most dev version describing
  the reorganisation.

No R code changes, no DESCRIPTION bump, no test changes.

## Test plan

- [x] \`pkgdown::clean_site() ; pkgdown::build_site_local()\` rebuilds clean
- [x] All 23 Methodology links resolve (manual click-through)
- [x] Common Pitfalls vignette knits in <= 30 s
- [x] \`R CMD check --no-manual\` Status: OK
- [x] Home page no longer links DOCS.md
- [x] README live-site line appears at the top of the rendered home page
- [ ] Reviewer eyeball: navbar shows [Get started] [Articles] [Methodology] [Reference] [Changelog]

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

Expected: returns the PR URL.  Leave the PR open for review; no auto-merge.

---

## Self-review (run by author before merging)

1. **Spec coverage check.**  Each item in `specs/2026-05-09-pkgdown-reorg-phase1-design.md`:

   - §5.1 navbar redesign → Task 2.
   - §5.2 Common Pitfalls vignette (5 sections, fixed template) → Tasks 3–8.
   - §5.3 single source of truth → Task 9.
   - §8 verification (8 steps) → Task 11 (Steps 1–6) + manual click-through.
   - §9 files-touched table — every file accounted for: `_pkgdown.yml` (Task 2), `vignettes/common-pitfalls.Rmd` (Tasks 3–8), `DOCS.md` (Task 9), `README.md` (Task 9), `NEWS.md` (Task 10).
   - §11 rollback — implicit; each commit is its own revertable unit (one commit per task).

2. **Placeholder scan.**  No `TBD` / `TODO` / "implement later" / "fill in details" in this plan.  The "(written in Task N)" lines inside the vignette skeleton are explicitly placeholders that subsequent tasks replace; they are not plan failures.

3. **Type / name consistency.**

   - Branch name: `docs/pkgdown-reorg-phase1` (Task 1, Task 12).
   - Vignette path: `vignettes/common-pitfalls.Rmd` and rendered `articles/common-pitfalls.html` (Task 2 navbar entry, Tasks 3–8 file edits, Task 11 verification).
   - Methodology dropdown component name: `methodology` (Task 2 YAML, Task 11 verification).
   - All 23 benchmark HTML basenames match between Task 1 verification loop and Task 2 YAML hrefs (verified character-by-character against the spec list).

No drift detected.
