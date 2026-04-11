# pigauto — documentation index

> **Live site**: <https://itchyshin.github.io/pigauto>
> **Source**: <https://github.com/itchyshin/pigauto>
> **One-line summary**: fill in missing species traits using a
> phylogenetic tree, cross-trait correlations, and optional
> environmental covariates. Handles continuous measurements, counts,
> binary, ordinal, and categorical variables in a single call.

This file is the single entry point for everything that ships with the
pigauto repository: tutorials, function reference, design notes,
changelog, and benchmark reports. Start here when you are trying to
find something and are not sure where to look.

The **live site** above is built automatically from this repo on every
push to `main` (GitHub Actions → `pkgdown`). It is the most convenient
way to browse tutorials, function reference, and the changelog as
clickable web pages. Everything linked from the live site also lives
in this repo, and the tables below give both the live URL (web-viewable)
and the in-repo source file (which you use to rebuild it) for each
artefact.

> **Note on in-repo HTML files.** The rendered HTML tutorials live at
> `inst/doc/*.html` and `docs/*.html`, but both directories are
> `.gitignore`d — they are *build outputs*, not committed source.
> Click the live-site links for the rendered page; click the "Source"
> column to see the script that generates it. To rebuild all HTMLs
> locally, see the "How this site is built" appendix at the bottom.

---

## 1. Start here

| What you want | Where to go |
|---|---|
| A one-page pitch and installation | [`README.md`](README.md) |
| First imputation, end-to-end | Live: [articles/getting-started.html](https://itchyshin.github.io/pigauto/articles/getting-started.html) · source: [`vignettes/getting-started.Rmd`](vignettes/getting-started.Rmd) |
| Why does pigauto work this way? | Live: [pigauto_intro.html](https://itchyshin.github.io/pigauto/pigauto_intro.html) · source: [`script/make_intro_html.R`](script/make_intro_html.R) |
| The full PCM workflow, mixed-type traits | Live: [pigauto_workflow_mixed.html](https://itchyshin.github.io/pigauto/pigauto_workflow_mixed.html) · source: [`script/make_workflow_mixed_html.R`](script/make_workflow_mixed_html.R) |
| How well does it actually work? | Section [3. Benchmarks](#3-benchmarks) below, or live: [scaling to 10,000 tips](https://itchyshin.github.io/pigauto/dev/scaling.html) / [AVONET missingness sweep](https://itchyshin.github.io/pigauto/dev/bench_avonet_missingness.html) |
| Package internals and gotchas | [`CLAUDE.md`](CLAUDE.md) |

---

## 2. Tutorials

The tutorials come in two shapes. **Knitted vignettes** (`*.Rmd` in
`vignettes/`) are built by `devtools::build_vignettes()` and become
`inst/doc/*.html`. **Static HTML walk-throughs** (`script/make_*_html.R`)
are hand-written R scripts that write pre-formatted HTML directly into
`inst/doc/` — they are faster to iterate on than Rmd because they do
not need to execute the code, but reader code snippets need to be
copy-pasted into an R session.

| Tutorial | Live link | Source (committed) | Rendered output (gitignored) | Covers |
|---|---|---|---|---|
| Getting started | [articles/getting-started.html](https://itchyshin.github.io/pigauto/articles/getting-started.html) | [`vignettes/getting-started.Rmd`](vignettes/getting-started.Rmd) | `inst/doc/getting-started.html` | First imputation on the bundled AVONET 300 dataset; how to read `pigauto_result` |
| Mixed-type traits | [articles/mixed-types.html](https://itchyshin.github.io/pigauto/articles/mixed-types.html) | [`vignettes/mixed-types.Rmd`](vignettes/mixed-types.Rmd) | `inst/doc/mixed-types.html` | Continuous + count + ordinal + binary + categorical in one model; round-trip encoding |
| Intro (architecture overview) | [pigauto_intro.html](https://itchyshin.github.io/pigauto/pigauto_intro.html) | [`script/make_intro_html.R`](script/make_intro_html.R) | `inst/doc/pigauto_intro.html` | The gated baseline + GNN ensemble, gate calibration, conformal intervals, trait-type handling |
| Mixed-type PCM workflow | [pigauto_workflow_mixed.html](https://itchyshin.github.io/pigauto/pigauto_workflow_mixed.html) | [`script/make_workflow_mixed_html.R`](script/make_workflow_mixed_html.R) | `inst/doc/pigauto_workflow_mixed.html` | Full multiple-imputation pipeline on mixed trait types, plus all three analysis paths (A/B/C) |
| Comparative study with covariates | [pigauto_walkthrough_covariates.html](https://itchyshin.github.io/pigauto/pigauto_walkthrough_covariates.html) | [`script/make_walkthrough_covariates_html.R`](script/make_walkthrough_covariates_html.R) | `inst/doc/pigauto_walkthrough_covariates.html` | Gloger's rule on 5,809 species: impute plumage lightness using phylogeny + environment, then phylogenetic regression with Rubin's rules |
| Multi-observation per species | [pigauto_walkthrough_multi_obs.html](https://itchyshin.github.io/pigauto/pigauto_walkthrough_multi_obs.html) | [`script/make_walkthrough_multi_obs_html.R`](script/make_walkthrough_multi_obs_html.R) | `inst/doc/pigauto_walkthrough_multi_obs.html` | CTmax-like multi-obs imputation with observation-level covariates; covariate-conditional predictions within species |

### The three analysis paths

The mixed-type workflow tutorial is organised around three ways of
turning M imputed datasets into pooled inference. This framing is
explicit in the tutorial and informal in the rest of the docs:

- **Path A — pigauto + glmmTMB + Rubin's rules** (frequentist). Fast
  iteration, regression-table output, works with any `coef()` +
  `vcov()` fit. Use `multi_impute()` → `with_imputations()` →
  `pool_mi()`.
- **Path B — pigauto + MCMCglmm + posterior concatenation** (Bayesian,
  pigauto as imputer). Use `multi_impute()` → `with_imputations()`
  with an `MCMCglmm()` call inside, then concatenate the `$Sol`
  chains. `pool_mi()` explicitly rejects `MCMCglmm` fits because
  Rubin's rules do not apply to posterior samples.
- **Path C — BACE integrated** (Bayesian, chained-equation MCMC).
  Use `BACE::bace()` + `BACE::pool_posteriors()` end-to-end. BACE is
  a separate in-tree package (`BACE/`); `pigauto::fit_baseline_bace()`
  wraps it as an alternative baseline but the integrated inference
  workflow lives in the BACE package itself.

---

## 3. Benchmarks

Two rendered HTML reports let you read pigauto's reported performance
without re-running anything locally. Both are physically hosted under
`pkgdown/assets/dev/`; the "Benchmarks" submenu on the live site's
Articles dropdown links straight to them.

| Report | Live link | Source / inputs | What it shows |
|---|---|---|---|
| **Scaling to 10,000 tips** | [dev/scaling.html](https://itchyshin.github.io/pigauto/dev/scaling.html) | [`script/make_scaling_html.R`](script/make_scaling_html.R) + `script/bench_scaling_v030.rds` + `script/bench_scaling_v031.rds` | Side-by-side v0.3.0 vs v0.3.1 wall-clock and peak memory at n ∈ {300, 1k, 2k, 3k, 5k, 7.5k, 10k}, plus an end-to-end validation run on the full AVONET3 + BirdTree dataset (9,993 species, 7 mixed-type traits). |
| **AVONET missingness sweep** | [dev/bench_avonet_missingness.html](https://itchyshin.github.io/pigauto/dev/bench_avonet_missingness.html) | [`script/make_avonet_missingness_html.R`](script/make_avonet_missingness_html.R) + `script/bench_avonet_missingness.rds` | Three methods (mean / mode, BM baseline, pigauto) × three MCAR missingness levels (20% / 50% / 80%) on the full 9,993-species bundled `avonet_full` dataset. Per-trait RMSE, Pearson r, accuracy. |
| **Continuous traits** | [dev/bench_continuous.html](https://itchyshin.github.io/pigauto/dev/bench_continuous.html) | `script/bench_continuous.R` + `script/make_bench_continuous_html.R` | 4 evolutionary models (BM, OU, regime shift, nonlinear) × 3 methods × 5 reps on 300-species simulated trees. Primary + secondary (missingness) sweep. |
| **Binary traits** | [dev/bench_binary.html](https://itchyshin.github.io/pigauto/dev/bench_binary.html) | `script/bench_binary.R` + `script/make_bench_binary_html.R` | Phylogenetic signal sweep (0.2–1.0) × class imbalance sweep. Mode vs phylo LP vs pigauto. |
| **Ordinal traits** | [dev/bench_ordinal.html](https://itchyshin.github.io/pigauto/dev/bench_ordinal.html) | `script/bench_ordinal.R` + `script/make_bench_ordinal_html.R` | Level count sweep (3–10) × signal sweep. Median vs BM vs pigauto. |
| **Count traits** | [dev/bench_count.html](https://itchyshin.github.io/pigauto/dev/bench_count.html) | `script/bench_count.R` + `script/make_bench_count_html.R` | Mean count sweep (5–500) × Poisson vs NegBin. Mean vs BM vs pigauto. |
| **Categorical traits** | [dev/bench_categorical.html](https://itchyshin.github.io/pigauto/dev/bench_categorical.html) | `script/bench_categorical.R` + `script/make_bench_categorical_html.R` | Category count sweep (K=3–12) × signal sweep. Mode vs phylo LP vs pigauto. |
| **Proportion traits** | [dev/bench_proportion.html](https://itchyshin.github.io/pigauto/dev/bench_proportion.html) | `script/bench_proportion.R` + `script/make_bench_proportion_html.R` | Signal sweep (0.2–1.0) × boundary density sweep. Mean vs BM vs pigauto. |
| **Zero-inflated counts** | [dev/bench_zi_count.html](https://itchyshin.github.io/pigauto/dev/bench_zi_count.html) | `script/bench_zi_count.R` + `script/make_bench_zi_count_html.R` | Zero fraction sweep (0.2–0.8) × mean non-zero count sweep. Mean vs LP+BM vs pigauto. |
| **Missingness mechanisms** | [dev/bench_missingness_mechanism.html](https://itchyshin.github.io/pigauto/dev/bench_missingness_mechanism.html) | `script/bench_missingness_mechanism.R` + `script/make_bench_missingness_mechanism_html.R` | MCAR vs MAR_trait vs MAR_phylo vs MNAR on mixed-type data. Cross-cutting validation of all types under realistic missingness. |
| **Tree uncertainty** | [dev/bench_tree_uncertainty.html](https://itchyshin.github.io/pigauto/dev/bench_tree_uncertainty.html) | `script/bench_tree_uncertainty.R` + `script/make_bench_tree_uncertainty_html.R` | `multi_impute_trees()` with 10 posterior trees vs single-tree MI. SE inflation (1.1–2.1×) and FMI rising with missingness. |
| **Environmental covariates (Delhey)** | [dev/bench_delhey.html](https://itchyshin.github.io/pigauto/dev/bench_delhey.html) | `script/bench_delhey.R` + `script/make_bench_delhey_html.R` | Real-data test: plumage lightness in 5,809 bird species with 4 environmental covariates (Delhey 2019). Covariates give ~zero lift when phylogenetic signal dominates. |
| **Covariate effectiveness** | [dev/bench_covariate_sim.html](https://itchyshin.github.io/pigauto/dev/bench_covariate_sim.html) | `script/bench_covariate_sim.R` + `script/make_bench_covariate_sim_html.R` | Simulated traits with varying phylogenetic signal (λ) and environmental effect (β). Covariates reduce RMSE by 8–15% when env effects are strong; ~0% when phylo signal dominates (gated safety). |
| **Multi-observation imputation** | [dev/bench_multi_obs.html](https://itchyshin.github.io/pigauto/dev/bench_multi_obs.html) | `script/bench_multi_obs.R` + `script/make_bench_multi_obs_html.R` | CTmax-like simulation: 200 species × variable obs/species × phylo signal × ARR sweep. species_mean vs pigauto (no cov) vs pigauto (with cov). Observation-level RMSE advantage with covariates when within-species effects are strong. |

For the authoritative numerical summary with discussion, see
[`README.md` → Benchmark results](README.md#benchmark-results). The
interactive `simulate_benchmark()` function in
[`R/simulate_benchmark.R`](R/simulate_benchmark.R) lets you regenerate
the synthetic simulation scenarios end-to-end.

---

## 4. Function reference

The full auto-generated reference lives at the pkgdown site
<https://itchyshin.github.io/pigauto/reference/>. A copy of every
exported function's `.Rd` page is also in `man/`.

| Group | Key functions | Purpose |
|---|---|---|
| One-call entry point | [`impute()`](man/impute.Rd) | Runs the full pipeline and returns a completed data frame |
| Multiple imputation | [`multi_impute()`](man/multi_impute.Rd), [`multi_impute_trees()`](man/multi_impute_trees.Rd), [`with_imputations()`](man/with_imputations.Rd), [`pool_mi()`](man/pool_mi.Rd) | Generate M datasets → fit M models → pool with Rubin's rules (including tree-uncertainty propagation) |
| Fine-grained pipeline | [`preprocess_traits()`](man/preprocess_traits.Rd), [`build_phylo_graph()`](man/build_phylo_graph.Rd), [`fit_baseline()`](man/fit_baseline.Rd), [`fit_pigauto()`](man/fit_pigauto.Rd), [`predict.pigauto_fit()`](man/predict.pigauto_fit.Rd), [`evaluate()`](man/evaluate.Rd) | The six stages `impute()` runs internally |
| Alternative baseline | [`fit_baseline_bace()`](man/fit_baseline_bace.Rd) | BACE (Bayesian chained-equation multiple imputation) as an alternative to the default phylogenetic BM baseline |
| Cross-validation and benchmarks | [`cross_validate()`](man/cross_validate.Rd), [`compare_methods()`](man/compare_methods.Rd), [`simulate_benchmark()`](man/simulate_benchmark.Rd), [`simulate_non_bm()`](man/simulate_non_bm.Rd) | Evaluation infrastructure |
| Reporting and plotting | [`pigauto_report()`](man/pigauto_report.Rd), [`plot.pigauto_fit()`](man/plot.pigauto_fit.Rd), [`plot.pigauto_pred()`](man/plot.pigauto_pred.Rd), [`plot.pigauto_benchmark()`](man/plot.pigauto_benchmark.Rd), [`summary.pigauto_fit()`](man/summary.pigauto_fit.Rd), [`calibration_df()`](man/calibration_df.Rd), [`confusion_matrix()`](man/confusion_matrix.Rd) | Interactive HTML reports + diagnostic plots |
| I/O | [`read_traits()`](man/read_traits.Rd), [`read_tree()`](man/read_tree.Rd), [`save_pigauto()`](man/save_pigauto.Rd), [`load_pigauto()`](man/load_pigauto.Rd) | Data loading and fit persistence |
| Bundled data | [`avonet300`](man/avonet300.Rd), [`tree300`](man/tree300.Rd), [`trees300`](man/trees300.Rd), [`avonet_full`](man/avonet_full.Rd), [`tree_full`](man/tree_full.Rd), [`delhey5809`](man/delhey5809.Rd), [`tree_delhey`](man/tree_delhey.Rd), [`ctmax_sim`](man/ctmax_sim.Rd) | 300-species and 9,993-species AVONET + 5,809-species Delhey plumage + simulated multi-obs CTmax + matching BirdTree phylogenies |

---

## 5. Design notes and architecture

For contributors, reviewers, and anyone trying to understand what
pigauto does under the hood:

- [`CLAUDE.md`](CLAUDE.md) — the canonical architecture reference.
  S3 class table, data flow, trait-type handling, gate initialisation,
  post-training calibration and conformal, multi-obs code path, and
  the gotchas that are easy to miss on a cold read of the source. Also
  documents the `BACE/` in-tree package and how it relates to pigauto.
- [`README.md`](README.md) — user-facing tour plus benchmark summary
  tables.
- [`NEWS.md`](NEWS.md) — per-version release notes.

Source-code map:

```
R/
├── impute.R                one-call entry point
├── multi_impute.R          M-dataset generator + pigauto_mi class
├── with_imputations.R      map a user fit function over M datasets
├── pool_mi.R               Rubin's rules + Barnard-Rubin df
├── preprocess_traits.R     trait-type detection + latent matrix
├── build_phylo_graph.R     adjacency + spectral coords
├── bm_internal.R           internal phylogenetic BM baseline
├── fit_baseline.R          BM + phylo label propagation
├── fit_baseline_bace.R     BACE wrapper (Suggests:)
├── fit_pigauto.R           training loop + post-training calibration
├── model_residual_dae.R    ResidualPhyloDAE torch::nn_module
├── predict_pigauto.R       MC-dropout sampling + conformal intervals
├── evaluate.R              per-trait metrics
├── cross_validate.R        k-fold CV driver
├── simulate_benchmark.R    synthetic benchmark scenarios
├── report.R                interactive HTML report builder
└── plot.R                  diagnostic plots
```

---

## 6. Developer reports

These HTML reports are dev artefacts living under `script/`. They do
not ship with the installed package (`^script$` is in `.Rbuildignore`)
but are committed to the repo and mirrored into `pkgdown/assets/dev/`
so the live site can serve them as rendered web pages. The two
benchmark reports (scaling, AVONET missingness) are also listed here
because they are technically developer artefacts, but the primary
place to find them is section 3 above.

| Report | Live link | Source |
|---|---|---|
| Test catalogue | [dev/tests_overview.html](https://itchyshin.github.io/pigauto/dev/tests_overview.html) | [`script/make_tests_html.R`](script/make_tests_html.R) |

---

## 7. Changelog

- [`NEWS.md`](NEWS.md) — release history, user-facing changes per version.
- Web version on the live site: <https://itchyshin.github.io/pigauto/news/index.html>.

---

## 8. Citation and license

- [`DESCRIPTION`](DESCRIPTION) — package metadata, authors, version.
- [`LICENSE`](LICENSE) — MIT.
- Canonical citation (BibTeX friendly):

  > Nakagawa S (2026). *pigauto: Phylogenetic Imputation via Graph
  > Autoencoder*. R package version 0.5.0.
  > <https://github.com/itchyshin/pigauto>.

---

## 9. How this site is built (appendix)

The live pkgdown site is produced by two things:

1. **`_pkgdown.yml`** at the repo root configures the site structure
   (navbar, reference groups, articles menu, template).
2. **`.github/workflows/pkgdown.yaml`** runs
   `pkgdown::build_site_github_pages()` on every push to `main` and
   deploys the `docs/` output via GitHub Pages' workflow-based
   deployment (`actions/upload-pages-artifact` +
   `actions/deploy-pages`). The site is served at
   <https://itchyshin.github.io/pigauto>.

### Static HTML tutorials: the dual-write pattern

The two `script/make_*_html.R` tutorial scripts (`make_intro_html.R`,
`make_workflow_mixed_html.R`) and the dev report builders
(`make_tests_html.R`, `make_scaling_html.R`,
`make_avonet_missingness_html.R`) each write their output to **two
places**:

- `inst/doc/xxx.html` (or `script/xxx.html` for dev reports) —
  the canonical location. For tutorials this is what ships with the
  installed package; readers open it via
  `browseURL(system.file("doc", "pigauto_workflow_mixed.html", package = "pigauto"))`.
- `pkgdown/assets/xxx.html` (or `pkgdown/assets/dev/xxx.html` for dev
  reports) — `pkgdown` copies `pkgdown/assets/*` verbatim into the
  `docs/` build output, so the same HTML becomes accessible on the
  live site at
  `https://itchyshin.github.io/pigauto/xxx.html`.

This keeps the installed-package copy and the web copy in lockstep:
one `Rscript script/make_workflow_mixed_html.R` rebuild touches both.

### Rebuilding the whole site locally

```r
# regenerate static HTML tutorials + dev reports
Rscript script/make_intro_html.R
Rscript script/make_workflow_mixed_html.R
Rscript script/make_walkthrough_covariates_html.R  # requires walkthrough_covariates.rds
Rscript script/make_tests_html.R
Rscript script/make_scaling_html.R                 # requires bench_scaling_v03*.rds
Rscript script/make_avonet_missingness_html.R      # requires bench_avonet_missingness.rds

# build vignettes (knitted from vignettes/*.Rmd)
devtools::build_vignettes()

# build the pkgdown site into docs/
pkgdown::build_site()
```

`pkgdown::build_site()` will pick up the vignettes from `vignettes/`,
the static HTMLs from `pkgdown/assets/`, the function reference from
`man/`, and `NEWS.md` for the changelog, and assemble them into
`docs/`. Open `docs/index.html` in a browser to preview before pushing.

### Artefact → source map

| Artefact | Source | Rebuild command |
|---|---|---|
| `inst/doc/getting-started.html` | `vignettes/getting-started.Rmd` | `devtools::build_vignettes()` |
| `inst/doc/mixed-types.html` | `vignettes/mixed-types.Rmd` | `devtools::build_vignettes()` |
| `inst/doc/pigauto_intro.html` | `script/make_intro_html.R` | `Rscript script/make_intro_html.R` |
| `inst/doc/pigauto_workflow_mixed.html` | `script/make_workflow_mixed_html.R` | `Rscript script/make_workflow_mixed_html.R` |
| `inst/doc/pigauto_walkthrough_covariates.html` | `script/make_walkthrough_covariates_html.R` | `Rscript script/make_walkthrough_covariates_html.R` |
| `script/tests_overview.html` | `script/make_tests_html.R` | `Rscript script/make_tests_html.R` |
| `script/scaling.html` | `script/make_scaling_html.R` | `Rscript script/make_scaling_html.R` |
| `script/bench_avonet_missingness.html` | `script/make_avonet_missingness_html.R` | `Rscript script/make_avonet_missingness_html.R` |
| `docs/` (pkgdown site) | `_pkgdown.yml` + everything above | `pkgdown::build_site()` locally, or push to `main` for CI build |
| `man/*.Rd` | roxygen blocks in `R/*.R` | `devtools::document()` |
