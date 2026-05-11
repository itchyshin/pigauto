# pkgdown Page Rethink - 2026-05-11

Scope: light audit of old static pkgdown pages and methodology pages. This is
not a benchmark rerun, and it does not make new performance claims. Evidence is
limited to files read in this pass and the two agent-team read-only audits.

## Static Pages

| Page | Current value | Main risk found | Disposition | Estimate |
|---|---|---|---|---:|
| `pkgdown/assets/pigauto_intro.html` | Useful high-level architecture sketch. | Generated from `README.md` / `CLAUDE.md` rather than a source-backed vignette (`script/make_intro_html.R:1-5`), and broad claims such as 95% coverage / 10,000 species need regime labels before public use (`pkgdown/assets/pigauto_intro.html:84-86`). | Keep hidden and archived. Rewrite as a short source-backed architecture vignette only after checking README, current gate wording, and benchmark scope. | 0.5 day |
| `pkgdown/assets/pigauto_workflow_mixed.html` | Valuable downstream-inference narrative: pigauto + Rubin, pigauto + MCMCglmm, BACE. | Tutorial scope is five exercised trait types while discussing all eight (`pkgdown/assets/pigauto_workflow_mixed.html:81-90`), had SE-era draw wording (`pkgdown/assets/pigauto_workflow_mixed.html:285-304`), and contains illustrative BACE output from an earlier batch run (`pkgdown/assets/pigauto_workflow_mixed.html:640-646`). | Keep hidden and archived. Split into a source-backed downstream-inference vignette; keep BACE as a separately scoped path, not a pigauto benchmark claim. | 1-2 days |
| `pkgdown/assets/pigauto_walkthrough_covariates.html` | Useful Delhey/Gloger workflow concept and covariate/regression distinction. | Version badge is v0.6.0 (`pkgdown/assets/pigauto_walkthrough_covariates.html:28`), the table shows little/no covariate lift in this dataset (`pkgdown/assets/pigauto_walkthrough_covariates.html:102-129`), and the page had broad "never hurt" / "always" claims (`pkgdown/assets/pigauto_walkthrough_covariates.html:131-138`, `pkgdown/assets/pigauto_walkthrough_covariates.html:286-293`). | Keep hidden and archived. Rewrite after choosing a current validation table and stating when covariates did or did not help. | 1 day |
| `pkgdown/assets/pigauto_walkthrough_multi_obs.html` | Useful concise demo of `species_col`, `ctmax_sim`, and observation-level covariates. | Version badge is v0.6.0 (`pkgdown/assets/pigauto_walkthrough_multi_obs.html:49`), generator hardcodes a local checkout path (`script/make_walkthrough_multi_obs_html.R:16-21`), and its BACE comparison had stale trait-type wording (`pkgdown/assets/pigauto_walkthrough_multi_obs.html:159-175`). | Keep hidden and archived. Rewrite as a small source-backed vignette if multi-observation tutorial coverage is still wanted. | 0.5-1 day |

Immediate mitigation applied in this pass:

- `_pkgdown.yml` keeps these pages out of the Articles dropdown and now points
  maintainers to this audit before re-adding them.
- The four tracked static HTML assets now show a visible archived-page banner.
- The tracked generators now carry archived-generator comments.
- The most obvious overclaims in those tracked assets were softened without
  pretending the pages are fully current.

## Methodology Pages

| Group | Decision | Evidence / risk | Estimate |
|---|---|---|---:|
| Active `_pkgdown.yml` Methodology menu | Keep shape, update rendered pages. | Active links are present in local `docs/dev`, but the ignored `docs/` pages lag the newer `pkgdown/assets/dev` copies. Example: continuous docs are older than the tracked asset (`docs/dev/bench_continuous.html:38`, `pkgdown/assets/dev/bench_continuous.html:38`). | 1-2h rebuild/check |
| Per-trait pages: continuous, binary, ordinal, count, categorical, proportion, zi_count | Keep, but label regimes tightly. | Useful evidence, mostly narrow simulated regimes such as `n = 300`, MCAR, limited reps (`script/bench_binary.md:3`, `script/bench_categorical.md:3`, `script/bench_zi_count.md:3`). | 2-3h |
| `dev/bench_avonet_missingness.html` | Update before relying on it. | Page-level summary says pigauto "matches or beats" and "never degrades" in places where table rows are more mixed (`pkgdown/assets/dev/bench_avonet_missingness.html:41`, `pkgdown/assets/dev/bench_avonet_missingness.html:73`, `pkgdown/assets/dev/bench_avonet_missingness.html:87`). | 1-2h |
| `dev/bench_tree_uncertainty.html` | Keep, update wording. | Navbar says 10 posterior trees, while validation-suite copy says 50 trees; the script is 10 trees x 5 imputations = 50 datasets (`_pkgdown.yml:86`, `script/bench_tree_uncertainty.md:3`, `pkgdown/assets/validation_suite.html:74`). | 30-60m |
| Delhey / covariate simulation / multi-observation pages | Keep with caveats. | Delhey shows essentially no covariate lift (`script/bench_delhey.md:63`); positive covariate/multi-observation results are simulated DGP-specific (`script/bench_covariate_sim.md:3`, `script/bench_multi_obs.md:8`). | 1-2h |
| Validation suite | Update hard; do not publish as-is. | Public/local paths and versions are inconsistent, and the validation asset duplicates rows while linking absent or hidden pages (`docs/index.html:443`, `docs/validation_suite.html:30`, `pkgdown/assets/validation_suite.html:30`, `pkgdown/assets/validation_suite.html:75`). | 0.5 day |
| `bench_multi_proportion.html` | Update, then add or keep. | It represents one advertised trait type but is commented out of the navbar; current asset commentary needs tighter wording because shown rows can be pigauto = baseline (`_pkgdown.yml:99`, `script/bench_multi_proportion.md:19`, `pkgdown/assets/dev/bench_multi_proportion.html:237`). | 1-2h |
| `bench_scaling_v090.html` | Hide/update before linking. | Validation suite claims max N = 10,000 / 42 stages, while the visible v0.9.0 script table observed here stops at 5,000; `docs/dev/scaling.html` is a separate v0.3.1 post-mortem (`pkgdown/assets/validation_suite.html:75`, `script/bench_scaling_v090.md:20`, `docs/dev/scaling.html:34`). | 1-2h |
| Phase 8 pages: signal, correlation, evo-model, clade-MAR | Keep hidden or publish as simulated dev evidence. | They are absent from the active navbar and are narrow synthetic regimes (`_pkgdown.yml:103`, `script/bench_signal_sweep.md:3`, `script/bench_evo_model_sweep.md:3`). | 2-4h |
| BACE head-to-head pages | Archive/hide until BACE actually runs. | Multiple "vs BACE" pages are pigauto-only because BACE was skipped (`script/bench_bace_avonet_head_to_head.md:3`, `script/bench_pantheria_bace_head_to_head.md:3`, `pkgdown/assets/dev/bench_avonet9993_bace_n3000.html:1`). | 0.5 day without rerun |
| Legacy `docs/dev/benchmark_report.html`, `benchmark_final_report.html`, `test_report.html`, `bench_tabpfn.html` | Archive/remove from public surface. | NEWS says these old pages were removed, but ignored local `docs/` still contains them and old nav links remain in generated pages (`docs/news/index.html:1078`, `docs/news/index.html:1156`, `docs/reference/fit_baseline_tabpfn.html:34`). | 30-60m |

## Next Patch Queue

1. Fix remaining roxygen warnings so `devtools::document()` is quiet enough to
   expose real doc regressions.
2. Decide whether `goodagents.md` and `AGENTS.md` should be tracked; if tracked,
   add the right build-ignore policy.
3. Do the deeper methodology pass: one table per active benchmark page with DGP,
   n, missingness, reps, compared methods, metric, and what claim is allowed.
