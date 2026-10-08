# Reader and public-site surface inventory

Scope: current checkout `release/cran-0.11-audit` compared with the public pkgdown deployment observed 2026-10-06. The tables below record the discrepancies found before this slice was edited; line references identify that initial source. Current implementation references are principally `R/multi_impute.R`, `R/mi_posterior.R`, and the defaults in `R/impute.R` / `R/fit_pigauto.R`. The 0.11.0 defaults are `gnn = FALSE`, eligible automatic posterior draws, separation residual prior, estimated lambda, and gate/floor off. Existing discrete-lambda and `sep_validation` results are intentionally not restated here.

## Disposition in the audit worktree

The source help, README, active vignettes, project guidance, and pkgdown navigation have been corrected in this worktree. In particular, the help source now identifies `"auto"` as the default, describes the active separation prior, limits the optional safety-floor blend to prediction diagnostics, and clarifies estimated discrete lambda. `devtools::document()` regenerated the affected Rd files. A local pkgdown site was built from the isolated pigauto 0.11.0 installation; the deployed public site still needs a post-merge check.

All 34 dated HTML benchmark pages and five supporting PNGs have been moved byte-for-byte from `pkgdown/assets/dev/` to the build-excluded `dev/archive/cran-011-public-pages/`. Their prior URLs are absent from the new site source and navigation. The exact files and SHA-256 values are in `retirement-manifest.md`; benchmark drivers and results remain under `script/`. Public URL retirement takes effect only after the new site is deployed. The excluded `vignettes/articles/simulation-study.Rmd` still records a historical GNN-on design and is not an active site article.

## Release-critical discrepancies

| Surface | Disposition | Evidence and action |
|---|---|---|
| `R/multi_impute.R` roxygen for `draws_method` details (around lines 265–285) | **Update** | The argument description correctly says `"auto"` is default (lines 38–44), and the function formal uses `c("auto", ...)` (line 415), but details still call `"conformal"` the default (line 272). Remove that contradiction and describe automatic eligibility/fallback consistently. `multi_impute()` resolves `"auto"` at lines 437–445. |
| `R/multi_impute.R` roxygen posterior model/prior equation (around lines 337–349) | **Update** | The prose gives `Sigma_E ~ IW(...)` as the active residual prior. Current `residual_prior` defaults to `"sep"` (`R/mi_posterior.R` around lines 943–961); its separated scale/correlation target is implemented around lines 406–418 and sampled through Metropolis around lines 717–720. Retain the inverse-Wishart formula only for the explicit legacy `"iw"` option. The phylogenetic `Sigma_W` inverse-Wishart description remains applicable. |
| `R/multi_impute.R` roxygen “Safety floor” section (near the end of the file) | **Update** | This section says `safety_floor = TRUE` is the default “since v0.9.1.9002”; current intended/release default is off. Also reconcile the claim that the three-way mean blend automatically propagates through draws with the current GNN-off default and call path. |
| `man/multi_impute.Rd` | **Replace by regeneration after roxygen fixes** | Generated help has the same stale default label (`draws_method = "conformal" (default)`, around line 298), stale active `Sigma_E` inverse-Wishart equation (around lines 361–362), and a safety-floor section claiming the default is TRUE (around lines 407–409). Live deployment reference page reproduces these errors verbatim. Do not hand-edit Rd independently. |
| `R/multi_impute.R` lambda prose (around lines 125–136) and generated help | **Update** | It says discrete traits stay at lambda 1 except an exact-route exception. That conflicts with the intended `discrete_lambda = "estimate"` default on the public default route. Change the wording to state the continuous-family `lambda_mode` scope and point to the `discrete_lambda` control without repeating its benchmark evidence. The live help page shows the stale lambda-1 statement. |
| `AGENTS.md` uncertainty section (around lines 429–432) | **Update** | Calls conformal draws the default and preferred method, omits posterior and auto resolution, and treats `multi_impute()` diagnostic paths as the current main route. Update this internal printed summary to match auto/posterior eligibility and the current analysis-aware route. The nearby claim of an unconditional “exactly ≥95%” conformal guarantee (line 426) should be qualified by exchangeability / calibration-match assumptions; README already warns about clade-shifted calibration. |
| `vignettes/multiple-imputation.Rmd` | **Update, then rebuild article** | The main routing table and auto-default text are current. Its short residual-prior note (lines 71–74) is historical and broadly correct; preserve that historical context. The rendered posterior-method account should also explain the current separation prior and not imply the old IW covariance formula is the active model. The site article URL could not be fetched through the web text cache, so live article contents remain unverified. |
| `vignettes/getting-started.Rmd` | **Keep with focused check** | It states the auto default and posterior eligibility at about line 413. Align the neighboring “supported pooling” wording with the narrow analysis-aware backend; rebuild and inspect after the roxygen/article edits. |
| `README.md` | **Update** | Its main description of auto routing is current (lines 63–68; defaults table around 158–169), but the quick-start does not state the new GNN-off default, estimated lambda, or gate/floor setting. Add a compact default statement / migration cue so a reader knows the default `impute()` fit is baseline-only and how to opt into GNN. Keep the existing caveat that default posterior pooling is restricted to all-continuous, single-row, no-covariate inputs. |
| `NEWS.md` | **Keep; targeted consistency edit only if needed** | Top entry records the estimated discrete lambda and gate/floor off defaults; subsequent entries document separation prior, auto posterior eligibility, and GNN off. These are historical release notes and their measured regimes should remain attached to the claims. Do not repeat or expand the existing discrete-lambda and `sep_validation` evidence in this audit. Older “conformal default” entries (e.g. around line 3018) are historical and should remain clearly under their old release heading. |

## Other reader surfaces

| Surface | Disposition | Evidence and action |
|---|---|---|
| `man/*.Rd` other than `multi_impute.Rd` | **Keep pending targeted scan after source changes** | Public help is source-linked in pkgdown. Run `devtools::document()` after the source edits and inspect changed Rd for any changed formal defaults; prioritize `impute.Rd`, `fit_pigauto.Rd`, `fit_baseline.Rd`, and `multi_impute_trees.Rd` for the five intended default changes. Do not infer the whole manual is aligned from `multi_impute` alone. |
| `vignettes/mixed-types.Rmd`, `common-pitfalls.Rmd`, `tree-uncertainty.Rmd`, `gnn-architecture.Rmd` | **Keep; targeted search/update if conflicting defaults appear** | These are active knitted articles listed in `_pkgdown.yml`. Search them for default GNN/gate/lambda claims during the final consistency pass; do not rewrite the broader methods explanation. |
| README and vignette URLs / examples | **Keep** | README links to the active getting-started, mixed-types, multiple-imputation, pitfalls, and tree-sensitivity articles. Current local source examples are coherent on the automatic posterior route; re-check their destinations after rebuild. |

## pkgdown navigation, static pages, and direct URLs

`_pkgdown.yml` names `https://itchyshin.github.io/pigauto` and exposes search. Before this edit, its navbar had six active knitted articles and twelve static benchmark links, and the source tree had 34 tracked `pkgdown/assets/dev/*.html` files. The 12 formerly navbar-linked pages are:

| Pages | Disposition | Evidence and action |
|---|---|---|
| `bench_continuous.html`, `bench_binary.html`, `bench_ordinal.html`, `bench_count.html`, `bench_categorical.html`, `bench_proportion.html`, `bench_zi_count.html`, `bench_avonet_missingness.html`, `bench_missingness_mechanism.html`, `bench_delhey.html`, `bench_multi_obs.html`, `bench_covariate_sim.html` | **Retired to archive** | The former navbar exposed these at direct `/dev/<name>.html` paths. Public `bench_continuous.html` explicitly reports five replicates, run date 2026-05-30, and commit `794537121b`; it labels a Rphylopars BM baseline and a “pigauto (BM + GNN)” arm. The archived bytes remain historical evidence, not results for 0.11 defaults. |
| `bench_amphibio.html`, `bench_fishbase.html`, `bench_avonet9993_bace.html`, `bench_avonet9993_bace_index.html`, `bench_avonet9993_bace_n3000.html`, `bench_avonet_full_local.html`, `bench_bace_avonet_head_to_head.html`, `bench_bien.html`, `bench_clade_missingness.html`, `bench_multi_proportion.html`, `bench_pantheria_bace_head_to_head.html`, `bench_pantheria_full.html`, `bench_scaling_v090.html`, `bench_signal_sweep.html`, `bench_tree_uncertainty.html`, `bench_correlation_sweep.html`, `bench_evo_model_sweep.html`, `calibration_grid.html`, `pantheria_summary.html`, `phase8_summary.html`, `scaling.html`, `tests_overview.html` | **Retired to archive** | These 22 static pages were not in the navbar but could remain reachable by direct URL when copied to the site. All were moved unchanged to the build-excluded archive; the manifest records each byte hash. |
| `bench_clade_missingness.png`, `bench_correlation_sweep.png`, `bench_evo_model_sweep.png`, `bench_signal_sweep.png`, `calibration_grid_plot.png` | **Retired with parent pages** | These five PNGs were moved unchanged to the same archive. |
| Four retired walkthroughs named in `_pkgdown.yml` comments (`pigauto_intro`, `pigauto_workflow_mixed`, `pigauto_walkthrough_covariates`, `pigauto_walkthrough_multi_obs`) | **Keep retired; no redirect currently evidenced** | The config says these HTML pages were retired to `dev/archive/` and should not return to the navbar. This checkout has no corresponding `dev/archive/` content, and a direct retired URL was not verified on the live host. If those old URLs matter, add explicit redirects to the closest active vignette during site build; keep any archived source/evidence in Git. Do not reanimate old HTML as current instructions. |
| Sitemap / search index / direct static-page inventory | **Rebuild and crawl** | No generated site tree or sitemap is present in this checkout (`pkgdown/` holds source assets only). The deployed home page exposes search, current article navigation, and the dated bench links. A full sitemap and all direct URLs were not enumerated here. After changes, inspect generated sitemap/search index, crawl each active nav target and retained direct page, and check retired URLs/redirects. |
| `pkgdown/clean-internal-pages.R`, `pkgdown/extra.css`, favicon and site assets | **Keep** | Build/configuration assets; no contradiction found in this bounded content scan. Verify the clean-page script does not silently remove newly retired pages or leave broken links during the full site build. |

The archived source files and corresponding `script/` drivers/results remain evidence-bearing snapshots. Their historical numbers were not edited.

## Pre-merge live deployment check (2026-10-06; superseded by post-merge proof)

Direct browser inspection reached the live home page and `reference/multi_impute.html`. Both display version 0.11.0.9002. The home page matches the current README's auto-default description. The live help page is generated from the stale source material above and visibly says conformal is the default, discrete traits stay at lambda 1, residual covariance has an IW prior, and the safety floor defaults TRUE. This is direct deployment evidence, not a cache inference. The live continuous benchmark is directly reachable and confirms the May 30 / five-replicate / commit `794537121b` metadata. The web text fetch could not access the live multiple-imputation article, so its deployed state is not claimed. Direct archived walkthrough URL lookup was also inconclusive; treat redirect/404 behavior as unverified.

## Local site validation and pre-merge deployment gate (historical)

`man/multi_impute.Rd` and `man/pool_mi.Rd` were regenerated from roxygen. The local help page visibly shows `draws_method = "auto"`, `gnn = FALSE`, the separation residual prior, estimated discrete lambda, and `safety_floor = FALSE` as the default. The first site render used the machine's previously installed pigauto, whose `fit_pigauto()` still defaulted to GNN on; it was stopped and discarded. The completed clean build used a separately installed 0.11.0 source copy with its version and `gnn = FALSE` formal verified before rendering.

After `pkgdown/clean-internal-pages.R`, the final local crawler checked 65 rendered HTML pages and 3,621 local links/assets/anchors. It found no broken local references, no indexed or served internal coordination pages, and none of the 34 retired `/dev/` pages. A planted retired search path and a planted broken link each made the crawler fail as intended. Historical NEWS text still names two old files; the crawler checks indexed page paths rather than mistaking that historical prose for active URLs. Desktop home, `multi_impute()` help, the multiple-imputation article, and the corrected `multi_impute_trees()` help were visually inspected; the `multi_impute()` page was also checked at a 390-pixel mobile viewport. No clipping beyond a horizontally scrollable code block was found.

Follow-up after checking Symbolizer's fitted-model scope: the README and getting-started article now distinguish supported imputation analyses from the broader set of model classes whose fixed effects can be extracted. A fresh full pkgdown rebuild, internal-page cleaning step, and crawl passed 65 rendered pages and 3,622 local references, with 34 retired direct pages absent. The one extra local reference is the getting-started article's link to the multiple-imputation article. The corrected text was checked in both rendered pages; the earlier desktop/mobile visual inspection was not repeated for this prose-only change.

The multiple-imputation article received the same scope correction: its posterior row now names linear analyses on the imputation scale, and its mixed-model examples explicitly separate fixed-effect extraction from inferential validation. A subsequent full rebuild, internal-page cleaning step, and crawl passed 65 rendered pages and 3,621 local references, with 34 retired direct pages absent. The new wording was checked in the rendered article. This prose-only follow-up did not repeat visual viewport inspection.

A final surface read caught three more stale statements before this rebuild. The generated tree help now makes GNN training conditional on `gnn = TRUE`, describes the GNN-off default draw path, and states that `multi_impute()` defaults to `"auto"`. The getting-started install block now makes the torch runtime conditional on opting into the GNN. The current NEWS entry states that drmTMB and gllvmTMB are optional without Suggests entries, places their real checks outside the built package, and distinguishes the earlier conformal default from the final 0.11.0 automatic default. The rebuilt tree help, getting-started article, and NEWS page contain those corrections.

### 2026-10-07 live-site correction

The earlier sentence that the current live help page is stale has been superseded by a fresh read of the homepage, `multi_impute()` reference, and multiple-imputation article: those pages report the corrected 0.11.0 behavior. Search for `bench_continuous` returned no suggestion, and direct requests to `/dev/bench_continuous.html` and `/articles/pigauto_workflow_mixed.html` returned pkgdown 404 pages. This establishes those sampled controls only. Sitemap access, every route in the 34-page retirement manifest, and a deployment commit/build identity remain unverified; G6 stays open.

The final reader check found a stale generated `vignettes/getting-started.R` companion and a NEWS sentence calling GNN-on the default. The companion now matches a fresh normalized `knitr::purl()` extraction from its R Markdown source, including the optional torch setup and explicit GNN opt-in examples. NEWS now calls GNN-on opt-in. A full local rebuild, internal-page cleaning step, and crawl passed 65 rendered pages and 3,622 local references, with 34 retired direct pages absent. The corrected phrases were checked in the rendered getting-started and NEWS pages.

At the time of the pre-merge local check, the public deployment had not yet been updated. The following post-merge verification supersedes that snapshot.

## Post-merge live deployment verification (2026-10-06)

Pages workflow run #570 succeeded for merge commit
`bb5835d1b214d783da7b0b414df99aa6ba926bc7`. The live homepage exposes Get started,
Articles, Reference, and Changelog, with no per-trait benchmark links. The current
getting-started default, separated-residual wording in the multiple-imputation article,
and `multi_impute()` reference are live. Direct checks returned 404 for all 34 retired
`/dev/` HTML pages and the four retired walkthrough articles. The sitemap has 64
locations and zero retired routes; the 574-path search index also has zero retired routes.
The exact URLs and statuses are retained in
`provenance/public-deployment-live-check-2026-10-06.tsv`.

## Historical live-site and candidate-source observations through 2026-10-07

The detailed 2026-10-07 observations from both sides of the merge are retained in the
branch history and summarized in `GATES.md`. They record successive deployed mismatches,
source fixes, and bounded checks. Those snapshots are not the current deployed-state
verdict: the 2026-10-08 cache-busted Chrome checks in `GATES.md` are the latest site
assessment and still leave the named tree/help corrections and full sitemap/retired-route
verification open. Candidate-source checks and the 2026-10-07 rebuild do not establish
deployed behavior or validate the final release artifact.
