# pkgdown reorganisation — Phase 1 design

> **Status**: Draft for review (2026-05-09).
> **Audience priority**: A > B > C
>   (A) New users like Bhavya — fastest path from landing to a correct
>   `impute()` call;
>   (B) Methodologically-skeptical readers — easy access to evidence;
>   (C) Active developers / contributors.
> **Phase 2** (math vignette conversion + extension) is a separate
>   spec, written after Phase 1 ships.

## 1. Motivation

The two open issues by @b1805 (#67 Mass "outliers", #68 Migration all
Resident) are the same root cause: the public docs make it easy for a
new user to call `impute()` on a fully-observed matrix, read
`result$prediction$imputed` as "the imputed values", and conclude
something is broken when nothing is. PR #69 fixed the function-level
docstring and the README Quick Start. The pkgdown site
(<https://itchyshin.github.io/pigauto>) still has three structural
problems that re-introduce the same confusion for any reader who
arrives via the site rather than the README:

1. The Articles dropdown in the navbar has 39 items (3 vignettes + 4
   walk-throughs + 14 type/missingness benches + 8 Phase-8 sweeps + 2
   developer reports), all under one nested menu separated by
   `text: "-------"` rows. New users see a wall.
2. There is no dedicated landing page for the most common pre-flight
   confusions (`prediction$imputed` vs `completed[imputed_mask]`,
   majority-class collapse on imbalanced K=3 ordinal traits, gate
   staying closed on low-signal traits, predictions exceeding observed
   range, signal-strength diagnosis).
3. Three documents describe pigauto in slightly different ways — the
   README, `DOCS.md` (282 lines), and the pkgdown home `description:`.
   They drift.

The methodology evidence (Phase 8 sims, `useful/MEMO_*.md`,
`useful/pigauto_math_description.html`) is currently invisible to a
site visitor. That is also a real problem but is addressed in Phase
2; Phase 1 does not introduce new methodology content.

## 2. Goals (in priority order)

1. **Cut the Articles dropdown to ≤ 8 first-class items.** Every
   non-tutorial bench HTML moves to a new top-level "Methodology"
   navbar dropdown.
2. **Ship a single Common Pitfalls / FAQ vignette** linked as the
   second item under Articles, so a confused user lands on it within
   one click.
3. **Make the README the single source of truth for the package
   index.** Delete `DOCS.md`. Remove the redundant pkgdown home link
   to it.
4. **No code changes**, no test changes, no DESCRIPTION bump. Phase 1
   is docs-only.

## 3. Non-goals

- Math vignette conversion + extension (Phase 2).
- Re-running any benchmarks. The 22 existing bench HTMLs are reused
  as-is; only their position in the navbar changes.
- Changing `home.description` in `_pkgdown.yml`.
- Reorganising the Reference section. The current reference grouping
  (one-call entry point, MI workflow, pipeline, active imputation,
  covariates, benchmarks, reporting/plotting, I/O, bundled data) is
  appropriate.
- Touching individual function docstrings beyond what PR #69 already
  shipped.
- Any changes to vignette content other than the new
  `vignettes/common-pitfalls.Rmd`.

## 4. Architecture

This is a configuration + content change. There is no new R code, no
new module, and no change to data flow inside the package. The two
units of change are independent:

- **Unit 1: navbar configuration.** A single file (`_pkgdown.yml`)
  rewrites its `navbar` and removes one entry from `home.links`. The
  unit is self-contained.
- **Unit 2: Common Pitfalls vignette.** A new file
  (`vignettes/common-pitfalls.Rmd`) plus its referenced position in
  the navbar (Unit 1).

Unit 2 depends on Unit 1 only because the navbar must list the new
vignette. The two units can be developed and reviewed together in a
single PR.

## 5. Components

### 5.1 Navbar redesign (Unit 1)

#### Before

```
navbar.left: [intro, articles, reference, news]

articles dropdown: 39 items
  - knitted vignettes (3)
  - walk-throughs (4)
  - benchmarks subsection (14)
  - phase 8 sweeps + head-to-heads (8)
  - developer reports (2)
```

#### After

```
navbar.left: [intro, articles, methodology, reference, news]

intro:        articles/getting-started.html  (unchanged href)

articles dropdown: 8 items
  - Get started                        articles/getting-started.html
  - Common pitfalls / FAQ              articles/common-pitfalls.html      [NEW]
  - Mixed-type traits                  articles/mixed-types.html
  - Propagating tree uncertainty       articles/tree-uncertainty.html
  -------- Walk-throughs --------
  - Architecture overview              pigauto_intro.html
  - Mixed-type workflow (Paths A/B/C)  pigauto_workflow_mixed.html
  - Comparative study with covariates  pigauto_walkthrough_covariates.html
  - Multi-observation per species      pigauto_walkthrough_multi_obs.html

methodology dropdown: 23 working items + Phase-2 placeholder
  -------- Per-trait benches (8) --------
  - Continuous (BM/OU/regime/nonlinear)  dev/bench_continuous.html
  - Binary (signal strength sweep)       dev/bench_binary.html
  - Ordinal (level count sweep)          dev/bench_ordinal.html
  - Count (Poisson / NegBin)             dev/bench_count.html
  - Categorical (K-level sweep)          dev/bench_categorical.html
  - Proportion (signal sweep)            dev/bench_proportion.html
  - Zero-inflated counts                 dev/bench_zi_count.html
  - Multi-proportion / compositional     dev/bench_multi_proportion.html
  -------- Cross-dataset benches (7) --------
  - AVONET full (n=9,993)                dev/bench_scaling_v090.html
  - AVONET missingness sweep             dev/bench_avonet_missingness.html
  - Missingness mechanisms (MCAR/MAR/MNAR) dev/bench_missingness_mechanism.html
  - Tree uncertainty (Rubin's rules)     dev/bench_tree_uncertainty.html
  - PanTHERIA mammals                    dev/pantheria_summary.html
  - Delhey 5,809-species birds           dev/bench_delhey.html
  - Multi-observation imputation         dev/bench_multi_obs.html
  -------- Phase 8 sims + head-to-heads (8) --------
  - Phase 8 summary                      dev/phase8_summary.html
  - Phylogenetic signal sweep            dev/bench_signal_sweep.html
  - Cross-trait correlation sweep        dev/bench_correlation_sweep.html
  - Evolutionary-model sweep             dev/bench_evo_model_sweep.html
  - Clade-correlated missingness         dev/bench_clade_missingness.html
  - Covariate effectiveness (sim)        dev/bench_covariate_sim.html
  - AVONET 300 head-to-head vs BACE      dev/bench_bace_avonet_head_to_head.html
  - PanTHERIA head-to-head vs BACE       dev/bench_pantheria_bace_head_to_head.html
  -------- Mathematical & algorithmic description --------
  - (Phase 2: vignettes/methodology.html)  [NOT WIRED IN PHASE 1]

reference: unchanged
news:      unchanged
```

The current Articles dropdown has 23 working benchmark links + 2
developer-internal artefacts ("Logo gallery", "Test catalogue") under
"Developer reports". Phase 1 moves all 23 benchmark links into
Methodology (organised into per-trait, cross-dataset, Phase 8
sub-sections) and drops the 2 developer-internal entries from the
navbar entirely. The developer artefacts remain accessible via the
build outputs but are no longer surfaced to public visitors.

The "Phase-2 placeholder" line is included as a comment in the YAML
so the next contributor understands where the math vignette will land
when Phase 2 ships. It is NOT a working link in Phase 1.

### 5.2 Common Pitfalls / FAQ vignette (Unit 2)

**File**: `vignettes/common-pitfalls.Rmd`. ~150–200 lines. Knits to
`articles/common-pitfalls.html`. Standard pkgdown vignette template.

**Structure** (one section per pitfall):

| Section | Pitfall (in user voice) | Source of evidence |
|---|---|---|
| 1 | "I called impute() and prediction$imputed looks like my input" | Issue #67; PR #69 docstring |
| 2 | "My ordinal trait predicted 100% majority class" | Issue #68; Phase H memo |
| 3 | "The gate stays closed and the GNN seems to do nothing" | CLAUDE.md safety floor section |
| 4 | "How do I know if my dataset has enough phylogenetic signal?" | `bench_signal_sweep.md` Phase 8 result |
| 5 | "Predictions are way bigger than anything I observed" | Phase G `clamp_outliers` memo + Casuarius case |

Each section follows a fixed template:

> **Symptom** — quoted user language or output snippet.
>
> **Why this happens** — 2–3 sentences of mechanism.
>
> **Diagnose** — 1–3 R commands the user can run on their result.
>
> **Fix** — concrete code change with explanation.
>
> **See also** — link to relevant `?function`, `useful/MEMO_*.md`, or
>   Methodology bench HTML.

All examples use bundled data (`avonet300`, `tree300`,
`ctmax_sim`/`tree300`). No new datasets, no network calls, no torch
calls in chunks marked `eval = TRUE`. Heavy chunks marked
`eval = FALSE` with output baked into the Rmd.

**Knit budget**: ≤ 30 s on a Mac M-series, since most chunks are
inert. The one chunk that calls `impute()` (section 1, 30-cell mask
demonstration) reuses the README's example with `epochs = 50L,
n_imputations = 1L` to keep the knit fast. If knit time exceeds 30 s
in practice, mark that chunk `eval = FALSE` and bake the output.

### 5.3 Single source of truth for the index (Unit 3)

- **Delete** `DOCS.md` from the repository root.
- **Remove** the `home.links` entry titled "Documentation index
  (DOCS.md)" in `_pkgdown.yml`.
- **Add** a single line near the top of `README.md`, immediately
  under the package tagline: `> Live documentation:
  <https://itchyshin.github.io/pigauto>`.

The pkgdown navbar plus the auto-generated Reference index together
serve the role `DOCS.md` was trying to play. A reader on GitHub sees
the README; a reader on the live site sees the navbar; both are
canonical and there is nothing to drift.

## 6. Data flow

Not applicable — no R code or runtime behaviour changes. Build
behaviour:

```
pkgdown::build_site_local()
  ├── reads _pkgdown.yml (changed)
  ├── knits vignettes/*.Rmd
  │     └── including new common-pitfalls.Rmd
  ├── writes docs/articles/*.html
  ├── writes docs/index.html (without the DOCS.md link)
  └── writes docs/dev/* (unchanged; built by separate scripts)
```

## 7. Error handling

- If `pkgdown::build_site_local()` fails because a navbar `href`
  points to a non-existent HTML, fix the path or remove the entry.
  The failure is loud and caught at build time.
- If `vignettes/common-pitfalls.Rmd` fails to knit (e.g. an example
  errors), fix the chunk before merge. The vignette is part of the
  R CMD check vignette-build step, so this also fails CI / `R CMD
  check`.

## 8. Testing / verification

| Step | Command | Expected |
|---|---|---|
| 1 | `devtools::document()` | No errors. (Should be a no-op since no R code touched.) |
| 2 | `devtools::check()` | Status: OK. |
| 3 | `pkgdown::build_site_local()` | Succeeds; no broken-link warnings. |
| 4 | Open `docs/index.html` in a browser | Navbar shows `[Get started] [Articles] [Methodology] [Reference] [Changelog]`. The `home.links` no longer mentions `DOCS.md`. |
| 5 | Click each Articles entry (8 items) | All resolve. |
| 6 | Click each Methodology entry (23 items) | All resolve. The Phase-2 placeholder is intentionally absent (commented in YAML). |
| 7 | Click "Common pitfalls / FAQ" | Renders cleanly. All five sections present. All internal `?function` links and external memo links resolve. |
| 8 | `git grep DOCS.md` | Returns no hits in tracked source. |

## 9. Files touched

| File | Change | Lines |
|---|---|---|
| `_pkgdown.yml` | Navbar rewrite; `home.links` entry removed | ~50 modified |
| `vignettes/common-pitfalls.Rmd` | NEW | ~200 added |
| `DOCS.md` | DELETED | 282 removed |
| `README.md` | One new line near the top with live-site link | 1 added |
| `NEWS.md` | New entry under upcoming version: "docs(pkgdown): reorganised Articles + new Methodology dropdown + Common Pitfalls vignette + DOCS.md retired in favour of README + live site as canonical." | ~6 added |

No DESCRIPTION bump. No `R/`, `man/`, `tests/`, or `data/` changes.

## 10. Risks and mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| A user has bookmarked `DOCS.md` on the GitHub repo. | Low (file existed for ~3 weeks) | Acceptable; README points at the live site; GitHub will 404 with the standard message. |
| The Methodology dropdown links 23 HTMLs that were previously rebuilt by separate scripts; some may not be present in `docs/dev/` on a fresh clone. | Medium | Pre-build verification step in §8 catches this. If any HTML is missing, regenerate it via the existing `script/make_bench_*_html.R` driver before the PR is merged. |
| The new vignette knit time pushes `R CMD check` over CRAN's runtime budget. | Low | Vignette is mostly inert chunks. Heavy chunks marked `eval = FALSE`. Knit budget ≤ 30 s. |
| An upstream contributor adds a new bench HTML and forgets to wire it into Methodology. | Medium | Convention documented in `_pkgdown.yml` comments; mention in CLAUDE.md "Repository layout" section. (Out of scope for Phase 1; can be added in a future small PR.) |

## 11. Rollback

If the merge produces broken navigation, revert the PR. Each unit is
independently revertable:

- Unit 1 (navbar) revert restores the 39-item dropdown.
- Unit 2 (Common Pitfalls vignette) revert deletes the file; navbar
  references it but pkgdown will warn about a missing target. Both
  units are normally reverted together.
- Unit 3 (DOCS.md deletion) revert restores `DOCS.md` and the
  `home.links` entry.

## 12. Phase 2 forward reference

The math vignette (item 3 in the original 6-piece bundle, plus item 4
the extension) is written in a separate spec
(`specs/YYYY-MM-DD-pkgdown-reorg-phase2-design.md`) after Phase 1
ships and we have feedback on the navbar reorganisation. The
Methodology dropdown's "Mathematical & algorithmic description"
placeholder lands as a working link in Phase 2.

---

## Self-review checklist (run by author before merge)

1. No `TBD` / `TODO` strings in this spec? — checked, none.
2. The "After" navbar tree in §5.1 lists every link present in the
   current `_pkgdown.yml` minus the explicitly-removed Logo gallery
   and Test catalogue entries? — checked.
3. The Common Pitfalls table in §5.2 references each issue / memo
   that motivates the section? — checked.
4. Files-touched table in §9 covers every file mentioned in §§5.1,
   5.2, 5.3? — checked (`_pkgdown.yml`, `vignettes/common-pitfalls.Rmd`,
   `DOCS.md`, `README.md`, `NEWS.md`).
5. Goals in §2 and out-of-scope in §3 do not contradict? — checked.
   Math vignette is explicitly Phase 2, not Phase 1.
6. Verification steps in §8 are executable as-listed without further
   setup? — checked.
