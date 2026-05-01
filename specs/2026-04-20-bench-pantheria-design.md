# PanTHERIA benchmark bring-up for pigauto

**Status:** Draft (awaiting user review)
**Author:** Claude Opus 4.7 (working with Shinichi Nakagawa)
**Date:** 2026-04-20

---

## Problem

pigauto's "general-purpose phylogenetic imputation" claim currently
rests on exactly **one** real dataset: AVONET (birds). Every other
real-data claim in the NEWS / validation suite is either a subset of
AVONET or a synthetic simulation. This is a single-point-of-failure.
If pigauto happens to work well on birds but not mammals, we don't know.

The roadmap (Part 1.A) flags **PanTHERIA** (Jones et al. 2009 *Ecology*)
as the HIGH-priority second real dataset: 5,416 mammals × ~55 traits,
with realistic missingness (30–60% per trait). Mammals differ from
birds in trait structure (continuous-heavy, different discrete types,
different ecological drivers) so it's a genuine generalisation test.

## Goal

Ship two reproducible bench scripts that answer:

1. **Does pigauto scale and work on a second real dataset?**
   Full-scale run (n = 5,416 mammals) with `pigauto_default` and
   `pigauto_em5` vs `mean_baseline`. Tests the general-purpose claim.
2. **Does pigauto match BACE on mammals, not just birds?**
   Subset run (n = 500) with `BACE::bace()` + OVR for head-to-head.
   Tests the cross-package parity claim on a second taxon.

Both benches follow the Phase 8 convention: RDS + MD + HTML, linked
from the pkgdown validation suite.

## Non-goals

- **No new R/ code.** Bench infrastructure only.
- **No PanTHERIA bundled as `data/`.** Fetched via a reproducible
  download script and cached in `script/data-cache/` (gitignored). The
  ESA archive at https://esapubs.org/archive/ecol/E090/184/ is stable;
  re-download is idempotent.
- **No multi-seed sweep** in the MVP (single seed = 2026). Multi-seed
  is a Phase PanTHERIA.1 follow-up.
- **No Open Tree / Upham pull.** We use the Fritz et al. 2009 mammal
  supertree (well-known pairing with PanTHERIA).
- **No full 5,416 × BACE run.** MCMCglmm at n = 5k takes hours. BACE
  comparison stays on the n = 500 subset for the MVP.

## Design decisions

### 1. Two bench scripts + aggregator integration

| Script | Purpose | Runtime (est.) |
|---|---|---|
| `script/bench_pantheria_full.R` | pigauto_default + pigauto_em5 + mean_baseline on n = 5,416 | ~45 min local CPU |
| `script/bench_pantheria_bace_head_to_head.R` | pigauto vs BACE on n = 500 subset | ~30 min (BACE dominates) |
| `script/make_bench_pantheria_full_html.R` | HTML render | <1 min |
| `script/make_bench_pantheria_bace_head_to_head_html.R` | HTML render | <1 min |

Each script writes `.rds` + `.md` and calls the HTML generator. Both
are linked from `make_validation_suite_html.R` via the same
RDS-exists → row, else pending pattern the Phase 8 hooks use.

### 2. Data acquisition (`script/fetch_pantheria_and_tree.R`)

One-shot helper that:
- Downloads PanTHERIA_1-0_WR93_Aug2008.txt from
  `https://esapubs.org/archive/ecol/E090/184/PanTHERIA_1-0_WR93_Aug2008.txt`
  to `script/data-cache/pantheria.txt` (skip if already present).
- Downloads the Fritz et al. 2009 mammal supertree from
  `https://datadryad.org/stash/downloads/file_stream/33013` (Dryad
  permanent URL) to `script/data-cache/fritz2009.tre`.
- If either URL is temporarily unreachable, prints a clear error with
  manual-download fallback instructions.

`script/data-cache/` is added to `.gitignore` (already partially
covered by gitignore rules; verify and extend if needed).

### 3. Trait selection for the MVP

PanTHERIA has ~55 traits; we use a canonical 8-trait subset that
captures mixed types and maps cleanly onto pigauto's type dispatcher:

| trait (PanTHERIA column) | pigauto type |
|---|---|
| `X5-1_AdultBodyMass_g` | continuous (log-transform) |
| `X13-1_AdultHeadBodyLen_mm` | continuous (log-transform) |
| `X9-1_GestationLen_d` | continuous (log-transform) |
| `X15-1_LitterSize` | count |
| `X17-1_MaxLongevity_m` | continuous (log-transform) |
| `X6-1_DietBreadth` | ordinal |
| `X12-1_HabitatBreadth` | ordinal |
| `X12-2_Terrestriality` | categorical (K = 3) |

Rationale: all 8 are well-populated (≥ 50% observed on most mammals),
mix types match what pigauto can handle, and the continuous + discrete
split tests both baselines + OVR categorical.

Species missing ALL 8 are dropped. Species missing 7/8 are kept.

### 4. Phylogeny alignment

Fritz 2009 has 5,020 tip labels. PanTHERIA has 5,416 species. We:
1. Normalise both to `Genus_species` form.
2. Take the intersection (typically ~4,100 matching species in
   published pairings; we'll report actual overlap).
3. Prune the tree to matched species via `ape::drop.tip()`.
4. Bench runs on the intersection at whatever n that turns out to be.

### 5. Missingness injection

MCAR at 30% per trait (matches the Phase 8 MVP convention). Seed =
2026. Hold-out per-trait; evaluate per-trait on held-out cells.

### 6. Methods compared (full-scale bench)

- `mean_baseline` — column mean for continuous/count, modal class for
  discrete.
- `pigauto_default` — `em_iterations = 0L`, `epochs = 500L`.
- `pigauto_em5` — `em_iterations = 5L`, `epochs = 500L`.

Expected wall time:
- `fit_baseline` at n ≈ 4,100: ~3–5 min (Rphylopars dominates).
- `fit_pigauto` at n ≈ 4,100, 500 epochs: ~15–25 min.
- Total per method: ~20–30 min.
- 2 pigauto methods + 1 mean = ~45 min total.

### 7. Methods compared (BACE head-to-head, n = 500 subset)

- `pigauto_default`
- `pigauto_em5`
- `bace_default` (OVR, with the same n = 500 split).

The subset is a stratified random 500 species from the matched set,
preserving the original trait-missingness pattern. `BACE::bace()` is
wrapped with `requireNamespace` + `tryCatch` (same pattern as the
AVONET head-to-head); if BACE isn't installed the bench completes
with a skip.

### 8. Metrics

Per trait, per method:
- Continuous / count: **RMSE**, **Pearson r**.
- Ordinal: **RMSE** (on integer class), **top-1 accuracy**.
- Categorical: **top-1 accuracy**, **macro log-loss**.

Headline TL;DR: one row per method × trait, plus a summary mean.

### 9. Output integration

- **RDS + MD**: per bench.
- **HTML**: per bench + update aggregate `phase8_summary.html` if we
  want (debatable — PanTHERIA is arguably its own bench family, not
  Phase 8). Decision: keep the Phase 8 aggregate unchanged and create
  a **new `pantheria_summary.html`** landing page.
- **Validation suite**: add two new rows below the Phase 8 block,
  labelled `PanTHERIA full (mammals)` and `PanTHERIA vs BACE`.
- **NEWS.md**: new `## PanTHERIA mammal benchmark (HIGH priority from
  roadmap)` entry above the Phase 8 section.

### 10. Where it runs

All local CPU. No Vulcan dependency. Scales to ~4,100 species which is
the biggest single pigauto run we've done locally (matches the v0.9.0
validation at n = 10,000 which was on Narval).

## API surface

No user-facing changes. Bench scripts are invoked:

```bash
Rscript script/fetch_pantheria_and_tree.R       # once
Rscript script/bench_pantheria_full.R
Rscript script/bench_pantheria_bace_head_to_head.R
Rscript script/make_bench_pantheria_full_html.R
Rscript script/make_bench_pantheria_bace_head_to_head_html.R
```

## Testing

- **No unit tests** (bench infrastructure).
- **Smoke check**: each bench runs to completion without error.
- **Existing suite stays green**: Phase PanTHERIA doesn't touch
  `R/` or `tests/testthat/`.
- **R CMD check**: 0/0/0 preserved (`script/` is `.Rbuildignore`'d).

## Documentation

1. `NEWS.md` entry:
   > **PanTHERIA mammal benchmark.** Two new scripts:
   > `bench_pantheria_full.R` on the full ~4,100-species match
   > between PanTHERIA traits and the Fritz 2009 supertree, and
   > `bench_pantheria_bace_head_to_head.R` on a 500-species subset
   > with `BACE::bace()` for cross-package parity testing on mammals.
2. Validation suite page: two new rows in the table.
3. Fetch script: `script/fetch_pantheria_and_tree.R` with clear
   download-fallback instructions in the header docstring.

## Backward compatibility

Fully additive. No API change, no default change. New files in
`script/` and `pkgdown/assets/dev/`; one NEWS bullet.

## Success criteria

- [ ] `Rscript script/fetch_pantheria_and_tree.R` acquires both
      datasets into `script/data-cache/` (or fails with a clear
      manual-fallback message).
- [ ] `Rscript script/bench_pantheria_full.R` completes in < 60 min
      and produces `.rds` + `.md`.
- [ ] `Rscript script/bench_pantheria_bace_head_to_head.R` completes
      (or cleanly skips BACE if missing).
- [ ] HTML generators produce readable reports.
- [ ] Validation suite regenerates with 2 new rows.
- [ ] pigauto headline numbers on PanTHERIA are reported honestly
      (whatever they turn out to be — the exercise isn't rigged for
      pigauto to win).
- [ ] R CMD check stays 0/0/0.

## Out of scope (deferred)

- **Multi-seed robustness** on PanTHERIA: Phase PanTHERIA.1.
- **Clade-correlated missingness** on PanTHERIA: Phase PanTHERIA.2.
- **Full n × BACE**: needs Vulcan. Separate PR when home-dir fix lands.
- **FishBase, TRY plant, AmphiBIO**: per roadmap Part 1.A, each deserves
  its own PR. PanTHERIA first; others follow.
