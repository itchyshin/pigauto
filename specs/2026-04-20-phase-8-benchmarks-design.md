# Phase 8: discriminative benchmark suite + AVONET head-to-head vs BACE

**Status:** Draft (awaiting user review)
**Author:** Claude Opus 4.7 (working with Shinichi Nakagawa)
**Date:** 2026-04-20

---

## Problem

Level-C phases 1–7 (plus B1/B2/B3 and Graph Transformer) have shipped, but
the claims in release notes — "pigauto matches or beats BACE on AVONET",
"Level-C baselines lift correlated-discrete accuracy", "scales to 10,000
species" — are drawn from ad-hoc scripts with inconsistent seeds, splits,
and metrics. A reader can't:

1. **Reproduce** the headline numbers from a single command.
2. **Disagree** — the experimental protocol is not pinned down.
3. **Extend** — adding a method or metric means rewriting scripts.

Phase 8 closes that gap by shipping a **discriminative benchmark suite**
(what the roadmap labelled as Phase 8, drafted 2026-04-15). This spec
covers the MVP; follow-ups (Phase 8.1–8.3) handle the other sweeps.

## Goal

Ship **two canonical benchmark scripts + one aggregate report** that
answer the two questions pigauto needs to answer honestly:

1. **Signal-strength behaviour.** Does pigauto discriminate over the
   baseline across the phylogenetic-signal range? Where does it break?
2. **BACE head-to-head on AVONET.** On identical splits / metrics / seed,
   who wins where, and by how much?

## Non-goals (deferred to Phase 8.x)

- **ρ sweep** (cross-trait correlation): Phase 8.1.
- **Evolutionary model sweep** (BM / OU / ACDC / rate-het): Phase 8.2.
- **Clade-correlated missingness**: Phase 8.3.
- **New bench infrastructure refactor** (the `R/sim_traits.R` exported
  API mentioned in the roadmap): deferred; this spec uses inline sim
  inside each bench script for now.
- **BACE parity on non-AVONET datasets** (PanTHERIA, FishBase, etc.):
  Phase 8.x or v0.11.0.

## Design decisions

### 1. Two bench scripts + one aggregator

| Script | Purpose | Runtime (est.) |
|---|---|---|
| `script/bench_signal_sweep.R` | Pagel's λ ∈ {0.1, 0.3, 0.5, 0.7, 0.9, 1.0} × mixed types (2 cont + 1 bin + 1 cat K=3) at n=300 | ~25 min CPU |
| `script/bench_bace_avonet_head_to_head.R` | Side-by-side pigauto + BACE on the bundled avonet300 + tree300 with identical seeds/splits | ~30 min CPU (BACE is slow) |
| `script/make_bench_phase8_aggregate_html.R` | Reads both RDS outputs, writes `phase8_summary.html` linked from the pkgdown validation suite | <1 min |

Each bench script writes `.rds` (machine-readable) + `.md` (human-readable)
and calls `make_bench_*_html.R` downstream.

### 2. Signal-strength sweep (`bench_signal_sweep.R`)

- **Data generator**: pigauto's existing `simulate_non_bm()` extended with
  Pagel's λ transformation applied to the BM covariance matrix:
  `V_lambda = lambda * V_BM + (1 - lambda) * I`. One helper `sim_lambda_mixed()`
  in-script (not exported) scales the tree covariance and draws:
  - 2 continuous traits (BM on rescaled covariance)
  - 1 binary trait (threshold on a latent continuous)
  - 1 categorical K=3 (argmax on K latent continuous, Pagel-scaled jointly)
- **Grid**:
  - λ ∈ {0.1, 0.3, 0.5, 0.7, 0.9, 1.0} (6 levels)
  - n_species = 300 (bundled tree300 size; cheap)
  - miss_frac = 0.30 (matches validation-suite default)
  - n_reps = 3 (seed ∈ {1, 2, 3})
- **Methods per cell**:
  - `mean_baseline`: column mean / mode (non-phylo reference)
  - `pigauto_LP`: pigauto with Level-C disabled (force `em_iterations = 0L`
    and old LP categorical path) — tests the legacy path
  - `pigauto_default`: pigauto v0.9.1.9000 with default args
  - `pigauto_em5`: pigauto with `em_iterations = 5L` (Phase 6) — tests
    the EM lift at each λ
- **Metrics** (per trait per cell):
  - Continuous: RMSE, Pearson r, 95% conformal coverage
  - Binary: accuracy, log-loss, Brier
  - Categorical: top-1 accuracy, macro log-loss
- **Expected result**: all methods converge at λ = 1.0 (pure BM — baseline
  is optimal). Methods separate at moderate λ (0.3–0.7). At λ = 0.1 all
  methods approach `mean_baseline` (no phylogenetic signal to exploit).
  `pigauto_em5` should show the biggest lift at moderate λ on the
  discrete traits.
- **Output**: `script/bench_signal_sweep.{rds,md}`, HTML via
  `make_bench_signal_sweep_html.R`.

### 3. AVONET head-to-head (`bench_bace_avonet_head_to_head.R`)

- **Data**: bundled `avonet300` + `tree300` (not the 9,993-species
  AVONET — the MVP uses the 300 subset for tractability; BACE at 9993
  takes too long to ship in-suite).
- **Splits**: `make_missing_splits(missing_frac = 0.30, seed = 2026, trait_map = ...)`
  — single seed for the MVP (Phase 8.x can add multi-seed).
- **Methods**:
  - `pigauto_default`: v0.9.1.9000 with `em_iterations = 0L` (fast path)
  - `pigauto_em5`: v0.9.1.9000 with `em_iterations = 5L`
  - `bace_default`: `BACE::bace()` with OVR enabled (matches the v0.9.0
    headline comparison setup)
- **Metrics per trait** (7 trait columns in avonet300):
  - Mass, Wing.Length, Beak.Length_Culmen, Tarsus.Length: RMSE, Pearson r
  - Trophic.Level, Primary.Lifestyle (K=3-ish): top-1 accuracy, log-loss
  - Migration (ordinal 3-class): accuracy, ordinal-aware metric if easy
- **BACE call**: wrap `BACE::bace()` with matching missingness patterns
  and extract per-trait predictions. A thin helper
  `call_bace_on_splits(pd, splits)` lives in the bench script (not
  exported from pigauto).
- **Output**: `script/bench_bace_avonet_head_to_head.{rds,md}` with a
  3-method-row, 7-trait-column table plus wall-clock per method.
- **Expected result**: pigauto ≥ BACE on all 7 traits when
  `em_iterations = 5L` (the v0.9.0 claim was Trophic.Level 77% vs 72%;
  we verify + extend with all traits).

### 4. Aggregate HTML report (`make_bench_phase8_aggregate_html.R`)

- Reads both RDS files, produces `phase8_summary.html` with:
  - A **TL;DR table** (one row per method, one column per metric group,
    summary mean across sweeps).
  - **Signal sweep plot** (λ on x-axis, metric on y-axis, one line per
    method; per-trait facets).
  - **AVONET head-to-head table** (7 traits × 3 methods with winner
    highlighted).
  - **Reproducibility block** (exact seeds, R session info, pigauto +
    BACE versions).
- Linked from the existing pkgdown validation suite page
  (`pkgdown/assets/validation_suite.html`) under a new "Phase 8" section.

### 5. Handling BACE being optional

- `BACE` is `Suggests:` in DESCRIPTION. The bench script checks
  `requireNamespace("BACE", quietly = TRUE)` and skips `bace_default`
  with a message when it's missing, so the script still produces a valid
  pigauto-only output.
- Same guard on the aggregator — if BACE wasn't run, render "BACE row
  unavailable — install BACE from in-tree source".

### 6. Where the scripts run

- Signal sweep (~25 min): local.
- AVONET head-to-head (~30 min): local.
- Both are well under the Vulcan threshold (2+ hr per task). No HPC
  needed for the Phase 8 MVP.
- Phase 8.x follow-ups (ρ sweep × evo-model sweep × clade-MAR) will
  multiply cost and probably need Vulcan — but that's a follow-up, not
  this spec.

### 7. No new `R/` code

Phase 8 is benchmark infrastructure. Everything lives under `script/`.
No changes to exported API, DESCRIPTION, or roxygen. NEWS.md gets a
single bullet linking to the aggregate report.

## API surface

No user-facing changes. Bench scripts are invoked via:

```bash
Rscript script/bench_signal_sweep.R
Rscript script/bench_bace_avonet_head_to_head.R
Rscript script/make_bench_phase8_aggregate_html.R
```

## Testing

- **No unit tests for bench scripts** (they're integration scripts).
- **Smoke check**: each bench script runs to completion without error on
  its smallest cell (e.g. λ = 0.5, n_reps = 1) during development.
- **Existing suite stays green**: Phase 8 doesn't touch `R/` or
  `tests/testthat/`.
- **R CMD check**: 0/0/0 preserved (script/ is in `.Rbuildignore`).

## Documentation

1. `NEWS.md` entry in the v0.9.1.9000 dev section:
   > **Phase 8 benchmark suite (MVP).** Two reproducible bench scripts:
   > `bench_signal_sweep.R` (Pagel's λ ∈ {0.1..1.0} × mixed types) and
   > `bench_bace_avonet_head_to_head.R` (pigauto vs BACE on
   > avonet300/tree300). Aggregate report at
   > `pkgdown/assets/dev/phase8_summary.html` summarises both.
   > Follow-ups (ρ sweep, evo-model sweep, clade-MAR) deferred to
   > Phase 8.x.
2. Validation suite page (`pkgdown/assets/validation_suite.html`):
   new "Phase 8: discrimination + BACE comparison" subsection linking
   to the aggregate report + the two per-bench HTMLs.
3. Each bench script gets a top-of-file docstring describing its
   protocol, output, and how to reproduce.

## Backward compatibility

Fully additive. No API change, no default change. Existing pkgdown site
gains one new subsection; nothing is removed or renamed.

## Success criteria

- [ ] `Rscript script/bench_signal_sweep.R` completes in <30 min local
      CPU and produces `.rds` + `.md`.
- [ ] `Rscript script/bench_bace_avonet_head_to_head.R` completes
      (or cleanly skips BACE if missing) and produces `.rds` + `.md`.
- [ ] `Rscript script/make_bench_phase8_aggregate_html.R` produces a
      `phase8_summary.html` that renders in a browser.
- [ ] Validation suite page links to the new report.
- [ ] `R CMD check` stays 0/0/0.
- [ ] `pigauto_em5` shows ≥2pp accuracy lift vs `pigauto_default` on
      at least one discrete trait at λ ∈ {0.3, 0.5} (smoke-level
      expected lift; actual magnitude determined by empirical run).
- [ ] AVONET vs BACE table reproduces the v0.9.0 headline claims
      (Trophic.Level ≥72%, Primary.Lifestyle ≥72%) with the new pipeline.

## Out of scope (deferred further)

- **ρ sweep** (cross-trait correlation): Phase 8.1. New script
  `bench_correlation_sweep.R`.
- **Evolutionary-model sweep**: Phase 8.2.
- **Clade-correlated missingness**: Phase 8.3.
- **Full AVONET 9,993 vs BACE**: likely needs Vulcan; deferred.
- **PanTHERIA, FishBase, AmphiBIO, EltonTraits, Tree-of-Sex**: each is
  a dataset bring-up task worth its own issue; out of Phase 8 MVP
  scope.
- **Exported `sim_traits()` API**: roadmap Part 2.B; deferred.
