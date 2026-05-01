# Phase 8 MVP Implementation Plan

> Bench infrastructure only — no `R/` code. Execute inline, commit per task.

**Goal:** Ship two reproducible bench scripts + one aggregate HTML report so pigauto's "discriminative across λ" and "beats BACE on AVONET" claims are auditable from a single command.

**Files created:**

| Path | Purpose |
|---|---|
| `script/bench_signal_sweep.R` | Pagel's λ sweep × mixed types × 4 methods |
| `script/make_bench_signal_sweep_html.R` | HTML generator |
| `script/bench_bace_avonet_head_to_head.R` | Side-by-side on avonet300/tree300 |
| `script/make_bench_bace_avonet_head_to_head_html.R` | HTML generator |
| `script/make_bench_phase8_aggregate_html.R` | Aggregate TL;DR report |
| `pkgdown/assets/validation_suite.html` | Add Phase 8 section (regenerated) |

## Task 1: `bench_signal_sweep.R`

- λ ∈ {0.1, 0.3, 0.5, 0.7, 0.9, 1.0}, n_species = 300, n_reps = 3, miss_frac = 0.30.
- Methods: `mean_baseline`, `pigauto_LP`, `pigauto_default`, `pigauto_em5`.
- Per-cell helper `sim_lambda_mixed(tree, lambda, seed)` generates 2 cont + 1 bin + 1 cat K=3.
- Per-trait metrics (continuous: RMSE/r/conformal cov; binary: acc/log-loss/Brier; categorical: top-1/macro-log-loss).
- Output: `script/bench_signal_sweep.{rds,md}`.

## Task 2: `make_bench_signal_sweep_html.R`

- Reads `.rds`, writes `phase8_signal_sweep.html` to both `script/` and `pkgdown/assets/dev/`.
- One line plot per method per trait, facets per trait type.

## Task 3: `bench_bace_avonet_head_to_head.R`

- Uses bundled `avonet300` + `tree300`.
- Splits: `make_missing_splits(missing_frac = 0.30, seed = 2026, trait_map = pd$trait_map)`.
- Methods: `pigauto_default`, `pigauto_em5`, `bace_default` (with `requireNamespace` guard).
- Thin helper `call_bace_on_splits()` inline in the bench script.
- Per-trait metrics for all 7 avonet300 traits.

## Task 4: `make_bench_bace_avonet_head_to_head_html.R`

- Reads `.rds`, renders 7-trait × 3-method table with winner-highlighted cells.
- Wall-clock column per method.

## Task 5: `make_bench_phase8_aggregate_html.R`

- Reads both bench outputs, writes `phase8_summary.html`.
- TL;DR table + embedded (or linked) signal-sweep plot + AVONET table + reproducibility block.

## Task 6: Link from validation suite

- Update `pkgdown/assets/validation_suite.html` (or the generator) to add a Phase 8 subsection.

## Task 7: NEWS + push + PR

- Append Phase 8 bullet to v0.9.1.9000 dev NEWS section.
- Verify R CMD check 0/0/0 (no R/ changes so should be fine).
- `gh pr create`.

## Runtime budget (local CPU)

- Task 1: ~25 min (4 methods × 6 λ × 3 reps × 1 min per fit, mostly overlapping).
- Task 3: ~30 min (BACE-default is the slow piece; pigauto is fast).
- Tasks 2, 4, 5, 6: <1 min each.

Total: ~55 min active wall-clock during development. Background-able.
