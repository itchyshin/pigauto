# PanTHERIA Benchmark Implementation Plan

> Bench infrastructure only. No `R/` code. Execute inline, commit per task.

**Goal:** Second real-data benchmark for pigauto: PanTHERIA mammals (4,100+ species × 8 mixed traits) head-to-head against BACE on a tractable subset.

**Files:**

| Path | Purpose |
|---|---|
| `script/fetch_pantheria_and_tree.R` | One-shot downloader, caches to `script/data-cache/` |
| `script/bench_pantheria_full.R` | Full-scale bench (pigauto + mean_baseline) |
| `script/make_bench_pantheria_full_html.R` | HTML render |
| `script/bench_pantheria_bace_head_to_head.R` | BACE parity on n=500 subset |
| `script/make_bench_pantheria_bace_head_to_head_html.R` | HTML render |
| `script/make_pantheria_summary_html.R` | Aggregate landing page |
| `.gitignore` | Add `script/data-cache/` if not already |
| `script/make_validation_suite_html.R` | 2 new rows for PanTHERIA |
| `NEWS.md` | PanTHERIA entry |

## Task 1: `fetch_pantheria_and_tree.R` + `.gitignore`

- Download `PanTHERIA_1-0_WR93_Aug2008.txt` from ESA archive.
- For the tree: try `rotl::tol_induced_subtree()` with PanTHERIA species; it returns a tree from Open Tree of Life.
- Cache both to `script/data-cache/` (gitignored).
- Print species-overlap stats (tree tips ∩ trait rows).

## Task 2: `bench_pantheria_full.R`

- Load cached data + tree.
- Align species (intersection), select 8 canonical traits, handle `-999` NA coding.
- MCAR 30% per trait, seed = 2026.
- 3 methods: mean_baseline, pigauto_default (em=0), pigauto_em5 (em=5).
- Output: RMSE, Pearson r per continuous/count; accuracy per discrete.
- Save `.rds` + `.md`.

## Task 3: `make_bench_pantheria_full_html.R`

- Per-trait metric table with winner highlighting.
- Wall-time column.

## Task 4: `bench_pantheria_bace_head_to_head.R`

- Stratified 500-species subset from the aligned set.
- pigauto_default + pigauto_em5 + `BACE::bace()` (with `requireNamespace` guard).
- Same splits across all three.
- Same metrics.

## Task 5: `make_bench_pantheria_bace_head_to_head_html.R`

- Winner-highlighted pivot table (method × trait × metric).

## Task 6: `make_pantheria_summary_html.R`

- Aggregate report: TL;DR table + sub-report links + reproducibility block.
- Parallel to `phase8_summary.html`.

## Task 7: Validation suite + NEWS

- Add 2 rows to `make_validation_suite_html.R` (full + BACE head-to-head).
- NEWS entry above Phase 7 section (or above Phase 8 if PR #38 merges first — rebase-safe).

## Task 8: Run everything locally

- Fetch.
- Full-scale bench (~45 min).
- BACE head-to-head (~30 min; BACE likely skipped in this env, but the script should run without erroring).
- Regenerate HTMLs + validation suite.
- Verify output renders.

## Task 9: R CMD check + push + Draft PR

- `devtools::check()` — should stay 0/0/0 (script/ is .Rbuildignore'd).
- `git push` + `gh pr create --draft` with the "⏸ hold until May 1" marker.
