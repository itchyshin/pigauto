# pigauto 0.11 pkgdown-011-repair — close report

Date: 2026-08-28  
Worktree: `/Users/z3437171/local-scratch/lanes/pigauto-011-ship`  
Branch: `codex/pigauto-0-11-trust-usability` @ `6fddd79` + uncommitted repair edits  
PR: https://github.com/itchyshin/pigauto/pull/174 — **still draft** (`gh pr view`: `isDraft: true`, MERGEABLE)  
Lane: pkgdown-011-repair (lease `cursor:pigauto-011-visual-audit` expanded to vignettes + README)

This run **replaces** the earlier visual-gate-only LOOP. Named G0 defects only. Florence/Tufte **gate**, not redesign. #175 parked. No `R/`, `tests/`, `NEWS.md`, `DESCRIPTION`, `man/` edits. No push / undraft / merge.

## Verdict: PASS (named defects)

Fail-to-ship items in the approved defaults are gone on the rebuilt site. CI green is still not this gate. #174 stays draft until Shinichi undrafts.

## What changed (source)

| # | Default | Change |
| --- | --- | --- |
| 1 | Drop Step 7 GNN-beats-BM table | `vignettes/getting-started.Rmd`: Step 7 now tells readers to read `r_cal`; high-signal morphometrics often close the gate. Fake RMSE comparison table removed. |
| 2 | Tombstone validation suite | `pkgdown/assets/validation_suite.html` replaced in place: withdrawn banner + `noindex`. Copied to `docs/validation_suite.html`. |
| 3 | Body link contrast | `pkgdown/extra.css`: `main` / `.contents` / `.page-header` links `#2c5e4f` (~7.4:1 on white). |
| 4 | Florence warning banners | README keeps GFM `[!WARNING]` (GitHub). Home sidebar `_pkgdown.yml` uses the same marker. `after_body` JS adds `pa-callout-warning` (cream `#fff8eb`, amber `#b45309`). |
| 5 | gnn-architecture claims / §5 | BACE 42%→72% attribution removed. Threshold-joint hands the matrix to the **in-house joint solver**, not `phylopars()`. §5 is three prose blocks (baseline SE, discrete score, conformal). |
| 6 | AGENTS.md / CLAUDE.md | Not touched. Leftover 42%/72% in `AGENTS.md` remains out of scope. |

## What was verified (Read the HTML, not only the Rmd)

- `pkgdown::build_site(preview = FALSE, devel = TRUE)` finished. Log: `/tmp/pigauto-pkgdown-repair/build.log`.
- `pkgdown::check_pkgdown()`: **No problems found.**
- `docs/articles/getting-started.html` Step 7: closed-gate prose; no `bm_rmse` / `gnn_rmse` table.
- `docs/articles/gnn-architecture.html`: “in-house joint solver”; §5 three blocks; no 42%/72%.
- `docs/validation_suite.html`: withdrawn + `noindex`; no pending-bench table.
- `docs/extra.css`: link colour `#2c5e4f`; `blockquote.pa-callout-warning` cream/amber.
- `docs/index.html`: `[!WARNING]` sits in `<blockquote>` plus the after_body script that paints `pa-callout-warning` on load.
- PR #174 still draft.

`R CMD check` was **not** run. This pass did not touch package code; `check_pkgdown()` was the named check.

## Rose fence

No #175 sentinel / Phase A/B numbers added. No new “GNN beats BM” public claim. `r_cal = 0` is described as a valid, expected outcome on high-signal morphometrics.

## Polish-later (not fail-to-ship)

- Warning banners need JS to strip the `[!WARNING]` token; static HTML still shows the marker until the script runs.
- Get started Visualise chunks remain `eval=FALSE` (G0: do not knit).
- Home still builds `AGENTS.md` / `CLAUDE.md` / `VALIDATION_LEDGER.md` / `goodagents.md` as extra pages (inventory leftover; out of this G0).
- `script/validation_suite.html` still has the old suite (not in lease).
- Full Methodology figure rewrite: out of scope.
- Overflow / mobile CSS already in `extra.css` from the parallel visual-audit pass; left in place.

## Melissa (plan vs actual)

See `docs/dev-log/plan-actual/2026-08-27-pkgdown-visual-gate.md`.
