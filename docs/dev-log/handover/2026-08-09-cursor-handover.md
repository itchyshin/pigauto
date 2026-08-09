# Session Handoff: pigauto BACE wrap lane (Cursor)

Meta: 2026-08-09 ~08:15 MDT · from Cursor Auto (this chat) · TARGET = Cursor  
Authoritative copy: this file. Chat is disposable.

You are **Cursor**, picking up a **new pigauto lane**. You inherit no chat history.

## Critical Context

1. **Remit is pigauto only.** Sister package BACE lives at `/Users/z3437171/Dropbox/Github Local/BACE` (`daniel1noble/BACE@ce8bc87`, in sync) and a stale nested clone `pigauto/BACE` (`de87d8c`, **86 commits behind**). **Do not modify either BACE tree.** AGENTS.md: do not treat in-tree `BACE/` as pigauto work. Do not vendor-sync unless Shinichi explicitly asks.
2. **Two lanes.** P0/#155 is still live (CI in progress). This handover starts the **BACE-wrap / re-bench** lane. Do not mix commits onto `fix/p0-review-blockers` / PR 155. Do not rebase onto `main` (43 behind; locked).

## Goals / mission

pigauto remains a mixed-type phylogenetic imputer (gated BM + GNN). BACE is the MCMCglmm sister. Dan Noble (`aut`; Shinichi `cre`) just fixed BACE MI propriety + priors + poisson points and published Study B (~10k sims). pigauto’s wrap still calls **`BACE::bace_imp()` only** (`R/fit_baseline_bace.R:92`), not `bace()` + `bace_final_imp()` — the path Dan validated. Any pigauto-vs-BACE sentence needs a re-bench against `@ce8bc87`.

## Plans / roadmap (beyond immediate steps)

- P1 claim-gate (DESCRIPTION still says joint-Σ + 95% coverage) — **not this lane**.
- Covariates into joint baseline (P1-8) — **not this lane**.
- Getting `fix/ci-install-libtorch` onto `main` — later G0.
- Restore EM (`max_iter>0`) — **forbidden** without new G0.

## What Was Accomplished (prior session, pigauto)

- P0 blockers + B1–B3 GNN train/cal on `fix/p0-review-blockers`: `58708ac`, hygiene `0ed6a73`, merge-from-parent `5800b01`.
- PR **#155** → base `fix/ci-install-libtorch` (not `main`): https://github.com/itchyshin/pigauto/pull/155
- Focused P0 tests green. Full testthat `FAIL 3 | WARN 289 | SKIP 32 | PASS 1661`. Local `--as-cran` 2E/2W/4N (none in P0 files).
- Rose READY narrow: `docs/dev-log/handover/2026-08-08-p0-rose.md`. Kill list still holds.
- Dan BACE synthesis: `docs/dev-log/handover/2026-08-09-bace-dan-progress.md` (HTML tables, not Slack).

## Current Working State

- **Working:** P0 code on `origin/fix/p0-review-blockers` @ `5800b01`. Standalone BACE @ `ce8bc87`.
- **In progress:** PR 155 R-CMD-check `ubuntu-latest (release)` IN_PROGRESS (run 31317655200). A merge agent may still be polling; **check before duplicating**. Merge 155 → parent only if CI green **or** red only for known non-P0 (PDF `Σ`, `jsonlite` in `test-gha-helpers`, missing Suggests). Stop if P0 test files fail.
- **Not working / blocked:** wrap ≠ Dan’s validated MI path; in-tree BACE stale; local BACE benches skipped / GHA snapshot 2026-05-13 pre-fix; DESCRIPTION overclaim (P1, other lane).

## Key Decisions & Rationale

- P0 only; no EM restore; no rebase onto main; B1–B3 kept.
- PR 155 base = parent `fix/ci-install-libtorch`, not `main`.
- BACE wrap/`bace_final_imp` / re-bench = **separate G0** (Shinichi confirmed).
- Do not ask Dan to re-run 10k. One wrap question if pinging (see Dan note).
- pigauto OVR ≠ BACE’s new multinomial default. Naming-only; do not copy BACE OVR-off into pigauto without a measurement.

## Landing State

Gate (`handoff_gate.sh`) **FAIL** — declared below. Other unpushed local branches (winbuilder, CRAN comments, experiments) are **out of remit**; do not touch.

| Artifact / branch | Committed | Pushed | PR | State |
|---|---|---|---|---|
| pigauto `fix/p0-review-blockers` `5800b01` | y | y | #155 OPEN → `fix/ci-install-libtorch` | **CARRIED-OVER** — wait CI; merge to parent per acceptance; **never merge to main** |
| pigauto dirty: `.gitignore` `AGENTS.md` `CLAUDE.md` `README.md` `_pkgdown.yml` `LOOP/checkpoint.md` | n | n | — | **CARRIED-OVER** — uinit + P1 banners + LOOP noise. **Do not stage.** |
| `dev/gnn_attribution_*` + `script/*gnn_attribution*` + `script/returned_gnn_attr/` | n | n | — | **CARRIED-OVER** — DRAC Fir ladder artifacts. **Do not stage.** |
| `docs/dev-log/handover/2026-08-08-*.md`, `2026-08-09-bace-dan-progress.md` | n (`docs/` gitignored) | n | — | **CARRIED-OVER** — local-only unless force-added. This handover file is force-added on `handover/2026-08-09-cursor`. |
| standalone BACE `@ce8bc87` | n/a | n/a | — | **PROTECTED** — do not edit |
| `pigauto/BACE` `@de87d8c` | n/a | n/a | — | **PROTECTED** — do not vendor-sync |

**Resume 155 (if still OPEN):**
```bash
cd "/Users/z3437171/Dropbox/Github Local/pigauto"
gh pr checks 155 --repo itchyshin/pigauto
gh pr view 155 --json state,mergeable,statusCheckRollup
# Merge to parent only after acceptance criteria; never to main:
# gh pr merge 155 --merge --repo itchyshin/pigauto
```

## Next Immediate Steps (OWED — BACE wrap lane only)

1. Run `tools/lane_preflight.sh` if present; `git status` + `gh pr view 155`. Classify items OWED/DONE/RETRACTED/PROTECTED. Do not redo P0. Do not finish 155 unless it is still OPEN and CI has concluded and no other agent already merged it.
2. Read: `AGENTS.md`, this file, `docs/dev-log/handover/2026-08-09-bace-dan-progress.md`, `R/fit_baseline_bace.R`, `docs/dev-log/coordination-board.md`.
3. **`/ultra-plan` G0 (STOP for approval)** — pigauto wrap: stay on cheap `bace_imp` chain averages **or** move to `bace()` + `bace_final_imp` to match Dan’s MI-proper path. Scope: `R/fit_baseline_bace.R` + tests + NEWS. Out of scope: editing BACE source, vendor-sync, DESCRIPTION 95%/joint-Σ rewrite, EM restore, DRAC, GNN.
4. After G0: install **standalone** BACE `@ce8bc87` (not in-tree). Re-bench wrap vs Dan’s object before any comparison claim. Existing scripts: `script/bench_avonet_bace.R`, `bench_bace_avonet_head_to_head.R`, `bench_pantheria_bace_head_to_head.R` (last local runs skipped).
5. Optional docs-only in same G0: pigauto OVR comment/roxygen no longer “the OVR strategy BACE uses” as current default. Footnote, not an estimator change.

## Blockers / Open Questions

- Wrap API choice needs Shinichi G0 (and optionally one Slack to Dan — wrap stay vs `bace_final_imp`).
- Installed BACE on this Mac may not be `@ce8bc87`. Verify `packageVersion` / git SHA before benches.
- #155 CI may still be red for Suggests-forced R CMD check. That does not block starting the wrap **plan**.

## Gotchas & Failed Approaches

- Slack 0.75–0.88 corr ≠ poisson Pearson (0.687/0.744). Use HTML Study B table.
- `bace_imp` already reinstated NAs; Dan’s imputed-as-observed bug was **`bace_final_imp`**. Switching wrap without understanding that double-counts the “fix.”
- BACE OVR default is now **FALSE**. pigauto OVR exists so Rphylopars stays full-rank. Different estimator.
- `docs/` is `.gitignore`d. Force-add any handover that must travel. Do not `git add -A`.
- Never `ScheduleWakeup` / archive this project. Never rebase P0 onto main.

## How to Resume (Cursor)

Working directory: `/Users/z3437171/Dropbox/Github Local/pigauto`  
Toolchain: R 4.6 + `devtools::load_all()` + `NOT_CRAN=true`. Torch required for GNN tests; wrap lane should not need GNN if scoped to `fit_baseline_bace`.  
Safe verify: `NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-fit-baseline-bace.R")'` (adjust if file name differs).  
Must not stage: uinit, README/_pkgdown banners, DRAC `script/`/`dev/` gnn attribution, other people’s unpushed branches.

Rehydrate:

```text
Read AGENTS.md and docs/dev-log/handover/2026-08-09-cursor-handover.md. Run the handover rehydration steps, reconcile them with the current git state, then continue only the OWED Next Immediate Steps.
```

Checkout for wrap work (do not use `fix/p0-review-blockers` for new commits):

```bash
git fetch origin
git checkout handover/2026-08-09-cursor
git pull
```

## Mission control

| | |
|---|---|
| Repo | `/Users/z3437171/Dropbox/Github Local/pigauto` |
| P0 branch / PR | `fix/p0-review-blockers` @ `5800b01` · #155 → `fix/ci-install-libtorch` · CI pending |
| Wrap handover branch | `handover/2026-08-09-cursor` |
| What shipped | P0 + B1–B3 + parent merge; Dan BACE note |
| This lane | pigauto wrap vs `bace_final_imp` + re-bench `@ce8bc87` · G0 first |
| Do not | edit BACE · vendor-sync · merge to main · restore EM |

> Related: `docs/dev-log/handover/2026-08-09-bace-dan-progress.md` · `docs/dev-log/handover/2026-08-08-p0-rose.md` · `docs/dev-log/coordination-board.md`
