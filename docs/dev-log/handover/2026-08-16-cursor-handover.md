# Session Handoff: pigauto Tier 0 CLOSED, P1 queue 8/9, Totoro campaign IN FLIGHT (Cursor)

Meta: 2026-08-16 ~05:40 MDT · from Claude (`AUTHOR = claude`) · `TARGET = cursor`
Authoritative copy: this file. Chat is disposable. You inherit **no** authoring chat.

You are **Cursor**, picking up **pigauto only**. A long overnight Claude lane closed Tier 0
and landed 8 PRs. **One campaign is still running on Totoro** — do not restart it. Do not
re-merge anything. Do not CRAN-submit.

## Critical Context

1. **`origin/main` = `2ec4af3`**, Version `0.10.0.9000`, Suggests `BACE`. Eight PRs merged
   overnight (#158–#165). CI green throughout, including post-merge pkgdown.
2. **A Totoro campaign is IN FLIGHT.** `~/pigauto_regime_map/` on
   `snakagaw@totoro.biology.ualberta.ca`, launched 2026-08-16 ~05:15, 5,400 jobs, 100
   workers, 0 failures at last check. **It is resumable and must not be relaunched from
   scratch** — `regime_cell.R` skips existing outputs, so a re-run of `run_campaign.sh`
   resumes rather than duplicates, but only if you leave `results/regime_map/` intact.
3. **Estimated 1.5 h, actually ~2.5 h** — logged as a D-139 overrun in
   `docs/dev-log/2026-08-16-campaign-go-no-go.md`. Not stopped: resumable, zero failures,
   103 workers ≤ the D-143 cap of 150, Totoro otherwise idle.
4. **The `main-direct` lane the preflight reports is NOT a foreign agent.** It is the
   overnight Claude lane's own PR merges, which commit under Shinichi's GitHub identity and
   therefore read as an unowned lane. Verified via `git log origin/main --format="%an"`.
   Do not go looking for a third party.
5. **PROTECTED / do not touch:** standalone BACE `/Users/z3437171/Dropbox/Github Local/BACE`
   `@ce8bc87`; in-tree `pigauto/BACE` `@de87d8c`; the 15 dirty uinit/GNN items on this
   checkout; the ~20 foreign unpushed branches; EM (`max_iter > 0`); DRM.jl. Never
   `git add -A`. Never archive this project.

## Goals / mission

pigauto remains a mixed-type phylogenetic imputer (gated BM/LP + GNN). Its five
architectural advantages, in order: (a) unified mixed-type imputation; (b) calibrated
conformal intervals; (c) multi-tree Rubin pooling; (d) multi-obs `obs_refine`;
(e) `suggest_next_observation()`. Any proposal that compromises one needs an explicit case.

The current program is **publication readiness**, mapped in
`docs/dev-log/2026-08-15-publication-readiness-roadmap.md` with a tick-list at
`docs/dev-log/2026-08-15-publication-checklist.md`. **Tier 0 is closed.** Tier 1 (GNN
evidence) is in flight. Tier 2 (latent multivariate) and Tier 3 (packaging) are open.

## Plans / roadmap (beyond the immediate steps)

- **Tier 1.2 regime map** — running now. Design: `2026-08-15-regime-map-design.md` (ADEMP;
  Morris 2019 + Williams 2024 audit). 108 cells, paired blend-vs-baseline estimand,
  gate-firing diagnostics.
- **Tier 1.3 coverage campaign** — designed and priced, NOT launched:
  `2026-08-15-coverage-campaign-design.md`. 24 cells, ~0.5 h at 100 workers. Shinichi
  approved launching both; only the regime map has actually been started.
- **Tier 2.1** — decide whether to scope the joint-Σ claim or restore the cell-refinement EM
  under a numerical gate. **Read the correction in the roadmap first**: Σ *is* estimated in
  closed form; only the refinement EM is off, because it diverged (0.53 vs 0.93 argmax).
- **Tier 2.2** — no known-Σ recovery simulation exists for the multivariate machinery.
- **P1-8** — covariates ignored by the joint/threshold-joint baseline. The last open Rose P1
  item, deliberately untouched: it is a real feature with a design choice in it (thread
  covariates through vs detect-and-message). **Shinichi's call, not a lane's.**
- Next CRAN cut must drop Suggests `BACE` or wait for BACE on CRAN (still 404).

## What Was Accomplished

Full narrative: `docs/dev-log/2026-08-16-overnight-report.md`. Summary:

| PR | Item |
|---|---|
| [#158](https://github.com/itchyshin/pigauto/pull/158) | P0 honesty fixes + two-layer #157 gate floor + 9-edit doc bundle (**closes #157**) |
| [#159](https://github.com/itchyshin/pigauto/pull/159) | pkgdown destination fix (post-merge fallout from #158) |
| [#160](https://github.com/itchyshin/pigauto/pull/160) | P1-7, P1-6 |
| [#161](https://github.com/itchyshin/pigauto/pull/161) | P1-13 + P1-10 fail-closed guards |
| [#162](https://github.com/itchyshin/pigauto/pull/162) | P1-5 joint-solver attribution |
| [#163](https://github.com/itchyshin/pigauto/pull/163) | P1-12 multi-obs phylo-signal gate |
| [#164](https://github.com/itchyshin/pigauto/pull/164) | P1-9 zi_count zeros in splits |
| [#165](https://github.com/itchyshin/pigauto/pull/165) | P1-11 documented (deliberately not fixed) |

**Three silent bugs, each proven against pre-fix code** — the discipline that matters here:
run the defect on the pre-fix commit before claiming a fix. #157 (surface mismatch),
P1-12 (λ was `NA`, now > 0.5), P1-9 (23 observed zeros, 0 ever held out).

**Bench re-run:** the held-out-context leak had **no detectable effect** on 6 of 8 trait
types (all < 1.15 MCSE, paired per-rep). Details + caveats:
`docs/dev-log/2026-08-16-bench-rerun-results.md`.

**Brain lessons filed** (D-37 approved): vault commit `0cf87d1`, 3 → `WHAT-WORKS`,
1 → `LESSONS`. The renice/fork lesson was already filed by the gllvmTMB session; not
duplicated.

## Current Working State

- **Working:** `origin/main` `2ec4af3`, all CI green, `--as-cran` 0 errors / 0 warnings /
  1 known note.
- **In progress (do not disturb):**
  - Totoro regime-map campaign (see Critical Context 2).
  - Two local bench retries on `.worktrees/bench-rerun` — `zi_count` (failed once on a PSOCK
    worker death) and `multi_proportion` (its stale `.rds` was moved to
    `/tmp/stale_mp_rds_backup.rds` so it cannot phantom-resume). If these are still running
    when you start, let them finish or kill them; nothing depends on them.
- **Not working / blocked:** CRAN submission (Suggests BACE, CRAN BACE 404). P1-8 pending a
  Shinichi decision. Coverage campaign approved but not launched.

## Key Decisions & Rationale

- **#157 lives on surface B.** The val-floor invariant is asserted on the surface a genuinely
  missing cell sees (val+test hidden), because that is what production delivers. Shinichi
  locked this.
- **Margin-based, not strict, BM floor.** `cal_min_rel_gain / 2` (1%). A strict check was
  tried in April 2026 and over-corrected on noise-level near-ties; the margin separates that
  regime from material losses.
- **P1-11 documented, not fixed.** Threading `draws_method` through `multi_impute_trees()`
  changes how a headline feature generates draws and embeds a default choice. Shinichi's.
- **The bench re-run relaxed the roadmap's blanket ban** on pre-fix numbers — cite the re-run
  values, noting the pre-fix ones agree within noise. It is **not** a claim the leak was
  harmless; power bounds effects only above ~0.01–0.02.
- **`r_cal = 0` remains a valid fallback.** Do not remove gate safety.

## Files Created / Modified

**Code on `main` (via the 8 PRs), diff `416561b..2ec4af3`:**
`R/fit_baseline.R` `R/fit_helpers.R` `R/fit_pigauto.R` `R/henderson_s_inv.R` `R/impute.R`
`R/mask_missing.R` `R/multi_impute.R` `R/multi_impute_trees.R` `R/phylo_signal.R`
`R/predict_pigauto.R` `R/preprocess_traits.R` · `DESCRIPTION` `NEWS.md` `README.md` ·
`man/{fit_baseline,fit_pigauto,impute,multi_impute,multi_impute_trees,predict.pigauto_fit,preprocess_traits}.Rd` ·
`vignettes/{gnn-architecture,getting-started,mixed-types}.Rmd` ·
`tests/testthat/test-{bm-internal,covariate-alignment,gbif-centroids,gnn-train-cal-symmetry,joint-baseline,joint-threshold-baseline,safety-floor,zi-count-conformal-mi,type-detection,phylo-signal-multiobs,zi-count-splits}.R`

**Docs on `handover/2026-08-09-cursor` (this branch):**
`docs/dev-log/2026-08-15-{brain-lesson-drafts,coverage-campaign-design,publication-checklist,publication-readiness-roadmap,regime-map-design}.md` ·
`docs/dev-log/2026-08-16-{bench-rerun-results,campaign-go-no-go,overnight-report}.md` ·
`docs/dev-log/after-task/2026-08-15-p0-land.md` ·
`docs/dev-log/plan-actual/2026-08-15-p0-land-reconcile.md` ·
`docs/dev-log/arc/2026-08-15-{p0-onto-main-arc-notes.md,p0-claim-gate-findings.md,floor-diagnostic.R,delta-surface-compare.R,gate-open-check-on-main.R,train-mask-heldout-experiment.R}` ·
`docs/dev-log/{coordination-board,phase-snapshot}.md` · this file.

**On Totoro** (`~/pigauto_regime_map/`): `regime_cell.R`, `run_campaign.sh`,
`regime_prerun_cell.R`, `results/`, `prerun/`.

## Landing State

`handoff_gate.sh` → **FAIL** (2026-08-16 ~05:40): 15 uncommitted on this branch + ~40
unpushed commits on ~20 other branches + a stale `index.lock`. All declared below; nothing
silently dropped.

| Artifact / branch | Committed | Pushed | PR | State |
|---|---|---|---|---|
| `origin/main` `2ec4af3` | y | y | #158–#165 all MERGED | **LANDED.** Do not re-merge. |
| `origin/handover/2026-08-09-cursor` (docs, this handover) | y | y | none — do **not** open a PR into main | **LANDED** process lane. |
| Totoro regime-map campaign | n/a | n/a | — | **IN FLIGHT.** Resume: see below. Do not relaunch from zero. |
| Coverage campaign | n/a | n/a | — | **CARRIED-OVER.** Approved but never started. Command below. |
| `.worktrees/bench-rerun` zi_count / multi_proportion retries | n | n | — | **CARRIED-OVER.** Local only, nothing depends on them. |
| 15 dirty uinit/GNN items on this checkout | n | n | — | **CARRIED-OVER — do not stage.** `.gitignore` `AGENTS.md` `CLAUDE.md` `README.md` `_pkgdown.yml` + `dev/gnn_attribution_*` + `script/*gnn_attribution*` + `script/returned_gnn_attr/` |
| ~40 unpushed commits on ~20 foreign branches | mixed | **n** | — | **CARRIED-OVER.** Not this lane. Do not push, rebase, or tidy. |
| Ghost `D`/`MM` entries in `git status` under `docs/` | — | — | — | **ARTEFACT, not deletions.** From a `GIT_INDEX_FILE` workaround around the stale lock. Files exist on disk and on origin — verified. One `git reset` clears them once the lock is gone. |
| stale `.git/index.lock` | — | — | — | **REPORT-ONLY.** Cursor gitWorker artefact, `lsof`-verified unheld. Harness blocks `rm`. Fix below. |

## Next Immediate Steps (OWED only)

Classify every item `OWED` · `DONE` · `RETRACTED` · `PROTECTED` against **current** git/gh
before acting. Re-merging PRs, CRAN submission, P1-8, EM restore, vendor-sync, DRM.jl, and
staging the dirty items are **not OWED**.

1. **Rehydrate + classify.** `~/shinichi-brain/tools/lane_preflight.sh "<repo>"`;
   `git fetch origin`; `git status`; `gh pr list --state merged --limit 10`; read the
   coordination board's **both** Active Lane Split rows. Expect `main-direct` in the census —
   that is the overnight lane's own merges (Critical Context 4), not a foreign agent.
2. **Check the Totoro campaign** — this is the one genuinely time-sensitive item:
   ```bash
   ssh snakagaw@totoro.biology.ualberta.ca 'cd ~/pigauto_regime_map; \
     echo "results=$(ls results/regime_map | wc -l)/5400"; \
     echo "workers=$(pgrep -fc regime_cell.R || echo 0)"; \
     echo "fails=$(grep -c "^FAIL " campaign_regime_map.log || echo 0)"'
   ```
   - Still running → leave it. Check again later.
   - Finished (workers 0, results ≈ 5400, fails 0) → **summarise it** into
     `docs/dev-log/` per the design's Reporting section, then launch the coverage campaign:
     `ssh snakagaw@totoro.biology.ualberta.ca 'cd ~/pigauto_regime_map && nohup bash run_campaign.sh coverage 100 > launch_coverage.log 2>&1 &'`
   - Failures > 0 → **stop and report**; do not analyse a partial grid as if complete.
3. **Finish the bench suite (optional, low value).** `zi_count` and `multi_proportion` are the
   only two of eight not re-run. Note `zi_count` needs a **third** run against today's `main`
   anyway, because P1-9 changed which cells are eligible for its splits.
4. **Do not start P1-8 or Tier 2** without Shinichi naming it.

## Blockers / Open Questions

- **P1-8** — the last open Rose P1 item; needs a design decision from Shinichi.
- **CRAN** — blocked while Suggests includes `BACE` and CRAN BACE is 404.
- **Issue [#135](https://github.com/itchyshin/pigauto/issues/135)** — OPEN; the regime map is
  what closes it.
- **Campaign overrun** — ~2.5 h actual vs 1.5 h estimated. Root cause recorded: the pre-run
  timed the `impute()` call only, not the whole job.
- The coordination board and this handover live on `handover/2026-08-09-cursor`, **not** on
  `origin/main`. A lane reading only `origin/main` cannot see them.

## Gotchas & Failed Approaches

- **Never `renice` a process that will later fork PSOCK workers.** Niceness is inherited at
  fork, so the workers miss `makePSOCKcluster()`'s connect window. This cost ~2 h across two
  sessions. Bound a cluster campaign by **worker count** (`MC_CORES`), never by priority.
- **`parLapply` logs nothing until every cell finishes.** A master at 0% CPU with a log
  frozen at "Dispatching cells" is **normal**, not hung. Healthy drivers were killed on this
  misread. `bench_proportion` legitimately took 103 min.
- **torch can hang the R process AFTER the work is written.** Drive chains off the driver's
  own completion marker in the log, not off process exit.
- **The bench drivers resume from their own output `.rds`.** A stale copy makes a "re-run"
  silently compare old numbers to themselves — it logs `Running 0 cells` and `0.0s`. Always
  verify non-zero cells and non-trivial wall time before comparing.
- **An experiment script's masking choice is a measurable assumption.** A val-only mask
  (leaving test cells pinned) flipped a conclusion's sign versus the package convention
  `c(val_idx, test_idx)`.
- `docs/` is `.gitignore`d — `git add -f` specific paths only. Never `git add -A`.
- Do not rebase `handover/2026-08-09-cursor` onto main.

## How to Resume (Cursor)

Working directory: `/Users/z3437171/Dropbox/Github Local/pigauto`
Branch: stay on `handover/2026-08-09-cursor` for docs. Code is on `origin/main` `2ec4af3`.
You inherit no chat, no extensions, no credentials from the authoring tool. Verify your own
toolchain before running anything.

Toolchain: R + `devtools` + `torch`. Env: `NOT_CRAN=true` for tests;
`MC_CORES=<n>` caps bench worker pools; `OPENBLAS_NUM_THREADS=1` on Totoro.

Safe verification (does not touch main, ~1 min):
```bash
cd "/Users/z3437171/Dropbox/Github Local/pigauto"
git fetch origin
git show origin/main:DESCRIPTION | grep -E 'Version|BACE'   # expect 0.10.0.9000, BACE
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-type-detection.R")'
```

Housekeeping Shinichi must do (the harness blocks it for agents) — once Cursor's git worker
is idle:
```bash
rm "/Users/z3437171/Dropbox/Github Local/pigauto/.git/index.lock"
git -C "/Users/z3437171/Dropbox/Github Local/pigauto" reset
```
The `reset` also clears the ghost `D`/`MM` entries described in Landing State.

Must not stage: `.gitignore` `AGENTS.md` `CLAUDE.md` `README.md` `_pkgdown.yml` ·
`dev/gnn_attribution_*` · `script/*gnn_attribution*` · `script/returned_gnn_attr/` ·
either BACE tree · any foreign branch.

Read order: `AGENTS.md` → this file → `docs/dev-log/2026-08-16-overnight-report.md` →
`docs/dev-log/coordination-board.md` (**both** Active Lane Split rows) →
`docs/dev-log/2026-08-16-campaign-go-no-go.md` → `docs/dev-log/2026-08-15-publication-checklist.md`.

Rehydrate (human starts a **fresh Cursor agent** in the repo and pastes):

```text
Read AGENTS.md and docs/dev-log/handover/2026-08-16-cursor-handover.md. Run the handover rehydration steps, reconcile them with the current git state, then continue only the OWED Next Immediate Steps.
```

## Mission control

| | |
|---|---|
| Repo | `/Users/z3437171/Dropbox/Github Local/pigauto` |
| main | `2ec4af3` · Version `0.10.0.9000` · Suggests `BACE` |
| Merged overnight | PRs #158–#165 (8) · CI green · `--as-cran` 0E/0W/1N |
| Docs branch | `handover/2026-08-09-cursor` (this handover + all dev-log) |
| In flight | Totoro regime-map campaign, 5400 jobs, resumable, ~2.5 h |
| Not launched | coverage campaign (approved, command in step 2) |
| Rose P1 | 8 of 9 addressed; **P1-8 open, needs Shinichi** |
| CRAN | pigauto 0.10.0 live · BACE **404** · do not submit |
| Do not | re-merge · CRAN-submit · relaunch the campaign from zero · start P1-8/Tier 2 unasked · edit BACE · `git add -A` · archive |

> Related: `docs/dev-log/2026-08-16-overnight-report.md` ·
> `docs/dev-log/2026-08-16-campaign-go-no-go.md` ·
> `docs/dev-log/2026-08-16-bench-rerun-results.md` ·
> `docs/dev-log/2026-08-15-publication-readiness-roadmap.md` ·
> `docs/dev-log/after-task/2026-08-15-p0-land.md` · `docs/dev-log/coordination-board.md`
