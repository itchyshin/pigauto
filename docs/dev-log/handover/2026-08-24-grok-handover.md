# Session Handoff: pigauto evidence programme → Grok

Meta: 2026-08-24 · from Codex · target Grok (Cursor-compatible handover)

## Critical Context

This is a three-stage evidence programme, not a feature-development lane. Do
not change pigauto defaults, GNN architecture, EM/refinement, or make CRAN or
public scientific claims. Stage A is complete only at a narrow scope: the
continuous, BM, lambda = 1, one-step recovery result. It does **not** support
a cross-type active-imputation claim.

Stage B's exact conditional solver is merged to `main`; the one-mask timing
smoke succeeded for default and opt-in exact pigauto. BACE was absent from
Totoro's isolated library, but the project record identifies its upstream
source as `daniel1noble/BACE`. The prior claim that no source was found was
too narrow and has been corrected locally, but the correction is not yet
committed or pushed.

## Goals / mission

Earn only the next claims that the retained evidence supports:

- Stage A: active sampling, measured as realised one-step recovery;
- Stage B: an honest post-merge competitiveness boundary for exact
  conditional prediction;
- Stage C: observed-cell real-data confirmation of Mondrian calibration.

Every campaign retains failures, reports family/regime separately, and needs
an independent claim review. No default changes follow automatically.

## Plans / roadmap

1. Land the carried-over BACE-source correction below.
2. Install and verify BACE in Totoro's *private* library, then obtain a
   bounded BACE-inclusive one-mask timing receipt.
3. Re-estimate full five-mask Stage B cost including BACE and request explicit
   approval before starting it (the base, non-BACE work was about 47 minutes).
4. Run the retained five-mask study only after that gate, then obtain an
   independent claim review.
5. Stage C remains a distinct later approval: the valid Tamia four-GPU ladder
   is operational evidence, not a calibration result.

## What Was Accomplished

- Merged `origin/main` into `codex/active-recovery-evidence` without rewriting
  branch history, bringing in merged PR #173's exact route.
- Verified `tests/testthat/test-exact-conditional.R` locally.
- Refreshed `/home/snakagaw/pigauto-stageb-timing-69f0fa5` on Totoro,
  reinstalled pigauto into its private `Rlib`, and checked that `impute()` has
  `predict_method`.
- Completed retained one-mask timing smoke `20260901` at
  `/home/snakagaw/pigauto-stageb-timing-69f0fa5/results/timing-smoke-20260901-exact-main-retry2`.
  Both continuous and mixed receipts exist. Default/exact wall times were:

| Regime | Default | Exact | Total without BACE |
|---|---:|---:|---:|
| Continuous core | 109.976 s | 172.388 s | 283.719 s |
| Mixed core | 110.074 s | 170.115 s | 281.163 s |

- Retained BACE as `BACE not installed`; `Rphylopars`, `phylolm`, and
  `missForest` completed. Do not use one-mask scores as performance evidence.
- Corrected Tamia operational ledger: job `424950` completed a valid four-H100
  bounded ladder in 24 seconds with four distinct physical GPU UUID/PCI pairs.
  This is not a CPU/GPU comparison or calibration evidence.
- Written and pushed Stage-B timing receipt and after-task report; the monitor
  was deleted because the smoke completed.

## Current Working State

- **Working:** nothing running on Totoro or Tamia from this lane.
- **In progress:** BACE upstream-source correction is present locally but
  uncommitted because the current sandbox cannot create Git index locks.
- **Blocked:** this session cannot access GitHub or Totoro: its existing SSH
  control socket returns `Operation not permitted`. A new Grok session must
  verify its own access rather than assuming this restriction applies there.
- **Scientific state:** no Stage-B competitiveness claim and no Stage-C
  calibration claim are earned.

## Key Decisions & Rationale

- Preserve BACE in the Stage-B protocol. It is an explicit comparator; do not
  silently drop it merely because the first private library lacked it.
- Use `daniel1noble/BACE`, identified by the routed project record, rather than
  the failed overly narrow `itchyshin` repository search.
- Use per-mask retained RDS receipts and separate run directories; never
  overwrite prior failed/smoke receipts.
- Keep `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` on Totoro. The full
  Stage-B campaign must receive a fresh estimate and approval after BACE is
  timed.

## Landing State

`handoff_gate.sh` was run and failed only because of the carried-over files
below plus unrelated unpushed branches elsewhere in the repository. Do not
touch those other branches.

| Artifact / branch | Committed | Pushed | PR | State |
|---|---|---|---|---|
| `codex/active-recovery-evidence` at `3f3ed40` | yes | yes | none | LANDED |
| `docs/dev-log/2026-08-23-stage-b-timing-smoke.md` | no | no | none | CARRIED-OVER |
| `docs/dev-log/after-task/2026-08-23-stage-b-timing-smoke.md` | no | no | none | CARRIED-OVER |
| `docs/dev-log/plan-actual/2026-08-18-pigauto-evidence.md` | no | no | none | CARRIED-OVER |
| this handover | no | no | none | CARRIED-OVER |

**Why carried over:** the current Codex sandbox cannot create
`.git/worktrees/active-recovery/index.lock` and cannot access the existing
Totoro control socket. The three carried-over edits correct the BACE-source
record from “no source” to the identified `daniel1noble/BACE` upstream. They
pass `python3 ~/shinichi-brain/tools/closeout.py check` for the existing
after-task report and `git diff --check`.

**Exact resume command after inspecting the working tree:**

```sh
git add -f docs/dev-log/2026-08-23-stage-b-timing-smoke.md \
  docs/dev-log/after-task/2026-08-23-stage-b-timing-smoke.md \
  docs/dev-log/plan-actual/2026-08-18-pigauto-evidence.md \
  docs/dev-log/handover/2026-08-24-grok-handover.md
git commit -m 'docs: identify BACE upstream source and hand over Stage B'
git push origin codex/active-recovery-evidence
```

## Files Created / Modified

Committed and pushed in this lane:

- `R/exact_conditional.R` and the merged exact-route plumbing from `main`
- `script/bench_external_continuous_core.R`
- `script/bench_external_mixed_core.R`
- `script/mondrian_confirmation/04_run_tamia_gpu_smoke.sbatch`
- `script/mondrian_confirmation/05_run_tamia_ladder_task.sh`
- `tests/testthat/test-exact-conditional.R`
- `docs/dev-log/2026-08-19-mondrian-scalability-correction.md`
- `docs/dev-log/2026-08-23-stage-b-timing-smoke.md`
- `docs/dev-log/after-task/2026-08-23-stage-b-timing-smoke.md`
- `docs/dev-log/plan-actual/2026-08-18-pigauto-evidence.md`

The last three paths plus this handover have further uncommitted BACE-source
corrections described in the Landing State ledger.

## Next Immediate Steps

Classify each item as `OWED`, `DONE`, `RETRACTED`, or `PROTECTED` after
rehydration; execute only `OWED`.

1. Run lane preflight, `git status --short --branch`, and compare this file
   with the current branch. Do not stage or alter the many unrelated unpushed
   branches reported by the gate.
2. Inspect the four carried-over files above. If their correction is still
   accurate, commit and push them on `codex/active-recovery-evidence` using
   explicit paths only.
3. Verify existing Totoro ControlMaster access first. Do not cause a Duo
   prompt or open a fresh unauthenticated login.
4. In the isolated Totoro checkout, install BACE only into its private library
   from `daniel1noble/BACE` and verify `library(BACE)` plus `MCMCglmm` load.
   Preserve `R_LIBS_USER=/home/snakagaw/pigauto-stageb-timing-69f0fa5/Rlib`.
5. Run a fresh, separately named BACE-inclusive one-mask timing receipt with
   the existing runners. Retain every error. Do not overwrite
   `timing-smoke-20260901-exact-main-retry2`.
6. Present the measured BACE-inclusive five-mask estimate and request explicit
   approval before launching the full campaign. A full run is expected to
   exceed 30 minutes.

## Blockers / Open Questions

- Does the receiving environment have working access to Totoro's existing
  `cm-snakagaw@totoro.biology.ualberta.ca:22` ControlMaster socket? If absent,
  ask Shinichi once to authenticate; never trigger Duo.
- BACE installation is still unverified. The upstream source is known, but its
  dependency/runtime and per-mask timing are unknown.
- Full Stage B requires a new explicit approval after the BACE timing receipt.
- Stage C requires separately scoped full-calibration approval; do not conflate
  it with Stage B.

## Gotchas & Failed Approaches

- A timing clone from `69f0fa5` lacked PR #173 and failed with
  `unused argument (predict_method = ...)`. Merge/refresh/reinstall before
  timing exact conditional code.
- First launch omitted `PIGAUTO_OUT_DIR`; second source-loaded via `devtools`,
  which was absent in the private library. Use the installed private package
  after checking `impute()` formals.
- Do not infer physical GPU distribution from task-local
  `CUDA_VISIBLE_DEVICES=0`. Validate UUID/PCI values. Jobs `419946` and
  `419948` were not valid four-GPU scaling; job `424950` is.
- Do not say BACE has no source. The `itchyshin` search was too narrow;
  project memory identifies `daniel1noble/BACE`.
- Do not use results from the one-mask Stage-B smoke as scientific
  competitiveness evidence.

## Mission-Control Summary

| Repo | Branch / main | CI | What shipped | Plan by leverage |
|---|---|---|---|---|
| pigauto | `codex/active-recovery-evidence` / `origin/main` | no PR open; check live CI before push | Stage-A campaign, Stage-B runners/timing receipt, Tamia ladder correction | land BACE correction → BACE timing → approval → Stage-B five masks → claim review |

## How to Resume

Working directory:

```sh
cd '/Users/z3437171/Dropbox/Github Local/pigauto/.worktrees/active-recovery'
```

Read `AGENTS.md`, this handover, the Stage-B timing receipt, and the
plan-versus-actual ledger before acting. Then run:

```sh
python3 ~/shinichi-brain/tools/route.py pigauto
bash ~/shinichi-brain/tools/lane_preflight.sh '/Users/z3437171/Dropbox/Github Local/pigauto' --file script/bench_external_continuous_core.R
git status --short --branch
```

Grok receives no inherited terminal, GitHub, or SSH credentials. It should
run the live R/remote work only if those capabilities are available; otherwise
it should preserve the carried-over state and report the access boundary.

Paste this to Grok:

```text
Read AGENTS.md and docs/dev-log/handover/2026-08-24-grok-handover.md. Run the handover rehydration steps, reconcile them with the current git state, then continue only the OWED Next Immediate Steps.
```
