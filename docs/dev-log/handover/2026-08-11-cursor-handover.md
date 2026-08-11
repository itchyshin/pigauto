# Session Handoff: pigauto wrap lane CLOSED pending CRAN (Cursor)

Meta: 2026-08-11 ~11:20 MDT · from Cursor (`AUTHOR = cursor`) · `TARGET = cursor`  
Authoritative copy: this file. Chat is disposable.

You are **Cursor**, picking up the **pigauto wrap lane**. You inherit no chat history. This is a dated state record, not an instruction to repeat restore or closeout.

**Remit is pigauto only.** DRM.jl has its own handover (`DRM.jl/docs/dev-log/handover/`). Do not start DRM from this file. If Shinichi wants a new lane, he names the repo.

## Critical Context

1. **Wrap lane is CLOSED pending CRAN.** Do **not** re-restore. Do **not** merge [PR #156](https://github.com/itchyshin/pigauto/pull/156) until the v0.10.0 CRAN cut has shipped **or** BACE is on CRAN **or** Shinichi explicitly says merge. Suggests BACE would block CRAN if merged now.
2. **Two pigauto lanes, one board.** Rehydrate via `docs/dev-log/coordination-board.md` **Active Lane Split**. A single `AGENTS.md` snapshot would orphan a sibling. P0/#155 is a sibling (merged to parent `fix/ci-install-libtorch`, **not** `main`). Never merge wrap or P0 to `main` until CRAN.
3. **PROTECTED trees.** Standalone BACE `/Users/z3437171/Dropbox/Github Local/BACE` @ `ce8bc87`. In-tree `pigauto/BACE` @ `de87d8c`. No vendor-sync. No EM restore (`max_iter>0`). Never `git add -A`. Never `ScheduleWakeup` / archive this project.

## Goals / mission

pigauto remains a mixed-type phylogenetic imputer (gated BM + GNN). The wrap lane’s locked goal was Option B-minus: keep the default `bace_imp` chain-average path bit-identical, and expose Dan’s Study-B object only via opt-in `final_imp = TRUE`. That implementation **already landed** on `feat/bace-wrap-restore`. This session’s mission for the next Cursor is **wait**, not rebuild.

## Plans / roadmap (beyond immediate steps)

- Merge #156 **only after** v0.10.0 CRAN cut ships or Shinichi names a new G0 that day. Until then wrap is not OWED.
- Getting `fix/ci-install-libtorch` (with merged #155) onto `main` — later G0, after CRAN. Never merge P0 to `main` until CRAN.
- DESCRIPTION 95% / joint-Σ claim-gate (P1) — **not this lane**.
- Restore EM (`max_iter>0`) — **forbidden** without a new G0.
- DRM.jl — **other repo**. Do not start from this handover.
- Slack to Dan — parked. Do not draft or send.

## What Was Accomplished

Wrap restore + closeout are **DONE** (2026-08-09). This 2026-08-11 file only records that closed state for a fresh Cursor.

- G0 locked: **Option B-minus**, `final_imp = FALSE` default, `n_final = 15L` (Shinichi AskQuestion pick). No Slack to Dan.
- Restore: `feat/bace-wrap-restore` @ `b180555` on origin, based on `origin/main` `bf46991` (still **Version 0.10.0**). Wrap branch Version `0.10.0.9000`.
- Draft PR: https://github.com/itchyshin/pigauto/pull/156 — `isDraft: true`, base `main`, **DO NOT MERGE**.
- S3: wrap-vs-wrap only. Simulated BM, n=100, 5 seeds: cover **0.672 → 0.940** (`final_imp`, `n_final=15`). **Not** pigauto-vs-BACE. No public comparison sentence.
- M1 full `--as-cran`: **2 WARN / 3 NOTE**, pre-existing (not wrap-caused). Focused tests: **FAIL 0 / SKIP 1 / PASS 18**.
- PR **#155 MERGED** (2026-08-09) into `fix/ci-install-libtorch`. Parent still has ubuntu-latest FAILURE history; do not treat that as a wrap task.
- After-task: `docs/dev-log/after-task/2026-08-09-bace-wrap-closeout.md`
- Melissa LIGHT: `docs/dev-log/plan-actual/2026-08-09-bace-wrap-closeout-reconcile.md`
- Plan: `docs/dev-log/handover/2026-08-09-bace-wrap-g0-ultra-plan.md`
- LOOP/: `GOAL.md`, `checkpoint.md`, `arcs.md`, `ultra-plan.md` on `handover/2026-08-09-cursor` (incl `0c99d37` closeout commit, pushed).

## Current Working State

- **Working:** Wrap code on `origin/feat/bace-wrap-restore` @ `b180555`. Docs/LOOP on `handover/2026-08-09-cursor` @ `0c99d37` (this handover commit will sit on top). `origin/main` @ `bf46991`, Version 0.10.0, cran-comments still a 0.10.0 resubmission. Standalone BACE @ `ce8bc87`.
- **In progress:** Nothing wrap-OWED. CRAN cut has **not** shipped (main still 0.10.0; #156 still draft).
- **Not working / blocked:** Merging #156 would put BACE in Suggests and undo the v0.10.0 CRAN surface (`b615579` deleted the wrap for that reason). Full `--as-cran` still 2 WARN / 3 NOTE. Dirty uinit + GNN artefacts on the handover checkout must stay unstaged.

## Key Decisions & Rationale

- **B-minus, not flip the default.** Default path stays chain averages. Proper-MI only when `final_imp = TRUE`. `n_final = 15L` is Shinichi’s pick (matches Study B). Do not flip to 50 without a new ask.
- **GitHub-dev only until CRAN.** Restore lives on `feat/bace-wrap-restore`, not `main`. Draft PR is a visible wait-gate, not a merge queue.
- **S3 is wrap-config vs wrap-config.** Coverage lift is internal. A public pigauto-vs-BACE claim is a Shinichi gate and was not measured.
- **NEWS must not claim an imputed-as-observed fix.** That defect lived in BACE’s `bace_final_imp()`, which pigauto never called until the opt-in path.
- **No vendor-sync. No EM. No Slack.** In-tree `BACE/` is not pigauto work.
- **Coordination board is pigauto-only.** Do not write DRM rows here.

## Landing State

`handoff_gate.sh` on `/Users/z3437171/Dropbox/Github Local/pigauto` **FAIL** (2026-08-11): 15 uncommitted on current branch + 44 unpushed commits on **other** branches. Declared below. Gate does not merge and does not push. Nothing below is silently dropped.

Verified just before write: `lane_preflight.sh` → PLATFORM cursor, ON BRANCH `handover/2026-08-09-cursor`, PR #156 open, no foreign claude/codex lane in last 12h (weak evidence). `gh pr view 156` → OPEN draft → `main`. `gh pr view 155` → MERGED. `origin/main` `bf46991` Version 0.10.0.

| Artifact / branch | Committed | Pushed | PR | State |
|---|---|---|---|---|
| pigauto `feat/bace-wrap-restore` `b180555` | y | y | [#156](https://github.com/itchyshin/pigauto/pull/156) OPEN draft → `main` | **LANDED** as GitHub-dev restore. Merge to `main` is **not OWED**. |
| pigauto `handover/2026-08-09-cursor` `0c99d37` + this handover commit | y (docs/LOOP) | y at `0c99d37`; push this commit | none (do not open a merge PR into `main`) | **LANDED** process lane. Dirty tree below is **not** part of this commit. |
| pigauto `origin/main` `bf46991` Version 0.10.0 | y | y | — | **LANDED** CRAN surface. Wrap deleted here on purpose (`b615579`). |
| pigauto `fix/p0-review-blockers` `5800b01` | y | y | [#155](https://github.com/itchyshin/pigauto/pull/155) **MERGED** → `fix/ci-install-libtorch` | **CARRIED-OVER sibling** — merged to parent, **not** `main`. Do not merge parent to `main` until CRAN. Resume: none for wrap. |
| Dirty on handover checkout: `.gitignore` `AGENTS.md` `CLAUDE.md` `README.md` `_pkgdown.yml` | n | n | — | **CARRIED-OVER** — uinit / banner noise. **Do not stage.** Why: not wrap; would entangle lanes. Resume: leave dirty. |
| Untracked `dev/gnn_attribution_*` + `script/*gnn_attribution*` + `script/returned_gnn_attr/` | n | n | — | **CARRIED-OVER** — Fir/Totoro GNN ladder artefacts. **Do not stage.** Why: other experiment. Resume: leave untracked. |
| Other local unpushed branches (gate: analysis/calibration-grid, chore/dedupe-agent-files, codex/winbuilder+CRAN-comment branches, experiment/*, feature/clade-mask, feature/curriculum, feature/swa, spec/vulcan-gpu-avonet9993, shannon-install, …) | mixed | **n** (44 commits) | not wrap | **CARRIED-OVER** — out of wrap remit. Why: foreign WIP; do not push or rebase them from this lane. Resume: do nothing unless Shinichi names that branch. |
| standalone BACE `@ce8bc87` | n/a | n/a | — | **PROTECTED** |
| `pigauto/BACE` `@de87d8c` | n/a | n/a | — | **PROTECTED** — no vendor-sync |
| DRM.jl (sibling repo) | n/a | n/a | own handover | **PROTECTED / other repo** — do not claim |

**Resume #156 (only after CRAN cut or explicit Shinichi merge G0):**

```bash
cd "/Users/z3437171/Dropbox/Github Local/pigauto"
gh pr view 156 --repo itchyshin/pigauto --json state,isDraft,baseRefName,mergeable
# Confirm origin/main Version is no longer the pre-cut 0.10.0 wait, or Shinichi said merge.
# Then, and only then, a NEW G0 that day. Do not `gh pr merge 156` from rehydration.
```

## Files Created / Modified

**This handover session (2026-08-11):**

- `docs/dev-log/handover/2026-08-11-cursor-handover.md` (this file; `git add -f` because `.gitignore` has `docs/`)
- `docs/dev-log/coordination-board.md` (wrap row → CLOSED / wait CRAN / PR 156; P0 sibling retained)
- `docs/dev-log/phase-snapshot.md` (pointer at the board’s Active Lane Split — not a single-lane orphan)

**Wrap restore already on `origin/feat/bace-wrap-restore` (`git diff --name-only origin/main...origin/feat/bace-wrap-restore`):**

- `DESCRIPTION` `NAMESPACE` `NEWS.md`
- `R/fit_baseline_bace.R` `R/ovr_categorical.R`
- `_pkgdown.yml` `man/fit_baseline_bace.Rd`
- `tests/testthat/test-fit-baseline-bace-final-imp.R` `tests/testthat/test-shipping-coverage.R`

**Closeout already on `handover/2026-08-09-cursor` @ `0c99d37`:**

- `LOOP/GOAL.md` `LOOP/arcs.md` `LOOP/checkpoint.md` `LOOP/ultra-plan.md`
- `docs/dev-log/after-task/2026-08-09-bace-wrap-closeout.md`
- `docs/dev-log/plan-actual/2026-08-09-bace-wrap-closeout-reconcile.md`
- `docs/dev-log/plan-actual/2026-08-09-bace-wrap-reconcile.md`
- `docs/dev-log/2026-08-09-bace-wrap-brain-write-proposal.md`
- `docs/dev-log/handover/2026-08-09-bace-wrap-g0-ultra-plan.md`
- prior `docs/dev-log/handover/2026-08-09-cursor-handover.md` (superseded for wrap OWED steps; keep as history)

**Must not stage (still dirty / untracked on this checkout):** `.gitignore` `AGENTS.md` `CLAUDE.md` `README.md` `_pkgdown.yml` `dev/gnn_attribution_*` `script/*gnn_attribution*` `script/returned_gnn_attr/`

## Next Immediate Steps (OWED only)

Classify every item `OWED` · `DONE` · `RETRACTED` · `PROTECTED` against **current** git/gh before acting. Execute only `OWED`. Restore, closeout, re-restore, and merge-#156-now are **not** OWED.

1. **Rehydrate.** `~/shinichi-brain/tools/lane_preflight.sh "/Users/z3437171/Dropbox/Github Local/pigauto"`. `git status`. `gh pr view 156 --repo itchyshin/pigauto`. Read the Active Lane Split. Do not redo wrap restore or closeout.
2. **If `origin/main` is still Version 0.10.0 and #156 is still draft: STOP wrap.** Do not merge. Optional only: check whether the CRAN cut shipped (`DESCRIPTION` Version / `cran-comments.md` / new `origin/main` commits).
3. **If and only if** the v0.10.0 CRAN cut has shipped **or** Shinichi explicitly says merge: then merge #156 as a **new G0 that day**. Until then wrap is not OWED.
4. **Do not start DRM.jl** from this pigauto handover. If Shinichi wants a new lane, he names the repo.

## Blockers / Open Questions

- **CRAN cut not shipped** (verified 2026-08-11: `origin/main` still 0.10.0). This is the merge blocker.
- M1c full check 2 WARN / 3 NOTE remains parked (pre-existing). Do not claim the ultra-plan “done” (D-43: two NOT-DONE rows).
- Board lives on `handover/2026-08-09-cursor`, not on `origin/main`. A lane that only reads `origin/main` cannot see it — use this branch or the force-added docs on it.

## Gotchas & Failed Approaches

- Do **not** tell the next session to re-restore or merge #156. That was the 2026-08-09 closeout; it is finished.
- Do not rebase `handover/2026-08-09-cursor` onto `main`. Leave dirty uinit / GNN unstaged.
- `docs/` is `.gitignore`d. Force-add handover paths only. Never `git add -A`.
- Slack 0.75–0.88 corr ≠ poisson Pearson. S3 numbers are wrap-vs-wrap coverage, not a pigauto-vs-BACE headline.
- `bace_imp` already reinstated NAs; Dan’s imputed-as-observed bug was `bace_final_imp`. Switching the default would not be “the same fix.”
- Never merge wrap or P0 to `main` until CRAN. #155 merged to **parent** only.
- Never `ScheduleWakeup` / archive this project (needs the exact phrase `Archive this project now.`).

## How to Resume (Cursor)

Working directory: `/Users/z3437171/Dropbox/Github Local/pigauto`  
Checkout: stay on `handover/2026-08-09-cursor` for docs; wrap **code** is the other worktree `/tmp/pigauto-bace-wrap-restore` @ `b180555` if you must inspect it. Do not invent a third restore branch.  
Toolchain: R + `devtools::load_all()`. Wrap verify (only if you must re-check, not as OWED work):

```bash
cd /tmp/pigauto-bace-wrap-restore   # or git fetch && checkout feat/bace-wrap-restore in a clean worktree
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-fit-baseline-bace-final-imp.R")'
```

Must not stage: uinit dirty files, GNN `dev/`/`script/` artefacts, either BACE tree, other people’s unpushed branches.

Rehydrate (paste into a **fresh** Cursor agent):

```text
Read AGENTS.md and docs/dev-log/handover/2026-08-11-cursor-handover.md. Run the handover rehydration steps, reconcile them with the current git state, then continue only the OWED Next Immediate Steps.
```

Read order: `AGENTS.md` → this file → `docs/dev-log/coordination-board.md` (both rows) → after-task + Melissa only if you need receipts → `gh pr view 156`.

## Mission control

| | |
|---|---|
| Repo | `/Users/z3437171/Dropbox/Github Local/pigauto` |
| Wrap code | `feat/bace-wrap-restore` @ `b180555` · draft [#156](https://github.com/itchyshin/pigauto/pull/156) · **DO NOT MERGE** |
| Docs / this handover | `handover/2026-08-09-cursor` |
| `origin/main` | `bf46991` · Version **0.10.0** · CRAN cut **not** shipped |
| P0 sibling | #155 **MERGED** → `fix/ci-install-libtorch` · not `main` |
| What shipped | B-minus wrap on GitHub-dev; S3 0.672→0.940 wrap-vs-wrap; focused 18/1/0 |
| Plan by leverage | **STOP wrap.** Optional CRAN-cut check. Merge #156 only after cut or explicit G0. |
| Do not | re-restore · merge #156 now · edit BACE · vendor-sync · EM · Slack · DRM.jl · `git add -A` · archive |

> Related: `docs/dev-log/coordination-board.md` · `docs/dev-log/after-task/2026-08-09-bace-wrap-closeout.md` · `docs/dev-log/plan-actual/2026-08-09-bace-wrap-closeout-reconcile.md` · `LOOP/checkpoint.md`
