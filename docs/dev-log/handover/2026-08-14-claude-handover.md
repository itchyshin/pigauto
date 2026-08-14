# Session Handoff: pigauto wrap LANDED; P0 not on main (Claude)

Meta: 2026-08-14 ~09:20 MDT · from Cursor (`AUTHOR = cursor`) · `TARGET = claude`  
Authoritative copy: this file. Chat is disposable. You inherit **no** authoring chat.

You are **Claude**, picking up **pigauto only**. Wrap merge G0 is **DONE**. Do not re-merge #156. Do not CRAN-submit from this `main`. Do not start DRM.jl from this file. Do not merge P0 onto `main` unless Shinichi locks that as a **new G0 that day**.

## Critical Context

1. **Wrap is on `origin/main`.** Merge commit `416561b043574174fcf7bf914bfed8b21bfa2146` (PR [#156](https://github.com/itchyshin/pigauto/pull/156) MERGED 2026-08-13T23:59:24Z). Version **0.10.0.9000**. Suggests `BACE`. `R/fit_baseline_bace.R` present. This is GitHub-dev, **not** a CRAN tarball.
2. **CRAN pigauto 0.10.0** published 2026-07-30. **CRAN BACE** still 404 (verified 2026-08-14). Next CRAN cut must drop Suggests BACE or wait until BACE is on CRAN.
3. **P0 sibling is CARRIED-OVER.** PR [#155](https://github.com/itchyshin/pigauto/pull/155) MERGED only to parent `fix/ci-install-libtorch` (`21d2ea6`), **not** `main`. Honesty fixes (per-tip BM SE, covariate alignment, zi_count conformal, train/cal symmetry) live on that parent only. Measured 2026-08-14: parent is **46 commits behind** `origin/main` and **9 commits ahead** of the merge-base (`c174c0f`).
4. **Shinichi memo verdict (2026-08-14): agree-with-corrections.** Next recommended G0 **if he names it**: land P0 on main, then claim-gate. Until he locks that, **STOP**. Do not invent the merge.
5. **PROTECTED trees.** Standalone BACE `/Users/z3437171/Dropbox/Github Local/BACE` @ `ce8bc87`. In-tree `pigauto/BACE` @ `de87d8c` (local main is 86 behind its own origin — do not “catch up”). No vendor-sync. No EM restore (`max_iter>0`). Never `git add -A`. Never `ScheduleWakeup` / archive this project.

## Goals / mission

pigauto remains a mixed-type phylogenetic imputer (gated BM + GNN). Architectural advantages, in order: (a) unified mixed-type imputation; (b) calibrated conformal intervals; (c) multi-tree Rubin pooling; (d) multi-obs `obs_refine`; (e) `suggest_next_observation()`.

Nothing wrap-OWED. Option B-minus is already on main: `final_imp = FALSE` default, `n_final = 15L`.

This session’s mission for Claude is **rehydrate, classify, then STOP** unless Shinichi explicitly locks a new G0.

## Plans / roadmap (beyond immediate steps)

- **If and only if** Shinichi locks “land P0 on main”: that is a **new G0 that day**. Rebase or merge `fix/ci-install-libtorch` (or `fix/p0-review-blockers`) onto current `origin/main` with care. Parent was already ~43 behind before wrap; measured **46 behind / 9 unique** after #156.
- After P0 is on main (only if that G0 happens): Rose claim-gate on DESCRIPTION / README / pkgdown / vignettes. P1 queue still open (see Blockers). GitHub [#135](https://github.com/itchyshin/pigauto/issues/135) still OPEN.
- Next CRAN cut: drop Suggests BACE or wait until BACE is on CRAN. **Do not CRAN-submit from this `main`.**
- Restore EM (`max_iter>0`) — **forbidden** without a dedicated numerical G0.
- DRM.jl — **other repo**. Do not start from this file.
- Slack to Dan — parked. Do not draft or send.
- Vendor-sync of in-tree `pigauto/BACE` — forbidden.

## What Was Accomplished

Wrap merge is **DONE** (2026-08-13, prior Cursor). This 2026-08-14 file records that landed state for a **fresh Claude** and declares everything the handoff gate still flags.

- G0 locked and executed: **Option B-minus** on `origin/main` via #156. `final_imp = FALSE` default, `n_final = 15L`.
- `origin/main` `416561b` · Version `0.10.0.9000` · Suggests `BACE` · `R/fit_baseline_bace.R` present.
- Post-merge CI **green** on `416561b`: R-CMD-check [31755795655](https://github.com/itchyshin/pigauto/actions/runs/31755795655) SUCCESS; pkgdown [31755795738](https://github.com/itchyshin/pigauto/actions/runs/31755795738) SUCCESS.
- CRAN pigauto 0.10.0 live (Published 2026-07-30). CRAN BACE 404.
- P0 #155 MERGED 2026-08-09T14:33:05Z to `fix/ci-install-libtorch` only. Not this G0. Not on `main`.
- Scoped Aug 2026 review was lanes **B+D+C**, not whole-repo. Rose P1 queue still open. #135 still OPEN.
- Memo verdict 2026-08-14: **agree-with-corrections**. Recommended next G0 *if named*: land P0, then claim-gate. Not locked.
- After-task (wrap merge): `docs/dev-log/after-task/2026-08-13-bace-wrap-merge.md`
- Prior Cursor handover (superseded for OWED steps; keep as history): `docs/dev-log/handover/2026-08-13-cursor-handover.md`
- P0 handover (sibling, still live): `docs/dev-log/handover/2026-08-08-p0-rose.md`
- Rose close (P1 queue): `docs/dev-log/handover/2026-08-08-rose-close.md`

## Current Working State

- **Working:** `origin/main` `416561b` Version 0.10.0.9000 with wrap. Post-merge CI green. Docs/LOOP on `handover/2026-08-09-cursor` (this commit will sit on top of `cf21e9b`).
- **In progress:** nothing wrap-OWED. P0 parent → main is **not** in progress unless Shinichi locks it.
- **Not working / blocked:** CRAN-submit from this `main` (Suggests BACE). P0 honesty fixes off `main` until a separate G0. Dirty uinit + GNN artefacts on this checkout must stay unstaged. 44 unpushed commits on other local branches — do not push them.

## Key Decisions & Rationale

- **Wrap on main is GitHub-dev, not a CRAN tarball.** `b615579` deleted the wrap to ship v0.10.0. #156 restored it after that cut. Suggests BACE would fail CRAN incoming while BACE is 404.
- **Do not re-merge #156.** It is MERGED. Re-restore / undraft / `gh pr merge 156` are retracted.
- **Do not merge P0 in this handover.** Shinichi has not locked that G0. “Merge or clean up” meant: wrap is already merged; clean the process lane for Claude.
- **B-minus, not flip the default.** Default BACE path stays chain averages. Proper-MI only when `final_imp = TRUE`. `n_final = 15L` is Shinichi’s pick.
- **S3 is wrap-config vs wrap-config.** Coverage 0.672 → 0.940 (`final_imp`, `n_final=15`) on simulated BM n=100, 5 seeds. **Not** a pigauto-vs-BACE headline.
- **Aug 2026 review was B+D+C only.** Do not generalise Rose findings to “the whole package” or “the GNN is broken.”
- **Coordination board is pigauto-only.** Do not write DRM rows here.
- **`r_cal = 0` remains a valid fallback.** Do not remove gate safety.

## Landing State

`handoff_gate.sh` on `/Users/z3437171/Dropbox/Github Local/pigauto` **FAIL** (2026-08-14 ~09:16 MDT): 15 uncommitted on current branch + 44 UNPUSHED on **other** branches. Declared below. Gate does not merge and does not push. Nothing below is silently dropped.

Verified just before write:

- `lane_preflight.sh` → PLATFORM cursor, ON BRANCH `handover/2026-08-09-cursor`, no foreign claude/codex lane in last 12h (weak evidence, D-87). Coordination board is in HEAD but not `origin/main` — push this docs branch.
- `gh pr view 156` → MERGED, `mergedAt: 2026-08-13T23:59:24Z`, merge commit `416561b`.
- `gh pr view 155` → MERGED to `fix/ci-install-libtorch`, `mergedAt: 2026-08-09T14:33:05Z`.
- `gh issue view 135` → OPEN.
- `origin/main` `416561b` Version 0.10.0.9000, Suggests BACE, `R/fit_baseline_bace.R` present.
- CI on `416561b`: R-CMD-check 31755795655 SUCCESS; pkgdown 31755795738 SUCCESS.
- `git rev-list --left-right --count origin/fix/ci-install-libtorch...origin/main` → `9	46`.

| Artifact / branch | Committed | Pushed | PR | State |
|---|---|---|---|---|
| pigauto `origin/main` `416561b` Version 0.10.0.9000 | y | y | [#156](https://github.com/itchyshin/pigauto/pull/156) **MERGED** | **LANDED** wrap. Do not re-merge. Not a CRAN tarball. |
| pigauto `feat/bace-wrap-restore` `a54e6a4` (NEWS) / merge `416561b` | y | y | #156 MERGED | **LANDED**. Historical restore branch. |
| pigauto `handover/2026-08-09-cursor` `cf21e9b` + this handover commit | y (docs/LOOP) | y at `cf21e9b`; push this commit | none (do **not** open a merge PR into `main`) | **LANDED** process lane. Dirty tree below is **not** part of this commit. |
| pigauto `fix/p0-review-blockers` `5800b01` → parent `fix/ci-install-libtorch` `21d2ea6` | y | y | [#155](https://github.com/itchyshin/pigauto/pull/155) **MERGED** → parent, **not** `main` | **CARRIED-OVER sibling.** Why: Shinichi has not locked “land P0 on main.” Honesty fixes are here only. Resume: only after a new G0 that day — see command below. |
| Dirty on handover checkout: `.gitignore` `AGENTS.md` `CLAUDE.md` `README.md` `_pkgdown.yml` | n | n | — | **CARRIED-OVER** — uinit / banner noise. **Do not stage.** Why: not wrap; would entangle lanes. Resume: leave dirty. |
| Untracked `dev/gnn_attribution_*` + `script/*gnn_attribution*` + `script/returned_gnn_attr/` | n | n | — | **CARRIED-OVER** — Fir/Totoro GNN ladder artefacts. **Do not stage.** Why: other experiment. Resume: leave untracked. |
| Other local unpushed branches (gate: 44 commits — `analysis/calibration-grid`, `chore/dedupe-agent-files`, `codex/fix-winbuilder-pagel-skip`, `codex/per-trait-v010-pages`, `codex/pigauto-final-comments`, `codex/pigauto-submit-record`, `codex/pkgdown-dev-status`, `codex/v010-mi-calibration`, `codex/v010-winbuilder-libtorch-skip`, `codex/v0100-cran-comments`, `codex/winbuilder-libtorch`, `experiment/covariate-honest-sim`, `experiment/coverage-simulation`, `experiment/overnight-sim-2026-05-19`, `experiment/test-pagel-lambda-bien`, `feature/clade-mask`, `feature/curriculum`, `feature/swa`, `fix_multi_impute_baseline_fix`, `pr-116`, `shannon-install`, `spec/vulcan-gpu-avonet9993`) | mixed | **n** (44 commits) | not this lane | **CARRIED-OVER.** Why: foreign WIP. Do not push, rebase, or “clean up” from this lane. Resume: do nothing unless Shinichi names that branch. |
| standalone BACE `@ce8bc87` | n/a | n/a | — | **PROTECTED** |
| `pigauto/BACE` `@de87d8c` (local main 86 behind its origin) | n/a | n/a | — | **PROTECTED** — no vendor-sync; do not fetch/rebase it |
| DRM.jl (sibling repo) | n/a | n/a | own handover | **PROTECTED / other repo** — do not claim |

**Resume P0 → main (ONLY after Shinichi locks a new G0 that day):**

```bash
cd "/Users/z3437171/Dropbox/Github Local/pigauto"
git fetch origin
# Confirm the lock in chat that day. Do not run this from rehydration.
# Parent is 46 commits behind origin/main (measured 2026-08-14).
# Preferred carrier: origin/fix/ci-install-libtorch (has #155 merge 21d2ea6).
# Alternate: origin/fix/p0-review-blockers (5800b01) if you must start from the PR head.
# Rebase or merge onto current origin/main with care. Do not force-push other people's branches.
# Do not CRAN-submit. Do not restore EM. Do not vendor-sync BACE.
```

## Files Created / Modified

**This handover session (2026-08-14):**

- `docs/dev-log/handover/2026-08-14-claude-handover.md` (this file; `git add -f` because `.gitignore` has `docs/`)
- `docs/dev-log/coordination-board.md` (wrap row + Current Rule + 2026-08-14 status; P0 sibling retained)
- `docs/dev-log/phase-snapshot.md` + `docs/dev-log/phase-snapshot-archive.md` (archive previous one-liner; new pointer at the board)
- `docs/dev-log/after-task/2026-08-14-cursor-to-claude-handover.md` (short closeout receipt)
- `LOOP/checkpoint.md` (banner + STATE/RESUME refreshed so stale “do not merge #156” is not the live instruction)

**Wrap already on `origin/main` `416561b` (do not re-land):**

- `DESCRIPTION` `NAMESPACE` `NEWS.md`
- `R/fit_baseline_bace.R` `R/ovr_categorical.R`
- `_pkgdown.yml` `man/fit_baseline_bace.Rd`
- `tests/testthat/test-fit-baseline-bace-final-imp.R` `tests/testthat/test-shipping-coverage.R`

**Must not stage (still dirty / untracked on this checkout):** `.gitignore` `AGENTS.md` `CLAUDE.md` `README.md` `_pkgdown.yml` `dev/gnn_attribution_*` `script/*gnn_attribution*` `script/returned_gnn_attr/` · either BACE tree · other people’s unpushed branches.

## Next Immediate Steps (OWED only)

Classify every item `OWED` · `DONE` · `RETRACTED` · `PROTECTED` against **current** git/gh before acting. Execute only `OWED`. Re-merge #156, CRAN-submit, P0-merge-now, DRM.jl, vendor-sync, EM, Slack, and staging dirty uinit/GNN are **not** OWED.

1. **Rehydrate + classify.** `~/shinichi-brain/tools/lane_preflight.sh "/Users/z3437171/Dropbox/Github Local/pigauto"`. `git fetch origin`. `git status`. `gh pr view 156` (expect MERGED). `gh pr view 155` (expect MERGED to parent, not main). Read the Active Lane Split (**both** rows). Classify this file’s items against current git/gh.
2. **If wrap is still landed and Shinichi has not named a new G0: STOP.** Do not re-merge #156. Do not CRAN-submit. Do not start P0 merge from rehydration.
3. **If and only if** Shinichi explicitly locks “land P0 on main”: that is a **NEW G0 that day**. Rebase or merge `fix/ci-install-libtorch` (or `fix/p0-review-blockers`) onto current `origin/main` with care (parent is 46 commits behind after #156). Then claim-gate. Until that lock, P0 merge is not OWED.
4. **Do not start DRM.jl** from this pigauto handover. If Shinichi wants a new lane, he names the repo.

## Blockers / Open Questions

- **P0 parent not on `main`.** Honesty fixes exist only on `fix/ci-install-libtorch`. Merge is a Shinichi G0, not a rehydration action.
- **CRAN-submit blocked** while Suggests includes BACE and CRAN BACE is 404.
- **Rose P1 queue still open** (from `2026-08-08-rose-close.md`; scoped B+D+C, not whole-repo): P1-5 joint-Σ claim-gate · P1-6 conformal coverage wording · P1-7 `$se` is three objects · P1-8 covariates ignored by joint/threshold-joint · P1-9 zi_count zeros excluded from val/test · P1-10 fail-closed encode gaps · P1-11 tree-MI has no conformal draws · P1-12 multi-obs phylo-signal gate no-op · P1-13 type detection / `read_traits` count mislabelling.
- GitHub [#135](https://github.com/itchyshin/pigauto/issues/135) OPEN (GNN evidence layer / gate-firing diagnostics).
- Board + this handover live on `handover/2026-08-09-cursor`, not on `origin/main`. A lane that only reads `origin/main` cannot see them — use this branch or the force-added docs on it.
- M1c full `--as-cran` 2 WARN / 3 NOTE remains pre-existing (data/bench rds; jsonlite in tests; version + BACE not in mainstream repos). Do not claim the old ultra-plan wholly done (D-43).

## Gotchas & Failed Approaches

- Do **not** tell the next session to re-restore or merge #156. That G0 is finished.
- Do **not** merge P0 because a handover said “clean up.” Clean up = process lane. P0 merge needs a new lock.
- Do not rebase `handover/2026-08-09-cursor` onto `main`. Leave dirty uinit / GNN unstaged.
- `docs/` is `.gitignore`d. Force-add handover / board / after-task / snapshot paths only. Never `git add -A`.
- Slack 0.75–0.88 corr ≠ poisson Pearson. S3 numbers are wrap-vs-wrap coverage, not a pigauto-vs-BACE headline.
- `bace_imp` already reinstated NAs; Dan’s imputed-as-observed bug was `bace_final_imp`. Switching the default would not be “the same fix.”
- Never `ScheduleWakeup` / archive this project (needs the exact phrase `Archive this project now.`).
- In-tree `pigauto/BACE` local main being 86 behind is **not** a pigauto task. Do not sync it.
- Claude plans / refactors / prose / logic+CI checks. Live `R CMD check` with torch, Totoro/DRAC benches, and CRAN tarball work belong on a machine that can run them (Codex or this host with R+torch). Do not claim a public result from chat memory.

## How to Resume (Claude)

Working directory: `/Users/z3437171/Dropbox/Github Local/pigauto`  
Checkout: stay on `handover/2026-08-09-cursor` for docs, **or** read the force-added docs on that branch. Wrap **code** is `origin/main` `416561b`. Do not invent a third restore branch. Do not rebase this docs branch onto main.

Toolchain: R + `devtools::load_all()` + `torch`. Safe verify (only if you must re-check wrap, not as OWED work):

```bash
cd "/Users/z3437171/Dropbox/Github Local/pigauto"
git fetch origin
git show origin/main:DESCRIPTION | rg 'Version|BACE'
# expect Version 0.10.0.9000 and Suggests BACE
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-fit-baseline-bace-final-imp.R")'
```

Must not stage: `.gitignore` `AGENTS.md` `CLAUDE.md` `README.md` `_pkgdown.yml` · `dev/gnn_attribution_*` · `script/*gnn_attribution*` · `script/returned_gnn_attr/` · either BACE tree · other people’s unpushed branches.

Rehydrate (human starts a **fresh Claude** session in the repo and pastes):

```text
Read AGENTS.md and docs/dev-log/handover/2026-08-14-claude-handover.md. Run the handover rehydration steps, reconcile them with the current git state, then continue only the OWED Next Immediate Steps.
```

Interactive equivalent: `claude "Rehydrate from docs/dev-log/handover/2026-08-14-claude-handover.md + the AGENTS.md snapshot, then continue with the Next Immediate Steps."`

Read order: `AGENTS.md` → this file → `docs/dev-log/coordination-board.md` (**both** Active Lane Split rows) → P0 sibling `docs/dev-log/handover/2026-08-08-p0-rose.md` → `gh pr view 156` + `gh pr view 155` + `git log -1 origin/main`.

## Mission control

| | |
|---|---|
| Repo | `/Users/z3437171/Dropbox/Github Local/pigauto` |
| Wrap code | `origin/main` `416561b` · PR [#156](https://github.com/itchyshin/pigauto/pull/156) **MERGED** · Version `0.10.0.9000` · Suggests `BACE` |
| Docs / this handover | `handover/2026-08-09-cursor` · `docs/dev-log/handover/2026-08-14-claude-handover.md` |
| CI on wrap merge | R-CMD-check 31755795655 SUCCESS · pkgdown 31755795738 SUCCESS |
| CRAN | pigauto 0.10.0 published 2026-07-30 · BACE **404** |
| P0 sibling | #155 **MERGED** → `fix/ci-install-libtorch` `21d2ea6` · **not** `main` · 46 behind / 9 unique |
| What shipped | B-minus wrap on GitHub-dev main; post-merge CI green |
| Plan by leverage | **STOP** unless Shinichi locks “land P0 on main” as a new G0. Then rebase/merge parent with care, then claim-gate. |
| Do not | re-merge #156 · CRAN-submit · merge P0 without a lock · edit BACE · vendor-sync · EM · Slack · DRM.jl · `git add -A` · archive |

> Related: `docs/dev-log/coordination-board.md` · `docs/dev-log/after-task/2026-08-13-bace-wrap-merge.md` · `docs/dev-log/handover/2026-08-13-cursor-handover.md` · `docs/dev-log/handover/2026-08-08-p0-rose.md` · `docs/dev-log/handover/2026-08-08-rose-close.md` · `LOOP/checkpoint.md`
