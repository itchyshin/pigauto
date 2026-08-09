# After-task — pigauto BACE wrap closeout

Date: 2026-08-09
Lane: wrap closeout via `/goal`
Restore branch: `feat/bace-wrap-restore` @ `b180555`
Handover checkout: `handover/2026-08-09-cursor` (docs + LOOP only)
Draft PR: https://github.com/itchyshin/pigauto/pull/156

## 1. Goal

Close out the GitHub-dev BACE wrap restore: re-verify focused wrap tests, write
after-task + Melissa LIGHT reconcile, open a **draft DO-NOT-MERGE** PR into
`main`, update the coordination board, then **stop until CRAN**. Do not merge.
Do not start DRM.jl, DESCRIPTION claim-gate, EM, Slack, vendor-sync, or a
public pigauto-vs-BACE sentence.

## 2. Implemented

Closeout-only. Restore code was already on origin; this arc did not re-restore
and did not edit `R/`.

- Re-ran focused wrap tests in `/tmp/pigauto-bace-wrap-restore` @ `b180555`.
  Fresh log: `[ FAIL 0 | WARN 0 | SKIP 1 | PASS 18 ]`.
- Opened draft PR [#156](https://github.com/itchyshin/pigauto/pull/156)
  (`feat/bace-wrap-restore` → `main`). Title and body say **DO NOT MERGE** until
  the v0.10.0 CRAN cut ships or BACE is on CRAN. Confirmed `isDraft: true`,
  `baseRefName: main`. Not merged.
- Wrote this after-task and Melissa LIGHT reconcile
  `docs/dev-log/plan-actual/2026-08-09-bace-wrap-closeout-reconcile.md`.
- Updated `LOOP/arcs.md`, `LOOP/checkpoint.md`, `LOOP/GOAL.md` (`n_final = 15L`
  is Shinichi's AskQuestion pick), and `docs/dev-log/coordination-board.md`.

Wrap API on the restore branch (unchanged this arc): `final_imp = FALSE`
default, `n_final = 15L`. NEWS already states this is not an imputed-as-observed
fix. `BACE` is in Suggests on this branch only.

## 3a. Decisions and Rejected Alternatives

- **Draft PR now, still do not merge.** Restore checkpoint said "no PR this
  turn." Closeout locked a draft DO-NOT-MERGE PR so the wait-for-CRAN gate is
  visible on GitHub. Rejected: mergeable ready-for-review PR (too easy to land
  by accident). Rejected: no PR at all (branch-only is easier to lose).
- **Lane-truth on handover; wrap code stays on restore.** LOOP / after-task /
  reconcile / board commit on `handover/2026-08-09-cursor`. Restore branch does
  not receive uinit or gnn artefacts. Rejected: committing closeout docs onto
  `feat/bace-wrap-restore` (would mix process files into the CRAN-wait code PR).
- **Do not re-restore.** Restore already reported green at `b180555`. Rejected:
  cherry-picking wrap files again.
- **D-43 honesty.** Two parked NOT-DONE rows remain (M1 full check still 2 WARN
  / 3 NOTE; wrap→main deferred). Do not claim the whole ultra-plan "done."
  Closeout of the GitHub-dev wrap lane is done; CRAN landing is not.
- **Slack / public BACE claim / vendor-sync / EM / DESCRIPTION 95% claim-gate /
  DRM.jl.** Left parked per the locked `/goal` paste. Rejected: "while we're
  here" starts.

## 4. Files Touched

This closeout arc (handover checkout):

- `LOOP/GOAL.md` (modified — `n_final` provenance + Slack parked wording)
- `LOOP/arcs.md` (modified — RESTORE + CLOSEOUT rows; landing text)
- `LOOP/checkpoint.md` (modified — closeout state + RESUME)
- `docs/dev-log/after-task/2026-08-09-bace-wrap-closeout.md` (new; `git add -f` because `.gitignore` has `docs/`)
- `docs/dev-log/plan-actual/2026-08-09-bace-wrap-closeout-reconcile.md` (new; same `-f`)
- `docs/dev-log/plan-actual/2026-08-09-bace-wrap-reconcile.md` (force-add; X1 was on disk, untracked)
- `docs/dev-log/2026-08-09-bace-wrap-brain-write-proposal.md` (force-add; X2 was on disk, untracked)
- `docs/dev-log/coordination-board.md` (modified — PR URL + wait-for-CRAN)

Restore branch (read + test only; no new commit):

- `/tmp/pigauto-bace-wrap-restore` @ `b180555` — `R/fit_baseline_bace.R`,
  tests, NEWS, DESCRIPTION Suggests `BACE`, NAMESPACE, `_pkgdown.yml`, T4 smoke,
  OVR footnote. Already committed earlier.

Left unstaged on handover (forbidden to stage):

- `.gitignore` `AGENTS.md` `CLAUDE.md` `README.md` `_pkgdown.yml`
- `dev/gnn_attribution_*` `script/*gnn_attribution*` `script/returned_gnn_attr/`

Not touched: either BACE tree, `main`, DRM.jl, EM, Slack, vault writes.

## 5. Checks Run

Focused wrap tests, restore worktree, 2026-08-09 this closeout (not memory):

```
cd /tmp/pigauto-bace-wrap-restore
NOT_CRAN=true Rscript -e 'devtools::load_all("."); testthat::test_file("tests/testthat/test-fit-baseline-bace-final-imp.R")'

ℹ Loading pigauto
SKIP: 'test-fit-baseline-bace-final-imp.R:127:3' ----------
Reason: installed BACE does export bace_final_imp(); nothing to assert

[ FAIL 0 | WARN 0 | SKIP 1 | PASS 18 ]
```

Worktree `git status -sb`: `## feat/bace-wrap-restore...origin/feat/bace-wrap-restore`
(clean). HEAD `b180555`. Standalone BACE still `ce8bc87d5b1dd4059bacfba23ab42c3f4dfe6080`.

PR: `gh pr view 156` → draft, base `main`, head `feat/bace-wrap-restore`,
title starts with `DO NOT MERGE`, state OPEN, not merged.

`python3 ~/shinichi-brain/tools/route.py pigauto` LOAD-FIRST loaded (Totoro/DRAC
compute default; recovery-to-truth; diff main before building; `r_cal = 0`
fallback; uncertainty first-class).

After-task structural check: `Rscript ~/shinichi-brain/tools/check-after-task.R`
on this file (run after write).

Not re-run this arc: full `devtools::check()`, full testthat suite, S3 re-bench.
Prior M1c receipt remains `/tmp/pigauto_m1c3/pigauto.Rcheck/00check.log`
(2 WARN / 3 NOTE, pre-existing).

## 6. Tests of the Tests

This closeout did not mutate production code, so there is no new red/green
cycle. The focused file's load-bearing guards, re-executed this turn:

- Formals test asserts `final_imp` defaults to `FALSE` and `n_final` to `15L`.
  A silent default flip would fail that test (it passed).
- Default-path same-seed reproducibility would fail if the `final_imp` branch
  leaked into `final_imp = FALSE` arithmetic (it passed, or would skip only if
  the BACE chain itself errored — it did not skip).
- The single SKIP is the negative branch: "error when BACE lacks
  `bace_final_imp`" cannot fire when installed BACE `@ce8bc87` does export it.
  That skip is expected, not a hole in this environment.

What the focused file still cannot catch: CRAN Suggests failure if this branch
were merged to `main`; a public pigauto-vs-BACE overclaim; full-suite
regressions outside `test-fit-baseline-bace-final-imp.R`.

## 7a. Issue Ledger

| ID | Issue | Status |
|---|---|---|
| M1c | Full `--as-cran` check 2 WARN / 3 NOTE (bench rds in `data/`, undeclared `jsonlite`, `.uinit`, top-level `results.tsv`/`run.log`, BACE not in mainstream repos) | Parked — pre-existing; wrap Rd `\usage` OK |
| LAND | wrap→main | Parked — wait for v0.10.0 CRAN cut / BACE on CRAN. Draft PR #156 is the pointer, not a merge offer |
| S3H | S3 harness lives in gitignored `docs/dev-log/`, not `script/` | Parked — needs a one-line scope OK to promote |
| X2S | Brain-write proposal still says `n_final = 15` is Ada's default; later AskQuestion made it Shinichi's pick | Residual in `docs/dev-log/2026-08-09-bace-wrap-brain-write-proposal.md`; vault still unwritten (D-37) |
| SLACK | Slack to Dan | Parked (`/goal`: "No Slack to Dan") |
| PUB | Public pigauto-vs-BACE sentence | Parked — S3 is wrap-config vs wrap-config only |

No new wrap bugs found this closeout.

## 8. Consistency Audit

- NEWS on restore (`git show b180555:NEWS.md`): states this is unrelated to the
  imputed-as-observed defect; "No pigauto-versus-BACE performance claim is made
  here." Matches the GOAL invariant.
- DESCRIPTION Suggests includes `BACE` on restore; `origin/main` `bf46991` does
  not. That is why merge-now is a CRAN blocker.
- NAMESPACE exports `fit_baseline_bace`; `_pkgdown.yml` lists it; T4 smoke
  restored in `test-shipping-coverage.R`. Diff vs `origin/main`: 9 files,
  +656 / −1.
- Handover dirty tree still holds uinit + gnn artefacts unstaged. Restore
  worktree clean.
- Neither BACE clone modified. In-tree `pigauto/BACE` not vendor-synced.
- Prior X1 reconcile
  (`docs/dev-log/plan-actual/2026-08-09-bace-wrap-reconcile.md`) still accurate
  for S0–X2 on the handover history; this closeout adds restore + PR + stop.
- Frozen `LOOP/ultra-plan.md` still says Ada's `n_final` default in one G0
  table; `LOOP/GOAL.md` and checkpoint now record Shinichi's AskQuestion pick.
  Frozen plan is not rewritten (immutability).

## 9. What Did Not Go Smoothly

- Handover `LOOP/checkpoint.md` still named handover HEAD `a5976b0` after
  restore-landing commit `88357ef` already existed. Closeout corrects the
  pointer.
- X2 brain proposal is stale on `n_final` provenance. Not silently rewritten
  into the vault; flagged here and in Melissa.
- `gh pr list --head feat/bace-wrap-restore` was empty before create; no
  duplicate PR. Draft create succeeded on the first try.
- Focused tests finished in ~5 s (tiny n=20, short MCMC). That is the smoke,
  not S3. Do not read speed as "BACE is cheap."

## 10. Known Residuals

- Whole ultra-plan is **not** complete under D-43: M1 full check not 0/0/0;
  wrap not on `main`. Both are parked/adaptive, stated honestly.
- Full testthat suite not claimed.
- S3 numbers remain regime-fenced (simulated BM, n=100, 3 traits, 5 seeds;
  1/5 `bace_final_imp` singularity). Not a pigauto-vs-BACE result.
- Brain `DECISIONS.md` wrap entry is still a proposal only (D-37).
- Handover branch remains dirty with uinit + gnn; 43-behind `main` history
  stays unrebased.
- Draft PR #156 is mergeable in GitHub's graph sense. That is not permission
  to merge.

## 11. Team Learning

Memory receipt: loaded `~/shinichi-brain/AGENTS.md` pointer via
`.cursor/rules/shinichi-brain.mdc`; `protocols/after-task.md` (12-section
Rose closeout); `agents/melissa.md` LIGHT reconcile; Cursor `/goal` adapter;
pigauto `AGENTS.md` / `CLAUDE.md` wrap invariants (add-a-parameter default-off;
NEWS must not claim imputed-as-observed; `r_cal = 0` fallback analogue;
never `git add -A`; never archive this project). `route.py pigauto` LOAD-FIRST
shaped compute (n/a this arc), "diff main before building" (restore already
from `origin/main` `bf46991`), and uncertainty-first (`se` propriety is why
`final_imp` exists). Brain MCP `search_notes` (`search_all_projects: true`)
returned BACE source / old head-to-head design notes, not a wrap-closeout
decision — consistent with Phase 0.25: no prior vault decision.

Golden Set: `tools/memory_regression.py` not run. No known-mistake class from
the Golden Set was in scope (this arc is process closeout + draft PR, not a
new estimator).

Durable lesson: a Suggests restore of an off-CRAN package is a **landing
gate**, not a code-quality gate. Focused tests can be green and the branch
still must not enter `main` until the CRAN surface is past. Draft + DO NOT
MERGE in the title is the mechanical reminder; merge remains a human gate.

## 12. Cross-Product Coverage

Cross-cutting flag this lane added: `final_imp` / `n_final` on
`fit_baseline_bace()`.

**Covers:** wrap `mu`/`se` derivation from chain datasets vs
`bace_final_imp()$all_datasets`; formals defaults (`FALSE` / `15L`); focused
test file; T4 shipping smoke; NEWS wording; OVR docs footnote; Suggests `BACE`
on `feat/bace-wrap-restore` only; draft PR #156 wait-for-CRAN pointer; S3
internal wrap-config coverage numbers (simulated BM, n=100, 5 seeds).

**does NOT cover:** `fit_pigauto()` gate calibration or `r_cal`; conformal
intervals; `multi_impute()` / Rubin pooling; mixed types beyond the wrap's
continuous+binary smoke (no count / ordinal / categorical / proportion /
zi_count / multi_proportion wrap cells); multi-obs / `obs_refine`; Level-C
joint / threshold / OVR estimator behaviour; `BACE::bace()` head-to-head vs
pigauto GNN; public pigauto-vs-BACE capability sentence; CRAN `main` tarball
(Suggests BACE would block); DESCRIPTION joint-Σ / 95% claim-gate; EM
(`max_iter > 0`); vendor-sync of in-tree `pigauto/BACE`; Slack to Dan; DRM.jl;
full `devtools::check()` 0/0/0; full testthat suite.
