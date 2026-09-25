# Plan vs actual: posterior multiple imputation (mi-posterior), 2026-09-24

Plan: `~/.claude/plans/read-agents-md-and-docs-dev-log-handover-curious-eagle.md` (GOAL block, slices
S0-S7, compute estimates, gates G1-G9, pre-authorisation envelope, and the overnight runbook appended
2026-09-24 evening). Worktree `pigauto-mi-posterior`, branch `arc/mi-posterior` (stacked on merged
#187), draft PR itchyshin/pigauto#189. Reconciled at HEAD `f03837b` (46 commits ahead of
`origin/main`), one local commit ahead of the PR's last pushed commit `3cd3139`.

## 1. Slices: planned vs delivered

All 46 commits on this branch carry `Co-Authored-By: Claude Opus 5.5`; Claude Code does not split git
authorship by sub-agent, so "who did it" below is read from the plan's own task assignment and from
the session's active-agent roster, not from `git blame`.

| Slice | Plan (owner, model) | Delivered | Evidence |
|---|---|---|---|
| S0 Orchestrator | session, Opus 5.5 | Worktree, branch, cherry-pick from `arc/mi-gls-attenuation`, `design.md`, `GATES.md`, pre-dispatch review | `abada3c`, `f14ffa8`, `.unlazy/mi-posterior/GATES.md` |
| S1 Sampler | ceiling child, Opus, high | `R/mi_posterior.R`, `test-mi-posterior.R`, plus mid-build additions not in the original design (collapsed Metropolis moves, automatic chain extension) | `d176e98`, `05438ea`, `74215f3`, design.md 2.5 |
| S2 Wiring | Sonnet, medium | `R/multi_impute.R`, provenance markers, NEWS | `f4d7498`, `0f6f2cf` |
| S3 Simulation v2 | Sonnet, medium | `script/mi_gls/{regimes.R,01_cell_v2.R,04_acceptance.R}`, later extended to regimes 25-40 | `a3b1d10`, `7a0f470`, `a595fa8`, `2b95ea9` |
| S4 Real data | Sonnet, medium | `script/mi_realdata/{01_run.R,02_summarise.R,03_acceptance.R}` | `2d8c686`, `041c2e0` |
| S5 Review | Sonnet, high, fresh | `review-design.md` (M1, pre-dispatch, verdict REVISE) and `review.md` (S5, verdict PROCEED pending final gates); plus an unplanned round-2 harness review and a claims-audit round, recorded only through their fix commits | `review-design.md`, `review.md`, `9597e18`, `79bf98e`, `8d2f612` |
| S6 Results doc | Haiku, low | `results.md`, `results_tables.md`, rewritten three times as evidence changed (round-1 fail, diagnosis, round-2 pass), plus an M2 statistical correction | `e19e5b4`, `e0e8951`, `3cd3139`, `f03837b` |
| S7 Reconcile | Melissa, Sonnet, low | This document | here |

Delivered but not named as a slice in the plan: `diagnosis.md` (807e7ef, dd7d3aa), the dedicated
investigation of the round-1 G6 failures. This was implied by the plan's verification loop but not
itemised, and it drove the CP2 scope changes (section 3).

## 2. Gates: planned threshold vs final threshold vs outcome

| Gate | Planned | Final threshold | Outcome | Changed by / when |
|---|---|---|---|---|
| G1 unit tests | pass, no count given in the compact plan | >= 30 expectations, 0 failures | MET (118 + 196 expectations) | raised at design review (Gauss/Rose finding R1), resolved `dcc1174`/`a7a24f4` |
| G2 exactness | n=40, K=2, 20k draws vs dense conditional | split into G2a (Sigma_E > 0) / G2b (Sigma_E ~ 0) per review R2 | MET, `EXACTNESS_OK` | design review, resolved before S1 finished |
| G3 recovery | n=1000, 3 seeds, lambda within 0.1 of truth | calibration check: n=1000, 4 settings, >=50 seeds each, 95% interval coverage in [0.88,1.00], \|bias\|<=0.03, REML agreement >=90% | MET, `RECOVERY_OK` (200 fits) | Shinichi at CP1, 2026-09-24: the original 0.1 tolerance is ~1.7 posterior SD, so honest fits would fail it |
| G4 convergence | R-hat<1.05, ESS>400 on 2 named regimes (smoke) | unchanged | MET, `CONVERGENCE_OK` | not changed (review B3 flagged this as a gap; not resolved, see section 5) |
| G5 package | suite green, check clean, solver untouched | unchanged, split into G5a/b/c | MET at commit `9059f87`; `R/`/`tests/` changed since, re-verification at final HEAD not yet re-run through the ledger | none |
| G6 simulation | 24 regimes x 200 reps; \|bias\|<=max(0.02,2.5 MCSE); coverage >= complete-0.05; absolute SE ratio in [0.90,1.15]; proper>improper SE ratio gated | 40 regimes (8,000 cells), only regimes 17-40 gated (1-16 reported stress test); SE ratio judged RELATIVE to complete data; proper-vs-improper reported, not gated | NOT MET: round 1 failed on 8 checks; round 2 fails on 3 rows only (twins 35, 36, 38, relative SE ratio under phylolm) | scope narrowed and SE rule changed by Shinichi at CP1 (design.md 5d) and CP2 (5e); M2 review corrected the noise estimate and confirmed the 3 failures are not Monte Carlo noise (`f03837b`) |
| G7 per-cell coverage | [0.92,0.98] exchangeable masks, all regimes | gated only on regimes 17-40 MCAR rows; 1-16 reported | MET, `CELL_COVERAGE_PASS` | scope narrowed with G6 at CP2 |
| G8 real data | all planned cells, model vs conformal coverage, 5% slope check | 5% slope check reported not gated (unchanged, orchestrator decision R3); acceptance tightened to fail closed after a bug let it print a false pass | PENDING: FishBase cell still running; 9 of 10 done, 2 of those non-converged | fail-closed fix `041c2e0`, 2026-09-24 repair round |
| G9 / M1 manual | S5 verdict PROCEED, numbers trace to files | M1 = design reviewed before dispatch | M1 MET (verdict REVISE then resolved); S5 (G9's basis) verdict PROCEED, pending final gate re-verification | design review Gauss+Rose |
| M2 | fresh stats+code review, verdict PROCEED | unchanged | PENDING in `GATES.md` ("EVIDENCE: pending"), even though its substantive finding (correlation-aware SE-ratio MCSE) already landed in `results.md` at `f03837b` | in progress |

## 3. Compute: planned vs actual

Plan (revised after the smoke, before any array): simulation about 50-80 core-hours, 1-2 h wall as a
fir array (24 regimes x 200 reps); real data about 5 core-hours (FishBase dominates).

Actual, by campaign:
- First campaign (commit `69670d4`): 4,800 cells (24 regimes x 200 reps), about 4.45 h wall on 130
  cores on Totoro. Failed G6 on 8 checks (`f12ed09`).
- G3 calibration re-run at the same commit: 200 fits, 10 cores, 7,566 s wall (`evidence/README.md`).
- Diagnosis of the round-1 failures: about 55 minutes, up to 60 cores (briefly 64 for about 3
  minutes), on Totoro (`diagnosis.md`).
- Round-2 campaign (commit `9597e18`): 3,244 cells (44 non-converged re-runs plus the 3,200 new
  in-model twin cells), about 5.9 h wall on 140 cores, under heavy Totoro contention.
- Real data: 10 cells on 10 cores planned; 9 finished, the FishBase cell has been running over 15 h.

Two incidents, both on Totoro:
- The first launch of the G3 calibration run exported no thread caps; load average reached about
  18,000 for about 10 minutes before the run was killed and relaunched with `OMP_NUM_THREADS=1`,
  `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1` (`evidence/README.md`).
- The account's Totoro use exceeded the 150-core cap (D-143) because other lanes (drmTMB,
  `exact_prerun_cell.R`) launched uncapped jobs concurrently, contributing to the round-2 campaign's
  contention.

The wall-time and core-count figures for the two full campaigns are not in a committed log file in
this worktree (Totoro's own `driver.log` is not copied into the repo); they are reported here from the
overnight runbook, not independently re-derived. The committed figures (G3 calibration, diagnosis) do
match. Total simulation compute is at minimum an order of magnitude above the planned 50-80
core-hours; the 1-2 h wall estimate held only for the first campaign's simulation stage in isolation,
not for the arc as a whole once diagnosis, round 2, and real data are counted.

## 4. Scope changes (all with Shinichi's approval)

1. **Collapsed Metropolis moves** added during the S1 build (design.md 2.5.1), because plain Gibbs
   mixed too slowly near the lambda boundary (bulk ESS 20-41 at n=1000 in early testing).
2. **G3 redefined as a calibration check** rather than a fixed per-fit tolerance (CP1, design.md 5d.2).
3. **G6 SE-ratio rule changed from absolute to relative** to the complete-data ratio under the same
   analysis model (CP1, design.md 5d.3).
4. **G6 proper-vs-improper SE ratio downgraded from gated to reported** after the direction reversed
   at campaign scale relative to the S1 fixture measurement (CP1, design.md 5d.4).
5. **In-model twin regimes 25-40 added**, doubling the simulation grid from 24 to 40 regimes (4,800 to
   8,000 cells), because regimes 1-16 (raw covariance of non-ultrametric trees) sit outside the
   sampler's model family (CP2, design.md 5e.1, driven by `diagnosis.md` finding H4).
6. **Automatic chain extension** added for non-converged fits, up to 4x the default chain length (CP2,
   design.md 5e.2), because round-1 non-convergence was shown to be an ESS shortfall only, not a
   sampler failure.
7. **G8 tightened to fail closed** after `03_acceptance.R` exited 0 and printed its pass token inside a
   header line even on failure (`041c2e0`, 2026-09-24 repair round).

None of these touch `R/joint_mvn_solver.R`, which stays untouched per G5b (`SOLVER_UNTOUCHED`).

## 5. Child agents and workflows vs the plan's cap

Plan: "<= 6 new children, <= 1 ceiling child" (S1-S7, S1 the ceiling). This session's own active-agent
roster (addressable via SendMessage) lists, beyond `main` and this task (`melissa-reconcile`): `s1-recon`,
`s1-sampler`, `s2-instrumentation`, `s3-sim-harness`, `s4-realdata-harness`, `mi-gls-builder`,
`gauss-stats-review`, `rose-rereview`, `m1-design-review`, `m2-audit`, `emmy-code-review`, `diag-bias`,
`s6-paper`, `s6b-mi-draws`, `s9-claim-gate`, `traceability-check`, plus `planner-A` and `planner-B`
(the D-280 plan-selection pair, which predate this arc's G0 and arguably should not count against it).
Excluding the two planners, that is 16 non-orchestrator agents beyond `main` and this reconcile, well
above the plan's cap of 6 new children plus 1 ceiling.

The cap was exceeded. Per this task's own brief, the reason is that Shinichi turned on ultracode
mid-arc; no committed file in this worktree records that decision or a revised cap, so this is
reported as told, not independently verified from the repo. One agent, `kohaku-r-install`, does not
map to any plan slice or to anything in the committed `docs/dev-log/mi-posterior/` record; its
connection to this arc, if any, is not evidenced here.

## 6. Still open

- **G6 relative SE-ratio rows** (twins 35, 36, 38 under phylolm): a decision is proposed in
  `results.md` (three options, none applied); the gate itself is unchanged pending Shinichi.
- **G8**: the FishBase real-data cell is still running (over 15 h); 2 of the 9 finished cells did not
  converge at default settings (PanTHERIA, seed 20260819, MCAR and structured masks).
- **M2**: its finding is already folded into `results.md` (`f03837b`), but `GATES.md` still shows
  `EVIDENCE: pending` for a formal PROCEED verdict.
- **Default decision**: the plan deferred changing `multi_impute()`'s default `draws_method` to
  Shinichi; nothing in this branch changes it, and no decision has been recorded either way.
- **PR #189 body is stale**: it was last updated at commit `3cd3139` (24-regime numbers, G6 "not met"
  on 8 checks); HEAD is one commit ahead (`f03837b`, round-2's 40-regime numbers and the M2 correction
  to 3 failing rows) and the body has not been refreshed.
- **After-task report**: `docs/dev-log/after-task/2026-09-24-mi-posterior.md` exists on disk but is an
  unfilled template; its auto-populated "Files Touched" section lists files from an unrelated task
  (brain-vault and cursor files), not this arc. It is also gitignored (`docs/` is ignored by default in
  this worktree; the mi-posterior docs above were force-added), so it is not committed.
- **Handover**: no `docs/dev-log/handover/2026-09-24-*mi-posterior*` file exists yet.
- **Lease release**: not checked here; out of this reconcile's scope.
