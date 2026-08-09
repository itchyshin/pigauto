# Ultra-plan — P0 review blockers (Phases 1–2)

Status: **awaiting G0** — do not implement until Shinichi approves

```
🎯 GOAL
Solo platform: Cursor
Deliverable: 4/4 P0 defects fixed + tests on current dirty tree (B1–B3 kept)
HEADLINE: silent wrong SE / cov mispair / zi_count MI / stale joint test
IN PARALLEL: S1 cov-align · S2 SE formula · S3 zi_count conformal MI · S4 joint_mvn_available test
DEFER: EM max_iter>0 · P1 claim/docs (joint Σ, 95% wording) · rebase onto main · GNN · DRAC
DISCIPLINE: verify=failing test then pass + focused testthat · compute=local · closure=4 P0s green
```

ARC PROGRAM: size · 2.5–4 h · Arc 0 = D-Blk1 cov align (or all four in one fan-out after G0)

Sweep: `docs/dev-log/handover/2026-08-08-ultra-plan-p0-fix-sweep.md` ([recon](7e6b882a-d634-4c50-909d-11b52624e503))  
Rose queue: `docs/dev-log/handover/2026-08-08-rose-close.md`

## DECISIONS LOCKED (Phase 0.4 defaults)
1. P0 only — not P1 DESCRIPTION/vignette claim rewrite  
2. Do **not** restore EM (`max_iter>0`)  
3. New branch `fix/p0-review-blockers` from **current dirty tree** (keep B1–B3); no rebase  

## SLICE TABLE

| Slice | Member | Model+effort | Bar | Time | Files | Dep |
| --- | --- | --- | --- | ---: | --- | --- |
| S0 RECON | recon | Composer low | Cursor Models | done | sweep receipt | — |
| S1 D-Blk1 | Emmy | Sonnet / Auto Cost | Other Models | 45–60 | `preprocess_traits.R` `impute.R` + shuffled-row×cov test | S0 |
| S2 B-Blk1 | Gauss/Emmy | Sonnet high | Other Models | 60–90 | `henderson_s_inv.R` `joint_mvn_solver.R` `joint_threshold_baseline.R` — SE = σ√(1−h_i) under max_iter=0; binary decode smoke | S0 |
| S3 C-Blk1 | Fisher | Sonnet high | Other Models | 45–60 | `fit_helpers.R` scores for zi_count; `multi_impute.R` draw SD; not E[X] SE | S0 |
| S4 B-Blk2 | Emmy | Sonnet medium | Other Models | 20–30 | `joint_mvn_baseline.R` `test-joint-baseline.R` | S0 |
| S5 VERIFY | recon | Composer low | Cursor Models | 20 | test_file the four test files; no invented paths | S1–S4 |
| S6 Rose | Rose | Auto Cost | Other Models | 15 | claim gate: do not say Level-C joint works | S5 |
| S7 Melissa | Melissa | Sonnet med | Other Models | 15 | plan-actual 2026-08-08-p0-fix.md | S6 |

PARALLEL: {S1,S2,S3,S4} after G0  
FAN-OUT BUDGET: checkpoint=`p0-fix-2026-08-08` · children ≤6 · scout=S0+S5 · build=4 · ceiling=0 (Rose on Sonnet/Auto unless claim fight)  
LUNA SUITABILITY: yes — S5 mechanical test run  
ULTRA EFFORT: no  
SEARCH: none  
ESTIMATE: ~2.5 h wall parallel · one `/goal` session · **START A FRESH TASK** (this planning chat is long)

## /goal paste (after G0)

```
/goal pigauto P0 review blockers (ultra-plan 2026-08-08)

READ: docs/dev-log/handover/2026-08-08-ultra-plan-p0-fix.md
      docs/dev-log/handover/2026-08-08-rose-close.md
      LOOP/GOAL.md if present — or scaffold LOOP/ for this lane

Branch: fix/p0-review-blockers from current dirty tree. Keep B1–B3. No rebase.
Execute S1–S4 (can parallel). Then S5–S7.
P0 only. Do not restore EM. Tests first.
Done: 4 P0s + NEWS + focused tests green.
```
