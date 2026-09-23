# Session Handoff (to Claude): MI under phylogenetic GLS

Meta: 2026-09-23, from Claude (Opus 5.5 session). Branch `arc/mi-gls-attenuation`,
worktree `pigauto-mi-gls`, pushed. No PR yet.

## Critical Context

1. pigauto's `multi_impute()` draws (conformal, either prediction route, and MC dropout)
   are not valid for a phylogenetic GLS analysis: PGLS slope bias -0.20 to -0.46 for the
   conformal draws, coverage near 0, in all 16 simulated regimes
   (`docs/dev-log/mi-gls/results.md`).
2. The internal prototype `draw_conditional_bm()` (`R/draws_conditional.R`), drawing
   missing cells jointly from their Gaussian conditional with an EM trait covariance, is
   unbiased at lambda 1 (|bias| <= 0.009, coverage 0.84-0.93 vs complete 0.88-0.98).
   At lambda 0.5 it is biased by +0.04 to +0.06 because it assumes lambda 1.
3. pigauto's joint-MVN plug-in covariance (`fit_mvn_bm_inhouse(max_iter = 0)` in
   `R/joint_mvn_solver.R`; `max_iter = 50` identical) shrinks cross-trait correlation under
   missingness (0.67 to 0.46 with 30% of both traits missing). That file belongs to the
   lambda lane in another Claude account; reported, not edited.

## Landing State

| Artifact / branch | Committed | Pushed | PR | State |
|---|---|---|---|---|
| `arc/mi-gls-attenuation` (2709e1a) | y | y | none | LANDED on branch; PR owed after Shinichi's decision below |

FINDING-OF-RECORD: conformal and MC-dropout MI draws are invalid for PGLS; proper conditional draws with an EM covariance fix it at lambda 1; the joint-MVN plug-in covariance shrinks cross-trait correlation under missingness  vault-note: [[pigauto-mondrian-realdata-and-mi-gls]] (proposed; draft on `arc/mondrian-realdata` at `docs/dev-log/mondrian-realdata/vault-note-draft.md`; not written without Shinichi's approval)

## Next Immediate Steps

1. Shinichi decides whether `draw_conditional_bm()` becomes a `multi_impute()` option
   (`draws_method = "conditional"`, continuous traits) and whether any default changes.
2. Make the draw lambda-aware: condition on the lambda-transformed tree, using the lambda
   estimate from `feat/joint-lambda-default` once it lands; rerun regimes 2, 4, 6, 8, 10,
   12, 14 and 16 in fast mode (`MI_GLS_FAST=1`, about 6 core-hours on fir).
3. Tell the lambda lane about the plug-in covariance shrinkage (finding 3), with the
   six-tree check in `docs/dev-log/mi-gls/RUNLOG.md`.

## Gotchas

- fir needs `module load cuda/12.6` plus `LD_LIBRARY_PATH=$EBROOTCUDA/targets/x86_64-linux/lib`;
  `/project` quota is full, use `/scratch` (purged after about 60 days; receipts are
  already copied into `script/mi_gls/returned*`).
- A fast-mode cell takes about 39 s at n = 1000 (dense EM); GNN-mode cells 5-10 min, a
  few n = 1000 both-missing cells exceed 1 h.

## How to Resume

```text
Read AGENTS.md and docs/dev-log/handover/2026-09-23-claude-handover-mi-gls.md on branch arc/mi-gls-attenuation, then continue only the steps Shinichi has approved.
```
