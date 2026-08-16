# Continuous-gap diagnosis — it's the in-house joint solver

2026-08-16 evening · Claude lane · follow-up to
`2026-08-16-external-comparison-results.md` (raw Rphylopars beats pigauto's continuous
output on AVONET300, 5 seeds). Script `dev/diagnose_continuous_gap.R` +
`.log`/`.rds` on `arc/bace-comparators`. **One mask (seed 2027 = that bench's rep 1) —
directions and rough magnitudes, not intervals.**

## Design

Four arms on the identical mask, scored with the bench's exact z-RMSE:
A = pigauto full 7-trait fit (blend AND decoded baseline); B = raw
`Rphylopars::phylopars()` on the 4 continuous traits, all observed cells; C = same as B
but with pigauto's val/test cells ALSO hidden (the calibration tax on the same solver);
D = pigauto on the 4 continuous traits only (isolates the mixed-type path).

## Result: one layer owns the gap

| layer (positive = hurts) | Mass | Beak | Tarsus | Wing |
|---|---:|---:|---:|---:|
| gate effect (blend − own baseline) | −0.075 | +0.001 | −0.010 | −0.002 |
| mixed-type path (7-trait − 4-trait baseline) | +0.002 | +0.046 | −0.424 | −0.118 |
| calibration tax (phylopars, fewer cells) | +0.046 | −0.001 | +0.200 | +0.072 |
| **solver wrapper (in-house − phylopars, comparable data)** | **+0.141** | **+0.360** | **+1.269** | **+0.278** |
| TOTAL (pigauto blend − raw phylopars) | +0.114 | +0.406 | +1.034 | +0.229 |

Three exonerations and one conviction:

1. **The GNN/gate is innocent** — the blend is at or slightly *better* than pigauto's own
   baseline on every trait (Mass −0.075). The safety machinery does its job on real data.
2. **The mixed-type path is innocent, even beneficial** — the 7-trait threshold-joint
   baseline beats the 4-trait continuous-only one on Tarsus (−0.424) and Wing: the
   discrete traits contribute cross-trait information. Unified mixed-type imputation is
   pulling its weight.
3. **The calibration tax is real but modest** (+0.05 to +0.20) — the price of holding out
   val/test cells from the baseline. Worth a look later (refit-on-all-observed at predict
   time), not the story.
4. **The in-house joint solver is the gap.** `fit_mvn_bm_inhouse()`
   (`R/joint_mvn_solver.R`, `max_iter = 0`: single-pass Henderson init, closed-form Σ,
   cross-trait EM refinement disabled 2026-05-17 after divergence) loses 0.14–1.27
   z-RMSE to phylopars' converged REML on comparable data. Verified in code: the
   threshold-joint path calls the in-house solver, NOT Rphylopars
   (`joint_threshold_baseline.R:318–322`) — CLAUDE.md/AGENTS.md still describe the joint
   path as "via `Rphylopars::phylopars()`", which is stale.

**This answers Tier 2.1's open question with a measurement.** "Is the closed-form Σ
estimate adequate?" — on AVONET300, no: the single-pass solver, not the neural network,
is pigauto's continuous-trait accuracy bottleneck. Both external losses (BACE, raw
phylopars) trace to the same layer.

## Caveats

- Single mask. The 5-seed bench pins the total; this decomposes one rep of it.
- Arm C emitted repeated `solve(): system is singular` warnings from phylopars
  (approximate solutions) — phylopars' numerics are not unconditionally robust on this
  data, which may be part of why the in-house solver exists. Any solver swap must watch
  for this.
- Tarsus's large numbers partly reflect this rep (bench mean Tarsus total ≈ 0.58; this
  rep's 1.03 is the high side).

## Recommended next slice (not started — Shinichi's call)

The in-house solver's output contract already matches phylopars
(`$anc_recon` / `$anc_var` / `$pars$phylocov` — by design, see the comment at
`joint_threshold_baseline.R:318`). So the cheap, reversible move per the "add a
parameter, default off" doctrine is an opt-in **`joint_solver = c("inhouse",
"rphylopars")`** on `fit_baseline()`: one dispatch line, no contract change, Rphylopars
already a hard requirement of the joint gate. Then re-run the 5-seed external bench with
the switch — expected recovery to ~arm-C level (0.37–0.87), i.e. most of the gap.
The deeper alternatives (restore EM under a known-Σ G0; improve the Henderson init)
remain on the Tier-2 ladder, but the switch gives users the accuracy now and turns the
EM question from a blocker into an optimisation.
