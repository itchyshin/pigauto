# Phase 6 + Phase 7 EM smoke benchmark

n = 150, correlated binary (shared BM liability) + 1 continuous.
Held out 30% of each binary trait for accuracy scoring.

| method | wall (s) | acc b1 | acc b2 |
|---|---:|---:|---:|
| plug_in | 1.32 | 0.511 | 0.556 |
| phase_6_diag | 1.28 | 0.511 | 0.556 |
| phase_7_offdiag | 1.17 | 0.511 | 0.556 |

**This is a functional smoke test, not a performance claim.** All three
paths produce the same held-out accuracy on this specific n=150 synthetic
— the refined liability posteriors collapse back to the same binary class
after threshold-at-0 decoding. A proper performance characterisation needs
(a) larger n, (b) stronger-correlated simulation, and (c) multi-seed runs
to quantify variance. See the Phase 6 / Phase 7 spec files for the target
regimes where lift is expected.

What the bench IS verifying:
- all three paths complete without error
- EM paths run to completion without phylopars singularities
- wall-time overhead of EM vs plug-in is modest at n=150 (under 2×)
- `em_state$em_offdiag` records which path fired

