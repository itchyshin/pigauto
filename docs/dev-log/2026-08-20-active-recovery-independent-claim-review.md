# Stage A: independent claim review

## Materials reviewed

This review is a separate read-only assessment of the completed Stage-A
campaign record in
`docs/dev-log/2026-08-18-active-recovery-results.md`, against the
pre-registered ADEMP design in
`docs/dev-log/2026-08-18-active-recovery-design.md`.  It assesses the wording
that the retained 8,000-replicate campaign can support; it does not rerun,
pool, or exclude any receipt.

## Registered gate

The pre-registration required each family-specific, lambda = 1 core
condition to have a paired 95% Monte-Carlo interval below zero for the
active-minus-random primary outcome, with no policy having more than a
one-percentage-point excess failure rate.  Negative contrasts favour active
selection.

## Findings

| Family | n = 100, lambda = 1 | n = 300, lambda = 1 | Gate interpretation |
|---|---:|---:|---|
| Continuous | -0.01223 (-0.02190, -0.00256) | -0.00431 (-0.00715, -0.00147) | Both paired intervals are below zero. |
| Binary | -0.00045 (-0.00219, 0.00129) | -0.00033 (-0.00130, 0.00064) | Both paired intervals include zero. |

All four policies recorded zero failures in every retained cell, so the
failure-rate safeguard does not alter the result.  It also cannot compensate
for the binary primary intervals including zero.

## Verdict

The broad Stage-A headline is **not earned**: this campaign does not establish
that `suggest_next_observation()` improves realised follow-up imputation
across the tested families.  The continuous BM, lambda = 1, one-step,
cell-level result is earned; the binary result is directionally favourable
but inconclusive.

The permitted wording is:

> In the preregistered continuous BM lambda = 1 setting, choosing the active
> recommendation reduced fixed-test normalized MSE versus a matched random
> choice at n = 100 and n = 300. The analogous binary contrasts were
> directionally favourable but inconclusive. This is a continuous-specific
> recovery result, not a cross-type active-imputation claim.

The campaign still does not support a benefit over uncertainty-only
acquisition, weak-signal benefit, mixed-type or species-level ranking,
batch acquisition, real-data benefit, GNN attribution, a default change, or
a general active-imputation claim.
