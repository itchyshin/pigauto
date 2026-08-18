# Stage A: active-recovery implementation and completed recovery campaign

## Status

This record preserves the 100-replicate-per-cell pilot and the completed,
pre-registered 1,000-replicate-per-cell campaign.  The broad headline
(`suggest_next_observation()` improves realised follow-up imputation across
the tested families) is not earned: both binary lambda=1 paired intervals
include zero.  A narrower continuous, lambda=1, one-step cell-level recovery
result is supported below.

The harness first completed a local 30-tip, five-epoch structural smoke for
both continuous and binary outcomes.  It fit initial data, selected
active/random/uncertainty candidates from the same candidate set, restored
exactly one true value, refit each policy with the same fitting seed, and
made finite predictions on protected fixed-test data.  Its small-validation
warnings were expected and non-evidential.

## Totoro pilot receipt

The Totoro campaign ran source SHA `beee5df360759f53cd22b595366436bca50845a5`
with 96 workers and one BLAS/OpenMP thread per worker.  It started at
2026-08-18T14:46:32Z and completed at 2026-08-18T14:53:23Z.  The retained
artefacts are under
`/home/snakagaw/pigauto-active-recovery-results/active-recovery-pilot-20260818-beee5df`.

The independent receipt audit found 800 receipts, 3,200 policy rows, exactly
four policy rows per replicate, one source SHA, zero failed fits, and the
expected treatment hashes: restored-policy data differed from initial data,
whereas no-acquisition data did not.  The launcher process group was empty
after completion.

Negative active-minus-random values favour active selection.  Intervals are
paired Monte-Carlo 95% intervals from 100 pilot replicates.

| n | lambda | family | active minus random (95% MC interval) |
|---:|---:|---|---:|
| 100 | 0.2 | binary | -0.00057 (-0.00246, 0.00132) |
| 300 | 0.2 | binary | -0.00057 (-0.00130, 0.00017) |
| 100 | 1.0 | binary | 0.00062 (-0.00523, 0.00647) |
| 300 | 1.0 | binary | -0.00195 (-0.00396, 0.00006) |
| 100 | 0.2 | continuous | 0.01018 (-0.01383, 0.03419) |
| 300 | 0.2 | continuous | 0.01109 (-0.00223, 0.02441) |
| 100 | 1.0 | continuous | -0.03542 (-0.05875, -0.01208) |
| 300 | 1.0 | continuous | -0.01055 (-0.02086, -0.00025) |

The lambda=1 continuous pilot cells are directionally compatible with active
selection improving fixed-test error.  The binary lambda=1 cells do not yet
give a reliably negative active-minus-random interval.  These pilot patterns
do not satisfy the pre-registered headline gate and must not be reported as
an active-imputation benefit claim.

## Full-campaign gate

For every cell, the pilot standard deviation implies the registered minimum
of 1,000 replicates per cell.  A complete campaign is therefore 8,000 total
replicates (7,200 additional after this pilot).  Mean four-policy
replicate time was 41.77 seconds; the observed cell medians ranged from
35.56 to 48.69 seconds.  At 96 workers with 20% scheduling/I/O headroom,
the projected cost is about 1.16 hours for all 8,000 replicates, or 1.04
hours for the additional 7,200.

## Full Totoro campaign

After explicit approval, the extension ran replicates 101--1,000 in every
cell at source SHA `ad3990c5bea62be7536d4e0ab2c81ea8f9ef686e`, retaining its
results separately at
`/home/snakagaw/pigauto-active-recovery-results/active-recovery-extension-20260818-ad3990c`.
It used 96 workers, with `OPENBLAS_NUM_THREADS=1` and the other BLAS/OpenMP
thread caps set to one.  It started at 2026-08-18T15:16:33Z and finished in
about 56 minutes, within the pilot-based 1.04-hour extension estimate.

An independent receipt audit passed before the combined estimator was read:
800 pilot plus 7,200 extension RDS receipts, 32,000 policy rows, exact
registered eight-cell / 1,000-replicate grid, four policy rows per replicate,
zero error rows, the two expected source SHAs, finite primary scores,
unchanged no-acquisition hashes, and changed hashes plus non-empty selected
species for all acquisition policies.  No worker remained after completion.
The retained audit and summary are `receipt_audit.rds` and
`combined_summary.rds` in the extension directory.

Negative active-minus-random values favour active selection.  Intervals are
paired Monte-Carlo 95% intervals from 1,000 replicates per cell.

| n | lambda | family | active minus random (95% MC interval) | active minus uncertainty (95% MC interval) |
|---:|---:|---|---:|---:|
| 100 | 0.2 | binary | 0.00016 (-0.00052, 0.00085) | 0.00004 (-0.00068, 0.00076) |
| 300 | 0.2 | binary | -0.00014 (-0.00044, 0.00017) | -0.00030 (-0.00062, 0.00003) |
| 100 | 1.0 | binary | -0.00045 (-0.00219, 0.00129) | -0.00042 (-0.00197, 0.00112) |
| 300 | 1.0 | binary | -0.00033 (-0.00130, 0.00064) | -0.00032 (-0.00119, 0.00056) |
| 100 | 0.2 | continuous | -0.00159 (-0.01010, 0.00692) | -0.00166 (-0.00934, 0.00601) |
| 300 | 0.2 | continuous | 0.00541 (-0.00013, 0.01096) | -0.00006 (-0.00435, 0.00423) |
| 100 | 1.0 | continuous | -0.01223 (-0.02190, -0.00256) | -0.00043 (-0.00681, 0.00594) |
| 300 | 1.0 | continuous | -0.00431 (-0.00715, -0.00147) | -0.00189 (-0.00402, 0.00025) |

All four policies had zero failures in every cell, so the no-more-than-one
percentage-point excess-failure safeguard passes.  It does not rescue the
headline gate: continuous lambda=1 passes against random at both sizes, but
binary lambda=1 is inconclusive at both sizes.  The independent claim review
therefore permits only this wording:

> In the preregistered continuous BM lambda=1 setting, choosing the active
> recommendation reduced fixed-test normalized MSE versus a matched random
> choice at n=100 and n=300.  The analogous binary contrasts were
> directionally favourable but inconclusive.  This is a continuous-specific
> recovery result, not a cross-type active-imputation claim.

The campaign does not support superiority over uncertainty-only acquisition,
weak-signal benefit, mixed-type or species-level ranking, batch acquisition,
real-data benefit, GNN attribution, or any default change.
