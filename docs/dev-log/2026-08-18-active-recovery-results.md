# Stage A: active-recovery implementation, smoke, and pilot receipt

## Status

This is a pre-registered 100-replicate-per-cell pilot, not a headline
performance claim.  It establishes treatment integrity, measured variance,
failure rates, and a costed full-campaign gate.

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

No full campaign runs unless explicitly approved after review of this pilot.
