# Stage A: active-imputation recovery pre-registration

## Aim and boundary

This study asks one narrow question: after one newly measured trait value, does
the cell selected by `suggest_next_observation()` improve realised fixed-test
imputation more than a random eligible acquisition?  It is not evidence about
batch acquisition, mixed-type global optimality, real-data benefit, or the GNN
component.

## ADEMP contract

The 8 cells are `n = 100, 300` crossed with `lambda = 1, 0.2` and family
`continuous, binary`.  Each replicate draws `tree ~ rtree(n)` and
`z ~ N(0, lambda R + (1-lambda) I)`, where `R` is the tree correlation
matrix.  The binary outcome thresholds `z` at its replicate median.

Species are partitioned once per replicate into 60% initially observed, 20%
eligible candidates, and 20% fixed test species.  Test species are never
eligible for selection.  Active selection filters the public cell ranking to
the common candidate set; random selection is uniform on that same set;
uncertainty-only selection uses decision-time predictive SE (continuous) or
binary predictive entropy.  Exactly one selected true value is restored, then
the ordinary single-observation `impute()` workflow is refit with the same
policy-independent fitting seed.

Continuous primary outcome is paired active-minus-random normalised MSE;
binary primary outcome is paired active-minus-random Brier score.  Negative
values favour active.  Secondary outcomes are active versus uncertainty,
no-acquisition loss, RMSE, log loss, accuracy, win rate, runtime, and retained
errors.  Every receipt records source SHA, selected cell, restored value,
initial/changed data receipts, and elapsed time.

## Replication and claim gate

Run 100 pilot replicates per cell first.  The eventual count per cell is
`max(1000, ceiling(pilot_sd^2 / 0.01^2))`; realised MCSEs are reported.  The
headline is supported only when both family-specific `lambda = 1` cells have
paired 95% Monte-Carlo intervals below zero and active has no more than a one
percentage-point excess failure rate.  Weak-signal results remain reported as
sensitivity evidence.

No full campaign runs until the pilot results, cost estimate, and integrity
review have been presented for approval.
