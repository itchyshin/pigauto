# Stage C: FishBase and PanTHERIA Mondrian pilot

This is a bounded execution and harness receipt, not calibration evidence.
It used source `9092224` (the runner explicitly caps Torch and OpenBLAS at
one compute thread), 300-tip deterministic subsets, one observed-only 20%
mask per dataset, seed `20260818`, and ten epochs. The two methods used the
identical mask within each data set. Retained Totoro receipts are under
`/home/snakagaw/pigauto-stagec-pilot-results/{fishbase,pantheria}300-pilot-threadcapped/`.

The initial attempt was stopped immediately because Torch ignored the
OpenBLAS-only cap and used multiple compute cores. It had already written
only its mask receipts. It is retained separately as a non-evidential
operational receipt. The rerun obeyed the Torch cap and completed both
methods for both datasets.

## Pilot outcome

All four fit receipts were `status = "ok"`, and the fail-closed summariser
wrote both paired tables. The small validation samples make many traits
incapable of supporting a 95% split-conformal claim; warnings were retained.

- FishBase: Mondrian changed interval width only for `Length` (82.99 split;
  95.32 Mondrian); all other pilot numeric trait widths and coverages were
  identical. `Weight` coverage was 0.30 on ten masked cells, so the pilot is
  explicitly not a calibration result.
- PanTHERIA: the two methods were identical for every reported numeric trait
  in this small pilot. Coverage ranged from 0.50 (`litter_size`) to 1.00 for
  two traits; this variation is exactly why it cannot support a scientific
  comparison.

## Gate

The harness, real-data preparation, observed-only mask protection, method
receipts, and one-thread control are operationally verified. A full Stage-C
campaign needs a new explicit compute approval after a design update that
increases the validation/calibration information and uses multiple retained
masks. It may not claim Mondrian improves calibration, alter the default, or
offer an unconditional coverage guarantee.
