# Stage B: corrected one-mask timing smoke

## Purpose and scope

This is a bounded operational receipt, not a claim-bearing external
comparison. It ran the locked first mask (`20260901`) of the two registered
AVONET regimes after `main`'s exact-conditional route was merged into the
evidence branch. It establishes that the two pigauto arms execute through the
registered runners and supplies a measured base runtime for later planning.
No one-mask score is interpreted here.

## Environment and provenance

The isolated Totoro checkout was
`/home/snakagaw/pigauto-stageb-timing-69f0fa5`, refreshed to source SHA
`d30b5e38872df1cac87e5dd9d529aac8644488c9`. The package was installed into
that checkout's private `Rlib`; `impute()` was checked to expose
`predict_method`. The run used `OPENBLAS_NUM_THREADS=1` and
`OMP_NUM_THREADS=1`, a single CPU process, and a two-hour `timeout` cap.

Results are retained at
`/home/snakagaw/pigauto-stageb-timing-69f0fa5/results/timing-smoke-20260901-exact-main-retry2`.

## Completion checks

Both regime receipts completed and wrote their aggregate RDS, markdown, and
per-mask receipt.

| Regime | pigauto default | pigauto exact | Total measured wall time |
|---|---:|---:|---:|
| Continuous core | 109.976 s | 172.388 s | 283.719 s |
| Mixed core | 110.074 s | 170.115 s | 281.163 s excluding unavailable BACE |

`Rphylopars`, `phylolm`, and `missForest` completed in both applicable
regimes. Both pigauto arms had no error rows. The log retained the expected
small-validation and numerical-solver warnings rather than suppressing them.

## Comparator availability boundary

`BACE` was unavailable in the isolated library and its rows are retained as
errors (`BACE not installed`) in both receipts. A read-only CRAN index query
found no BACE package. An initial public GitHub search under `itchyshin` found
no source, but the routed project record identifies the upstream package as
`daniel1noble/BACE`. Therefore this smoke is sufficient to time the
implemented arms but is not a complete Stage-B comparator run. Install and
verify that upstream source in the private library before a claim-bearing
campaign; do not silently remove BACE from the protocol.

## What this permits and does not permit

At the observed base timings, five masks in each of the continuous and mixed
regimes would take about 47 minutes sequentially before BACE. That is only a
planning estimate, not approval for a full campaign: BACE timing remains
unknown and its source/installation is unresolved. This receipt makes no
competitiveness, parity, default-change, mixed-type, or exact-solver
performance claim.
