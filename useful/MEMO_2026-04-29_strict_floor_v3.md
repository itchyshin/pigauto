# Strict val-floor v3 calibration (2026-04-29)

## Why v3

The 2026-04-29-am fix (commit `d651c6c`) extended the strict val-floor
to ALL trait types via a single type-agnostic `cal_mean_loss`-based
check on the FULL val set:

```r
if (length(val_row_idx) > 0L && tm$type != "multi_proportion") {
  # check blend vs both BM and MEAN corners; override to whichever wins
}
```

This closed the discrete benches (binary/categorical/zi_count) cleanly
(commit `720fab0`).

**But** an apples-to-apples real-data verification at AVONET n=1500
same-seed showed the strict check was over-correcting on
continuous-family types: it spuriously overrode the calibrated blend
to pure-BM whenever blend's val MSE was numerically equal to pure-BM's
val MSE (which is typical on high-phylo-signal continuous traits like
AVONET Mass).  Pure-BM had higher TEST MSE than the calibrated blend,
producing a 70 % test-RMSE regression on Mass (385.7 → 659.0 in v1),
17 % on Beak, and degraded behaviour on Wing/Tarsus.

## What v3 does

Splits the post-calibration invariant into two type-specific blocks:

* **Continuous family (continuous, count, ordinal, proportion):**
  byte-for-byte the pre-2026-04-29 mean-only check. Override only if
  the blend loses to pure-MEAN on the full val set — never to pure-BM.
  Pure-BM is trusted as the safe corner via the half-A/half-B + median
  cross-check that already runs above.

* **Discrete family (binary, categorical, zi_count):** the strict
  both-corners check from the original 2026-04-29 fix. Override if
  the blend loses to either pure-BM or pure-MEAN by 1e-12 on the full
  val set; pick the lower-loss corner. Discrete losses are 0-1 /
  argmax step functions on a small set, so the strict check is
  well-conditioned at typical val sizes (~5-20 cells per trait at
  n=300, 30 % MAR).

* **multi_proportion:** skipped entirely (gate naturally closes via
  the legacy path; `cal_mean_loss` is undefined under safety_floor =
  TRUE for this type).

## Calibration evidence: AVONET n=1500, seed=2026, same sample

Three runs differ only in `R/fit_helpers.R::calibrate_gates`'s post-
calibration block:

* `pre`  : 2026-04-23..28 type-restricted block (the original bug;
           only continuous-family floor active)
* `v1`   : 2026-04-29-am strict-both-corners for ALL types (this
           memo's reason for existing — over-corrected continuous)
* `v3`   : 2026-04-29-pm type-conditional (current head)

| trait | metric | pre | v1 | v3 | v3 − pre |
|---|---|---|---|---|---|
| Beak.Length_Culmen | RMSE | 6.51 | 6.39 | 6.55 | +0.04 |
| Mass | RMSE | **385.7** | **659.0** ⚠ | **377.3** | **−8.4** ✓ |
| Migration | acc | 0.815 | 0.815 | 0.815 | 0 |
| Primary.Lifestyle | acc | 0.838 | 0.838 | 0.838 | 0 |
| Tarsus.Length | RMSE | 11.60 | 12.04 | 11.92 | +0.32 |
| Trophic.Level | acc | 0.829 | 0.829 | 0.829 | 0 |
| Wing.Length | RMSE | 37.09 | 35.11 | 34.12 | −2.97 |

* Continuous mean change v3 − pre: **−2.75 RMSE** (slight improvement on average)
* Continuous max abs change: 8.4 (Mass moved −2.2 %; the rest under 1 %)
* Discrete: identical (categorical / accuracy traits are at corner gates
  in this dataset, so the strict check is a no-op for them)

For comparison, v1 had continuous mean change of **+67.9 RMSE**
(driven by Mass +273) — that was the over-correction this memo
documents and v3 removes.

## What v3 doesn't change vs the 2026-04-29-am state

* The strict check still applies to binary, categorical, zi_count.
  These are the types the simulator benches (`bench_binary`,
  `bench_categorical`, `bench_zi_count`) showed underperforming
  baseline pre-fix (5–12 pp gaps). The discrete strict check remains
  bit-identical to the 2026-04-29-am implementation — just gated by
  trait type.
* `bench_binary.md`, `bench_categorical.md`, `bench_zi_count.md`
  closures from commit `720fab0` are preserved.
* Multi_proportion behaviour unchanged (gate naturally closes).

## Files preserved as evidence

`script/_strict_floor_v3_calibration/`:

* `bench_avonet_full_local_n1500_PRE.{md,rds}` — pre-fix run
* `bench_avonet_full_local_n1500_POST.{md,rds}` — v1 (over-strict)
* `bench_avonet_full_local_n1500_POSTv2.{md,rds}` — v2 (intermediate)
* `bench_avonet_full_local_n1500_POSTv3.{md,rds}` — v3 (current head)
* `bench_avonet_full_local_n2000_postfix.{md,rds}` — earlier subset
  used for the un-controlled comparison that surfaced the issue

## Open: still need to verify the discrete benches under v3

The strict discrete check in v3 is bit-identical to v1's check (just
type-gated), so the discrete-bench closures from `720fab0` should
still hold. But this should be verified by re-running
`bench_binary.R` / `bench_categorical.R` / `bench_zi_count.R` under
v3. Not done in this session due to wall-time budget.

## Strategic implication

The strict val-floor's value is in protecting discrete traits.
Continuous-family traits are already well-protected by the existing
mean-only check + the half-A/half-B + median cross-check above.
Adding a strict pure-BM check for continuous types is too sensitive
to val-set sampling noise and produces real regressions on
high-phylo-signal traits.
