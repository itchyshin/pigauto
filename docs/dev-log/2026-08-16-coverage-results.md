# Conformal coverage campaign — results

Complete 2026-08-16 ~08:15 MDT on Totoro. **1,920 jobs, ok=1920 skip=0 fail=0**, wall 2152 s.
Design: `2026-08-15-coverage-campaign-design.md`. Pipeline `main` @ `c655d75`.
Full table: `2026-08-16-coverage-summary-table.md`.

## Headline — the intervals undercover at n = 100

Nominal level 0.95. A cell is "adequate" iff its upper 95% MCSE bound reaches 0.95.

| n | coverage (range over family × λ × split) | verdict |
|---:|---|---|
| **100** | **0.853 – 0.861** | **INADEQUATE — 6/6 cells fail** |
| 300 | 0.953 – 0.957 | adequate, 6/6 |
| 1000 | 0.948 – 0.963 | adequate, 6/6 |

The split is clean: **every n = 100 cell fails and every n ≥ 300 cell passes**, across both
DGP families, both λ values, and both `conformal_split_val` settings.

## `conformal_split_val` does not fix it (S1)

| n | split = FALSE | split = TRUE |
|---:|---:|---:|
| 100 | 0.853 – 0.861 | 0.854 – 0.861 |
| 300 | 0.953 – 0.957 | 0.953 – 0.957 |
| 1000 | 0.948 – 0.952 | 0.954 – 0.963 |

At n = 100 the two settings are indistinguishable. The documented rationale for
`conformal_split_val` — that reusing validation cells for both gate calibration and the
conformal quantile causes undercoverage "most visible at small n_val" — is **not** what is
driving the n = 100 failure, because removing that reuse does not repair it.

## The mechanism is arithmetic, not statistical

Split conformal takes the `⌈(1−α)(n_val+1)⌉`-th order statistic of the calibration residuals.
At 95% that needs `⌈0.95 × (n_val+1)⌉ ≤ n_val`, i.e. **n_val ≥ 19**. At n = 100 with the
default splits, n_val per trait is roughly 10 — the campaign's own fits logged
"Small validation set … (n = 10)" repeatedly. With 10 calibration residuals the largest
achievable coverage is `n_val/(n_val+1) ≈ 0.909`, **so 95% is unreachable before any noise
is considered.** Observed 0.85–0.86 sits below that ceiling, consistent with quantile noise
on top of a hard arithmetic bound.

This also retires a loose end: the only prior number on record, **0.884–0.887** from the
v0.9.3 ablation, was dismissed as leak-tainted. It was not wrong. It was measuring this.

## S2 — the accuracy cost of `split_val = TRUE` is small

| family | n | mean cost | MCSE | % of baseline loss |
|---|---:|---:|---:|---:|
| F1 | 100 | +0.00088 | 0.0086 | +0.29% |
| F2 | 100 | +0.00106 | 0.0187 | +0.21% |
| F1 | 300 | −0.00066 | 0.0041 | −0.24% |
| F2 | 300 | +0.00159 | 0.0102 | +0.34% |
| F1 | 1000 | +0.00316 | 0.0040 | +1.29% |
| F2 | 1000 | +0.00349 | 0.0097 | +0.86% |

**0.2 – 1.3%**, and at n = 300 for F1 it is negative. The roxygen for
`conformal_split_val` says enabling it "regresses the AVONET300 / OVR-categorical / BIEN
safety-floor smoke benches by 2-26%". On simulated data at these sizes the cost is an order
of magnitude smaller. Those bench numbers may be real for those specific fixtures, but the
"2-26%" figure should not be read as the general cost of the option.

## What this licenses, and what it obliges

**Obliges.** `DESCRIPTION`, `README`, and `gnn-architecture.Rmd` §7 present ≥95% marginal
coverage as the conformal guarantee. That guarantee **cannot hold at n = 100 with default
splits** — it is arithmetically unreachable, not merely unmet. This is the same class of
defect as #157 and the "Exact under BM" SE claim: a documented property the shipped
configuration does not deliver. It needs either a documented minimum n, an automatic guard,
or both.

Concrete options (not chosen here — Shinichi's call):
1. **Warn or error** when `n_val < 19` for a trait, naming the achievable ceiling.
2. **Raise `val_frac` automatically** at small n so n_val clears 19.
3. **Document a minimum n** for the conformal claim and scope the guarantee to it.

**Licenses.** At n ≥ 300 the intervals are well calibrated (0.948–0.963, six for six). That
is a real, defensible, quantified claim the package could not previously make.

## Regime

Simulated BM/λ DGPs (F1 linear, F2 nonlinear), continuous traits only, MCAR at m = 0.30,
single-obs, 100 reps (n = 100/300) and 40 (n = 1000), one machine. No real-data anchor —
AVONET300 conformal coverage remains a Tier-3 item.

## CORRECTION (2026-08-16, same day, after Shinichi asked "didn't we check coverage before?")

**This document's original scope claim was wrong, and the error was mine.** It read the
n >= 300 simulated cells (0.948-0.963) as "well calibrated at n >= 300". That generalises a
**MCAR simulation with equal per-trait observation counts** to the package's behaviour at
large n. The repo's own real-data benches contradict it, and predate this campaign:

| bench | n | conformal coverage |
|---|---:|---|
| `bench_fishbase.md` | 10,654 | 0.893-0.955 |
| `bench_bien.md` | 4,745 | 0.896-0.977 |
| `bench_pantheria_full.md` | 4,027 | 0.868-0.950 |
| `bench_avonet_full_local.md` | 1,500 | 0.913-0.958 |
| `bench_pantheria_bace_head_to_head.md` | 500 | **0.653**-0.943 |

**These are two different defects and must not be merged.**

1. **n = 100, simulated (this campaign).** Arithmetic. Split conformal needs
   n_val >= 19 for 95%; n = 100 gives n_val ~ 10, ceiling n_val/(n_val+1) ~ 0.909.
   Real, cheap to document, and unremarkable -- ordinary small-sample behaviour.
2. **Large n, real data (the benches, long-standing).** **NOT** the ceiling. fishbase
   `Length` has 3,103 test cells, so the order-statistic ceiling is ~0.998; observed
   coverage is 0.906. The arithmetic explanation cannot apply. The leading candidate is
   the one assumption split conformal actually requires: **exchangeability between the
   calibration and test cells.** Real missingness is not MCAR -- species missing a trait
   differ systematically from those that have it -- so validation residuals are not a
   valid reference distribution for the missing cells. That is a methodological
   limitation of applying split conformal to phylogenetic imputation, not a sample-size
   artifact, and it is the half that bears on publication.

**Untested and directly relevant:** `script/bench_missingness_mechanism.R` varies
MCAR/MAR/MNAR and has not been run against this question. It is the experiment that
would confirm or kill the exchangeability explanation.

**Do not repeat the original error:** this campaign simulated MCAR only. No coverage
claim derived from it transfers to real data.
