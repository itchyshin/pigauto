# Mondrian real-data confirmation: M2 re-derivation and decision-rule audit

Auditor: independent statistical review (Gauss lens), read-only on all files except this one.
Repo: pigauto-mondrian-realdata, branch arc/mondrian-realdata.

## M2 RE-DERIVATION

Target: PanTHERIA, structured arm, trait `gestation_d`, far stratum, seeds 20260818/19/20.

Method: read `split.rds` and `mondrian.rds` directly from each of
`script/mondrian_confirmation/returned/pantheria-structured-m2026081{8,9}/` and
`.../pantheria-structured-m20260820/`, filtered `$cells` to
`trait == "gestation_d" & stratum == "far"`, computed per-mask coverage as
`mean(truth >= lo & truth <= hi)`, and pulled `n_far` from
`mondrian.rds$mondrian$gestation_d$n_far`. Pooled coverage across the 3 masks as
`sum(cov_i * n_test_i) / sum(n_test_i)` (equivalent to total-covered / total-n_test),
pooled `n_far` as the sum across masks, then applied the MCSE formula stated in the
pre-registration (`00-preregistration.md`, "Uncertainty" section) and in `results.md`'s
header: per-method `sqrt(coverage*(1-coverage)/n_test + alpha*(1-alpha)/(n_s+2))`, paired
MCSE `sqrt(mcse_mondrian^2 + mcse_split^2)`. None of `02_summarise_masked_confirmation.R`
or `12_build_results_doc.R` was sourced or imported; this is a from-scratch computation
against the raw receipts.

Per-mask figures (my own code):

| seed | n_test | cov_split | cov_mondrian | n_far |
|---|---|---|---|---|
| 20260818 | 164 | 0.975610 | 0.975610 | 68 |
| 20260819 | 148 | 0.932432 | 0.986486 | 68 |
| 20260820 | 128 | 0.968750 | 0.976562 | 60 |

Pooled (my own code):

| quantity | my value | rounded (4dp) |
|---|---|---|
| n_test_split = n_test_mondrian | 440 | 440 |
| n_far (pooled) | 196 | 196 |
| coverage_split | 0.959091 | 0.9591 |
| coverage_mondrian | 0.979545 | 0.9795 |
| paired_diff (coverage_gain) | 0.020455 | 0.0205 |
| mcse_split | 0.018140 | -- |
| mcse_mondrian | 0.016895 | -- |
| mcse (paired) | 0.024789 | 0.0248 |

Regenerated `docs/dev-log/mondrian-realdata/results_table.csv` with:
`OPENBLAS_NUM_THREADS=1 Rscript script/mondrian_confirmation/12_build_results_doc.R`
(writes only `results.md` / `results_table.csv`, confirmed via `git status` before/after —
both files are `.gitignore`d, not tracked, so nothing else in the repo changed). Note: a
concurrent process in this shared session deposited `fishbase-structured-m20260818/`
receipts between my first read of `results.md` and this regeneration, which added
FishBase rows to Table 1/2 — irrelevant to the PanTHERIA row under audit, which is
unaffected (confirmed byte-identical before/after for that row).

Row from the regenerated `results_table.csv`
(`pantheria,structured,gestation_d,far`):

```
coverage_split=0.959090909090909, coverage_mondrian=0.979545454545454,
coverage_gain=0.0204545454545454, mcse=0.0247892391117775,
n_test_split=440, n_test_mondrian=440, n_val=196, n_near=197, n_far=196
```

**Agreement**: matches my independent computation to 4 decimals on every reported
quantity (coverage_split 0.9591, coverage_mondrian 0.9795, coverage_gain 0.0205,
mcse 0.0248), and to 6 decimals on the raw values I printed. No discrepancy found.

## RULE AUDIT

Reviewed `script/mondrian_confirmation/08_apply_decision_rule.R` and
`script/mondrian_confirmation/12_build_results_doc.R` against
`00-preregistration.md` and its two amendments.

### 1. Condition 1 — structured arm, far stratum, cond1_eligible, median gain >= 0, min far cov >= 0.90

**DEVIATION** — `08_apply_decision_rule.R:76-83`.

Pre-registration's opening clause: "Flip the default ... only if all of the following
hold **on every dataset** where Mondrian activates for at least one trait, with Holm
adjustment across datasets for **the one-sided tests**." Read together with condition 2
being the only condition phrased as a significance test (the only one needing Holm
correction), the natural reading is: conditions 1 and 3 are point-estimate thresholds
evaluated **per dataset**, each of which must individually pass; condition 2 is the one
p-value-based test, Holm-adjusted **across** datasets. `12_build_results_doc.R`'s Table 2
instantiates exactly this per-dataset reading (`build_table2()`, one row per dataset,
`median_far_gain_structured` computed per dataset at line 230).

`08_apply_decision_rule.R` does not split condition 1 by dataset. `far_struct` (line 76)
pools eligible structured-far rows from **all** datasets together, and
`stats::median(far_struct$coverage_gain)` (line 81) is one global median across every
dataset's traits combined. The `all(far_struct$coverage_mondrian >= 0.90)` clause
(line 82) is unaffected by this pooling (`all()` is associative across a partition), so
only the median-gain half of condition 1 is at risk.

Demonstrated with the current (partial, live) receipts: eligible structured-far traits
pooled globally give median gain 0.0155, versus PanTHERIA-only median gain 0.0110 and
FishBase-only median gain 0.0164 computed per dataset — the three numbers differ
materially, confirming the two aggregation rules are not interchangeable. In the current
data neither per-dataset value is negative, so this particular divergence would not have
flipped today's verdict, but the mechanism is real: a dataset with few, strongly positive
traits can pull a weak or negative dataset's pooled median across zero, which the
pre-registration's "on every dataset" language was written to prevent.

### 2. Condition 2 — near stratum non-inferiority, 2pp margin, Holm across datasets

**PASS** (with one adjacent latent-default issue, see item 5) —
`08_apply_decision_rule.R:87-103`, `one_sided_binom_z()` at lines 51-59.

Per-dataset split (`per_dataset <- split(near, near$dataset)`, line 90) correctly
implements "Holm adjustment across datasets for the one-sided tests." p-value convention,
confirmed: `one_sided_binom_z()` tests H0: `p_mondrian - p_split <= -0.02` (Mondrian is
inferior by more than the margin) vs H1: `p_mondrian - p_split > -0.02` (non-inferior);
it returns `1 - pnorm(z)`, the upper-tail p under H0. **Small p rejects inferiority, i.e.
small p = non-inferior**, and the rule correctly requires `all(adj <= 0.05)` to pass —
this is the standard, correct one-sided non-inferiority convention. Arm handling matches
Amendment 1 ("condition 2 ... use every arm that exists"): no arm filter is applied, so
whatever arms exist in the near stratum for a dataset (structured, mcar, or both) all
feed the pooled per-dataset coverage.

### 3. Condition 3 — near stratum, paired median half-width ratio <= 1.10

**DEVIATION**, same root cause as item 1 — `08_apply_decision_rule.R:106`.

`stats::median(near$width_ratio) <= 1.10` pools `near` (all datasets, all arms) into one
global median, not one check per dataset. `12_build_results_doc.R`'s `build_table2()`
again computes this per dataset (`near_width_ratio`, line 264, `stats::median(t1_near$width_ratio)`
inside the per-dataset loop). Demonstrated on current data: global pooled median width
ratio is 0.8247, versus per-dataset values 0.8078 (avonet), 0.7766 (fishbase), 0.8340
(pantheria) — again materially different numbers from the same underlying rows,
confirming the same inconsistency as condition 1.

### 4. Fallback traits count as no evidence, never a pass

**PASS** — `08_apply_decision_rule.R:69`, `active <- long[!isTRUE_vec(long$fallback), ]`.
Fallback rows are dropped before any of the three conditions are computed, so a
fallback trait can never contribute a passing data point. A dataset where every trait
fell back contributes zero rows to `active` and is therefore implicitly excluded from
`far_struct`, `near`, and the per-dataset Holm split — which also correctly implements
"datasets where Mondrian never activates are excluded from every activating dataset"
without needing an explicit filter for it.

### 5. Condition 2's default when there is no near-stratum evidence at all

**DEVIATION (latent, currently non-exploitable)** — `08_apply_decision_rule.R:88` vs
`:106`.

`cond2 <- TRUE` is the initial value (line 88), overwritten only `if (nrow(near))`
(lines 89-103). If `near` has zero rows, `cond2` stays `TRUE` — a vacuous pass with no
evidence. This is the opposite default from `cond3` on line 106
(`if (!nrow(near)) FALSE else ...`), which fails closed in exactly the same
zero-near-rows scenario, and it is the opposite of the explicit "no evidence, never a
pass" principle the pre-registration states for fallback traits. As written today this
cannot by itself flip the verdict, because `pass <- cond1 && cond2 && cond3` and `cond3`
is `FALSE` whenever `near` is empty — `cond3`'s fail-closed default currently masks
`cond2`'s fail-open one. But the two conditions should not rely on each other's default
to stay safe; a future edit to `cond3` (e.g. changing its empty-`near` branch, or adding
an arm/dataset filter to `near` upstream of both conditions) could silently make an
empty-evidence dataset pass condition 2. This is worth a one-line fix
(`cond2 <- FALSE` as the initial value) independent of anything else in this audit.

### 6. Way a flip could pass on partial data

**GAP, not literally a pre-registration deviation, but directly on point** —
`08_apply_decision_rule.R` (whole file) and `12_build_results_doc.R:189`.

Neither script checks `n_masks` against the pre-registered count (3 for PanTHERIA and
AVONET, 1 for FishBase) before treating a dataset's pooled row as evidence.
`12_build_results_doc.R` computes and reports `n_masks = length(unique(g$seed))`
(line 189) purely descriptively; nothing downstream, including
`08_apply_decision_rule.R`, reads or gates on that column. If, say, one of PanTHERIA's
three structured-arm masks failed or was never produced, `results_table.csv` would
silently carry `n_masks = 2` for the affected trait/stratum rows (narrower evidence, no
explicit flag), and `08_apply_decision_rule.R` would treat those rows identically to a
complete n_masks = 3 row — nothing stops the decision rule from computing a verdict off
partial mask coverage. Missing FishBase specifically is handled correctly and is not a
gap: Amendment 1 explicitly authorizes evaluating condition 1 "on PanTHERIA, and on
FishBase if it runs," and conditions 2/3 "use every arm that exists" — a dataset or arm
with zero rows is correctly and intentionally excluded rather than counted as a failure.
The gap is specifically the absence of a **completeness gate** (e.g., asserting
`n_masks == 3` for PanTHERIA/AVONET rows, or `n_masks == 1` for FishBase rows, before
`apply_decision_rule()` runs) — combined with item 1/3's global pooling across datasets,
a dataset that completed fewer masks than its registered count, or fewer traits than its
siblings, has no floor stopping it from moving a globally pooled median.

### Other checks (no deviation found)

- Amendment 2's 5% real-missingness eligibility filter: `12_build_results_doc.R:41`
  (`COND1_MIN_REAL_MISSING <- 0.05`) and `:199`
  (`cond1_eligible = ... real_missing_frac >= COND1_MIN_REAL_MISSING`), consumed correctly
  by `08_apply_decision_rule.R:77-79`. `real_missing_frac` is computed from
  `mask_receipt.rds$truth` (`is.na` before masking), matching the amendment's stated
  definition. PASS.
- REGISTRY (`12_build_results_doc.R:47-53`) matches the pre-registration/Amendment 1
  mask-and-arm plan exactly: PanTHERIA mcar+structured x 3 seeds, AVONET mcar-only x 3
  seeds, FishBase mcar+structured x 1 seed. PASS.
- Between-mask SD: `NA` when fewer than 2 masks, `sd()` over up to 3 otherwise
  (`12_build_results_doc.R`, `sd_or_na()`), matching the results.md header's stated
  convention. PASS.
- Paired MCSE formula (`sqrt(mcse_mondrian^2 + mcse_split^2)`, each term the per-method
  analytic formula applied to pooled coverage/n): confirmed against raw receipts in the
  M2 re-derivation above. PASS.

## VERDICT

M2 (PanTHERIA structured, gestation_d, far) agrees with `results_table.csv` to 4 decimal
places on every reported quantity; no discrepancy.

Rule audit: **3 deviations** flagged —
(1) condition 1's median-gain pools traits across all datasets instead of per dataset,
(2) condition 3's median width ratio pools across all datasets instead of per dataset,
(3) condition 2's zero-evidence default (`cond2 <- TRUE`) is fail-open, inconsistent with
condition 3's fail-closed default in the same situation (currently masked by cond3, not
independently safe) —
plus one non-deviation **gap** (no `n_masks` completeness gate before the decision rule
runs, which compounds directly with deviations 1 and 2: an incompletely-run dataset has
no floor stopping it from moving a globally pooled statistic).
