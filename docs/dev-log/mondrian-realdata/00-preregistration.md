# Mondrian real-data confirmation: pre-registration

Written 2026-09-23, before any 0.11.0 result exists. Branch `arc/mondrian-realdata`.
Nothing below may change after the first campaign receipt is written, except by an
appended, dated amendment that says why.

## Question

On real trait databases, does `conformal_method = "mondrian"` improve interval coverage
where missing cells actually sit, without degrading the cells where split conformal is
already valid?

## Why the design has two mask arms

A random mask of observed cells makes test cells exchangeable with the validation cells,
so split conformal is already valid there and cannot fail in the way Mondrian repairs.
That arm is kept only as a no-harm control. The structured arm places test cells where
genuinely missing cells sit.

- **MCAR arm.** Per trait, 20% of observed cells masked uniformly at random.
- **Structured arm.** Per trait, fit a logistic regression of "cell is missing in the
  real data" on the first k phylogenetic eigenvectors (k as chosen by
  `build_phylo_graph()`), then mask 20% of observed cells with probability proportional
  to the fitted propensity. Test cells then resemble the real missing cells in tree
  position.

Masks per arm: 3 for PanTHERIA (n = 4,027) and AVONET full (n = 1,500); 1 for FishBase
(n = 10,484) if Shinichi approves the GPU campaign after its pre-run. Seeds 20260818,
20260819, 20260820. Split and Mondrian are fitted on the identical masked data with the
identical seed. pigauto 0.11.0 plus the S2 instrumentation, 500 epochs, one imputation.

## Per trait, per stratum, per mask, per method

Coverage, mean and median half-width, Winkler interval score, n_test, and the stratum
of each test cell (near or far, using the validation median locality). From the fit:
n_val, n_near, n_far, fallback. Test cells in a trait that fell back carry no Mondrian
evidence.

## Uncertainty

Per trait and stratum, the Monte Carlo SE combines the binomial term on n_test with the
conditional-coverage SD of the calibration order statistic, approximately
sqrt(alpha (1 - alpha) / (n_s + 2)) with n_s the stratum's validation count. The
between-mask SD is reported beside it (2 df). FishBase is descriptive only.

## Decision rule for the default

Flip the default to `"mondrian"` only if all of the following hold on every dataset
where Mondrian activates for at least one trait, with Holm adjustment across datasets
for the one-sided tests:

1. Structured arm, far stratum: the median across traits of the paired coverage gain
   (Mondrian minus split, same mask, same cells) is at least 0, and no activated trait
   has Mondrian far-stratum coverage below 0.90.
2. Structured and MCAR arms, near stratum: Mondrian coverage is at least split
   coverage minus 2 percentage points (one-sided non-inferiority).
3. Near stratum: the paired ratio of median half-widths, Mondrian over split, is at
   most 1.10. Far-stratum widening is the intended effect and is reported, not capped.
   (Shinichi chose this reading on 2026-09-23, before any result.)

Traits that fall back count as no evidence, never as a pass. If the rule fails, the
default stays `"split"` and the fallback message names the realised stratum sizes.

## What this study cannot show

It cannot establish an unconditional 95% guarantee. Coverage on masked observed cells
is a proxy for coverage on the cells users impute, and the structured arm is a better
proxy than the random arm, not the real thing.

## Amendment 1 (2026-09-23, before any campaign receipt)

AVONET is almost fully observed (at most 20 missing cells per trait in the 9,993-species
bundle), so the structured arm has no real missingness pattern to fit a propensity to.
AVONET therefore runs the MCAR arm only, as a no-harm control, on the same 1,500-species
subset the committed bench used (`set.seed(2026); sample(rownames, 1500)`). The
decision rule's condition 1 (structured arm, far stratum) is evaluated on PanTHERIA, and
on FishBase if it runs; condition 2 (near-stratum non-inferiority) and condition 3 (near
width ratio) use every arm that exists. Campaign receipts: PanTHERIA 2 arms x 3 masks x
2 methods = 12; AVONET 1 arm x 3 masks x 2 methods = 6; total 18 before FishBase.

## Amendment 2 (2026-09-23, before any structured-arm receipt was read)

The structured arm needs real missingness to fit a propensity to. A trait whose real
missing fraction is below 5% gives a propensity that is flat or degenerate, so its
structured mask is effectively random. Decision-rule condition 1 (structured arm, far
stratum) is therefore evaluated only on traits with at least 5% real missingness. Those
traits are all PanTHERIA traits (13.8% to 74.9%) and FishBase Weight, DepthRangeDeep and
Troph (54% to 82%). FishBase Length (2.4%), BodyShapeI (0.4%) and Vulnerability (0%)
are reported, but in the MCAR-like condition 2 and 3 tables only. The trigger was
glm non-convergence and 0/1 fitted probabilities in the FishBase propensity fit, seen in
the run log; the per-trait propensity quantiles are recorded in each mask receipt.
