## Mondrian conformal: real-data confirmation, default decision, paper section, MI draw scale

Keeps the default `conformal_method = "split"` after a pre-registered real-data test of
`"mondrian"`. Also records stratum sizes in Mondrian fits, gives Mondrian fits per-cell
MI draw scales, and adds the paper's uncertainty section.

### What was tested

Originally observed cells were masked, so there was truth to score against, in three
databases under a pre-registered design (`docs/dev-log/mondrian-realdata/00-preregistration.md`,
with two amendments written before the outcomes they govern).

| Database | n | Mask arms | Masks |
|---|---:|---|---:|
| PanTHERIA | 4,027 | random and phylogenetically structured | 3 each |
| AVONET | 1,500 | random only (almost no real missingness) | 3 |
| FishBase | 10,484 | structured | 1 (descriptive) |

### Result

Pooled over traits, Mondrian raised far-stratum coverage from about 0.92-0.93 to
0.94-0.96. It narrowed near-stratum intervals by about a fifth, and near coverage fell
from 0.97-0.98 to 0.95-0.97. The pre-registered rule keeps split. Conditions 1 (far
gain) and 3 (near width) passed on every dataset. Condition 2 (near non-inferiority
within 2 points) failed for AVONET, and FishBase could not show it from one mask.
Per-trait tables are in `docs/dev-log/mondrian-realdata/results.md`.

### Also in this PR

- `compute_conformal_scores()` stores `n_val`, `n_near` and `n_far`, and a verbose
  fallback message names the realised stratum sizes.
- `mondrian_cell_scores()` is extracted from `predict()` with its output unchanged.
  `multi_impute(draws_method = "conformal")` uses per-cell Mondrian scores for Mondrian
  fits.
- Paper section 8 with verified references. One citation is corrected: Bostrom and
  Johansson 2020.
- MI finding: conformal draws, split or Mondrian, halved a phylogenetic GLS slope in
  simulation, while OLS on the same draws was unbiased. This is pre-existing on the
  default split path. NEWS carries a caveat, and the investigation is on
  `arc/mi-gls-attenuation`.

### Checks

- `devtools::test()`: FAIL 0, PASS 2509, SKIP 8, on 9f8f2b8. The mondrian filter was
  re-run after the final code change.
- Every gate in the acceptance ledger was re-verified with `gate-check --reverify`.
- Reviews are in `docs/dev-log/review/`:
  - method audit: hand re-derivation matched to 4 decimals; three rule-code deviations
    fixed, verdict unchanged;
  - traceability: 0 mismatches;
  - claim gate: 1 blocking and 14 required items, all addressed.
- Compute ran on Totoro, DRAC fir and kohaku. No GitHub Actions compute was used.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
