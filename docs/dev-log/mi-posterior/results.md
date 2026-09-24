# Posterior multiple imputation: results

2026-09-24. Branch `arc/mi-posterior`, stacked on PR #187. Status: **simulation complete; G6 fails on 8
checks (diagnosis running); real data 9 of 10 cells done (FishBase running); not ready to merge.**

Every number below comes from a committed file:
- `sim_summary.csv` and `cell_coverage.csv` (written by `script/mi_gls/03_summarise_v2.R`);
- `results_tables.md` (written by `script/mi_gls/06_results_tables.R` from those two files);
- `evidence/` for G3, the smoke run and the calibration check.

Campaign code: commit 69670d4, a frozen `git archive` on Totoro. Later commits change only messages,
provenance markers, docs and tests. Proof that the sampler is unchanged:
- `.mip_fit()` output is bitwise identical to 69670d4 on a test fixture;
- `review.md` records the function-by-function comparison.

## What was built

`multi_impute(traits, tree, draws_method = "posterior")` for continuous traits. It draws the missing
cells jointly from vec(Y) ~ N(1 mu', Sigma_P %x% R + Sigma_E %x% I_n), with full Sigma_P and
Sigma_E and R = cov2cor(vcv(tree)).
- Sigma_P, Sigma_E (and so each trait's lambda) and mu are redrawn by Gibbs data augmentation, with
  collapsed Metropolis moves, so the m completed datasets carry parameter uncertainty.
- Per-cell 95% predictive intervals come from 1,000 kept posterior-predictive draws.
- The GNN is not used.
- Design: `design.md`. Review: `review.md`.

## Simulation (24 regimes x 200 replicates; `results_tables.md`)

Design:
- two traits, n = 300 or 1000, 30% missing (x only, or both traits);
- MCAR or clade-biased (MAR_phylo) masks;
- regimes 1 to 16: one lambda (1 or 0.5) with trait correlation 0.7;
- regimes 17 to 24: two lambdas (0.3/0.9 or 0.7/0.7) with different phylogenetic and residual
  correlations.

Downstream analyses are `nlme::gls(corBrownian)` and `phylolm(model = "lambda")`, compared with the
complete data in the same replicate.

Convergence: 4,756 of 4,800 fits converged (99.1%).

### Headline 1: per-cell predictive coverage is about 95%

| Mask | Posterior 95% interval coverage | Conformal (`impute(gnn = FALSE)`) coverage | Posterior width / conformal width |
|---|---|---|---|
| MCAR (20 regime x trait rows) | 0.933 to 0.957 | 0.954 to 0.970 | 0.57 to 0.92 |
| Clade-biased (20 rows, reported only) | 0.932 to 0.949 | 0.936 to 0.959 | 0.60 to 0.93 |

Every MCAR row is inside the G7 band [0.92, 0.98]. The posterior intervals are 8% to 43% narrower
than conformal and still cover about 95%.

### Headline 2: the downstream slope

**Paired bias (MI pooled slope minus complete-data slope, same replicate):**
- gls: -0.013 to +0.017;
- phylolm: -0.026 to +0.001;
- for comparison, conformal MI draws gave -0.20 to -0.46 in regimes 1 to 16 of the earlier sweep
  (same data-generating seeds; `arc/mi-gls-attenuation`, `docs/dev-log/mi-gls/results.md`).

**Coverage of the pooled 95% interval:**
- 0.77 to 0.975, against 0.725 to 0.97 for complete data under the same model;
- within 0.05 of complete data in 47 of 48 regime x analysis rows;
- conformal MI draws covered 0 to 17%.

**SE ratio under phylolm, relative to complete data:** 0.88 to 1.12, inside [0.90, 1.15] in 22 of 24
regimes. The Monte Carlo SE of this ratio is about 0.07.

**Proper vs plug-in** (reported, not gated). Holding Sigma_P and Sigma_E at their posterior means:
- lowers the SE ratio in all 32 both-missing pairs;
- lowers coverage (plug-in 0.77 to 0.94, against proper 0.875 to 0.975).

So propagating parameter uncertainty matters downstream.

### Gate verdicts (pre-registered; `.unlazy/mi-posterior/GATES.md`; rules revised only at CP1, before any result)

| Gate | Result | Detail |
|---|---|---|
| G1 unit tests | met | 29 tests, 200 expectations in the two posterior test files, 0 failures |
| G2 exactness | met | `EXACTNESS_OK` (fixed-parameter draws match the dense conditional) |
| G3 calibration | met | `RECOVERY_OK`: 200 fits, 95% interval coverage of lambda and rho_P 0.92 to 0.98 |
| G4 smoke convergence | met | `CONVERGENCE_OK` |
| G5a suite / G5b solver / G5c check | met | 2,840 pass, 0 fail; `R/joint_mvn_solver.R` untouched; R CMD check 0 errors 0 warnings |
| **G6 simulation acceptance** | **not met** | 8 failing checks, listed below |
| **G7 per-cell coverage** | **not met** | Every coverage value passes; the only failures are the shared non-convergence rule in regimes 5, 21, 23 |
| G8 real data | pending | 9 of 10 cells done; FishBase running |

Ledger run (`gate-check.mjs`, worktree root, 1 h timeout, 2026-09-24):
- G1 to G5c: PASS (9 met, counting M1);
- G6, G7: FAIL;
- M2: pending.

That run also recorded a false G8 PASS. `03_acceptance.R` exited 0 on failure and printed its token
inside a header line; gate-check matches the token as a substring and requires exit 0. This was fixed
in 041c2e0 (exit 1, token only on the pass line) and G8 was reset to unmet. Every other gate script
was checked for the same two faults; none had them.

The 8 G6 failures:
1. Non-converged fits above 2% in three n = 300 regimes: regime 5 (2.5%), 21 (5.5%), 23 (8.0%).
2. Paired bias beyond max(0.02, 2.5 MCSE) under phylolm in regime 1 (-0.026) and regime 9 (-0.021).
   Both have lambda = 1 and n = 300.
3. Coverage shortfall 0.0515 (limit 0.05) under phylolm in regime 3 (lambda = 1, n = 1000, x only):
   0.854 against complete 0.905.
4. Relative SE ratio just below 0.90 under phylolm in regime 8 (0.879) and regime 15 (0.897). Both
   are within about 1.5 Monte Carlo SE of the floor.

Pattern: every lambda = 1 regime has a small negative paired bias (phylolm -0.008 to -0.026; gls
-0.002 to -0.013). The plug-in arm shows the same bias, so the cause is not parameter uncertainty. The
cause is not yet known. A diagnosis is running (`diagnosis.md` when it lands):
- an oracle arm with the true covariances;
- a prior-sensitivity arm;
- longer chains for the non-converged fits.

## Real data (preview, 9 of 10 cells; G8 pending)

Data: PanTHERIA (4,027 species, 4 continuous traits; 3 MCAR and 3 clade-structured masks) and AVONET
(1,500 species, 4 traits; 3 MCAR masks). Source: `script/mi_realdata/`, receipts on Totoro.
- **Convergence:** 7 of 9 cells converged. Two PanTHERIA cells have min ESS 368 and 366 (R-hat at
  most 1.01).
- **Per-trait coverage:** model-based 0.87 to 0.96 across 36 trait-cells, against split conformal
  0.90 to 0.98 and Mondrian 0.91 to 0.99. So on real data the model-based intervals cover somewhat
  less than conformal.
- **Slopes:**
  - AVONET: pooled slopes within about 3% of the complete-row reference.
  - PanTHERIA: some cells differ by up to 3 reference SEs (for example longevity ~ body mass: 0.148
    against 0.180).
  - Likely reason: the PanTHERIA columns are already log values and pigauto's default
    `log_transform = TRUE` logs them again. The imputation model is then linear on a log(log) scale
    while the pre-registered analysis is on the log scale, which is outside the stated congeniality
    scope.
  - A `log_transform = FALSE` sensitivity run is the natural follow-up.
- These numbers are a preview and are not yet the G8 report.

## What this does not cover

- Discrete, mixed and multi-observation traits.
- GNN blending.
- Covariates in the imputation model.
- Analyses with nonlinear terms or external covariates.
- More than two traits in simulation (real data has 4 or 5).
- Trees with more than about 10,000 tips. FishBase, at 10,484 tips and 5 traits, takes about 1 s
  per sweep, so a default fit takes hours.

The default draws method is unchanged (Shinichi's decision).
