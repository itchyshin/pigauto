# Arc C decision map: the AVONET300 continuous-trait gap to raw Rphylopars

2026-09-19, Claude Code, branch `arc/avonet-liability-gap` off main ebbf63e. Started on Shinichi's "merge and start the AVONET liability gap arc".

## Destination

pigauto's default baseline, with the GNN off and the mixed-type path intact, matches or beats raw Rphylopars on AVONET300's four continuous traits within Monte Carlo error over 20 seeds, without regressing on the simulated DGPs, and the with/without-GNN campaign tables are re-stated for that default.

## Decisions so far

- **The gap is the in-house single-pass joint solver, not the liability step.** Recorded in `docs/dev-log/2026-08-16-continuous-gap-diagnosis.md` (one mask, four-way decomposition: mixed-type path beneficial, calibration cost modest, in-house `fit_mvn_bm_inhouse()` closed-form Sigma owns 0.14 to 1.27 z-RMSE per trait) and confirmed on 5 masks by PR #169 (`joint_solver = "rphylopars"`). That note lives on the old handover branch, not on main's docs; this arc brings its conclusion forward. The campaign results note (2026-09-19) said "all four traits go through the threshold-joint fit, which is why"; that inference was too strong and is corrected here: the 08-16 decomposition and this arc's `cont_only_pure` arm both show the mixed-type path is not the cause.
- **`sigma_method = "fisher_ml"` does not help on AVONET300** (one seed so far): it is bit-identical to single-pass, meaning its `optim()` non-convergence fallback fired. To be confirmed over 20 seeds.
- **`joint_solver = "rphylopars"` closes the gap and, with the mixed-type path, beats raw Rphylopars** (one seed: 0.721 vs 0.751 mean z-RMSE; Tarsus 0.59 vs 0.64; Wing 0.27 vs 0.38). To be confirmed over 20 seeds and on the simulated DGPs.
- `lambda_mode = "bayes"` alone moves the in-house solver by about 1%; `pheno_correlated` makes no difference to raw Rphylopars on this data.

## Not yet specified

- Whether to change the default `joint_solver` to `"rphylopars"`. Costs: Rphylopars is `Suggests`, so a default cannot depend on it without a fallback; the OVR categorical path makes K phylopars calls (120 s vs 1.4 s per fit at n = 300); phylopars emits `solve(): system is singular` warnings on some AVONET masks. Benefit: about 30% better continuous-trait accuracy on real data, and the campaign's BACE and Rphylopars comparisons become fair to pigauto's mixed-type design. Usability is the principle that does not bend (D-139): a 100x slower default with a soft dependency needs Shinichi's word.
- Whether to repair the in-house solver instead (a converged REML or ML Sigma with cross-trait EM; the 2026-05-17 divergence that disabled EM refinement would need a root cause). That is a methods slice with its own recovery study, not this arc.
- Whether the campaign's GNN-on arms should be re-run with the better solver (the tax finding is independent of the solver, but absolute numbers change).

## Out of scope

- Multi-observation, multi-proportion, zi_count, covariates, MNAR: unchanged from the campaign's exclusions.
- BACE tuning.
- pkgdown (pre-existing red on main since 2026-08-28).

## Run in flight

`script/campaign_solver_cell.R` on Totoro: 200 cells (avonet x 20 seeds; bm_mixed, ou_mixed, bace_dgp x n {100, 300, 1000} x 20 seeds), same masks as the campaign, 11 arms, GNN off only. Results land in `~/gnn-off/results_solver/`; aggregate and write-up follow in this note's companion `2026-09-19-avonet-gap-results.md`.
