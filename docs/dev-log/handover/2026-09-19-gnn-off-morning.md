# Morning report: `gnn = FALSE` (arc A done; arc B pre-run done, waiting for you)

Written 2026-09-18 late evening, Claude Code session, for Shinichi to read first thing.

## Two decisions for you

1. **Merge PR #180?** https://github.com/itchyshin/pigauto/pull/180 (draft). `feat/gnn-off` off `origin/main` 7f2aff3 (v0.11.0), 11 commits, last 31c5358. Full suite green, `R CMD check` E0 W0 N0, ledger ALL MET, Fable adversarial pass plus a three-reviewer completion panel all DONE. Nothing in it changes GNN-on predictions (bit-identical to main on CPU). I did not merge; that is yours.
2. **Submit the campaign?** The pre-run says the locked design (7 arms, 3 DGPs, n 100/300/1000, 20 seeds, plus AVONET300) costs about **2 to 3 h wall on Totoro at 36 concurrent cells**, dominated by BACE (23 min per n = 1000 cell). Two cheaper variants: BACE only at n <= 300 (about 1 h), or 10 BACE seeds at n = 1000. Numbers and the command are in `docs/dev-log/arc/2026-09-18-campaign-prerun.md`. Say "full", "BACE <= 300", or "10 seeds" and it runs; it is over the 30-minute line so it waits for you (D-139).

## What exists now

- `impute(traits, tree, gnn = FALSE)`: the phylogenetic baseline through the whole pipeline with **no torch call at all** (no device probe, no tensors, no GPU): completed data, conformal intervals, MI datasets, `multi_impute()`, `multi_impute_trees()`, `cross_validate()`, `simulate_benchmark()`, `evaluate()`, `summary()`, `print()`, `pigauto_report()`, `save_pigauto()`/`load_pigauto()`.
- Speed (AVONET300, n = 300): GNN-off 8 s, pure GNN-off 1.4 s, GNN-on 114 s at 300 epochs. On Totoro at n = 1000: GNN-off 136 s (that is the gate calibration, not the baseline), pure GNN-off 2 s, GNN-on 523 s, BACE 1,381 s, Rphylopars 2 s.
- Two arms are available with existing flags: the default `gnn = FALSE` keeps the GNN arm's safety floor and phylo-signal gate (so with-vs-without differs only in the GNN term); `gnn = FALSE, safety_floor = FALSE, phylo_signal_gate = FALSE` is the pure traditional-statistics arm for the Rphylopars / BACE comparison.
- `fit$baseline_full`: a second baseline fit with no held-out cells, used only for production predictions; every scorer stays on the held-out `fit$baseline`. This removes the held-out-cell cost that the 2026-08-16 note suspected. `fit_baseline()` now records which dispatch produced each trait (`path`).
- Docs: NEWS entry, AGENTS.md "GNN-off arm", roxygen on every changed export. Tests: `tests/testthat/test-gnn-off.R`, 15 tests. Reports: `docs/dev-log/after-task/2026-09-18-gnn-off.md` (12 sections), `docs/dev-log/plan-actual/2026-09-18-gnn-off.md` (Melissa), `docs/dev-log/arc/2026-09-18-gnn-off-verify.md` (AVONET smoke).

## What the pre-run hints (2 seeds, BM-correct mixed DGP; the campaign decides)

- GNN-off and raw Rphylopars are within 0.005 z-RMSE of each other at n = 100 and n = 1000, as designed.
- GNN-on is a little worse than GNN-off on continuous traits (0.320 vs 0.296 at n = 100; 0.099 vs 0.090 at n = 1000); predicting the same GNN-on fit from the tax-free baseline closes about half the gap at n = 100 and all of it at n = 1000. That is the tax mechanism, now measurable.
- BACE at 50k iterations trails every pigauto arm on continuous traits and matches on discrete ones.
- Coverage of the GNN-off production interval: 0.88 at n = 100, 0.95 at n = 1000.

## Residuals I want you to know about

- Conformal intervals centred on `baseline_full` have no formal split-conformal guarantee (scores come from the held-out fit); conservative under BM-correct only. The campaign measures coverage on user-masked cells.
- Not exercised under `gnn = FALSE`: multi-observation data, multi_proportion and zi_count traits, `conformal_method = "bootstrap"`, `simulate_benchmark()` / `compare_methods()` through the switch. Listed in the report's section 12.
- The default GNN-off arm spends its time in `calibrate_gates()` (cv folds over a simplex grid); at n = 1000 that is 136 s of a 138 s fit. A cheaper calibration for the GNN-off case would be a later slice, not this one.

## Housekeeping

- Worktree: `/Users/z3437171/Dropbox/Github Local/pigauto-gnn-off`. Your main checkout on `handover/2026-08-09-cursor` is 129 commits behind `origin/main`; nothing there was touched.
- Totoro: branch installed in a separate library `~/R/lib-gnn-off` (the #175 lane's pigauto install is untouched); runner and results under `~/gnn-off/`.
- Lane lease `claude:pigauto-gnn-off:94331` released at the end of the session.
- Brain: nothing written to the vault. A proposed AGENT_LOG entry and a proposed decision record are staged at `docs/dev-log/arc/2026-09-18-brain-draft.md` for your approval.

Next action for you: reply "merge" and one of "full" / "BACE <= 300" / "10 seeds".
