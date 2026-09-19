# Brain draft (staged, not written to the vault; pigauto rule: propose, never write)

## Proposed `memory/AGENT_LOG.md` entry

- 2026-09-18 · Claude Code · pigauto · **`gnn = FALSE` full-pipeline arm** (ultra-plan; PR #180, draft). The phylogenetic baseline now runs through the whole pipeline with zero torch calls; `baseline_full` gives production predictions without the held-out-cell cost of gate calibration while every scorer stays on the held-out baseline; `fit_baseline()` records its dispatch `path`. Fable plan review changed two design defaults before code was written (keep the safety machinery in the GNN-off arm; two baselines rather than one refit), and a three-reviewer D-43 panel closed it. AVONET300: GNN-off 8 s vs GNN-on 114 s at 300 epochs. Campaign pre-run on Totoro (2 seeds, n 100 and 1000, seven arms incl. BACE at 50k iterations and raw Rphylopars): GNN-off matches raw Rphylopars within 0.005 z-RMSE under a BM-correct DGP; the GNN-on gap to GNN-off is mostly the tax. Full campaign (about 2 to 3 h Totoro) waits for approval. Reports: `docs/dev-log/after-task/2026-09-18-gnn-off.md`, `docs/dev-log/arc/2026-09-18-campaign-prerun.md`.

## Proposed decision record (number to be allocated by committing)

**D-xxx: pigauto `gnn = FALSE` semantics.** (a) The GNN-off arm keeps `safety_floor` and `phylo_signal_gate` exactly as the GNN arm, so with-vs-without comparisons differ only in the GNN term; the pure traditional-statistics arm is `gnn = FALSE, safety_floor = FALSE, phylo_signal_gate = FALSE`. (b) Two baselines: `fit$baseline` (held-out) feeds every scorer; `fit$baseline_full` (no held-out cells) is read only by production `predict()`. Routing a scorer to `baseline_full` is test-cell leakage. (c) The GNN-off path makes no torch call; the preflight device probe is skipped. Origin: Shinichi 2026-09-18 ("GNN off must be super fast, no GPU or ML"; "fair comparison with Rphylopars and BACE"); Rose (Fable) plan review supplied (a) and (b).

## Proposed `WHAT-WORKS` bullet

- A frozen slot contract, committed before fan-out, let three Sonnet builders edit one S3 object concurrently with no merge conflicts (pigauto gnn-off, 2026-09-18); the two mid-run contract corrections were cheap because they were explicit. Pair it with gate scripts that get their own adversarial pass: three of six gates were wrong before the code was.
