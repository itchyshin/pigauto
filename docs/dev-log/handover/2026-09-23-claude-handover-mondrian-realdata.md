# Session Handoff (to Claude): Mondrian real-data confirmation closed; MI under GLS in flight

Meta: 2026-09-23, from Claude (Opus 5.5 session; Fable 5.1 ceiling reviewer). Two lanes,
two branches, two worktrees. Platform-neutral; any tool resumes from this file.

## Critical Context

1. **The Mondrian default decision is made: keep `"split"`.** The pre-registered rule
   (`docs/dev-log/mondrian-realdata/00-preregistration.md`) failed on condition 2 only
   (near-stratum non-inferiority within 2 points): AVONET failed it, FishBase could not
   show it from one mask, PanTHERIA passed. Conditions 1 and 3 passed on every dataset.
   NEWS says so. Do not flip the default without a new pre-registered rule and
   Shinichi's decision.
2. **The failing condition penalises removal of over-coverage.** Pooled over traits,
   Mondrian lifted far-stratum coverage from 0.92-0.93 to 0.94-0.96 and brought
   near-stratum over-coverage (0.97-0.98) down to 0.95-0.97. Whether to write a new rule
   (for example "near coverage at or above nominal minus a margin") is Shinichi's call;
   it would need a new confirmation, not a re-reading of these data.
3. **A separate, pre-existing MI problem was found and is its own lane.**
   `multi_impute(draws_method = "conformal")` halves a phylogenetic GLS slope in
   simulation (split and Mondrian alike); OLS is fine; an oracle proper imputation is
   unbiased. Default draws method unchanged; NEWS carries a caveat.

## What Was Accomplished

- `arc/mondrian-realdata`: two-arm real-data campaign (PanTHERIA 12 receipts, AVONET 6,
  FishBase 2), results table regenerated from receipts (ROWS_MATCH), decision rule
  (per dataset, fail-closed, mask-completeness gate), NEWS, paper section 8 with
  Table S-UQ2 and verified references, MI memo with a 500-rep simulation and two
  diagnostics, stratum sizes in Mondrian fits and in the fallback message.
- Reviews in `docs/dev-log/review/` and `docs/dev-log/mondrian-realdata/`: method audit
  (hand re-derivation matched to 4 decimals), traceability (0 mismatches), claim gate
  (Fable; 1 blocking, 14 required, all addressed).
- kohaku now has user-space R 4.5.3 with torch on CUDA 12.8 in
  `~/micromamba/envs/r45` (record: `docs/dev-log/mondrian-realdata/kohaku-r-install.md`).

## Current Working State

- Mondrian lane: complete; PR as recorded in the Landing State table.
- MI-GLS lane (`arc/mi-gls-attenuation`, worktree `pigauto-mi-gls`): prototype
  `draw_conditional_bm()` in `R/draws_conditional.R` (39 tests pass), 16-regime sweep
  in `script/mi_gls/`. Sweep running on DRAC fir as array `61137481` (1,920 cells;
  240 done at 13:10 MDT); results land in `/scratch/snakagaw/pigauto-mi-gls/results/`
  (scratch is purged after about 60 days: copy back promptly). Smoke: new
  conditional-BM draws gave PGLS slopes 0.643 and 0.625 against complete-data 0.649 and
  0.690; conformal draws 0.33-0.49.

## Landing State

| Artifact / branch | Committed | Pushed | PR | State |
|---|---|---|---|---|
| `arc/mondrian-realdata` | y | y | draft PR https://github.com/itchyshin/pigauto/pull/188 | LANDED as draft; not merged; merge is Shinichi's call |
| `arc/mi-gls-attenuation` (67a8df8) | y | n | none | CARRIED-OVER: sweep running on fir |
| Checkout `handover/2026-08-09-cursor` (18 dirty files) | n | n | none | PROTECTED, another lane's |

## Next Immediate Steps

1. MI-GLS: when array `61137481` finishes, `rsync` the results into
   `script/mi_gls/returned/`, run `Rscript script/mi_gls/02_summarise.R`, and write
   `docs/dev-log/mi-gls/results.md` with regime tables (bias, SE ratio, coverage, MCSE,
   OLS and PGLS). Check `sacct -j 61137481` for OUT_OF_MEMORY or TIMEOUT first.
2. MI-GLS: from those tables, recommend whether `draw_conditional_bm()` should become a
   `draws_method` option (a `multi_impute()` change; the lambda lane in another Claude
   account owns the baseline files, so read their branch `feat/joint-lambda-default`
   before touching the baseline). Default changes need Shinichi.
3. Mondrian: Shinichi decides whether to pre-register a new rule; nothing else owed.

## Gotchas

- fir: `/project` quota is full, use `/scratch`; cu126 libtorch needs
  `module load cuda/12.6` and `LD_LIBRARY_PATH=$EBROOTCUDA/targets/x86_64-linux/lib`
  even on CPU nodes; build packages on `$SLURM_TMPDIR` (a Lustre glitch broke one
  install on scratch); size `--mem` from a measured peak (PanTHERIA n = 4,027 peaks at
  16.4 GB).
- Totoro was at its 150-core cap from another lane's `script/campaign_s...` job today.
- kohaku: use `OPENBLAS_NUM_THREADS=8` for dense phases at n about 10,000 (one thread ran
  for 2 h without reaching training); cap 8 of 32 vCPU; a neighbour holds about 50 GB VRAM.
- `R --vanilla` drops `~/R/lib` on Totoro; set `R_LIBS=<lib>:$HOME/R/lib`.

## How to Resume

```text
Read AGENTS.md and docs/dev-log/handover/2026-09-23-claude-handover-mondrian-realdata.md on branch arc/mondrian-realdata. Run the rehydration steps, reconcile with git and fir job 61137481, then continue only the OWED MI-GLS steps.
```
