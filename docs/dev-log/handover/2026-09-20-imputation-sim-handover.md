# Handover: the four-arm imputation simulation, mid-campaign

2026-09-20, from the Claude Code lane that built and launched it. Branch `arc/imputation-sim` off
`origin/main` fd0b513, 23 commits, worktree
`/Users/z3437171/Dropbox/Github Local/pigauto-imputation-sim`.

**Read `.unlazy/imputation-sim/PROGRESS.md` first.** It is the live record: goal, completion
criteria, what is done, every output location, the compute rules, the findings not to re-derive, the
pending submissions, and the overnight authority Shinichi granted. This file is the orientation; that
file is the state.

## Where it stands

The campaign is **running**, not finished.

| stream | machine | state at handover |
|---|---|---|
| core slice, fast arms | Totoro | 2,922 of 3,600 replicate-jobs; n = 100 and n = 300 complete, n = 1000 filling |
| core BACE, n = 100 and 300 | nibi | running |
| core BACE, n = 1000 | fir | running, 2 h 46 m per replicate |
| factorial BACE, n = 1000 half A | fir | running |
| factorial BACE, n = 1000 half B | nibi | running, two replicates per task |
| factorial BACE, n = 100 | nibi | running |
| factorial fast arms, AVONET300, covariate sensitivity | not yet submitted | wait for Totoro to free |

rorqual is **retired** from this campaign. Its project allocation is at its file-count quota, so
result writes fail, and its nodes are slow enough that a replicate which takes 2 h 46 m on fir hit a
3 h 30 m wall there. Nine salvaged cells remain in its results directory and must be included when
pooling.

## What is already delivered

- Pre-run note with measured walls and a re-derived budget: `docs/dev-log/arc/2026-09-20-simulation-prerun.md`.
- Methods write-up for the BACE paper: `docs/dev-log/arc/2026-09-20-simulation-methods-bace.md`.
- The four-arm article, drafted but deliberately unlisted and unpublished:
  `vignettes/articles/simulation-study.Rmd`.
- The private results board: https://claude.ai/artifact/FkA8scunNMgVaH2791fFfx, currently marked
  partial. Republish the same scratchpad file path to keep that URL.
- Lane after-task report: `docs/dev-log/after-task/2026-09-20-imputation-sim-lane.md`.

## The one command that refreshes everything

```
bash script/campaign_sim_pool.sh core /tmp/pig_pool2
Rscript script/campaign_gnn_off_aggregate.R --dir /tmp/pig_pool2/core --out /tmp/pig_pool2/agg --reference floor
bash script/campaign_sim_page_build.sh /tmp/pig_pool2/agg <scratchpad>/sim-results.html "<status line>"
```

then republish that html to the artifact URL. Run the same three for `factorial`.

## Traps that already cost time

- **Pool per host, never flatten.** Totoro and the clusters write the same filename for the same
  cell, carrying different arms. `rsync --ignore-existing` silently discarded 1,034 Bayesian cells.
- **Never wrap MCMCglmm in `mclapply`.** It segfaults. Parallelism is one process per cell.
- **BACE wants one core and 32 GB.** Measured CPU efficiency was 24.9% of four cores, and 200 tasks
  were killed at 16 GB.
- **BACE's `runs` are imputation iterations, not chains.** Gelman-Rubin across them is meaningless;
  its own `assess_convergence()` needs at least three and the campaign uses five.
- **DRAC's 1,000-job cap counts array tasks**, so bundle replicates per task rather than submitting
  one task per replicate.
- Connect only through the existing ControlMaster sockets, or Duo fires.

## What remains

Finish the compute; pool across machines; aggregate; run gates G10, G11, G14 against the pooled tree;
run AVONET300 and the covariate sensitivity on Totoro once it frees; refresh the board; fill the
results sections of the methods note and the article; Melissa's plan-versus-actual reconcile; push
the branch and open a draft PR.

## Two gates that need Shinichi

- **Publication.** The article stays unlisted and unpublished until he has read the board. He has
  already approved everything short of publishing and merging.
- **Merge.** No agent merges a PR here.

## Residual unknowns

G6, the coverage gate for arms 1 and 2, has not returned a number. Whether `runs = 5` converges at
n = 1000 was never measured; only n = 100 was probed. The half-B task shape, two replicates per task,
has not been timed end to end.
