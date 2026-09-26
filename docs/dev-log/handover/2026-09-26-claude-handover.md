# Session Handoff: extend the freq-vs-BACE Rubin campaign to discrete traits

Meta: 2026-09-26 · from Claude Code (Opus 5.5) to a fresh Claude Code session · lane `claude:pigauto:rubin`
(released; re-claim it) · worktree `/Users/z3437171/Dropbox/Github Local/pigauto-rubin-freq-bace` · branch
`arc/rubin-freq-bace` (pushed, no PR)

You are Claude, picking up the pigauto Rubin lane. The continuous campaign is finished and published. Shinichi
wants the discrete traits (binary, ordinal, categorical) **added to the same campaign**: the same 18 cells, the
same seeds and therefore the same simulated datasets, so discrete and continuous results pair dataset by dataset.
Your job is to plan that extension and bring Shinichi a pre-run plan with a measured time estimate. It will
exceed the 3-hour line, so do not launch the campaign without his approval (D-287).

## Mission control

| Repo / worktree | Branch | State | What shipped | Next, by leverage |
|---|---|---|---|---|
| pigauto (worktree above) | `arc/rubin-freq-bace`, pushed, no PR | LANDED on the branch; not merged | Continuous campaign: 6,600 fits; report, summary and accuracy pages | 1. Discrete extension (this handover). 2. Shinichi decides whether the branch gets a PR. |

## Critical context

1. **All continuous results are kept, but not the imputed datasets.** Per-fit `rds` files hold the summary
   scores (`cells`, `estimands`), the true data (`truth`, which includes `bin`, `ord`, `cat3`) and the `mask`.
   They do not hold the 20 imputations, so the discrete traits cannot be scored from them. Copies:
   - cluster originals: nibi `~/projects/def-snakagaw/snakagaw/pigauto_rubin/results/`, fir
     `~/pigauto_rubin/results/` (BACE) and `~/pigauto_rubin_f3/results/freq` (frequentist of record),
     rorqual `~/projects/def-snakagaw/snakagaw/pigauto_rubin/results/`, Totoro `~/pigauto_rubin/results/`;
   - Mac pool `~/pigauto_rubin_pool` (556 MB, outside Dropbox; `agg/report.json` feeds every page). It is not
     backed up anywhere else on the Mac.
2. **The datasets can be rebuilt exactly.** `rubin_cell.R` builds each dataset with `make_cell()` from
   (n, λ, ρ, seed) under `RNGkind("L'Ecuyer-CMRG")`. Each arm has its own seed offset (`arm_seed()`: freq A 101,
   freq B 202, BACE 303, chained BACE 404). The `truth` and `mask` in each stored `rds` let you confirm a
   rebuilt dataset is identical before reusing anything.
3. **The frequentist continuous results can be reused unchanged.** Only the new discrete arms and BACE need
   running. BACE must be re-run because its imputations were not saved. With the same seed it should
   reproduce the stored continuous rows. Test that on a few datasets first: if it holds, the re-run both
   confirms and extends the continuous results; if it does not (cross-host BLAS, library versions), say so
   and treat the re-run as a replicate.
4. **The continuous campaign was the frequentist method's home ground**, and it did not test BACE's claimed
   advantage (one model for mixed trait types). Simulation v1 measured discrete per-value accuracy only (see
   Q1), without multiple imputation, Rubin pooling or a downstream estimate.
5. **Hold Dan.** Shinichi's standing instruction: do not contact Dan Noble. The pages are private.

## What was accomplished (continuous campaign, 2026-09-24 to 26)

- 18 cells (n 100/300/1000 × λ 0.3/0.7/1 × ρ 0/0.5), M = 20. Arms: freq A (proper: parametric bootstrap
  then joint conditional draw), freq B (improper: one fit), BACE as shipped, chained BACE, BACE + residual
  (negative control), and complete data. BACE 3,000 fits (100 per cell at n = 1000, 200 elsewhere),
  frequentist 3,600; 0 missing.
- Verdict: BACE as shipped under-covers at λ 0.3/0.7 (about 0.91 per value; downstream correlation 0.90 at
  n = 300 and 1000) because its 20 final runs start from one shared dataset. Chaining the final runs fixes
  it: chained BACE matches freq A on coverage and accuracy at n ≥ 300. At n = 100 it is less accurate and pulls
  the slope toward zero (about −0.02 against 0.001 for freq A).
- Pages (private, owned by Shinichi):
  - report https://claude.ai/artifact/PFkoRFtTjox4tndEbuBtPQ
  - summary https://claude.ai/artifact/9vNrmz2PUfjW9vZwuKmvFc
  - accuracy https://claude.ai/artifact/CGGmi4Km6unGPhxEiv94qk
  - BACE settings pre-run https://claude.ai/artifact/WMHxQ8idafD5sLAVEJtg9Q
- Full record, every fix and residual: `docs/dev-log/after-task/2026-09-25-rubin-campaign.md`.

## Prior discrete work to reuse (simulation v1, lane `arc/imputation-sim`, draft PR #184)

- Frequentist discrete arm: **castor Mk hidden-state prediction**, `run_freq()` in
  `script/campaign_gnn_off_lib.R`. It uses `castor::hsp_mk_model` with equal rates for binary and
  categorical and stepwise (SUEDE) for ordinal, and it keeps the tip likelihoods as class probabilities.
- v1 scored discrete per-value accuracy only (summary `docs/dev-log/arc/2026-09-23-simulation-v1-summary.md`,
  board https://claude.ai/artifact/M5HtGRnNGfwsK2Se4gMX24 v10):
  - at λ = 0.3, BACE led by 11 to 17 accuracy points (0.587 against 0.422 at n = 1000), and every other arm
    was at or below the mode floor;
  - at λ = 1, castor and pigauto led (0.96 to 0.97 against BACE's 0.85).
- v1's BACE ran before the chaining fix and the input cleaning, so its discrete numbers may move.
- Known failure: castor errors on a trait with only one observed class (common at λ = 1 with fixed
  thresholds); v1 scored those at the floor.
- AGENT-INFERRED, untested: castor sees one trait at a time. At low signal the tree carries little
  information, so BACE's use of the correlated continuous traits may be what gives it the lead.

## Current working state

- Working: all runners, tests and build scripts. The lane suite `testthat::test_dir("script/tests-rubin")`
  gives FAIL 0, PASS 234. All 16 lane gates are met (re-approved and re-run 2026-09-26).
- In progress: nothing. No jobs on any cluster.
- Blocked: nothing in this lane.

## Key decisions and rationale (still binding)

- Never edit `R/` or `BACE/` in this lane; BACE is used as installed. New code goes in `script/rubin_*` and
  `script/tests-rubin/`.
- BACE campaign settings: runs 5, nitt 50,000, burnin 10,000, thin 25 (Shinichi; BACE's own convergence check
  does not measure mixing).
- Coverage target is the population value. The correlation estimand is GLS-whitened on Fisher's z. PGLS is
  scored with the fast eigenbasis REML estimator (`est_pgls_slope_fast`), which is more exact than
  `nlme::gls` when λ is at a boundary.
- Compute rules:
  - D-287: estimate first; over 3 h needs a plan, a pre-run with its results, and approval.
  - D-143: Totoro ≤150 cores for this user, across all lanes.
  - D-64: attach only via the `~/.ssh/cm-*` sockets; never trigger Duo.
  - Never compute on a DRAC login node.
  - Proven hosts: fir, rorqual, nibi, Totoro. narval's CPUs cannot run the fir-built library.
- `ScheduleWakeup` is forbidden in this project (see `CLAUDE.md`); use Monitor or background Bash.

## Landing state

| Artifact / branch | Committed | Pushed | PR | State |
|---|---|---|---|---|
| pigauto `arc/rubin-freq-bace` (tip includes this handover) | y | y | none | LANDED on branch |
| `.unlazy/rubin-freq-bace/gates/*` (16 gates) | untracked by design | n/a | n/a | ALL MET |
| `.unlazy/imputation-sim/gates/*` (5 ledgers unmet) | n/a | n/a | n/a | PROTECTED: another lane's (`arc/imputation-sim`) |
| "89 UNPUSHED on other branch(es)" in the gate output | n/a | n/a | n/a | PROTECTED: other branches in the shared repo |
| Result copies (clusters and Mac pool, above) | n/a | n/a | n/a | data of record; never delete |

FINDING-OF-RECORD: BACE as shipped gives improper multiple imputation (all final runs start from one shared
dataset); chaining fixes coverage; chained BACE ties proper frequentist imputation at n ≥ 300 on continuous
traits. vault-note: NOT WRITTEN. Brain writes need Shinichi's explicit approval; the canonical source is this
branch (`docs/dev-log/after-task/2026-09-25-rubin-campaign.md`). Also propose a brain note that castor is
pigauto's frequentist discrete arm: semantic search did not surface it, and only a repo search found it.

## Files (this lane, `git diff --stat 3f28e26..HEAD`)

- Runners: `script/rubin_lib.R`, `rubin_freq.R`, `rubin_bace.R`, `rubin_cell.R`, `rubin_checks.R`,
  `rubin_convcheck.R`, `rubin_prerun_summary.R`
- Cluster tooling: `script/rubin_campaign_nibi.sh` (grid or `RETRY_FILE`, `TAG`, task cap, `CONFIRM=yes`),
  `rubin_prerun_nibi.sh`, `rubin_prerun_totoro.sh`, `rubin_totoro.sh`, `rubin_status.sh`, `rubin_retry_prep.R`,
  `rubin_pool.sh`, `rubin_missing.R`
- Aggregation and pages: `script/rubin_campaign_aggregate.R`, `rubin_report_build.sh`, and the templates
  `rubin_{report,summary,accuracy}_template.html`
- Tests: `script/tests-rubin/test-*.R` (13 files)
- Docs (force-added, because `docs/` and `.unlazy/` are gitignored): `docs/dev-log/arc/2026-09-2{4,5}-rubin-*`,
  `docs/dev-log/after-task/2026-09-2{4,5}-rubin-*.md`, `docs/dev-log/plan-actual/2026-09-24-rubin-steps1-2.md`,
  `.unlazy/rubin-freq-bace/{GOAL,arcs,checkpoint}.md`, and this file.

## Next immediate steps (OWED)

1. Run `~/shinichi-brain/tools/lane_preflight.sh "/Users/z3437171/Dropbox/Github Local/pigauto-rubin-freq-bace"`,
   re-claim the lease, and classify this handover's items OWED / DONE / RETRACTED / PROTECTED.
2. Open an ultra-plan for the discrete extension. Before spending compute, say back to Shinichi in one
   sentence what you will build.
3. Settle Q1 to Q3 below with Shinichi, with your recommendation drafted for each.
4. Mac smoke (under 3 h):
   - pick one stored dataset (n = 100) and confirm the rebuilt `truth` and `mask` match the stored ones;
   - re-run BACE on it and check the continuous rows reproduce;
   - run the castor arms;
   - save the imputed discrete columns for all M;
   - score the chosen discrete estimands.
5. Pre-run plan with a measured time estimate, then STOP for approval (D-287).

## Open questions for Shinichi (draft a recommendation for each)

**Q1. Frequentist discrete arm.** Suggested: castor, as in v1, for all three types.
- Plug-in (fixed fitted rates) is the freq B analogue.
- A parametric bootstrap of the Mk rates is the freq A analogue: simulate states on the tree at the fitted
  rates, refit, recompute the tip probabilities, then draw.
- Draws from castor's per-tip marginal probabilities ignore dependence between tips; check whether a joint
  draw is needed.
- Option for binary only: a phylogenetic probit GLMM (`glmmTMB`, installed) inside chained equations, which
  would also use the continuous traits. Needs a feasibility check that the random effect for species with a
  missing response can be recovered.
- There is no off-the-shelf phylogenetic frequentist method that imputes ordinal or categorical traits jointly
  with continuous ones. Report that as a finding.

**Q2. What to score.**
- Per value: accuracy of the pooled class, Brier score and calibration of the pooled class probabilities.
- Coverage analogue: how often the true class falls in the prediction set built from the 20 draws.
- Downstream, with Rubin pooling: for example the PGLS slope of c1 on `bin`, or a phylogenetic logistic
  regression `bin ~ c1` (`phylolm::phyloglm`). Neither has ρ as its true value, so compute the truth,
  analytically on the liability scale or as a large Monte Carlo mean of the complete-data estimate.

**Q3. Scale.**
- Suggested: n = 100 and 300 first (all λ and ρ, 200 datasets per cell, the same seeds as the continuous
  campaign), with n = 1000 decided after.
- castor is cheap. BACE dominates the cost: about 6 h per fit and 19 GB at n = 1000 on fir. Time n = 100
  and 300 from the smoke.

## Gotchas and failed approaches (do not repeat)

- **Save the imputations this time**, at least the discrete columns for all M.
- `RNGkind("L'Ecuyer-CMRG")` leaks between test files unless it is restored.
- BACE stops when a discrete level is unobserved (common at λ = 1). `fit_bace_mi()` drops unused levels and
  leaves out a discrete trait with fewer than 2 observed classes. For discrete scoring, decide how a dropped
  trait or level counts.
- Seeded BACE failures on fixed code are deterministic; retry only after a code change. MCMCglmm segfaults are
  not deterministic; rerun those.
- Request cluster memory from measured peak plus margin (40 GB blocked fir scheduling; 24 GB ran).
- Monitor queue counts: `rubin_status.sh` prints `queued=` twice on the fir line. Count distinct seeds with
  `rubin_missing.R`, not host totals.
- `git add --force` for `docs/` and `.unlazy/`, as a separate command.
- The destructive-command hook blocks `rm -rf`; move files to the scratchpad instead.
- Shipped prose: no em dashes; run `python3 ~/shinichi-brain/tools/slop_check.py <absolute path>`.
- Commit attribution line: `Co-Authored-By: Claude Opus 5.5 (1M context) <noreply@anthropic.com>`.

## How to resume

Environment: R with BACE (installed), MCMCglmm, Rphylopars, castor, glmmTMB, phylolm, nlme, ape, testthat.
- Safe check from the worktree: `Rscript -e 'testthat::test_dir("script/tests-rubin")'` (about 5 minutes;
  expect FAIL 0).
- Cluster status: `bash script/rubin_status.sh`.
- Never stage `.unlazy/imputation-sim/` or files outside this lane.

Read in order:
1. this file;
2. `docs/dev-log/arc/2026-09-23-simulation-v1-summary.md`;
3. `docs/dev-log/after-task/2026-09-25-rubin-campaign.md`;
4. `docs/dev-log/arc/2026-09-24-rubin-freq-bace-plan.md`;
5. `script/rubin_cell.R` and `script/campaign_gnn_off_lib.R` (`run_freq()` for castor);
6. the repo `AGENTS.md`.

```text
Read AGENTS.md and docs/dev-log/handover/2026-09-26-claude-handover.md. Run the handover rehydration steps, reconcile them with the current git state, then continue only the OWED Next Immediate Steps.
```
