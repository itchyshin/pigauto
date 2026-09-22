## 1. Goal

Run the ADEMP phylogenetic-imputation simulation (frequentist stack, BACE, pigauto GNN off/on) to completion on the corrected design; deliver the private results Artifact, the BACE-paper methods note and the unlisted pkgdown article.

## 2. Implemented

Compute, on the corrected design (Pagel lambda on the tree correlation with no double-counted residual,
fixed population thresholds, an always-observed MAR driver with its own 0.35 loading):

- Core slice, 18 cells x 200 replicates, every arm: frequentist stack at Rphylopars' `model = "BM"`
  default (`freq`) and at `model = "lambda"` (`freq_lambda`, the fifth arm Shinichi added on 2026-09-21),
  BACE (100 replicates), pigauto GNN off with the in-house and the Rphylopars joint solver, pigauto GNN on,
  mean/mode floor. Totoro.
- Factorial, 56 cells x 200 replicates, every fast arm complete; BACE at 100 replicates (n = 100) and 30
  (n = 1000, Shinichi 2026-09-21), 98.4% complete with the clade n = 1000 tail recomputing. Totoro (fast
  arms; the wave the permission classifier had refused three times overnight was launched 2026-09-21 09:19
  and landed 2026-09-22 09:02), nibi and fir (BACE).
- AVONET300 case study, 20 seeds, every arm. Covariate sensitivity on the 18 core cells (ncov = 2),
  3,584 of 3,600 on fir plus the freq_lambda wave on Totoro.
- Aggregation with paired MCSE, cross-host dedupe on (filename, arm set), a missing-fraction key, coverage
  as its own metric row, pooled ECE with a bootstrap MCSE, and (this session) a divergence rule that scores
  finite-but-absurd fits at the floor and counts them.

Deliverables:

- (a) Private results Artifact, https://claude.ai/artifact/M5HtGRnNGfwsK2Se4gMX24, v2: Frequentist vs
  BACE (Dan's view), All arms, Factorial, AVONET, Cost and failures, Decisions (six, each with my reading).
- (b) `docs/dev-log/arc/2026-09-20-simulation-methods-bace.md`: ADEMP methods for the BACE paper plus
  Results for the core slice and the factorial, every number with its regime and MCSE.
- (c) `vignettes/articles/simulation-study.Rmd`: the four-arm accuracy and coverage report, pre-rendered
  from committed csv under `script/campaign_sim_results/`, `.Rbuildignore`d, unlisted, unpublished.

## 3a. Decisions and Rejected Alternatives

- **Report the frequentist stack twice** (`freq` at the BM default, `freq_lambda` with the signal estimated)
  rather than disclose the default as a caveat. Shinichi first chose "BM only and disclose", then reversed
  on 2026-09-21 after the single-seed probe. The difference (0.20 to 0.25 z-RMSE at lambda = 0.3) changes
  the BACE paper's headline, so a caveat would have hidden a finding.
- **Divergence rule** (2026-09-22): z-RMSE above 3x the floor or interval score above 1e3 is a failure that
  did not throw; scored at the floor, interval dropped, counted. Rejected: trimming or winsorising (hides
  the rate), or dropping the replicates (violates "never dropped"). The 1e3 threshold is two orders of
  magnitude above the largest sane value observed (12) and is stated wherever the rule is used.
- **BACE at 30 replicates for n = 1000 in the factorial** (Shinichi, 2026-09-21), after a successful fit
  measured at about 3 h and over 16 GB and a 72% singular-failure rate at lambda = 1. MCSE widens by about
  1.8x there and every table says so. Rejected: 100 replicates (about 28,000 core-hours against an
  approved 11,564).
- **BACE convergence judged by BACE's own `assess_convergence()` verdict and median ESS**, not Gelman-Rubin
  across `runs` (which are sequential imputation iterations, not chains). `runs = 5`, `n_final = 20`.
- **Failures scored at the floor and reported beside every metric**, including BACE's 16 to 39% error rate
  at lambda = 1. Rejected: reporting BACE on its successful fits only.
- **In-house joint solver stays pigauto's default.** The Rphylopars solver diverges at high signal and
  small n in simulation (34 divergent replicates, z-RMSE 1.654 vs 0.487 at lambda = 1, n = 100) yet is the
  most accurate arm on AVONET; that split is recorded, not resolved here (a package change is fenced).
- **Per-host pool subdirectories with recursive read**, after `rsync --ignore-existing` silently dropped
  1,034 BACE cells that shared filenames with fast-arm cells.

## 4. Files Touched

Worktree `../pigauto-imputation-sim`, branch `arc/imputation-sim` (from `origin/main`). No edits to `R/`,
`BACE/`, PR #175 files, or the dirty `handover/2026-08-09-cursor` checkout.

- `script/campaign_gnn_off_lib.R` (DGP, mechanisms, scoring, `freq_lambda` dispatch, BACE diagnostics)
- `script/campaign_gnn_off_aggregate.R` (recursive read, cross-host dedupe, `frac` key, coverage rows,
  pooled ECE, divergence rule + `n_divergent`)
- `script/campaign_sim_design.R`, `script/campaign_sim_cell.R`, `script/campaign_sim_checks.R`
- `script/campaign_sim_totoro.sh`, `script/campaign_sim_nibi_array.sh`, `script/campaign_sim_pool.sh`
- `script/campaign_sim_page_data.R`, `script/campaign_sim_page_build.sh`, `script/campaign_sim_page.template.html`
- `script/campaign_sim_results/{summary,paired,failures,ece,per_trait}.csv`
- `docs/dev-log/arc/2026-09-20-simulation-prerun.md`
- `docs/dev-log/arc/2026-09-20-simulation-methods-bace.md`
- `docs/dev-log/arc/2026-09-22-imputation-sim-after-task.md` (this report)
- `vignettes/articles/simulation-study.Rmd`, `.Rbuildignore` (`^vignettes/articles$`)
- `.unlazy/imputation-sim/{GOAL,GATES,PROGRESS}.md`, `.unlazy/imputation-sim/gates/leaf-*.md` (on disk only,
  git-ignored by design)
- Scratchpad (not in the repo): `sim-results.html` (the Artifact source), per-replicate extracts and the
  table scripts that cross-check the aggregator.

## 5. Checks Run

Gates (`script/campaign_sim_checks.R`, evidence in `.unlazy/imputation-sim/gates/leaf-*.md`):

- G14 cross-host reproducibility: `--dir totoro_fl,nibi`, `totoro_ff,fir`, `nibi,fir`: "5 overlapping
  cell(s) checked, truth/mask/freq identical / G14 PASS" x3 (the freq clause is vacuous on these pairs;
  truth and mask are the evidence).
- G10 core completeness: "18 design cells, 10023 rds present, 24 replicate-arms missing", all 24 BACE,
  every fast arm 200/200. FAIL pending the fir recovery array 60949560.
- G11 factorial completeness: "56 design cells, 39254 rds present, 57 replicate-arms missing", all 57
  BACE (56 in the four Brownian clade n = 1000 cells), every fast arm 200/200 in all 56 cells. FAIL
  pending fir arrays 60938890 (12 h) and the n100B tail.
- G13a superlatives: grep over the methods note and the article, 0 hits, PASS.
- G13b slop check: methods note 1.3 per 1000 words, article 3.0 per 1000, 0 em dashes in both,
  FINDINGS 0, PASS.
- G13d `.Rbuildignore` carries `^vignettes/articles$`, PASS.
- G12 results csv: pending the re-aggregation with the divergence rule (running at the time of writing;
  to be filled in before this report is closed).

Other checks:

- Aggregator divergence rule smoked on the one cell known to diverge: "divergence rule: 5 arm-replicates
  (9 arm-replicate-traits) scored at the floor"; `failures.csv` gained `n_divergent` per arm; that cell's
  `freq_lambda` mean fell from 1e15 to 0.181 and nothing else moved.
- Aggregator output cross-checked against an independent per-replicate read of the rds: core
  `freq_lambda` at n = 1000 0.70292 (aggregator) vs 0.7029 (direct); factorial stratum means agree to
  three decimals once both apply the same floor rule.
- `rmarkdown::render` of the article on the committed csv: exit 0; the rendered z-RMSE table matches the
  direct read to three decimals.
- Results board: `node --check` on the inline script, PASS; every `getElementById` target present.
- `agent_mention_check.py --text` on the PR body: "Agent @-mention check (G-20) — clean".
- Palette for the board validated with the dataviz validator in light and dark mode (light mode
  "relief required", satisfied by direct labels and full tables).

## 6. Tests of the Tests

- The divergence rule was tested on the cell where divergence was known to exist from the raw rds
  (OU, lambda = 1, MCAR 10%, n = 1000): it caught exactly the 3 `freq_lambda` and 2 Rphylopars-solver
  replicates identified beforehand by the per-replicate read, and no others. On the full pools its
  counts (15, 34, 2, 1) match the independent read.
- G14 was shown to be partly vacuous rather than assumed to be strong: the frequentist clause cannot fire
  when one host holds only BACE, so the report says the evidence is truth and mask, not predictions.
- The aggregator race was caught because a second, independent path (direct rds read) disagreed with the
  first; the disagreement was traced to a partial pool rather than papered over, and the arithmetic that
  reconciles them (0.703 = mean of 0.878, 0.766, 0.464) is in the record.
- The failure counts in the report were taken from the `failed` flags inside the rds, not from log lines;
  the two agree (160, 160, 157 + 51), which is the check that the runner records what it logs.
- The BACE clade result was deliberately split into failed and fitted replicates before being read,
  because a floored-failure mean and a fitted-accuracy mean would have told different stories; they did.
- Negative control not run: the gates were not fed a deliberately broken pool to confirm they fail
  (G10 and G11 did fail on the real BACE gaps, which is weaker evidence than a planted defect).

## 7a. Issue Ledger

Fixed this arc:
- BACE `n_final` set by a chain-length formula to 400 (about 80x cost) -> env `PIG_BACE_NFINAL`, 20.
- `runs = 2` below `assess_convergence()`'s `min_iterations = 3`, every cell "not converged" -> `runs = 5`.
- Gelman-Rubin misapplied to BACE's sequential runs -> BACE's own verdict + ESS.
- macro-F1 deflated by classes absent from the masked truth -> averaged over present classes.
- MAR driver degenerate at rho = 0 -> its own 0.35 loading with a PD check.
- MCMCglmm SIGSEGV under `mclapply` -> one process per cell.
- `rsync --ignore-existing` dropping 1,034 BACE cells -> per-host subdirectories.
- Coverage never emitted as a metric row (board's coverage panel empty on every version) -> emitted.
- Cross-host double counting (650 to 1,047 duplicates) -> `(filename, arm-set)` dedupe.
- MCAR 0.10 and 0.30 sharing a key -> `frac` in every key.
- Duplicated BACE wave (wrong env var), stale array left running, rorqual inode exhaustion, fir OOM at
  16 GB (-> 1 core, 32 GB), nibi/fir TIMEOUTs after `runs` 2 -> 5 (-> longer `--time`, resume-skip).
- The factorial fast-arm wave never having run (three overnight refusals) -> launched 2026-09-21 09:19.
- Divergent fits dominating stratum means -> divergence rule in the aggregator (1ec3e98).
- Over-broad coverage claim in the methods note -> scoped to the two arms it compares (bee0e87).
- My own misdiagnosis of an "aggregator bug" -> it was a race with my rsync; corrected in the record.
- Mac rsync 2.6.9 lacking `--info=stats1`, hidden by a grep -> flag removed.
- A stuck covsens array (0 running, "JobArrayTaskLimit") -> scancel + resubmit with resume-skip.

Deferred (out of scope, recorded):
- pigauto `gnn_on` dimension error on a monomorphic one-hot at lambda = 1 (160 replicates in the factorial,
  32/29/16 in the core): a package robustness defect, floored here, fixed separately (test file
  `tests/testthat/test-monomorphic-discrete.R` exists on the main checkout's dirty tree).
- Rphylopars `model = "lambda"` type error after a singular solve (51 replicates, one OU cell): upstream.
- Rphylopars joint-solver instability at high signal: the arm-3 default-solver decision is a package change.
- Szymek has not yet signed off on the corrected Pagel-lambda parameterisation.

## 8. Consistency Audit

Same-class sweeps after each defect found:

- Filename-collision class (after `rsync --ignore-existing` dropped BACE cells): every pool now uses
  per-host subdirectories, the aggregator reads recursively and dedupes on (filename, arm set), and the
  article, the board and the methods note all read the deduped aggregate. Checked that AVONET and covsens
  pools follow the same layout.
- Partial-pool class (after the rsync race): every aggregation launched after a pull was started only
  after the pull's file count matched the source; the covsens waiter pulls only when both hosts are at
  3,600.
- Divergence class: the rule is applied in the aggregator for both stages, and the per-replicate scripts
  apply the identical rule, so no reported number depends on which path produced it. Checked that the
  discrete metrics show no divergence (max accuracy 1, Brier in [0, 1]) so the rule's continuous-only
  scope loses nothing.
- Arm-parity class (after `freq_lambda` was added): the fifth arm was run on every stage where the
  others were (core, factorial, AVONET, covsens), and the discrete-path identity between `freq` and
  `freq_lambda` was checked in the tables (0.5 of a point) as an internal control.
- Over-broad-claim class (after the coverage sentence): grepped both shipping documents for "no arm",
  "every arm", "always", "never" and re-scoped each to the arms and cells it was measured on.
- Stale-count class (after the plan's "40 cells"): checked the design table, G11, the ledger, the article
  and the methods note all say 56; the plan file itself is left as the historical record.
- Clock class: PROGRESS.md timestamps from 2026-09-22 morning were written from Totoro's clock and my
  own labels drifted about an hour ahead of the Mac's; the entries are ordered correctly and the labels
  are not load-bearing, so they are left as written with this note.

## 9. What Did Not Go Smoothly

- The factorial fast arms were never run overnight: the permission classifier refused the Totoro wave three
  times and the progress record logged it each time, but no wave landed until I found the 3,945 BACE-only
  files on 2026-09-21 morning. That cost about 24 h on the critical path.
- The in-session watchers die with every session gap; the compute never did. Six re-arms.
- I misdiagnosed a race with my own rsync as an aggregator bug and said so to Shinichi before checking; the
  correction was made within the hour and the arithmetic (0.703 = mean of the three per-lambda values) is in
  the record.
- The pre-run's `--time` values were sized before `runs` went from 2 to 5, so the first nibi and fir arrays
  timed out at scale; the clade n = 1000 BACE cells then timed out again at 5 h and needed a 12 h array.
- rorqual ran out of inodes (500K/500K on /project) and was retired; fir OOM-killed 200 tasks at 16 GB.
- Two BACE-only clusters made G14's "frequentist arm identical" clause vacuous; the evidence there is truth
  and mask bit-identical across three machines, and the report says exactly that.
- The board could not be looked at before publishing: the in-app browser is not signed in and I will not
  sign it in. The JS was parse-checked and every element id verified instead.

## 10. Known Residuals

- BACE is 96% complete in the core (576 of 600) and 98.4% in the factorial (3,493 of 3,550). Recovery
  arrays are running on fir; every BACE figure in the deliverables is on the replicates present, and its
  MCSE says so. G10 and G11 will pass only when they drain.
- Covariate sensitivity (S6d) is at 3,584 of 3,600 plus the freq_lambda wave; it is not yet in any
  deliverable.
- G12 has not yet been run on the aggregate produced with the divergence rule (re-aggregation in
  progress at the time of writing).
- The results board has not been looked at by the agent: the in-app browser is not signed in and
  signing in is out of bounds. Its script parses and every referenced element exists; layout and
  rendering are unverified until Shinichi opens it.
- G13c (Shinichi has read the board and recorded the six publication decisions) is open; the pkgdown
  article stays unpublished until it is closed.
- pigauto `gnn_on` errors on a monomorphic one-hot at lambda = 1 (160 factorial and 77 core replicates,
  floored). A package defect, out of scope here, with a test file already on the main checkout's dirty tree.
- Rphylopars at `model = "lambda"` fails with a type error after a singular solve in one OU cell (51
  replicates) and diverges numerically on 15 others; upstream behaviour, floored and counted.
- The Rphylopars joint solver's instability at high signal versus its accuracy on AVONET is recorded and
  not resolved; the default-solver question is a package change.
- Szymek has not signed off on the corrected Pagel-lambda parameterisation.
- G14's frequentist clause was never exercised, because no host pair ran the frequentist arm twice.
- The board's Decisions tab describes the article as "currently rendering the core slice"; it now also
  carries the factorial. To be refreshed at v3 with covsens.

## 11. Team Learning

Memory receipt: loaded `route.py pigauto` LOAD-FIRST (compute default, recovery-to-truth, prediction-path
audit lane, r_cal = 0 fallback), the repo `AGENTS.md`/`CLAUDE.md` "how to be useful here" standards (regime
on every number, no superlatives, smallest defensible claim), the brain's D-139 estimate-before-you-run,
D-143/D-254 core caps, D-64 no-Duo sockets, and `70-missing-data-simulation-design` (coverage MCSE binds
replicates). The ones that shaped the work: regime-on-every-number (every table here carries n, lambda,
mechanism, arm, replicates, MCSE), the failure-partition rule from *A simulator can fail in your package's
favour* (the BACE clade result was split into failures vs fitted accuracy before it was read), and the
pre-run discipline (the 30-min Totoro pre-run caught `n_final = 400`, `runs = 2`, and the misapplied
Gelman-Rubin before any budget was spent).

Golden Set: not in scope (no known-mistake class from `tools/memory_regression.py` applies to a simulation
campaign).

Durable lessons (candidates for WHAT-WORKS / LESSONS, to be filed with approval):
1. A fit that returns a finite absurd number is a failure that did not throw; give the aggregator a
   divergence rule and a counter, or one replicate rewrites a stratum. (WHAT-WORKS)
2. Partition a "collapse" by mechanism before reading it: failures-floored and fitted-accuracy are two
   results, and they pointed in opposite directions here (inaccurate at low signal, fragile at high).
3. Never aggregate a pool while an rsync into it is still running; a partial pool produces confident wrong
   numbers with the right schema.
4. On a shared-filename pool, dedupe on (filename, arm set), never on filename.
5. A permission refusal is a silent stall: if a launch is refused, the record must say the wave did NOT run
   in the completion summary, not only in the log line, or the next reader assumes it did.
6. Size Slurm `--time` after the last cost-changing parameter change, not before.

## 12. Cross-Product Coverage

Cross-cutting things this arc touched, and what they do and do NOT cover.

- **Evolutionary model of the frequentist stack (`model = "BM"` vs `"lambda"`)**: covers ✓ the continuous
  family via Rphylopars in the core, the factorial, AVONET and covsens; does NOT cover ✗ the discrete path
  (castor Mk is shared, so the two arms are identical there by construction, and the tables confirm it to
  0.5 of a point), nor Rphylopars' OU/EB/kappa/delta models.
- **Divergence rule**: covers ✓ zRMSE, coverage, width and interval_score in the aggregator for every arm
  and stage; does NOT cover ✗ accuracy, Brier, macro-F1 or ECE (no divergence was observed on the discrete
  metrics, and the rule would need a different criterion there); does NOT cover ✗ the per-cell rds
  themselves, which keep the raw values.
- **Failure-at-floor scoring**: covers ✓ errored and divergent fits in every reported mean; does NOT cover ✗
  a BACE fit that converged badly without erroring (its `converged` verdict and ESS are stored per cell and
  reported as a rate, not folded into the means).
- **Cross-host reproducibility (L'Ecuyer-CMRG)**: covers ✓ truth and mask on Totoro, nibi, fir; does NOT
  cover ✗ the frequentist arm's predictions across hosts (no host pair ran that arm twice), nor torch
  non-determinism in the GNN arm (single host).
- **Pagel-lambda DGP correction**: covers ✓ BM cells (G4 recovery within 0.05 at n = 1000); does NOT cover
  ✗ OU cells, which are exempt from G4 by design because the OU correlation is standardised before mixing.
- **Conformal intervals**: covers ✓ MCAR and MAR (0.90 to 0.95); does NOT cover ✗ clade-biased missingness,
  where exchangeability fails and coverage drops to 0.86, exactly as the design predicted.
- **Design stages**: core, factorial and AVONET complete for every fast arm; does NOT cover ✗ multi-obs,
  multi-proportion, zero-inflated traits, MNAR, trees above 1,000 tips, or the TDIP/GAIN arm (fenced).
