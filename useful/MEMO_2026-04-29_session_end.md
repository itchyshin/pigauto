# Session-end synthesis: 2026-04-29

This memo summarises the day's work end-to-end so you have everything
in one place at 5 am.

## Branch state

`experiment/gnn-earnings-sim`, unpushed.  Commits added today (newest
first):

```
a4511d9 docs: bisect Phase 6 Migration regression to commit a541dbd
74af0a5 sim: AmphiBIO discrete bench script (blocked by Rphylopars bug)
0b3deed docs+sim: verify v3 strict val-floor preserves discrete bench closures
1ac34b1 fix(safety-floor): split val-floor into discrete-strict + continuous-legacy (v3)
39cdcfc docs(claude.md): add 'How to be useful here' guardrails
720fab0 docs+sim: refresh discrete benches post strict-val-floor fix
ba018dc docs+sim: AE-attribution smoke + 2026-04-29 strategic memo
d651c6c fix(safety-floor): extend strict val-floor to ALL trait types
4e33153 docs: fix two pre-existing R CMD check Rd warnings
5620512 test/lint: post-Opus cleanup
fb4461e fix: address four Opus adversarial-review correctness bugs
```

R CMD check: 0 errors / 0 warnings / 0 notes.
Full `devtools::test()`: FAIL 0 / SKIP 2 / PASS 1081.

## What this session shipped (in order of consequence)

### 1. Strict val-floor fix v3 (commit 1ac34b1, refined from d651c6c)

The morning fix (`d651c6c`) extended the post-calibration "strict
val-floor guarantee" check to all trait types via a single
`cal_mean_loss`-based block.  Apples-to-apples real-data
verification at AVONET n=1500 same-seed surfaced an over-correction:
on high-phylo-signal continuous traits the strict check spuriously
overrode the calibrated blend to pure-BM (Mass test RMSE went
386 -> 659 = +71 % regression).

v3 (commit `1ac34b1`) splits the post-calibration block by trait
family:

* Continuous family (continuous, count, ordinal, proportion):
  byte-for-byte the pre-2026-04-29 mean-only check.  Override only
  if blend loses to pure-MEAN.  Pure-BM is trusted via the
  half-A/half-B + median cross-check above.
* Discrete family (binary, categorical, zi_count): strict
  both-corners check.  Override if blend loses to either pure-BM or
  pure-MEAN; pick the lower-loss corner.
* multi_proportion: skipped (gate naturally closes via legacy path).

Verification (apples-to-apples, AVONET n=1500 seed=2026):

| trait | metric | pre-fix | v1 (overstrict) | v3 (current) | v3 - pre |
|---|---|---|---|---|---|
| Mass | RMSE | 386 | 659 | 377 | -8 |
| Beak | RMSE | 6.51 | 6.39 | 6.55 | +0.04 |
| Tarsus | RMSE | 11.60 | 12.04 | 11.92 | +0.32 |
| Wing | RMSE | 37.09 | 35.11 | 34.12 | -2.97 |
| categorical traits | acc | unchanged in all versions |

Continuous mean change v3 vs pre = -2.75 RMSE units (slight
improvement, within MC-dropout noise).  Discrete benches under v3
preserve the closures from commit 720fab0 (mean acc gap binary
-0.005, categorical -0.004; zi_count RMSE ratio 0.986).

Calibration evidence preserved at
`script/_strict_floor_v3_calibration/`.  Memos at
`useful/MEMO_2026-04-29_strict_floor_and_AE.md`,
`useful/MEMO_2026-04-29_discrete_bench_reruns.md`,
`useful/MEMO_2026-04-29_strict_floor_v3.md`.

### 2. CLAUDE.md guardrails (commit 39cdcfc)

Eight standards I committed to following after you pushed back on
sloppy "replace BM with phylolm" framing:

* Distinguish what was measured from what was claimed.
* Pigauto's actual architectural advantages (mixed-type unification,
  conformal intervals, multi-tree workflow, multi-obs covariate
  refinement) -- recommendations that compromise these need
  explicit case for the trade-off.
* BM is a kernel, not a method (touches Phase 2/3/6 + liability +
  GNN gate).
* Avoid superlatives without evidence.
* Scope changes precisely (file paths, in/out scope, expected effort).
* Distinguish pre-existing limitations from new bugs.
* Compare against existing benches before claiming novelty.
* Answer "what could we do more" with two or three defensible items,
  not ten.

### 3. Phase 6 Migration regression bisected (commit a4511d9)

Opus E flagged a Migration regression on AVONET 300 Phase 6 between
Apr 16 (lift +0.011 vs LP) and Apr 29 (lift -0.085 vs LP).  Bisected
via `git show <commit>:<file>` checkouts of the five baseline-relevant
R/ files at each candidate commit.

Result: regression appears at exactly commit `a541dbd`
("B3: full threshold-model ordinal baseline").  Mechanism: B3
swapped the ordinal baseline from LP to a threshold-joint phylopars
fit with interval-censored Gaussian liability.  AVONET's Migration
has only K=3 levels (resident / partial / full migrant), and at K=3
the threshold-joint estimate cannot resolve the two thresholds well
-- predictions are less accurate than LP.

Note: Phase 6's bench calls `fit_baseline()` directly, bypassing
`calibrate_gates`.  In the full `impute()` pipeline the strict
val-floor would catch this IF the safety floor's corner set
included LP for ordinal traits.  It currently does not (strict
floor uses BM and MEAN corners; LP is implicit via the BM-corner
for binary/categorical, but ordinal goes through a separate path).

Proposed surgical fix (not implemented this session): one-line guard
in `R/joint_threshold_baseline.R::build_liability_matrix` to skip
threshold-joint for ordinals with K <= 3.  Memo at
`useful/MEMO_2026-04-29_phase6_migration_bisect.md`.

### 4. AE-attribution gate-zeroing experiment (committed yesterday's
session, ba018dc, but underlies today's BM-vs-phylolm discussion)

27-cell smoke (single-obs n=200 continuous DGPs) showed:

* GNN silent in 13/27 cells (gate closes when GNN doesn't help).
* When GNN engages, biggest gains on linear DGPs (counterintuitive).
* Pigauto loses to phylolm-lambda BLUP in 6/9 single-obs cells by
  5-8 % due to BM (lambda=1) vs phylolm-fitted-lambda baseline.

This experiment provoked my morning recommendation to "replace BM
with phylolm-lambda" -- which you correctly pushed back on as
sloppy because:

* Phylolm is mono-type continuous; pigauto's USP is mixed types.
* BM is a kernel, not a method.
* Single-obs continuous is the worst regime for pigauto's GNN.
* Even within that regime, the right fix is adding lambda-fitting
  to pigauto's existing internal BM machinery, not external
  replacement.

The CLAUDE.md guardrails (commit 39cdcfc) document the lesson.

## What this session did NOT settle

### A. Real-data discrete verification of strict-floor v3

The strict-discrete check fires on binary, categorical, and zi_count
traits.  AVONET's discrete columns are too phylo-conserved (gate
already at corner) to exercise the strict check.  AmphiBIO has
better candidates (Fos/Ter/Aqu/Arb habitat indicators) but triggers
a pre-existing Rphylopars internal type-mismatch bug in
`fit_joint_threshold_baseline` -> `Rphylopars::phylopars()`.

Status: simulator benches confirm the strict-discrete check works
correctly (commit 0b3deed).  Real-data verification on AmphiBIO
binary is blocked by the upstream bug.  Other real-data candidates:

* BIEN with growth_form column (tree/shrub/herb/grass) -- different
  kingdom, weak phylo signal, probably exercises the strict check.
* PanTHERIA mammal diet categorical.
* LepTraits butterfly host plant.

### B. AVONET-missingness sweep

Skipped tonight.  Original full bench (n=9993, three missingness
levels) takes 4+ hours.  At n=1500 it would still be ~1 hour and
mostly confirm what we already know (AVONET categoricals stay at
corner gates regardless of missingness fraction; continuous family
behaviour is bounded by the n=1500 30 % MAR run we already did).

### C. Push to remote

11 commits unpushed on `experiment/gnn-earnings-sim`.  Not pushed
because you have not explicitly authorised it.  Also, none of the
commits are urgent (no production blocker; the discrete-trait
regression they fix was a 5-12 pp accuracy gap on the package's own
simulator benches, not a real-user-facing failure).

### D. Step 2 -- peer comparison on real mixed-type data

The bigger strategic experiment I proposed earlier (pigauto vs
phylolm + phyloglm + ape::ace + Rphylopars + BACE on AVONET or
similar).  Not started; ~1 week of focused work scoped tightly.

## What to consider when you return

In rough priority order, all defensible per the CLAUDE.md
guardrails:

1. **Push the 11 commits** to remote on `experiment/gnn-earnings-sim`
   (or merge to main if the strict val-floor work meets your bar).
   Cheap, recoverable.

2. **Fix the Phase 6 K-aware ordinal guard.**  One-line change in
   `R/joint_threshold_baseline.R`; closes the Migration regression
   without touching the threshold-joint path for K >= 4 ordinals.
   ~30 min.

3. **Fix the Rphylopars-AmphiBIO bug** that blocks AmphiBIO discrete
   benching.  Worth ~1-2 hours of debugging (the C++ trace points
   to `tp()` inside `estim_pars`).  Once fixed, run the AmphiBIO
   discrete bench (script committed at 74af0a5) for real-data
   discrete verification.

4. **BIEN growth_form bench** as an alternative real-data discrete
   verification.  ~1 day.

5. **Step 2 peer comparison** on mixed-type real data.  ~1 week.
   The strategic experiment that grounds the "best of best in
   mixed-type unified imputation" claim or refutes it.

6. **The lambda-fit experiment** I proposed (more carefully) at the
   end of yesterday's pushback discussion.  Add ML lambda-fit to
   `R/bm_internal.R::bm_impute_col` and re-run AE-attribution.
   ~1 day.  Defer until Step 2 finishes -- it's downstream of
   knowing whether pigauto wins peers in its target regime.

What I am NOT recommending (based on tonight's evidence):

* Self-supervised pretraining (Item #8 from the prior strategic
  review).  AE-attribution showed the GNN's bottleneck is engagement,
  not feature quality; pretraining would not address this.
* Architectural replacement of any pigauto kernel.  The session's
  most concrete lesson was that even well-intentioned fixes to
  internal kernels can produce real-data regressions (the v3 vs v1
  contrast).  Architecture changes deserve explicit user oversight.

## Files / locations of evidence

* `useful/MEMO_2026-04-29_strict_floor_and_AE.md` -- AE-attribution
  + morning strict-floor fix synthesis
* `useful/MEMO_2026-04-29_discrete_bench_reruns.md` -- three-way
  comparison (Apr 17 / Apr 28 / Apr 29) on discrete simulator benches
* `useful/MEMO_2026-04-29_strict_floor_v3.md` -- v1 vs v3 calibration
  on AVONET n=1500
* `useful/MEMO_2026-04-29_phase6_migration_bisect.md` -- bisect of
  the Apr 16 -> Apr 29 Migration regression to commit a541dbd
* `script/_strict_floor_v3_calibration/` -- raw bench outputs from
  the four AVONET n=1500 calibration runs (pre / v1 / v2 / v3)
* `script/bench_amphibio_discrete.R` -- discrete real-data test
  script (preserved for the day the Rphylopars bug is resolved)
* `script/ae_attribution_smoke.R` -- gate-zeroing AE-attribution
  experiment (yesterday)
* `script/smoke_strict_floor.R` -- focused smoke on the worst
  pre-fix discrete cell
* `/tmp/phase6_bisect2_results.csv` -- raw bisect output
