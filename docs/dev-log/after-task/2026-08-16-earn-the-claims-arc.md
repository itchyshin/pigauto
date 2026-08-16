# After-task report — "Earn the claims" arc (λ attribution · conformal honesty · comparators · Tier-2)

Date: 2026-08-16 · Lane: Claude (solo) · Plan approved by Shinichi (all four slices;
conformal fix contingent on diagnosis — both conditions honoured).

## 1. What was asked

Plan a "next big big arc to improve the existing functionalities," then execute it.
Four slices locked: (A) λ-attribution campaign; (B) conformal coverage diagnosis →
contingent fix; (C) working external comparisons; (D) Tier-2 functionality fixes.

## 2. What was delivered

- **A (headline):** 1,920-job Totoro campaign, 0 failures.
  `2026-08-16-lambda-attribution-results.md`. The regime map's low-λ GNN gains are
  **entirely** reproduced by fitting λ in the baseline (closure 1.00–1.17 at λ=0.2);
  the GNN's F2-nonlinear lift **survives** a λ-fitted baseline (5.9–8.5 MCSE);
  `lambda_mode="bayes"` eliminates the ML boundary collapse (P(λ̂<.05): 0.47→0.00).
- **B1:** 200-job mechanism-coverage campaign, 0 failures.
  `2026-08-16-mechanism-coverage-results.md`. Exchangeability failure CONFIRMED on the
  genuinely-missing surface (MAR_phylo −3.4 pp, does not wash out with n); MCAR control
  healthy — the pipeline is fine, the distribution shift is real.
- **B2 (contingency fired):** `conformal_method = "mondrian"` implemented
  ([PR #168](https://github.com/itchyshin/pigauto/pull/168)), verified on a paired
  200-job re-run: at n=1000 all MAR/MNAR arms recover to within 3×MCSE of nominal,
  MCAR untouched, width +2.3%. The first verification run caught a design flaw (13-cell
  strata pinned at the 13/14 ceiling) → stratum floor raised to 19.
- **C:** First working pigauto-vs-BACE head-to-head in the repo's history
  (`2026-08-16-external-comparison-results.md`): the historical "BACE skipped" was an
  API-signature bug, not a missing install. BACE wins continuous (up to −49% RMSE, with
  SHORT chains); pigauto wins categorical (+30 pp). `lambda_mode="bayes"` closes part of
  the continuous gap at a categorical price (force_per_column drops the joint baseline).
  Rphylopars/phylolm standalone bench built (`script/bench_external_comparators.R`,
  arc/bace-comparators).
- **D:** [PR #167](https://github.com/itchyshin/pigauto/pull/167): λ-under-covariates
  warning; P1-8 joint-covariates warning; n_val<19 ceiling warning (including the
  previously-silent 10–18 band); `conformal_split_val` exposed in `impute()`; NEWS for
  the undocumented cv/bayes modes. Full suite FAIL 0 / PASS 2017.

## 3. What was verified, and how

Every campaign: per-job RDS, failure counts printed, MCSEs on every aggregate,
pre-registered decision rules written before results. Code: targeted tests → full suite →
PR CI; `--as-cran` gates both PRs before merge (running at close). B2 verified against a
pre-registered acceptance target on paired seeds. The B2 flaw was caught by the
verification campaign itself — the campaign-first discipline paid for itself same-day.

## 4. What is NOT covered (regime honesty)

- λ-attribution: simulated, MCAR, continuous-only, no-covariate, n≤300, single-obs.
  λ never reaches discrete paths. λ̂ diagnostics approximate (recomputed).
- Mechanism coverage: F1 λ=1 simulated; magnitudes smaller than the worst real-data
  cases; discrete types unmeasured. Mondrian: no real-data re-run yet; MNAR only
  partially repaired; inactive below ~38 val cells/trait (by arithmetic necessity).
- BACE comparison: ONE dataset, ONE seed, short chains — a direction, not a publishable
  head-to-head. The λ-modes appendix likewise single-seed.
- The two per-type benches never re-run post-#158: proportion refinements pending from
  the previous arc are unaffected by this arc.

## 5. Surprises

1. The n=300 Mondrian failure mode was the arc's OWN documented ceiling arithmetic —
   filed in the morning as a warning fix, bit the fix in the afternoon.
2. BACE was installed all along; the skip string conflated two causes and hid a
   signature bug for months. Corollary: pigauto itself was not installed in the local
   R 4.6 library — "bench skipped" locally had the same shape.
3. P1-9, P1-12, P1-13 were already fixed on main — the checklist was stale; the arc's
   Slice D shrank accordingly rather than re-doing done work.
4. D-143 breach mid-arc: three lanes (two foreign) ran campaigns on Totoro under one
   account, 192 cores total, each lane individually compliant. Lesson drafted for the
   brain: `2026-08-16-lesson-compute-lane-preflight.md` (compute-host preflight, not just
   repo preflight). My campaign finished before the throttle landed; breach ~10 min.

## 6. Decisions taken (and whose they were)

- Mondrian over weighted conformal (mine, from B1's MAR_phylo ranking; recorded).
- B2 floor 19 (mine, forced by measurement).
- P1-11 threading + zi_count intervals DEFERRED with reasons (mine): the conformal
  sampler lives in `multi_impute()`'s caller (refactor, doc already scoped honestly);
  zi_count needs interval design for E[X]=P·E, not a naive ±q.
- Paper framing: NOT taken — Shinichi's, and the evidence is now assembled for it.

## 7. Cost

Totoro: ~2,440 jobs across 4 campaign waves, ~1 h wall total, 0 failed jobs. Local: two
bench runs, full suite ×2, 2 PR builds. Token economy: 2 Explore scouts + 4 Sonnet
implementation agents; Fable/Opus kept for design, gates, and verification review.

## 8. State of the tree

- `origin/main` = `1beda8e` (#167 + #168 MERGED; CI green ×3 platforms; --as-cran 0E/0W/1N each).
- `handover/2026-08-09-cursor`: +8 doc commits (designs pre-registered before results;
  results; corrections; checklist ticks; this report).
- Worktrees added: `.worktrees/arc-tier2` (PR #167), `.worktrees/arc-mondrian`
  (PR #168), `.worktrees/arc-bace` (arc/bace-comparators, pushed `da78100`). All three prunable.
- 15 carried-over files: untouched throughout (verified at every commit).
- Totoro: `~/pigauto_regime_map/` gains lambda_attr (1,920 rds), mech_cov (200),
  mech_cov_mondrian (200 valid + 120 archived bad), summaries, lib_mondrian.

## 9. Open items for the next session

1. ~~Merge #167 and #168~~ — **DONE** (main `1beda8e`; both --as-cran 0E/0W/1N; the #168
   CI failure was a platform-dependent hard-coded test reference, rewritten deterministic).
2. ~~Commit + push arc/bace-comparators~~ — **DONE** (`da78100`).
3. **Diagnose the continuous gap** (NEW, front page): raw Rphylopars beats pigauto's
   continuous output on AVONET300 (5 seeds, C2 bench). Testable order: calibration tax
   (~49% vs 70% visible cells), which baseline path actually fired, λ=1 default.
4. Multi-seed external comparison (the single-seed BACE result begs it) — Totoro-sized.
5. Per-type λ dispatch (λ per-column for continuous + joint for discrete) — removes the
   bayes categorical price; natural next code slice.
6. Real-data Mondrian re-run (fishbase/pantheria).
7. P1-11 threading; zi_count interval design; P1-6/7 docs consolidation.
8. Paper framing decision (Shinichi) — the evidence set is now: F2-nonlinear lift on
   λ-fitted baselines + safety-gated architecture + honest coverage story + one working
   external comparison.

## 10. Claims discipline check

No claim in this arc's docs exceeds its regime; every table carries MCSEs or a
single-seed disclaimer; the two "wins" (λ-attribution answer, Mondrian repair) are
each bounded by an explicit failure or limitation in the same document.

## 11. Handover

Done in the same close: ⚡ UPDATE block atop `handover/2026-08-16-claude-handover.md`,
phase-snapshot replaced (old entry archived verbatim), coordination-board Status entry
added. START HERE for the next session = the handover's UPDATE block.
