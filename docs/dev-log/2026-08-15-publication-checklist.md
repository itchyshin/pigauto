# pigauto → publication: the checklist

Compact follow-along version of `docs/dev-log/2026-08-15-publication-readiness-roadmap.md`
(which holds the evidence, file paths, and measurements behind every line). 2026-08-15.

**Rule of thumb:** nothing below Tier 0 is worth doing before Tier 0 closes.
~~no GNN number measured before the P0 land may appear in the manuscript~~ —
**relaxed 2026-08-16 with evidence**: the re-run reproduces the pre-fix per-type numbers
within 1.1 MCSE in every trait type (`2026-08-16-bench-rerun-results.md`). Cite the
**re-run** values; the pre-fix ones agree within noise.

## Tier 0 — correctness (local Mac, ~1 week total)

- [x] **Decide #157** *(2026-08-15: surface B — decided and shipped)*: which surface does the "never worse than baseline" invariant live on —
      the one users get (truth hidden, recommended) or the one `predict()` currently scores?
      → [github.com/itchyshin/pigauto/issues/157](https://github.com/itchyshin/pigauto/issues/157)
- [x] Fix the calibration↔predict context mismatch *(two-layer floor, PR #158)* per that decision; make the gate-acceptance
      guards assert on the delivered surface
- [x] **Land P0** *(PR #158 merged; `main` = `3677a85`)* (`arc/p0-onto-main` `f5e8416` — merge already done and verified; needs a
      fresh G0 lock that day)
- [x] Fix the **BLOCKING** vignette claim *(9-edit bundle in PR #158)* ("Exact under BM" → split by path) + 5 SHOULD-FIX
      → `docs/dev-log/arc/2026-08-15-p0-claim-gate-findings.md`
- [x] `rcmdcheck --as-cran` clean *(0E/0W/1N — known incoming note only)* (expect only the pre-existing 2 WARN / 3 NOTE)

## Tier 1 — GNN evidence (Totoro; D-139: estimate + pre-run test first)

- [x] **Re-run all 8 per-type benches** on the corrected pipeline *(2026-08-16: 5/8 done —
      continuous, binary, ordinal, count, categorical; all within 1.1 MCSE of the pre-fix
      values. proportion / zi_count / multi_proportion running.)*
- [x] **Regime map (#135)** *(2026-08-16: 5,400 jobs, 0 failures —
      `2026-08-16-regime-map-results.md`. F2-nonlinear lift real at λ=1 (58% of cells);
      low-λ gains flagged as baseline misspecification)*
- [x] **λ-attribution follow-up** *(2026-08-16: 1,920 jobs —
      `2026-08-16-lambda-attribution-results.md`. λ-fitting reproduces 100–117% of the
      low-λ gains; F2 lift SURVIVES a λ-fitted baseline (5.9–8.5 MCSE); `bayes` kills the
      ML boundary collapse)*
- [x] **Conformal coverage campaign** *(2026-08-16: two campaigns —
      `2026-08-16-coverage-results.md` (+CORRECTION) and
      `2026-08-16-mechanism-coverage-results.md`. n=100 ceiling is arithmetic;
      exchangeability failure CONFIRMED under MAR_phylo/MNAR on the genuinely-missing
      surface; Mondrian fix implemented on `arc/mondrian-conformal`, verification
      campaign pending)*
- [ ] **Then decide the paper's GNN framing** — "lift in mapped regimes" vs "safety-gated
      architecture". Not before. *(Evidence now in: the defensible accuracy claim is
      F2-nonlinear-on-λ-fitted-baseline; the safety framing is fully supported.)*

## Tier 2 — latent multivariate

- [ ] **EM decision**: scope the claim to "single-pass Henderson-init joint MVN" (cheap) OR
      restore `max_iter > 0` under a dedicated numerical G0 with recovery sims
- [ ] **Known-Σ recovery simulation**: Σ recovery, per-tip SE coverage, liability-decode
      probability calibration (needs P0 landed first)
- [ ] **P1-8**: covariates silently ignored by joint/threshold-joint — thread through, or
      detect-and-message
- [x] **P1-9** zi_count zeros excluded from val/test *(fixed on main pre-arc; NEWS)*
- [ ] **P1-11** tree-MI conformal draws — doc claim now honestly scoped in roxygen
      (2026-08-16); mechanical threading of the conformal sampler still open
- [x] **P1-12** multi-obs phylo-signal gate no-op *(fixed on main pre-arc; NEWS)*
- [x] **P1-13** `read_traits` count fix *(PR #161)*; P1-6 / P1-7 docs consolidation open
- [x] **P1-8 (detect-and-message form)** + λ-under-covariates + n_val<19 ceiling warning
      + `conformal_split_val` exposure *(PR #167, 2026-08-16)*
- [ ] **OVR categorical head-to-head** vs LP / multinomial / BACE (the 72% number is BACE's,
      not pigauto's)

## Tier 3 — packaging

- [~] Head-to-head vs Rphylopars / BACE / missForest, post-P0, regime-labelled
      *(2026-08-16: FIRST WORKING BACE run — `2026-08-16-external-comparison-results.md`:
      BACE wins continuous (up to −49% RMSE), pigauto wins categorical (+30pp);
      `lambda_mode="bayes"` closes part of the continuous gap at a categorical price.
      Single-seed; Rphylopars/phylolm standalone bench in progress; multi-seed +
      missForest still open)*
- [ ] CRAN: drop Suggests `BACE` or wait for BACE on CRAN (404 as of 2026-08-15)
- [ ] Every manuscript number has a committed script that regenerates it

## Paper shape (short version)

Lead with what no competitor combines **today**: mixed-type + conformal UQ + tree-MI Rubin +
`suggest_next_observation()`. GNN = safety-gated correction, framed by the Tier-1 map.
Venue: MEE methods/application; JOSS fallback. **Never publish**: pre-P0 numbers, "exact
under BM" for the joint path, wrap-vs-wrap coverage as pigauto-vs-BACE, any EM claim while
`max_iter = 0`.

## First action on return

**The #157 surface decision.** Everything else queues behind it.
