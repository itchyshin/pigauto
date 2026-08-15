# pigauto → publication: the checklist

Compact follow-along version of `docs/dev-log/2026-08-15-publication-readiness-roadmap.md`
(which holds the evidence, file paths, and measurements behind every line). 2026-08-15.

**Rule of thumb:** nothing below Tier 0 is worth doing before Tier 0 closes, and no GNN
number measured before the P0 land may appear in the manuscript.

## Tier 0 — correctness (local Mac, ~1 week total)

- [ ] **Decide #157**: which surface does the "never worse than baseline" invariant live on —
      the one users get (truth hidden, recommended) or the one `predict()` currently scores?
      → [github.com/itchyshin/pigauto/issues/157](https://github.com/itchyshin/pigauto/issues/157)
- [ ] Fix the calibration↔predict context mismatch per that decision; make the gate-acceptance
      guards assert on the delivered surface
- [ ] **Land P0** (`arc/p0-onto-main` `f5e8416` — merge already done and verified; needs a
      fresh G0 lock that day)
- [ ] Fix the **BLOCKING** vignette claim ("Exact under BM" → split by path) + 5 SHOULD-FIX
      → `docs/dev-log/arc/2026-08-15-p0-claim-gate-findings.md`
- [ ] `rcmdcheck --as-cran` clean (expect only the pre-existing 2 WARN / 3 NOTE)

## Tier 1 — GNN evidence (Totoro; D-139: estimate + pre-run test first)

- [ ] **Re-run all 8 per-type benches** on the corrected pipeline — every existing `.rds` is
      leak-tainted (May 30 / Jun 11 / Apr 28, pre-dating the Aug 8–9 fix)
- [ ] **Regime map (#135)**: GNN value vs signal × n × missingness × type + gate-firing
      diagnostics, multi-seed
- [ ] **Conformal coverage campaign**: coverage vs n per type (only number on record is
      0.884–0.887, leak-tainted, tiny n_val)
- [ ] **Then decide the paper's GNN framing** — "lift in mapped regimes" vs "safety-gated
      architecture". Not before.

## Tier 2 — latent multivariate

- [ ] **EM decision**: scope the claim to "single-pass Henderson-init joint MVN" (cheap) OR
      restore `max_iter > 0` under a dedicated numerical G0 with recovery sims
- [ ] **Known-Σ recovery simulation**: Σ recovery, per-tip SE coverage, liability-decode
      probability calibration (needs P0 landed first)
- [ ] **P1-8**: covariates silently ignored by joint/threshold-joint — thread through, or
      detect-and-message
- [ ] **P1-9** zi_count zeros excluded from val/test (blocks any zi_count metric)
- [ ] **P1-11** tree-MI conformal draws (or scope the tree-uncertainty claim)
- [ ] **P1-12** multi-obs phylo-signal gate no-op
- [ ] **P1-6 / P1-7 / P1-13** docs consolidation + `read_traits` count fix (small)
- [ ] **OVR categorical head-to-head** vs LP / multinomial / BACE (the 72% number is BACE's,
      not pigauto's)

## Tier 3 — packaging

- [ ] Head-to-head vs Rphylopars / BACE / missForest, post-P0, regime-labelled
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
