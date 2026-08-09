# Brain-write PROPOSAL — pigauto BACE wrap decision

**Status: DRAFT / NOT WRITTEN.** Per D-37 nothing goes into `~/shinichi-brain` without
explicit approval. This file stages the entry; Shinichi says yes or edits, then it lands.

Target: `~/shinichi-brain/memory/DECISIONS.md` (new entry). The Phase 0.25 sweep confirmed
`DECISIONS.md` currently has **zero** `bace` hits, so this is a genuine gap, not a duplicate.

---

## Proposed entry

> **D-xx (2026-08-09) — pigauto's BACE wrapper gets an opt-in proper-MI path; the cheap
> chain-average path stays the default.**
>
> `fit_baseline_bace()` gains `final_imp = FALSE` (default) and `n_final = 15L`. With
> `final_imp = TRUE` it appends `BACE::bace_final_imp()` to the `bace_imp()` chain it already
> runs and builds `mu` / `se` from those `n_final` independent draws.
>
> **Why.** The chain datasets are successive sweeps of one chained-equations chain —
> autocorrelated, still carrying convergence transient — so the between-dataset SD pigauto
> was reporting as `se` had no coverage property. Measured on simulated BM (n=100, 3 traits,
> 5 seeds): empirical 95% coverage **0.672** on the default path vs **0.940** with
> `final_imp = TRUE`, against a 0.95 nominal. RMSE also improved slightly (0.227 → 0.216).
> Cost is ~3.9x runtime.
>
> **What this is not.** It is not a fix for the imputed-as-observed defect Dan Noble repaired
> upstream in August 2026. That defect lived in BACE's `bace_final_imp()`, which pigauto had
> never called, so no pigauto output was ever affected by it. Framing this change as a bug
> fix would double-count Dan's work.
>
> **`n_final = 15` is Ada's default, not Shinichi's choice.** It matches Dan's Study B.
> BACE's own default is 50 and gives better tails. Revisit deliberately, not by drift.
>
> **Known caveat.** The final phase refits MCMCglmm per draw and is less robust than the
> chain: 1 of 5 simulated seeds hit "Mixed model equations singular". The wrapper errors with
> attribution rather than silently returning chain averages.
>
> **Unresolved at time of writing.** `origin/main` deleted `R/fit_baseline_bace.R` in
> `b615579` (2026-07-10, v0.10.0 CRAN release surface). The wrap exists only on
> `handover/2026-08-09-cursor`. Whether it should exist at all is an open question that this
> decision does not settle.

---

## Also worth a brain note (separate, smaller)

> **pigauto's OVR categorical path is not BACE's.** pigauto runs K threshold-joint
> Rphylopars fits so phylopars stays full-rank. BACE's `ovr_categorical` default became
> `FALSE` (true multinomial) in 2026-08 because OVR binaries hit quasi-separation in its
> one-observation-per-species MCMCglmm regime. Shared name, different estimator, different
> motivation. Do not port BACE's default into pigauto without a pigauto-side measurement.

---

## Not proposed

No entry claiming a pigauto-vs-BACE performance result. No such comparison was run — S3
compared two pigauto wrapper configurations against each other, not pigauto against BACE.
