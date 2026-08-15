# Claim-gate audit: pigauto `origin/main` (416561b) vs the four unmerged P0 honesty fixes on `fix/ci-install-libtorch`

Audited read-only via `git show origin/main:<path>`. Working tree was not read (it is on a
different branch). The four P0 fixes are taken from commit `58708ac` ("Fix P0 review blockers
and GNN train/cal symmetry") on `origin/fix/ci-install-libtorch`, whose `NEWS.md` diff and
`R/henderson_s_inv.R` / `R/preprocess_traits.R` / `R/fit_pigauto.R` diffs I read directly to
confirm what each fix actually changes, rather than trusting the fix names alone.

**What the four fixes actually do (grounded in the diff, for reference):**

1. **Per-tip BM SE** — under the *shipped default* `max_iter = 0` (cross-trait EM not restored),
   the Henderson-init SE used by the K≥2 joint-MVN and threshold-joint baselines (`fit_baseline()`
   dispatch: fires whenever `Rphylopars` is installed and ≥2 BM-eligible latent columns exist —
   this is the **default** path, not an opt-in flag) recycled **one empirical SD for every missing
   tip**, independent of phylogenetic position. The fix replaces it with the real per-tip
   conditional variance from the same sparse Cholesky already used for the mean. This SE also
   feeds `decode_binary_liability()`'s `Φ(μ/√(1+σ²))`, so pre-fix it degrades binary/ordinal
   *predicted probabilities*, not only reported interval widths.
2. **Covariate alignment** — pre-fix, `preprocess_traits()` only checked `nrow(covariates) == n_obs`,
   not species/observation identity. A covariate table in a different row order than `X_scaled`
   (after tree-tip reordering) is silently mispaired with no error.
3. **zi_count conformal MI** — pre-fix, `multi_impute(draws_method = "conformal")` skipped
   held-out non-zero scores for `zi_count` and drew from `pred$se` of `E[X]` (which conflates gate
   uncertainty with count SD). Post-fix it scores the magnitude latent directly.
4. **GNN train/cal/predict symmetry** — pre-fix, the DAE's training-time encoder input used `t_X`
   (full truth) rather than `t_X_eval` (held-out cells replaced by baseline) as its base tensor, so
   val/test truth was visible as *context* (not as a loss target) during training. Gate calibration
   and conformal scoring also ran a single forward pass rather than the same iterative
   `refine_steps` loop `predict()` uses. Both are fixed. This means any benchmark number that
   depended on a full `fit_pigauto()` GNN training run pre-fix was generated on a training loop with
   this held-out-context leak; direction and magnitude of its effect on reported numbers is unknown
   without re-running.

---

## DESCRIPTION

| Claim | Location | Why it is a problem | Severity | Suggested wording |
|---|---|---|---|---|
| "Provides conformal prediction intervals for continuous, count, and ordinal traits" | DESCRIPTION, `Description:` field, line ~16 | Under-claim, not over-claim: `R/fit_helpers.R::compute_conformal_scores()` on `origin/main` dispatches on `tm$type %in% c("continuous","count","ordinal","proportion")` — `proportion` traits **do** get conformal intervals, but DESCRIPTION's list omits them. Not related to the P0 fixes; a pre-existing listing gap. `zi_count` is correctly omitted (code explicitly skips it via `next`), consistent with fix #3 above not yet existing on main. | NOTE | "Provides conformal prediction intervals for continuous, count, ordinal, and proportion traits" |

Nothing else in DESCRIPTION is contradicted by the four fixes; the "Stochastic graph-network and
posterior-tree completions are prediction diagnostics rather than validated inferential
imputations" sentence is already appropriately hedged and does not depend on any of them.

## README.md

| Claim | Location | Why it is a problem | Severity | Suggested wording |
|---|---|---|---|---|
| "Covariates must be fully observed. Numeric columns are z-scored; factor/ordered columns are one-hot encoded automatically." (no mention of row/species-identity alignment) | README.md, "## Using environmental covariates", lines 87–106 | Absence, not falsehood: on `origin/main`, `preprocess_traits()` matches `covariates` to traits by row count only, not by name/species identity (fix #2 is unmerged). The README's own worked example is safe (covariates and traits are both subset from the same `ctmax_sim` object, preserving row order), but the surrounding prose gives no warning that a covariate table in a different order than the trait table is silently mispaired with no error — a plausible mistake for a user building covariates from an external table (e.g. joined/sorted by species name). This is a "lacking stated scope" gap, not a false statement. | SHOULD-FIX | Add one sentence: "Covariate rows are currently matched to trait rows by position, not by species name — supply covariates already in the same row order as `traits` (e.g. by subsetting the same source table), or verify alignment yourself before calling `impute()`." (Revise once the covariate-alignment fix lands, since named alignment will then be enforced with an error on mismatch.) |

Everything else in README.md is already carefully hedged: the top-of-file WARNING box ("no cell's
interval coverage is certified, and covariance routes have focused-test evidence only"), the
"Caveats from multi-seed evidence" section (explicitly scoped by dataset/N/seed count), and the
"Benchmarks" section's refusal to print headline numbers are all consistent with what the P0 fixes
would change. No BLOCKING items found here.

## _pkgdown.yml

Clean. The `home:` sidebar status text and `description:` field duplicate the same hedged language
already in README.md's WARNING box ("Point estimates are the supported claim; no cell's interval
coverage is certified, and covariance routes have focused-test evidence only") and are not
contradicted by any of the four fixes. The `reference:` section's descriptions for
`multi_impute`/`multi_impute_trees` ("These paths failed or were outside the downstream inferential
gate and must not be used as analysis-aware multiple imputations") and for
`multi_impute_analysis` are consistent with NEWS.md's dated campaign results. No table entries.

## NEWS.md

Audited the current-development section (`# pigauto 0.10.0.9000`) and the most recent release
section (`# pigauto 0.10.0`) in full; spot-checked the remaining ~2,600 lines of historical
per-version entries for the specific terms the fixes touch (Henderson, joint MVN, threshold-joint,
per-tip, zi_count, covariate). I did **not** exhaustively re-read all historical entries — they are
dated record-of-what-shipped-then, and re-litigating every past release note is out of scope for a
current-claims audit. Flagging this as a scope limitation on this section rather than a clean bill.

| Claim | Location | Why it is a problem | Severity | Suggested wording |
|---|---|---|---|---|
| "a 60-replicate ablation ... confirmed `pigauto_OFF` beats both `pigauto_ON` and `pigauto_ON_L0` ... z-RMSE 1.038 vs 1.056 vs 1.057; conformal coverage stable at 0.884–0.887" (also mirrored in `gnn-architecture.Rmd`, see below) | NEWS.md, `v0.9.3` "within-row cross-trait attention" entry | This ablation depended on a full `fit_pigauto()` GNN training run, which pre-fix (fix #4) trained with held-out val/test truth visible as encoder context. The number is not necessarily wrong — the direction of any leakage effect on a null-result ablation (GNN doesn't help either way) is unclear and could plausibly be conservative — but the entry states the result as settled without disclosing that the training pipeline it ran on has since been found to leak held-out context. Same caveat applies to any other dated entry whose measurement pipeline included `fit_pigauto()` GNN training (I did not enumerate all of them). | SHOULD-FIX | Add a forward note once P0 lands: "Numbers in this and earlier entries were measured before the 2026-08 GNN train/cal-symmetry fix (held-out context leak in training); not re-verified against the corrected pipeline." |

Nothing else in the current dev/most-recent sections is contradicted — the BACE `final_imp`
restoration entry and the OVR-vs-BACE documentation-only entry both explicitly disclaim any
performance comparison ("No pigauto-versus-BACE performance claim is made here").

## vignettes/getting-started.Rmd

| Claim | Location | Why it is a problem | Severity | Suggested wording |
|---|---|---|---|---|
| "### Assembling covariates from public data" worked example: `clim <- pull_worldclim_per_species(...)`; `result <- impute(my_traits, my_tree, covariates = clim)` with no rowname/species-identity step shown | getting-started.Rmd, lines 467–501 | Same gap as the README covariate section, more concretely exposed here: this example builds a covariate table from an *external* source (GBIF/WorldClim) rather than subsetting the trait object itself, so there is no guarantee the helper's output row order matches `my_traits`' row order (e.g. a species dropped for zero GBIF occurrences would silently shift the row alignment for all species after it, while `nrow` might still coincidentally match if another species also drops out). I have not read `pull_gbif_centroids()`/`pull_worldclim_per_species()`'s source to confirm whether they guarantee row-for-row correspondence with the `species=` argument order — flagging this as an inference from the documented interface, not a confirmed defect. | SHOULD-FIX | Add: "Both helpers return one row per requested species in the order given to `species=`, with no rows dropped for species that returned zero occurrences [if true — verify against the helper's own contract]. `covariates` is currently matched to `traits` by row position; confirm this ordering holds for your species list before calling `impute()`." |

Nothing else here is contradicted by the fixes — `suggest_next_observation()`'s variance-reduction
formula section, the phylo-signal-gate section, and the tips/next-steps section are all
independent of the four fixes.

## vignettes/mixed-types.Rmd

| Claim | Location | Why it is a problem | Severity | Suggested wording |
|---|---|---|---|---|
| "The SE matrix provides type-appropriate uncertainty" | mixed-types.Rmd, line 144, immediately before `head(pred$se)` | Vague but still an affirmative accuracy claim, and the synthetic example it sits above (`mass` continuous + `clutch` count = 2 BM-eligible latent columns) is exactly the shape that triggers the K≥2 joint-MVN baseline when `Rphylopars` is installed — the default path carrying the per-tip-SE bug (fix #1). Less specific than the `gnn-architecture.Rmd` table claim below, so lower severity, but it is the same underlying defect reaching a second, more example-driven audience. | SHOULD-FIX | "The SE matrix provides type-appropriate uncertainty (see the [uncertainty-mechanisms table](gnn-architecture.html#5-uncertainty-quantification) for what each column's SE does and does not guarantee)." |

## vignettes/gnn-architecture.Rmd

This is the file the P0 fixes bite hardest, because it is the one place the package states its
SE/coverage claims precisely enough to be falsifiable.

| Claim | Location | Why it is a problem | Severity | Suggested wording |
|---|---|---|---|---|
| `pred$se` (cont./count/ordinal/prop.) — "Conditional-MVN SD from the BM baseline, delta-method back-transformed. **Exact under BM, model-dependent.**" | gnn-architecture.Rmd, §5 "Uncertainty quantification" table, line 349 | Unconditionally false for the default multi-trait path under current `origin/main`. §2.3 of this same vignette states the joint-MVN baseline is used "when `Rphylopars` is installed and ≥ 2 continuous-like latent columns exist" — i.e. by default for any mixed dataset with 2+ BM-eligible traits and the (Suggests-only, commonly installed) `Rphylopars` package. Under the shipped default `max_iter = 0`, that path's SE (pre-fix #1) is a single empirical SD recycled across every missing tip — not phylogenetically informed at all, and specifically *not* "exact under BM" by any reading, since it ignores each tip's actual conditional variance. The table gives no scope carve-out between the single-column `bm_internal.R` path (where "exact under BM" is a defensible claim) and the Henderson-init joint/threshold-joint path (where, pre-fix, it is not). This also isn't only a display/interval-width issue: §2.4 shows this same SE feeding `decode_binary_liability()`'s `Φ(μ/√(1+σ²))`, so the defect degrades *predicted binary/ordinal probabilities* wherever the joint/threshold-joint path fires, not merely `pred$se`'s reported width. | BLOCKING | Split the row: "`pred$se` (cont./count/ordinal/prop., **single-column BM path**) — exact under BM. `pred$se` (cont./count/ordinal/prop., **joint-MVN/threshold-joint path**, fires by default when ≥2 BM-eligible columns + Rphylopars) — conditional-MVN SD from a Henderson-init approximation; per-tip accuracy has [not yet been / now been, once #1 lands] verified against the closed-form GLS SE." |
| "$\hat{y}_i \pm s_t$ ... giving $\ge 95\%$ marginal coverage on the original scale." | gnn-architecture.Rmd, §7, line 402 | Not false on its own terms (split-conformal marginal coverage under exchangeability is a real, correctly-stated guarantee), but it drops the caveat the vignette itself states three sections earlier in §5's table ("... under the split-conformal exchangeability and fixed calibration-procedure assumptions; not an MI draw distribution"). Restating the headline number without the caveat in the "worked equation" summary risks a reader quoting §7 in isolation as an unconditional claim. Independent of the P0 fixes, but worth fixing alongside them since this section is being touched anyway. | NOTE | "giving ≥95% marginal coverage on the original scale, under the split-conformal exchangeability assumptions stated in §5" |
| "a 60-replicate ablation ... z-RMSE 1.038 vs 1.056 vs 1.057; conformal coverage stable at 0.884–0.887" | gnn-architecture.Rmd, §3.5, lines 300–307 | Same issue as the NEWS.md entry above (duplicated text) — measured on a `fit_pigauto()` GNN training pipeline that pre-fix (#4) leaked held-out truth into training-time encoder context. Direction/magnitude of any effect on this specific null-result ablation is unknown. | SHOULD-FIX | Add: "(measured before the 2026-08 GNN train/cal-symmetry fix; not yet re-verified against the corrected training pipeline)" |
| "33.7% RMSE lift on simulated correlated BM data" (joint MVN baseline bench) | gnn-architecture.Rmd, §2.3, line ~124 | This bench (`script/bench_joint_baseline.R`) compares the baseline only (no GNN training loop involved), so it is **not** touched by fix #4 (GNN symmetry). It *is* potentially touched by fix #1 (per-tip SE) only if the bench's metric depends on SE rather than point mu — RMSE on point predictions (mu) is unaffected by the SE bug, since the fix only changes the `se` output, not `mu`. I checked this and believe it is **not** a defect, but flagging the reasoning explicitly since it sits one section above the item that is. | NOTE (no action needed) | — |

## vignettes/tree-uncertainty.Rmd

Clean. This vignette is already maximally hedged ("What this does not establish": "Variation across
these datasets is not a calibrated standard error, confidence interval, or fraction of missing
information") and none of its claims depend on the specific mechanics the four fixes change — it
never states an SE/coverage number, only that the tool is descriptive-only.

## vignettes/common-pitfalls.Rmd

Clean. None of its five FAQ entries (imputed-vs-observed-cells confusion, ordinal majority-class
collapse, closed-gate behaviour, phylogenetic-signal diagnosis, mass outlier clamping) touch SE,
covariates, zi_count, or the GNN train/cal loop in a way the P0 fixes would change.

---

## Verdict

**Is `main`'s public surface honest as it stands?** Mostly yes, and unusually so for a package at
this stage — the top-level WARNING box, the refusal to print headline benchmark numbers in the
README, and the three-way SE/conformal/MI-draws distinction in `gnn-architecture.Rmd` §5 are all
evidence of a documentation culture that already anticipates most of what an adversarial audit
would ask for. The one clear exception is the `gnn-architecture.Rmd` §5 table's unqualified "Exact
under BM, model-dependent" claim for `pred$se`, which is currently false for the default multi-trait
baseline path (joint-MVN/threshold-joint, active whenever `Rphylopars` is installed and ≥2
BM-eligible columns exist) and also degrades binary/ordinal decoded probabilities, not just
reported SE width. That is the one BLOCKING item.

**What specifically must change if the P0 fixes land:**

1. **Once fix #1 (per-tip BM SE) lands**, the `gnn-architecture.Rmd` §5 table's "Exact under BM"
   row becomes accurate for the joint/threshold-joint path too — but only if it stays scoped to
   the specific claim the fix delivers (per-tip conditional variance from the Henderson Cholesky),
   not silently extended to imply the un-restored cross-trait EM (`max_iter > 0`) is now default.
   Update the row (or split it) rather than leaving the old text as retroactively "now true."
2. **Once fix #2 (covariate alignment) lands**, the README and getting-started.Rmd covariate
   sections' current silence on row-identity matching stops being a silent-failure risk (mismatches
   will error instead of mispairing) — at that point the SHOULD-FIX suggestions above should be
   revised from "warn the user to check manually" to "named alignment is now enforced; unnamed
   covariates still follow input-row order," matching the fix's actual NEWS.md language.
3. **Once fix #3 (zi_count conformal MI) lands**, no current doc claim needs correcting — DESCRIPTION
   and the vignettes are already silent on zi_count conformal coverage, which was correct given the
   pre-fix code path skips it. Fix #3 is an opportunity to *add* a claim (zi_count conformal MI is
   now supported), not to remove a wrong one.
4. **Once fix #4 (GNN train/cal symmetry) lands**, the historical benchmark numbers flagged
   SHOULD-FIX above (the v0.9.3 within-row-attention ablation, duplicated in NEWS.md and
   gnn-architecture.Rmd §3.5) should be re-run and either confirmed or revised — they were measured
   on a training pipeline with a now-documented held-out-context leak, and no doc currently
   discloses that provenance caveat.

Overall: one BLOCKING defect (the SE-exactness table row), five SHOULD-FIX items (covariate-alignment
silence in two files, the vague SE claim in mixed-types.Rmd, and the two pre-fix-pipeline benchmark
citations in NEWS.md/gnn-architecture.Rmd), and two NOTE-level items (DESCRIPTION's conformal-type
list, and the §7/§5 coverage-caveat inconsistency). `_pkgdown.yml`, `tree-uncertainty.Rmd`, and
`common-pitfalls.Rmd` are clean.

**Uncertainty in this audit:** the getting-started.Rmd covariate-helpers finding rests on an
inference about `pull_gbif_centroids()`/`pull_worldclim_per_species()`'s row-correspondence
contract that I did not verify by reading their source (out of scope for a doc-claims audit of the
four named files/vignette set) — flagged as SHOULD-FIX with that caveat rather than BLOCKING. The
NEWS.md historical-entries scope was spot-checked, not exhaustively read line-by-line across all
~2,600 lines of past-release record; I would not be surprised if a fuller pass turned up one or two
more instances of the same pre-fix-pipeline benchmark-citation pattern.
