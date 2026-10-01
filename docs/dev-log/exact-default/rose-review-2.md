# Rose review 2: feat/exact-default (predict_method = "auto"), commits b29ac47..13bd3e6

Reviewer: Rose (fresh adversarial pass, D-43 own-the-verifier). Base: origin/main eccb298. Read-only review:
no code edited, no full suite or R CMD check run. Three small probes were run from the worktree with
`devtools::load_all()` (scripts in the session scratchpad, results quoted below).

## Verdicts

| item | verdict |
|---|---|
| 1. Split correctness (route half vs calibration half) | NOT-DONE |
| 2. Route reuse in multi_impute / multi_impute_trees / cross_validate / predict | DONE (minor notes) |
| 3. NEWS claims vs round-4 evidence | NOT-DONE (numbers trace; two claims overreach) |
| 4. User-facing trip hazards | NOT-DONE (stale docs; one silent regression) |

## Blocking findings (ranked)

### B1. The split is done on latent cells, not on trait rows, so multi-column traits leak and shrink

`.pigauto_split_route_score()` (R/fit_baseline.R:1483-1524) samples individual latent cells from
`val_idx[val_col %in% tm$latent_cols]`. But `make_missing_splits()` (R/mask_missing.R:176-186, via
`expand_trait_idx_to_latent()`) holds out whole trait rows, so a categorical row contributes K cells and
a zi_count row contributes 2 (gate and magnitude, the magnitude NA for zeros). Consequences:

- Route choice (`.pigauto_route_val_loss()`, R/fit_baseline.R:1395-1408) scores a categorical row if ANY
  of its K cells is in the route half. Gate calibration (`calibrate_gates()`, R/fit_helpers.R:260) keeps
  a row only if its FIRST latent column is in the score half. The same row therefore often informs both.
- Measured on the bundled AVONET 300 (`impute(..., gnn = FALSE, seed = 1)`, 3 s probe): Trophic.Level
  (K = 4, 29 held-out rows) calibrates on 10 rows, 8 of which also chose the route; Primary.Lifestyle
  (K = 5, 43 rows) calibrates on 22 rows, 20 of which also chose the route. A synthetic zi_count with 20
  rows: 12 calibration rows, 7 shared with the route half. Single-column traits are correctly disjoint
  (Beak.Length_Culmen 25 / 24, overlap 0).
- The 38-cell threshold counts cells, so a K = 4 categorical splits at 10 rows. On AVONET 300 this raised
  a new user-visible warning, "Small validation set for 1 trait(s): Trophic.Level (n=10)", for a trait
  with 29 held-out rows that main would calibrate on in full.
- zi_count conformal scores use the magnitude column only (R/fit_helpers.R:858), and that column is NA
  for observed zeros, yet those NA cells count toward the 38. The "at least 19 conformal cells" rationale
  in the code comment (R/fit_baseline.R:1496-1502) does not hold for zi_count.
- multi_proportion is never route-chosen (`.pigauto_route_val_loss()` returns n = 0, R/fit_baseline.R:1418-1421),
  yet its K-cell rows are still split, so gate calibration loses rows for no decision at all.

This also falsifies the NEWS sentence "so no cell is used twice" at the level that matters (the trait
row). Tests (g) and the 38-cell test (tests/testthat/test-route-choice.R:253-296, 432-450) check cell
disjointness only, which is why they pass. Fix: split trait rows (the `val_idx_trait` space, or unique
rows per trait), count rows (and for zi_count, non-NA magnitude rows) against the threshold, and skip
multi_proportion. Add a test on row-level disjointness with a categorical and a zi_count trait.

### B2. "auto" splits even when there is nothing to choose, which silently degrades conformal intervals

The split runs before, and independently of, the route comparison (R/fit_baseline.R:1603-1607). When
both candidate fits are identical for a trait (a single trait, a trait the exact route cannot reach,
exact fallback because K < 2, label-propagation-only discrete traits), half the validation cells are
still removed from calibration. Probe: AVONET 300, `Mass` only, 60 NA, `gnn = FALSE`, seed 3. Baseline
`mu` identical to `predict_method = "per_column"` (max abs difference 0), yet `score_val_n` = 30 of 60 and
the conformal score changes from 1.180 (per_column) to 1.929 (auto), a 63% wider interval from one
draw. Every single-trait user and every trait the joint path never touches pays this under the new
default with no benefit. None of the round-4 cells exercise it (all are 7-trait mixed designs). Fix: only
split a trait when the two candidate fits differ on its columns (or re-merge the halves after a tie),
and add a regression test that `auto` and `per_column` give identical `conformal_scores` on a single
trait.

### B3. NEWS overclaims on the real-data result

The numbers trace exactly (below), but the sentence "On 13 real-data cases, z-RMSE falls in 12 (... by
0.1 to 2.5% elsewhere)" omits that the real-data runs are 5 seeds, continuous traits only, GNN off
(`lambda_datasets_auto4_report.txt:3`), with no Monte Carlo error reported. Changes of 0.1 to 0.4%
(LepTraits 300, AmphiBIO 300, GlobTherm 300, BIEN 300, LepTraits 2000) cannot be called falls at 5 seeds.
Required wording: state the regime (5 seeds, continuous traits, without the GNN) and say that 5 cases
improve clearly (AVONET x3, PanTHERIA x2, by 17 to 51%), the rest change by less than 3% either way.

## Non-blocking findings (ranked)

N1. Stale roxygen and man. `fit_baseline()` `predict_method` still says the split is at "10 or more
validation cells" and that below 10 the route "falls back to exact" with all cells kept for calibration
(R/fit_baseline.R:152-163, man/fit_baseline.Rd:150-151; also the internal comment at R/fit_baseline.R:1474-1477).
Current behaviour: split at 38, and below 38 the SAME cells both choose and calibrate (chooser defaults
to exact only under 5). `impute()` and `fit_pigauto()` `lambda_mode` text says `predict_method = "exact"
(the default)` (R/impute.R:71, same block in R/fit_pigauto.R). `multi_impute()` `lambda_mode` still says
discrete traits stay at lambda = 1 (R/multi_impute.R:74-75), true only under per_column. README.md:157-158
and vignettes/getting-started.Rmd:473-477 still present "exact" as an opt-in technical control without
saying the default now chooses per trait.

N2. The default user path is unmeasured. Every round-4 cell and real-data case is `gnn = FALSE`, while
`impute()` defaults to `gnn = TRUE`. Under the GNN the split halves the cells the discrete gate floor
(2 / 1 absolute cells, R/fit_helpers.R:389-393) and the gate grid search see, which is the safety the
project instructions warn about. NEWS is honest that the evidence is GNN-off; a GNN-on check on a few
core cells is advisable before release.

N3. Coverage context. At lambda = 1, n = 100 the conformal coverage is already well under nominal on
main (c1 0.889, cnt 0.868) and auto lowers it further (c1 0.882, cnt 0.851). NEWS reports only the
changes; the absolute level should be stated, since a reader will assume near 0.95. Also, "coverage
rises by up to 0.010 at lambda 0.3 and 0.7" is the c1/c2 mean; the count trait falls at lambda 0.7,
n = 100 (-0.003).

N4. Mixed aggregation in the accuracy sentence. "Rises by up to 0.050" is per trait averaged over rho
(ord, lambda 0.7, n 1000: +0.0496); "the largest fall, 0.009 for cat3" is a single rho = 0.5 cell. Both
are defensible, but say which. The cat3 falls are small but consistently negative at lambda 1 and at
lambda 0.7, n 300 (rho-averaged -0.006 to -0.008), so "within Monte Carlo error" is true per cell but a
consistent sign across cells deserves a word.

N5. Runtime claim. "1.4 to 1.5 times as long as a single fit" comes from two unreplicated
`fit_baseline()` timings (S5b-route-choice-report.md:225-231, 1.42x and 1.50x). It is the baseline cost,
not `impute()` wall time. Reword to "the baseline fit takes about 1.4 to 1.5 times as long (two timings)".
`multi_impute_trees()` pays this per tree.

N6. Warning under the default in the production refit. `.fit_baseline_route()` always passes
`predict_method_explicit = TRUE` (R/fit_baseline.R, the `fits[[r]] <- .fit_baseline_core(...)` loop), so
under `gnn = FALSE` a default user whose replayed "exact" route falls back in the `splits = NULL` refit
gets a `warning()` every time rather than the once-per-session `message()` NEWS promises. Rare, but it
contradicts the stated contract.

N7. No guard that `baseline$score_val_idx` is a subset of `splits$val_idx`
(`.pigauto_calibration_val_idx()`, R/fit_pigauto.R, near line 1764). A user passing a baseline fitted
with other splits into `fit_pigauto()` would calibrate on cells that may be training cells. An
`intersect()` with a warning would close it.

N8. Evidence file missing. S5-benchmark-round4.md:48 cites "R CMD check in S5d-check.log"; that file is
0 bytes. The suite log (FAIL 0 / PASS 2724) is present. Do not cite the check until the running check
lands.

N9. Multi-obs: the split is at observation level, so two observations of the same species can land in
different halves and the shared species-level `mu` means the route half is not fully independent of the
score half. Weaker than B1; worth a sentence in the docs.

N10. `multi_impute_trees()` shared-GNN path seeds each tree's split with `seed + t`
(R/multi_impute_trees.R:583), so tree 1 does not reproduce the reference fit's own split under
`gnn = TRUE`. Harmless, but surprising if someone compares them.

## Item 2 detail (DONE)

- `impute()` passes its "auto" baseline into `fit_pigauto()`, which uses `baseline$score_val_idx` in
  both the GNN-on (R/fit_pigauto.R:1306-1307) and GNN-off (R/fit_pigauto.R:616-617) calibration paths.
- The `gnn = FALSE` production refit in both `impute()` (R/impute.R:551-567) and `fit_pigauto()`
  (R/fit_pigauto.R:566-582) replays `baseline$predict_method_by_trait` through `predict_route`.
- `predict.pigauto_fit()` never refits the baseline; it uses the stored one, so a saved fit reproduces
  its routes.
- `multi_impute()` goes through `impute()` (R/multi_impute.R:224, 259) and inherits the default.
- `multi_impute_trees()` shared-GNN path: `gnn = FALSE` replays the reference routes per tree; `gnn =
  TRUE` re-decides per tree on `splits_ref`, which is the documented intent. The explicit request is
  carried through `dots`, since `model_config` does not store `predict_method` itself.
- `cross_validate()` calls `fit_pigauto()` per fold with fold splits; auto runs per fold. B1 applies here
  too (fold val sets hold whole categorical rows).

## Item 3 trace (numbers recomputed from the CSVs)

Simulation (main vs auto4, merged on all design keys, rho averaged):
- z-RMSE, mean over c1, c2, cnt, prp: lambda 0.3: -3.1 / -4.6 / -6.6%; lambda 0.7: -6.1 / -8.5 / -9.4%;
  lambda 1: -1.2 / -0.3 / +0.4%. Matches round 4 and NEWS.
- Coverage, mean of c1 and c2: +0.003 / +0.004 / +0.009; +0.008 / +0.005 / +0.010; -0.007 / -0.003 /
  +0.002. Matches.
- Discrete accuracy, mean of bin, cat3, ord: +0.009 / +0.014 / +0.022; +0.012 / +0.025 / +0.032;
  -0.002 / -0.002 / -0.001. Matches.
- 10 of 36 rho-averaged trait-cells below main on coverage; cnt -0.017 and prp -0.013 at lambda 1,
  n 100. Matches.
- cat3 largest per-cell fall -0.0094 (z -0.70) and -0.0093 (z -1.47). Matches "at most -0.009, |z| < 1.5".

Real data (column `estimate`, main vs auto4): AVONET 300 bundled -50.5%, AVONET 300 -35.8%, AVONET 2000
-33.0%, PanTHERIA 300 -17.4%, PanTHERIA 2000 -16.8%, AmphiBIO 2000 -2.5%, BIEN 2000 -0.9%, BIEN 300
-0.4%, LepTraits 2000 -0.4%, AmphiBIO 300 -0.2%, GlobTherm 300 -0.2%, LepTraits 300 -0.1%, GlobTherm 1969
+0.8%. All match the round-4 table and NEWS; the objection is scope (B3), not arithmetic.

## What this review does not cover

It does not verify the exact-conditional algebra in R/exact_conditional.R or R/joint_mvn_solver.R (the
first Rose review and S1/S4 reports cover that), does not run the suite or R CMD check, and does not
measure GNN-on behaviour.
