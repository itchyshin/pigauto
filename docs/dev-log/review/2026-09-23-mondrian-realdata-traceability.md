# Mondrian real-data traceability review

Reviewer: independent traceability check (read-only except this file).
Scope: NEWS.md sections "Mondrian conformal on real data: the default stays `"split"`"
and "Change: Mondrian-scaled conformal MI draws, with a caveat for phylogenetic GLS";
`useful/paper_section_draft.md` section 8.3 (including Table S-UQ2).

Method: all numbers were recomputed independently in fresh R scripts against
`docs/dev-log/mondrian-realdata/results_table.csv`,
`script/mondrian_confirmation/returned/**/mask_receipt.rds`,
`script/mondrian_confirmation/returned/mi-sim/rep_*.rds` (all 500 files), and
`docs/dev-log/mondrian-realdata/mi_se_summary.rds` (this file is plain-text markdown
despite its `.rds` extension). No repo script (`0*_*.R`, `1*_*.R`) was sourced. Table
S-UQ2 was recomputed per `results.md`'s stated pooling rule: coverage = n_test-weighted
mean over traits within dataset x arm x stratum; width ratio = median over traits of
`width_ratio`.

## Table

| claim | file:line | stated value | recomputed value | source | verdict |
|---|---|---|---|---|---|
| PanTHERIA n | NEWS.md:11 | n = 4,027 | 4,027 (nrow of `truth` in `pantheria-mcar-m20260818/mask_receipt.rds` and `pantheria-structured-m20260818/mask_receipt.rds`) | mask_receipt.rds | MATCH |
| PanTHERIA arms/masks | NEWS.md:11-12 | random and structured masks, 3 each | 6 receipt dirs: `pantheria-mcar-m2026081{8,9,20}`, `pantheria-structured-m2026081{8,9,20}` | `script/mondrian_confirmation/returned/` listing | MATCH |
| AVONET n | NEWS.md:12 | n = 1,500 | 1,500 (nrow of `truth` in `avonet-mcar-m20260818/mask_receipt.rds`) | mask_receipt.rds | MATCH |
| AVONET arm/masks | NEWS.md:12 | random mask, 3 | 3 receipt dirs: `avonet-mcar-m2026081{8,9,20}`; no `avonet-structured-*` dir exists | directory listing | MATCH |
| FishBase n | NEWS.md:13 | n = 10,484 | 10,484 (nrow of `truth` in `fishbase-structured-m20260818/mask_receipt.rds`) | mask_receipt.rds | MATCH |
| FishBase arm/masks | NEWS.md:13 | structured mask, 1 | 1 receipt dir: `fishbase-structured-m20260818`; `results.md` explicitly logs `fishbase-mcar-m20260818: not run` | directory listing + results.md | MATCH |
| Rule verdict | NEWS.md:14 | KEEP_SPLIT / default stays "split" | Per-trait near-stratum diffs (coverage_mondrian - coverage_split) fall below the pre-registered -2pp margin for multiple AVONET traits (Beak.Length_Culmen -0.0231, Mass -0.0341, Wing.Length -0.0278) and multiple FishBase traits (DepthRangeDeep -0.0266, Vulnerability -0.0251); condition 2 of the pre-registered rule therefore fails for those datasets, so the rule as written (00-preregistration.md, "flip only if all conditions hold on every dataset") cannot flip the default. Independently corroborated by a prior in-repo audit (`docs/dev-log/mondrian-realdata/review-m2-rule-audit.md`), which re-derived one Table-1 row from raw receipts to 4-6 decimals and found no deviation that flips the verdict on the current data. | results_table.csv (own aggregation) + review-m2-rule-audit.md | MATCH (qualitative; exact per-dataset p-values in Table 2 of results.md use a script formula not independently re-derived here, but the underlying per-trait/aggregate evidence supports the stated verdict) |
| far coverage split range | NEWS.md:16-17 | about 0.92-0.93 | split far coverage across the 4 dataset x arm cells: 0.919, 0.929, 0.930, 0.931 (range 0.919-0.931) | results_table.csv, own n_test-weighted aggregation | MATCH |
| far coverage mondrian range | NEWS.md:17 | about 0.94-0.96 | mondrian far coverage: 0.940, 0.949, 0.957, 0.964 | results_table.csv, own aggregation | MATCH |
| near-stratum narrowing | NEWS.md:18 | about 20% | near width_ratio (median over traits): 0.78, 0.81, 0.83, 0.84 i.e. 16-22% narrowing | results_table.csv, own aggregation | MATCH |
| near coverage split range | NEWS.md:18-19 | 0.97-0.98 | split near coverage: 0.966, 0.974, 0.978, 0.984 | results_table.csv, own aggregation | MATCH (approximate, as worded) |
| near coverage mondrian range | NEWS.md:19 | 0.95-0.97 | mondrian near coverage: 0.950, 0.958, 0.961, 0.968 | results_table.csv, own aggregation | MATCH |
| near non-inferiority margin | NEWS.md:20-21 | no more than 2 points below split's | matches 00-preregistration.md condition 2 verbatim ("at least split coverage minus 2 percentage points") | 00-preregistration.md | MATCH |
| MI sim reps/n | NEWS.md:36-37 | 500 replicates, n = 1000 | 500 rep files found in `mi-sim/`; all have `$n == 1000` | rep_1.rds .. rep_500.rds (own loop over all 500) | MATCH |
| MI GLS slopes | NEWS.md:37-38 | 0.35 and 0.32 against a true 0.70 | true_beta = 0.70 (all 500 reps); mean split estimate = 0.352 (bias -0.348); mean mondrian estimate = 0.320 (bias -0.380) | own aggregation over all 500 rep_*.rds `$split$estimate` / `$mondrian$estimate` | MATCH |
| structured mask description | paper_section_draft.md:329-341 | masking + propensity design; PanTHERIA both arms/3 masks, AVONET random-only, FishBase structured/1 mask | matches directory listing and 00-preregistration.md/Amendments | own check | MATCH |
| PanTHERIA n (paper) | paper_section_draft.md:338 | 4,027 mammals | 4,027 | mask_receipt.rds | MATCH |
| AVONET n (paper) | paper_section_draft.md:339 | 1,500 birds | 1,500 | mask_receipt.rds | MATCH |
| FishBase n (paper) | paper_section_draft.md:340 | 10,484 fishes | 10,484 | mask_receipt.rds | MATCH |
| mask fraction | paper_section_draft.md:334, 336 | 20% of observed cells (both arms) | 0.20 (avonet mcar), 0.20 (fishbase structured), 0.2002 (pantheria mcar and structured) | own computation, mask TRUE / observed-cell count in each mask_receipt.rds | MATCH |
| Table S-UQ2, PanTHERIA structured far | paper_section_draft.md:348 | split 0.919, Mondrian 0.940, width ratio 1.18 | 0.9189, 0.9401, 1.1775 | results_table.csv, own n_test-weighted coverage + median width_ratio | MATCH |
| Table S-UQ2, PanTHERIA structured near | paper_section_draft.md:349 | split 0.974, Mondrian 0.958, width ratio 0.84 | 0.9736, 0.9583, 0.8382 | results_table.csv, own aggregation | MATCH |
| Table S-UQ2, PanTHERIA random far | paper_section_draft.md:350 | split 0.931, Mondrian 0.957, width ratio 1.23 | 0.9308, 0.9565, 1.2318 (dataset=pantheria, arm=mcar) | results_table.csv, own aggregation | MATCH |
| Table S-UQ2, PanTHERIA random near | paper_section_draft.md:351 | split 0.978, Mondrian 0.968, width ratio 0.83 | 0.9777, 0.9679, 0.8298 | results_table.csv, own aggregation | MATCH |
| Table S-UQ2, FishBase structured far | paper_section_draft.md:352 | split 0.930, Mondrian 0.949, width ratio 1.14 | 0.9296, 0.9487, 1.1385 | results_table.csv, own aggregation | MATCH |
| Table S-UQ2, FishBase structured near | paper_section_draft.md:353 | split 0.966, Mondrian 0.950, width ratio 0.78 | 0.9662, 0.9501, 0.7766 | results_table.csv, own aggregation | MATCH |
| Table S-UQ2, AVONET random far | paper_section_draft.md:354 | split 0.929, Mondrian 0.964, width ratio 1.67 | 0.9291, 0.9637, 1.6657 (dataset=avonet, arm=mcar) | results_table.csv, own aggregation | MATCH |
| Table S-UQ2, AVONET random near | paper_section_draft.md:355 | split 0.984, Mondrian 0.961, width ratio 0.81 | 0.9844, 0.9610, 0.8078 | results_table.csv, own aggregation | MATCH |
| far coverage rise | paper_section_draft.md:359 | rises to 0.94-0.96 | mondrian far coverage across 4 rows: 0.940, 0.949, 0.957, 0.964 | Table S-UQ2 recomputation above | MATCH |
| near narrowing | paper_section_draft.md:360 | about a fifth | near width_ratio: 0.78, 0.81, 0.83, 0.84 -> 16-22% narrowing | Table S-UQ2 recomputation above | MATCH |
| near coverage floor | paper_section_draft.md:360-361 | stays at or above 0.95 | mondrian near coverage values 0.950, 0.958, 0.961, 0.968, all >= 0.95 | Table S-UQ2 recomputation above | MATCH |
| non-inferiority margin (paper) | paper_section_draft.md:361-362 | no more than two percentage points below split's | matches 00-preregistration.md condition 2 | 00-preregistration.md | MATCH |
| condition not met for AVONET, FishBase | paper_section_draft.md:362-364 | condition was not met for AVONET and FishBase | per-trait near-stratum diffs below -2pp for AVONET (Beak.Length_Culmen, Mass, Wing.Length) and FishBase (DepthRangeDeep, Vulnerability); PanTHERIA's per-arm aggregate near diffs (-0.0098 mcar, -0.0153 structured) stay within the 2pp margin | results_table.csv, own per-trait diff computation | MATCH (qualitative; same caveat as the NEWS verdict row above) |
| "Mondrian activated for every continuous trait" | paper_section_draft.md:341 | all traits activated (no fallback) | `fallback` column in results_table.csv is FALSE for all 38 rows across all 3 datasets (8 avonet, 10 fishbase, 20 pantheria) | results_table.csv | MATCH |

## Sample-size cross-check (mask_receipt.rds truth row counts vs quoted n)

| dataset | quoted n | mask_receipt.rds truth nrow | verdict |
|---|---|---|---|
| PanTHERIA | 4,027 | 4,027 | MATCH |
| AVONET | 1,500 | 1,500 | MATCH |
| FishBase | 10,484 | 10,484 | MATCH |

## MI simulation rep-file cross-check (500 replicates, n = 1000)

All 500 files `script/mondrian_confirmation/returned/mi-sim/rep_1.rds` .. `rep_500.rds`
exist and load. Every rep has `n = 1000`, `epochs = 500`, `m = 20`, `true_beta = 0.7`,
`split$failed = FALSE`, `mondrian$failed = FALSE` (0 failures in either arm across 500
reps). Mean `n_missing_x` = 300.9 (range 267-340), consistent with the ~30% MAR_phylo
target described in `useful/mondrian-mi-se-justification.md`. Aggregating the 500 reps
independently reproduced every number in `docs/dev-log/mondrian-realdata/mi_se_summary.rds`
and the results table in `useful/mondrian-mi-se-justification.md` ("## Result" section)
to 3-4 decimal places: split mean SE 0.0376, empirical SD 0.0550, SE ratio 0.6829 (stated
0.683); Mondrian mean SE 0.0372, empirical SD 0.0608, SE ratio 0.6111 (stated 0.611);
reference-arm bias +0.0008 (stated +0.001); split per-cell draw coverage 0.8456 (stated
0.846); Mondrian far-stratum draw coverage 0.8680 (stated 0.868), near-stratum 0.8959
(stated 0.896).

## Notes and limitations

- `results.md`'s Table 2 (per-dataset decision-rule statistics: `near_noninferiority_p`,
  `min_mondrian_far_cov`, etc.) is produced by
  `script/mondrian_confirmation/08_apply_decision_rule.R` /
  `12_build_results_doc.R`. Per this task's instruction not to source repo scripts, those
  exact per-dataset test statistics were not independently re-derived here; instead the
  qualitative claims that depend on them ("KEEP_SPLIT", "condition not met for AVONET
  and FishBase") were checked against a from-scratch reading of the raw per-trait
  coverage numbers in `results_table.csv`, which supports the same conclusion. A
  separate, more thorough prior audit of that decision-rule script
  (`docs/dev-log/mondrian-realdata/review-m2-rule-audit.md`) independently re-derived one
  Table-1 row from the raw `split.rds`/`mondrian.rds` receipts (matching to 4-6 decimals)
  and reviewed the rule logic itself, finding three script deviations from the
  pre-registration's per-dataset aggregation intent, none of which currently flip the
  verdict.
- `docs/dev-log/mondrian-realdata/mi_se_summary.rds` has a `.rds` extension but is a
  plain-text Markdown file, not a serialized R object; `readRDS()` fails on it. It was
  read as text and its numbers were independently reproduced from the raw
  `mi-sim/rep_*.rds` files rather than taken on faith.
- The NEWS.md coverage-range claims ("about 0.92-0.93", "0.97-0.98", "0.95-0.97") are
  explicitly hedged as approximate and use only 2 significant figures; the recomputed
  4-cell ranges (e.g. near split 0.966-0.984) span slightly beyond the literal 2-sig-fig
  wording at the edges (0.984 rounds to 0.98, 0.966 rounds to 0.97), which is consistent
  with "about" and was scored MATCH.

## VERDICT

0 mismatches found. All checked numbers in NEWS.md's two new sections and
paper_section_draft.md section 8.3 (including Table S-UQ2) trace to
`docs/dev-log/mondrian-realdata/results_table.csv`, `mask_receipt.rds` truth-row counts,
`mi-sim/rep_*.rds`, `mi_se_summary.rds`, and `useful/mondrian-mi-se-justification.md`,
and were independently recomputed to match.
