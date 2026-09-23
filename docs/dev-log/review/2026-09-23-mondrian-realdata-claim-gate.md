# Mondrian real-data confirmation: public-claim gate (Rose)

Reviewer: adversarial claim-vs-evidence gate, read-only except this file.
Worktree `pigauto-mondrian-realdata`, branch `arc/mondrian-realdata`.
HEAD reviewed: `ef88f13` (HEAD moved twice during the review, from `f682073` through
`5709a7f` to `ef88f13`; every line number below is against `ef88f13`). Working tree clean
apart from one untracked log (`returned/fishbase-structured-m20260818/run-1thread.log`).

Files read: `00-preregistration.md`, `RUNLOG.md`, `results.md`, `results_table.csv`,
`review-m2-rule-audit.md`, `review/2026-09-23-mondrian-realdata-traceability.md`,
`NEWS.md:1-43`, `useful/paper_section_draft.md:240-373`,
`useful/mondrian-mi-se-justification.md` (whole), the after-task report, the harness
scripts `01`, `02`, `06`, `08`, `12`, `13`, `R/fit_helpers.R:895-980`,
`R/multi_impute.R` and `R/predict_pigauto.R` diffs against `origin/main`, and every
receipt under `script/mondrian_confirmation/returned/` (opened in R: field names,
trait lists, stratum sizes, fallback flags, mask columns, propensity summaries, file
mtimes). I re-ran `08_apply_decision_rule.R` on the committed CSV and re-derived the
pooled table, the per-dataset conditions, and the condition-2 p-values in my own R.

Numbers: I agree with the traceability review that every number in NEWS and Table S-UQ2
traces to `results_table.csv`. The problems below are about wording, timeline, and
provenance, not arithmetic.

## Timeline evidence (question 1)

| event | time (MDT) | source |
|---|---|---|
| pre-registration committed | 05:31:30 | `b36a90b` |
| decision-rule script with the condition-2 statistic (unpaired pooled-binomial z, Holm) | 05:59:55 | `7af133a` |
| first mask receipt written (pantheria-mcar-20260818) | 06:01:28 | file mtime |
| Amendment 1 committed | 06:01:35 | `e2bfb33` |
| first method (outcome) receipt (pantheria-mcar-20260818 mondrian) | 06:46:41 | file mtime; RUNLOG:28 records that its stratum sizes were read |
| FishBase structured mask receipt with propensity summary | 06:49:15 | file mtime |
| AVONET method receipts (3 masks, both methods) on fir | 06:55 to 07:04 | file mtimes |
| Amendment 2 committed | 07:12:46 | `b0c8a88` |
| results-doc generator and rule on the CSV (Amendment 2 filter only; statistic unchanged) | 07:26:35 | `11a179a`, verified by `git diff 7af133a 11a179a` |
| first structured-arm method receipt (pantheria-structured-20260818 mondrian) | 07:46:15 | file mtime |
| PanTHERIA and AVONET receipts committed | 08:59:34 | `25c7f15` |
| results, NEWS, paper 8.3 committed; verdict KEEP_SPLIT | 11:50:22 | `f682073`, results.md says `Source SHA: b4d34e7` |
| decision-rule script rewritten (per-dataset conditions 1 and 3, fail-closed cond 2, mask-completeness gate) | 11:54:54 | `5709a7f`, after every result had been read |

Findings from the table:

- Both amendments precede every outcome they govern. Amendment 1 (AVONET) precedes any
  AVONET receipt by 48 minutes. Amendment 2 (condition 1 eligibility) precedes any
  structured-arm outcome by 34 minutes.
- Two headers are inexact. Amendment 1 says "before any campaign receipt": the first mask
  receipt was written 7 seconds before the commit (no outcome existed). Amendment 2 says
  "before any structured-arm receipt was read": the FishBase structured mask receipt
  (06:49:15) and its run-log propensity diagnostics were read, and that reading is the
  amendment's own stated trigger. The claim that holds is "before any structured-arm
  outcome existed".
- Amendment 2 did not move the verdict. Condition 1 passes with and without the 5%
  filter (FishBase, all five traits: median far gain 0.0164, minimum Mondrian far
  coverage 0.9235; eligible three: 0.0164 and 0.9630; PanTHERIA unchanged at 0.0110 and
  0.9217). Nothing in `results.md`, NEWS or the paper says so.
- The rule script was changed after the data were read (`5709a7f`). The change aligns
  the code with the pre-registration's "on every dataset" wording (the M2 audit found the
  earlier version pooled conditions 1 and 3 across datasets), adds a fail-closed default
  for condition 2, and adds a mask-completeness gate. I ran both versions: verdict
  KEEP_SPLIT under both; condition 1 and 3 pass per dataset and pooled; condition 2 fails
  either way. This is benign, but it is a post-result edit of the analysis script and it
  is disclosed only inside `review-m2-rule-audit.md`. `results.md` still cites the
  pre-edit source SHA.
- The condition-2 statistic was never written into the pre-registration; it was fixed in
  code at 05:59, before any receipt, and never changed. It is conservative in two ways
  not stated anywhere: it treats the paired Mondrian and split coverages as independent
  binomials (the same cells are scored twice, so the true paired SE is smaller), and it
  omits the calibration-order-statistic term the pre-registration's Uncertainty section
  puts in the MCSE. Both lower the power to demonstrate non-inferiority, so both push
  toward KEEP_SPLIT. This should be stated, not fixed now.
- FishBase is pre-registered as "descriptive only" (Uncertainty section) yet enters the
  condition-2 Holm family as a decision dataset. The verdict does not depend on it:
  AVONET alone fails (raw one-sided p 0.74; observed pooled near gap -0.0234, beyond the
  margin). FishBase's failure is a non-demonstration, not an exceedance: its observed
  pooled near gap is -0.0161, inside the 2 pp margin, with raw p 0.22 and Holm p 0.43
  from one mask.

## BLOCKING

B1. `NEWS.md:25-27` claims "a verbose fit names the trait and its stratum sizes when
Mondrian falls back", and `00-preregistration.md:62` promised "the fallback message names
the realised stratum sizes". The code does not do this in the case that matters.
`R/fit_helpers.R:952-957` resets `n_near` and `n_far` to `NA_integer_` whenever
`fallback` is TRUE, and only then (`:960-965`) formats the message, so a trait that falls
back because a stratum has fewer than 19 residuals prints
`n_near=NA, n_far=NA (floor 19)`. The realised sizes are printed only when they were
never computed. `tests/testthat/test-mondrian-conformal.R:170-172` asserts the NA, so
the test enforces the gap. Two acceptable fixes: (a) keep the realised `n_near`/`n_far`
in the message (and preferably in `mondrian_info`) while still resetting the scores and
threshold, and change the test to expect the realised sizes; or (b) rewrite NEWS:25-27
to "a verbose fit names the trait when Mondrian falls back; the realised stratum sizes
are recorded only for traits that activated" and add a one-line note under
`00-preregistration.md:62` saying the promised message was not delivered in this arc.
Either resolves the public claim; (a) is the smaller honest fix.

## REQUIRED

R1. `NEWS.md:22-23`: "`"mondrian"` remains available and opt-in; it is the better choice
when missing species sit in poorly sampled clades." The pre-registered rule decided the
default and failed; no criterion for "better" was pre-registered, and six of nineteen
near-stratum rows drop below 0.95 under Mondrian (AVONET Mass 0.947, FishBase Troph
0.924 both methods, Vulnerability 0.948, PanTHERIA gestation_d 0.948 random, body_mass_g
0.939 and gestation_d 0.949 structured). Replace with:
"`"mondrian"` remains available and opt-in. In the far stratum it raised coverage in
all 18 trait-by-arm rows (pooled +0.02 to +0.035) at the cost of wider intervals there,
and it lowered near-stratum coverage by 1 to 2 points, in six trait rows to below 0.95.
Users whose missing cells sit in undersampled clades may prefer that trade; the study
measured it on masked observed cells, not on real missing cells."

R2. `useful/paper_section_draft.md:364-366`: "We report Mondrian as the recommended
option when missing species are concentrated in poorly sampled clades." Replace with:
"Mondrian remains an opt-in. The far-stratum gains under the structured mask (condition 1
passed on both databases that ran it) are the evidence for choosing it when missing
species are concentrated in poorly sampled clades; the price is wider far intervals and a
1 to 2 point drop in near-stratum coverage, which fell below 0.95 for six of nineteen
trait rows."

R3. `useful/paper_section_draft.md:362-364` and `NEWS.md:19-22`: "That condition was not
met for AVONET and FishBase, where removing split's over-coverage lowered near coverage
by more than that margin" (paper) and "which the removal of split's near over-coverage
did not meet for AVONET and FishBase" (NEWS). False for FishBase: its pooled near gap is
-0.016, inside the margin; it failed because non-inferiority could not be demonstrated
from one mask (raw p 0.22, Holm-adjusted 0.43). Paper replacement: "That condition was not
met. For AVONET the pooled near-stratum drop (2.3 points) exceeded the margin. For
FishBase the drop (1.6 points) was inside the margin but, with a single mask, could not
be shown non-inferior (one-sided p = 0.22, Holm-adjusted 0.43). PanTHERIA passed (drop
1.2 points, adjusted p = 0.025). The AVONET result alone decides the verdict, so the
default remains split." NEWS replacement: "which failed for AVONET (pooled near coverage
2.3 points below split's) and could not be demonstrated for FishBase from its single
mask (1.6 points below, p = 0.22); PanTHERIA passed. The default remains `"split"`."

R4. `useful/paper_section_draft.md:357-361`: "The split quantile undercovers in the far
stratum and overcovers in the near stratum in every database, including under the random
mask. Mondrian moves both towards nominal: ... near-stratum intervals narrow by about a
fifth while their coverage stays at or above 0.95." These hold only for the eight pooled
rows of Table S-UQ2. Per trait: split far coverage is 0.959 (PanTHERIA gestation_d,
structured) and 0.954 (FishBase Troph), so it does not undercover there; Mondrian moves
those two rows and PanTHERIA litter_size (random, 0.949 to 0.975) away from nominal;
Troph near is 0.924 under both methods; near coverage is below 0.95 in six trait rows
(list in R1); PanTHERIA head_body_length_mm near widened (ratio 1.12). Insert "pooled
over traits" after "in every database" and after "stays at or above 0.95", and add one
sentence: "Per trait the picture is noisier: six of nineteen near rows fall below 0.95
under Mondrian, and in three far rows split was already at or above nominal."

R5. `useful/paper_section_draft.md:277-278`: "the coverage guarantee holds per stratum
rather than only on average." An implied guarantee. Replace with: "the marginal guarantee
then holds within each stratum, provided the cells to be imputed are exchangeable with
the validation cells of their stratum; the structured mask approximates that condition
for real missing cells but cannot establish it (Section 8.3)." This matches
`00-preregistration.md:64-68`.

R6. `useful/paper_section_draft.md:341` ("Mondrian activated for every continuous trait
in all three databases") is true (fallback FALSE on all 38 rows; receipts confirm) but
silent on scope. The masks also covered ordinal traits (PanTHERIA diet_breadth,
habitat_breadth; AVONET Migration), Mondrian activated on them, and
`01_run_masked_confirmation.R:153` scores only `is.numeric()` truth, so they carry no
coverage evidence. Add after line 341: "Count and ordinal traits also received Mondrian
intervals but were not scored; the evidence covers the continuous and count traits
listed in Table S-UQ2 only." Make the same scope note in
`docs/dev-log/after-task/2026-09-23-mondrian-realdata.md:112` ("Covers:
single-observation continuous-family traits" should read "continuous and one count
trait (litter_size); ordinal traits were masked but not scored").

R7. `00-preregistration.md:88-89` says FishBase BodyShapeI "are reported, but in the
MCAR-like condition 2 and 3 tables only." BodyShapeI is a factor, receives no conformal
interval, and appears in no table. Append a dated note under Amendment 2 (do not edit the
amendment body): "Note (2026-09-23, after results): BodyShapeI is categorical and has no
conformal interval, so it appears in no table; the sentence above should have named
Length and Vulnerability only."

R8. Amendment headers. `00-preregistration.md:70`: change "before any campaign receipt"
to "before any campaign outcome; the first mask receipt was written at 06:01:28 and this
amendment committed at 06:01:35, and it concerns AVONET, whose first receipt is 06:49".
`00-preregistration.md:81`: change "before any structured-arm receipt was read" to
"before any structured-arm outcome existed (earliest 07:46); the trigger was the FishBase
structured mask receipt's propensity summary, read at about 06:49". Add to
`results.md` (Method notes) the robustness sentence from the timeline section above:
condition 1 passes with and without Amendment 2, with the numbers.

R9. Post-result rule-script change. Add to `results.md` Method notes and to the
after-task section 3a: "The decision-rule script was rewritten at `5709a7f` (11:54,
after all results were read) to evaluate conditions 1 and 3 per dataset as the
pre-registration says, to fail closed on condition 2 with no evidence, and to gate on
mask completeness. Verdict under the earlier pooled version and the current version:
KEEP_SPLIT in both; conditions 1 and 3 pass under both readings and condition 2 fails
under both." Then regenerate `results.md` at HEAD so `Source SHA` is not `b4d34e7`, or
add "rule script: `5709a7f`" beside it.

R10. Condition-2 statistic. Add to `results.md` Method notes: "The one-sided
non-inferiority test (08_apply_decision_rule.R, `one_sided_binom_z`, fixed at `7af133a`
before any receipt) uses an unpaired pooled-binomial SE on n_test-weighted coverage per
dataset and omits the calibration term of the pre-registered MCSE; both make it
conservative for demonstrating non-inferiority. FishBase, pre-registered as descriptive,
is nevertheless in the Holm family; the verdict is unchanged with it removed, because
AVONET fails alone."

R11. Table S-UQ2 (`paper_section_draft.md:343-355`): FishBase rows rest on one mask and
the pre-registration calls FishBase descriptive. Add "(one mask; descriptive)" to the
two FishBase rows or to the caption, and change the caption's "median width ratio" to
"median over traits of the per-trait width ratio" (coverage is test-cell weighted; the
width column is not).

R12. MI memo provenance. `useful/mondrian-mi-se-justification.md:198-212` reports
numbers the committed diagnostic does not produce: the x-only oracle residual SD 0.325,
the joint conditional residual SD 0.239, `predict_method = "exact"` 0.289, the
independent-noise-around-oracle slope 0.617, and "10 draws" for the oracle, whereas
`13_mi_gls_attenuation_diag.R:58` draws five. No output log of the script is committed.
Either extend the script to print each quoted quantity and commit its log, or mark those
five figures "(from an unscripted follow-up on the same tree; not reproducible from the
committed script)". A re-run of the committed script is recorded at the end of this
file.

R13. `NEWS.md:39-41`: "The draw spread is calibrated on point-prediction error, which is
larger than the proper conditional spread, and the default prediction route does not use
the other observed traits." Stated as fact; the evidence is one tree, n = 400, MCAR, 150
epochs. Replace with: "A one-tree diagnostic points to two causes: the draw spread is
calibrated on point-prediction error, which was 1.3 times the proper conditional spread
there, and the default prediction route does not use the other observed traits." The
memo's Scope paragraph (`:229-233`) already says this correctly; NEWS should match it.

R14. After-task report is stale against HEAD. `2026-09-23-mondrian-realdata.md:4-5`,
`:31`, `:49`, `:57`, `:64`, `:100` still say DRAFT or PENDING for items that now exist
(results table with FishBase, decision, NEWS entry, M2 re-derivation, review panel), and
`:27-28` ("each committed before the data it governs were read") needs the R8 wording.
Fill the PENDING sections from HEAD before the PR.

## SUGGESTIONS

S1. `RUNLOG.md` ends with FishBase split "now running"; add a MEASURED line for the
split finish (split.rds mtime 11:45:54) and the elapsed time, and a line for `5709a7f`.

S2. Pooling across masks before taking the median across traits (results.md Method
notes) was not specified in the pre-registration; it was fixed in `12_build_results_doc.R`
at 07:26, before any structured outcome. Say so in one sentence.

S3. `NEWS.md:16` "in every dataset and arm": insert "pooled over traits" so the
sentence and Table S-UQ2 make the same claim.

S4. Consider reporting the paired condition-2 test (difference of correlated indicators
on the same cells) beside the pre-registered one as a sensitivity line, clearly labelled
post hoc; it can only make non-inferiority easier to show, so it costs nothing in
credibility and answers the obvious referee question.

S5. `results.md` Table 1 has 38 rows over 19 trait-by-arm pairs; the paper text should
say 19 pairs and 38 stratum rows when it counts.

## VERDICT

NOT READY. One blocking item (a public NEWS claim and a pre-registration promise that the
code contradicts and a test enforces), fourteen required wording and provenance items,
five suggestions. The KEEP_SPLIT verdict itself is honest and robust: it survives both
versions of the rule script, the removal of Amendment 2, and the removal of FishBase. The
overreach is in the sentences that follow the verdict, which recommend the option the
rule declined to make default without saying what that recommendation rests on.

## Appendix: diagnostic re-run

I re-ran the committed `13_mi_gls_attenuation_diag.R` at HEAD (`OPENBLAS_NUM_THREADS=1`,
local Mac, about 3 minutes; log in my scratchpad, not committed). Memo figure, re-run
figure: imputed-x to y correlation 0.81, 0.812; true 0.80, 0.798; single-imputation GLS
slope 0.65, 0.66; per-draw OLS slopes 0.81 to 0.86, 0.81 to 0.86; OLS on truth 0.85,
0.851; per-draw GLS slopes 0.33 to 0.38, 0.32 to 0.37; oracle GLS slopes 0.62 to 0.67
with mean 0.655 over "10 draws", 0.64 to 0.66 with mean 0.651 over the script's 5 draws;
GLS on complete truth 0.638, 0.638; oracle conditional SD 0.218, 0.218; pigauto draw SD
0.293, 0.307; point residual SD 0.325, 0.331. The memo's corrected oracle numbers
reproduce. The scripted figures move by up to 0.015 between runs (torch is not
bit-reproducible here), so the memo should say "about" or commit a log. The five
figures named in R12 (0.239, 0.289, 0.617, the x-only oracle 0.325, and "10 draws") are
not produced by the committed script; R12 stands.
