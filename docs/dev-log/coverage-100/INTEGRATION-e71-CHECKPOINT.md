# Integration qualification checkpoint

The campaign remains active. Source `e71e5ddbf37d5bcd1770fc02e7be7ff48290c5fe` combines the Wave 2 tests with the six draft repairs #284, #285, #286, #289, #290 and #291. This is a frozen qualification snapshot on `codex/coverage-100-integration`, not a main-branch result or the final campaign candidate.

The full profile covers **12,679 of 13,003 locations (97.5083%)**, leaving **324 uncovered locations**. All four neural methods, including their 253 tracked locations, are covered. Collection with per-test records used one CPU thread and took 17 minutes 40 seconds, within the 30-minute limit.

Conservative Git-hunk accounting reconciles the denominator exactly:

| Prior catalog disposition | Locations |
|---|---:|
| Unchanged mapped locations covered by named test sets | 12,660 |
| Unchanged mapped locations still uncovered | 324 |
| Changed or removed locations held for separate repair approval | 18 |
| Total prior catalog | 13,002 |
| Added current locations, all covered | 19 |
| Current denominator: 13,002 minus 18 plus 19 | 13,003 |

The original 1,904 misses reconcile separately: 1,573 are covered by named test sets, 324 remain uncovered, and seven changed locations remain unresolved pending repair approval. No original entry disappears. Among unchanged mapped locations, no previously covered baseline line lost coverage. The remaining physical gaps are prediction/baselines 141, inputs/evaluation 78, GNN/fitting 61, MI 40 and presentation four.

Of the 12,679 covered current locations, 12,667 have a single named-test witness. Twelve require collective test sets because different tests exercise overlapping expressions on the same line. Collective sets must not be described as individual test witnesses. The raw counter RDS is retained on Totoro; its SHA-256, the CSV hash and the exact extractor hash are recorded in `named-test-provenance.json`.

The strict verifier remains unchanged and fails closed on the five modified runtime files. The separate reconciliation report maps only lines outside zero-context Git diff hunks, verifies exact text, and leaves changed entries unresolved. It always reports `complete: false`. Independent Sol review reproduced the counts and checked all mapped rows and catalog identities. It also tested the mapper against 200 seeded multi-edit Git-diff fixtures. Whole-file deletion currently fails closed rather than producing a report; no file in this snapshot was deleted.

The exact-source package check passes with zero errors, warnings or notes, including rebuilt vignettes; the PDF manual was disabled. Its installed suite reports 4,886 passes, 192 warnings and 34 skips. The instrumented installed suite reports 4,889 passes, 192 warnings and 33 skips. Both condition inventories retain individual locations and reasons. The package-check skips comprise five missing optional fixtures, one missing HTML Tidy tool, three dependency-absence controls when those dependencies are present, and 25 source-only assertions. All Suggests dependencies were available.

The main source-only gate passes 192 assertions. NEWS adds three independently verified clean assertions from the saved supplementary test results. The supplementary whole-file wrapper reported other test warnings and failed its strict zero-warning condition; its failed log is retained and is not described as a clean whole-file run.

Independent test and tool reviews found no blocker in the inspected coverage increment. The mutation control now first proves that its unmodified baseline passes, then requires the specific gate-aggregation assertion to fail after a mean-to-sum perturbation, with no unrelated failures. Restoration runs through `on.exit`. Sol independently reran that control.

Further real CPU probes found two training issues outside this snapshot. Draft #292 clones the validation-best tensors so later optimizer steps cannot change the checkpoint; its regression failed before the fix and passes 11 assertions afterward. Draft #293 addresses exact-zero gates during iterative prediction and rejects nonfinite active contributions. Its first one-forward proposal was rejected by independent review after an eight-step covariate fixture still returned NaN. The revised head passes fresh bounded review; combined verification remains pending. These findings show why line coverage alone does not establish correctness.

The next candidate must include the reviewed training repairs and any further test increments. G7 remains unmet. Final campaign acceptance is pending. No PR has been merged, no public issue reply has been posted, and no CRAN release action has been taken.
