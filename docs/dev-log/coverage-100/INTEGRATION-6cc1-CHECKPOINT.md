# Coverage qualification checkpoint

Frozen source `6cc1ffa73b6e9fffa9fe203bbfe2a0715eb9156a` covers **12,868 of 13,016 executable locations (98.862938%)**, leaving **148 uncovered**. The preceding `0aad929c263419ab238df326951d1f2ab186651a` snapshot covered 12,843 of 13,019 locations, with 176 uncovered. Both evidence bundles are retained beside this report. These are combined integration results; no percentage is attributed to main or the coverage-only PR head.

The current snapshot contains coverage tests through `9164473d` and separate repair PRs #284, #285, #286, #289, #290, #291, #292, #293, #294, #295 and #296. Later tests, the benchmark diagnostic repair #297 and redundant-guard cleanup #298 are outside this measured source.

| Prior catalog disposition | Locations |
|---|---:|
| Unchanged mapped locations covered by named test sets | 12,818 |
| Unchanged mapped locations uncovered | 148 |
| Changed or removed locations pending approval | 36 |
| Prior catalog total | 13,002 |
| New current locations, all covered | 50 |
| Current denominator: 13,002 minus 36 plus 50 | 13,016 |

The original 1,904 misses reconcile as **1,744 covered by named test sets, 148 uncovered and 12 changed locations pending approval**. No original entry disappears. All 253 tracked locations in the four neural methods are covered. Of current covered locations, 12,856 have single named-test witnesses and 12 require collective witnesses. No previously covered unchanged location regressed, across all mapped files.

Independent Sol review recomputed the artifact hashes, checked every catalog identity, mapped Git line and counter value, and verified the new-location complement. It independently parsed the raw logs and inspected the source, package-check and CPU receipts. The raw counter RDS remains remote; its byte hash is recorded in the receipt and was not independently downloaded for review. The strict verifier still rejects full source identity for the eleven modified runtime files. Conservative reconciliation preserves those changed entries as unresolved and always reports incomplete.

The instrumented suite passes **5,301 assertions**, with 192 inventoried warnings and 33 classified skips. The exact-source package check passes with zero errors, warnings or notes, including rebuilt vignettes and excluding the PDF manual. Its installed suite passes **5,298 assertions**, with 192 warnings and 34 skips. The extra installed-package skip is the NEWS source test; its three assertions execute in the instrumented run. The separate source-only gate passes 192 assertions without conditions. Missing optional fixtures and HTML Tidy remain disclosed; all condition rows reconcile without an unclassified entry.

Collection took **17 minutes 38 seconds**, using one CPU thread and torch 0.17.0, within its 30-minute cap. The runtime receipt records actual CPU tensor execution and unavailable CUDA/MPS. Execution of the whole campaign has reached approximately seven hours. Recent profiles reconcile roughly 20 to 30 additional missing locations per cycle, but the remainder increasingly consists of numerical failure paths and contract-dependent guards. Further completion time is uncertain; the next exact-source verification remains bounded to 30 minutes per run.

The current remainder is prediction/baselines 89, inputs/evaluation 15, MI 23 and GNN/fitting 21. Each missing line is retained in the residual CSV. Provisional source classifications are not proofs of unreachability. The original-baseline remainder is 160 locations when the 12 changed entries pending approval are included. G7 remains unmet; final G8 acceptance awaits the final candidate.

A later combined candidate, `a89152a36e72fa1877d4c9587364515121c2434d`, has passed 1,560 focused assertions across its 20 changed test files and its package check. Its instrumented coverage collection is running separately. The benchmark repair preserves a caught pipeline failure and attempted-replicate count when no evaluations survive. Generation errors still propagate directly; partial summaries still require a separate failure ledger for scientific comparisons. The cleanup in #298 has passed 42 public-path preservation assertions before and after removal, but has no measured coverage credit yet.

Targeted-run reporting was corrected: testthat's `nb` includes condition records. The related #295 run records 165 passes and 17 warnings, rather than 182 passing assertions. Broad passing-assertion totals without separate retained pass counts were removed from the descriptions of #292 and #293. Full-suite PASS totals in these profile reports come directly from the final testthat Rout summary and are unaffected.

Coverage does not establish calibration, GPU correctness or singleton resolution. No PR has been merged, no public issue update has been posted and no CRAN release action has been taken.
