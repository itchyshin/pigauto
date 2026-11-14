# Coverage qualification checkpoint

Frozen candidate `a89152a36e72fa1877d4c9587364515121c2434d` covers **12,892 of 13,028 executable locations (98.956095%)**, leaving **136 uncovered**. These are integrated candidate results, not results on main or the coverage-only PR head.

The original 1,904 misses reconcile as 1,756 covered by named test sets, 136 uncovered and 12 changed locations pending separate repair approval. The prior 13,002-location catalog reconciles as 12,830 named-covered, 136 uncovered and 36 changed locations. All 62 new current locations are covered. Thus 13,002 minus 36 plus 62 equals the current denominator of 13,028. No entry disappears and no unchanged previously covered location regresses. All 253 tracked neural-method locations are covered. Current covered locations have 12,880 single-test witnesses and 12 collective witnesses.

The instrumented suite passes 5,352 assertions with 192 inventoried warnings and 33 classified skips. The installed suite passes 5,349 assertions with 192 warnings and 34 skips. The extra installed NEWS skip accounts for three source-only assertions executed by the instrumented run. Package check reports zero errors, warnings and notes, with rebuilt vignettes and no PDF manual. The separate source gate passes 192 assertions without conditions. Collection took 17 minutes 35 seconds on one CPU thread; the runtime receipt records torch 0.17.0, CPU tensor execution, unavailable CUDA/MPS, and available zip, unzip and Pandoc.

Independent Sol review checked all catalog keys, mapped source text and hunks, coverage counters, witnesses, new-location complement and raw condition logs. The raw RDS hash is recorded remotely rather than independently downloaded. Strict verification continues to reject modified source identity; conservative reconciliation holds changed locations unresolved. G7 remains unmet with 148 original entries unresolved, and final G8 awaits a final candidate.

This source includes tests through 4dcd0f40 and benchmark diagnostic repair #297. The later ambient-seed test, redundant-guard cleanup #298 and private conditional-BM guard repair #300 are outside this measured source. The later candidate 3898a77b1a953341a414b377c16a3e0e4e95ba5c passes 1,634 focused assertions across 22 changed files, with full verification running separately.

No coverage exclusion was added. Line coverage does not establish statistical calibration, GPU correctness, inferential validity or resolution of the singleton limitation. Draft PRs remain unmerged; no issue update or CRAN release action was taken.
