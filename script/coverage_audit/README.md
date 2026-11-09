# Coverage evidence for a frozen development candidate

The baseline XML came from main run 42537328573 at source 033649b9. It records 10,845 covered of 12,749 tracked lines. The 1,904-row inventory is an initial routing map; function ranges come from R's retained parse references, while suggested test files are heuristic hints, not coverage evidence. Classification stays provisional until reviewed.

The original XML omits GraphTransformerBlock and records only ordinary helpers in model_residual_dae.R, leaving both factories' embedded R6 methods untracked. The shipped helper-coverage-torch.R exposes the four actual method closures on both generators to covr's tracing primitives during the full test suite. It preserves existing counters and registers restoration of all patched functions. Normal tests do not register instrumentation. trace_torch.R independently inventories the four methods in the returned coverage object. A complete claim must include all their pinned tracked locations, including unexecuted branches. R/data.R contains six character-literal dataset documentation anchors and no calls or function definitions.

Qualified method names keep Cobertura summaries specific to each source. The verifier reads authoritative class-level line records rather than repeated method summaries. The private covr registration API is guarded and checked with a small positive/negative-control probe: an ordinary function's counters survive registration, both modules acquire counters, an unexecuted method stays at zero, and original method bodies are restored afterwards. Unsupported covr shapes fail instead of skipping instrumentation.

Run full collection on Totoro, against an exported committed tree with a SOURCE_SHA marker and all Suggests installed. Set NOT_CRAN=true explicitly; otherwise supported tests guarded by skip_on_cran() are skipped. Record every environment flag and classify remaining source-layout, cached-data, optional-tool and unavailable-dependency skips. The first pilot omitted NOT_CRAN and is timing/debugging evidence only.

```sh
Rscript script/coverage_audit/collect.R SOURCE_ROOT OUTPUT_DIR FULL_SHA
Rscript script/coverage_audit/check.R SOURCE_ROOT CHECK_DIR FULL_SHA
python3 -B script/coverage_audit/audit.py inventory --baseline docs/dev-log/coverage-100/baseline.xml --out docs/dev-log/coverage-100/inventory.csv
python3 -B script/coverage_audit/audit.py selftest
python3 -B script/coverage_audit/audit.py verify --candidate OUTPUT_DIR/cobertura.xml --manifest OUTPUT_DIR/manifest.json --expected-sha FULL_SHA --git-root . --out OUTPUT_DIR/audit --allow-partial
```

A partial report may exit zero only with explicit --allow-partial; it keeps complete=false and lists missing lines, missing source classes, regressions and source-identity changes. To prove completion, use --require-complete instead. That mode requires Git binding, all executable source classes, unchanged baseline R source identities, matching coverage configuration, a successful collection manifest and every candidate line covered. No source-line mapping for repairs is implemented: changed R files receive no baseline credit, so repair candidates remain separate.

The collection runner retains coverage RDS/XML, line-level CSV, source/test/configuration hashes, installed packages, session records and process output. The package checker disables only the PDF manual, records notes separately, and rejects errors or warnings. These are development checks, not a CRAN release gate or exact-artifact release evidence.

Inventory warnings and skips from the retained installed-package output:

```sh
python3 script/coverage_audit/conditions.py --selftest
python3 script/coverage_audit/conditions.py /path/to/testthat.Rout /path/to/conditions
```

The parser reconciles every individual warning and skip with the final testthat counts. Its categories describe the fixture or prerequisite; they do not establish that a skipped path was exercised. Any unclassified condition remains available for review in the CSV. The audit summary's `uncovered` field is the unique tracked-line deficit; `failures` also counts overlapping provenance and neural-method diagnostics.
