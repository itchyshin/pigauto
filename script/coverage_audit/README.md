# Coverage evidence for a frozen development candidate

The baseline XML came from main run 42537328573 at source 033649b9. It records 10,845 covered of 12,749 tracked lines. The 1,904-row inventory is an initial routing map; function ranges come from R's retained parse references, while suggested test files are heuristic hints, not coverage evidence. Classification stays provisional until reviewed.

The XML omits the graph-transformer factory's embedded R6 methods. `trace_torch.R` exposes the actual generator to covr's R6 traversal and runs its existing tests. It does not replace mathematical code or alter the installed package. `collect.R` retains the primary and supplemental coverage objects and combines their disjoint source files. A 100% claim must include this added denominator. The only non-executable source record is data.R: six character-literal dataset documentation anchors, with no calls or function definitions.

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
