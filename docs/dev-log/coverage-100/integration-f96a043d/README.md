# Frozen f96a coverage evidence

Source: `f96a043db32ebf795b34eee348aa997109d92a9c`. These are combined-candidate results, separate from main and the coverage-only PR head. The full raw RDS remains on Totoro at `~/sim-archive/pigauto/coverage100-20261109/integration-f96a043d-coverage/coverage.rds`, SHA-256 `20d098bd52c54b67f478b620c0e2a3b4b5eb3fb418465e97230263967b89f858`.

The provenance field `source_manifest_sha256` binds this directory's coverage `manifest.json`. The separately named `source-manifest.json` records the source-only test gate. Compressed files retain the complete line ledger, named evidence and raw test output. Dependency versions omit machine library paths.

The strict verifier remains partial. All 98 uncovered current locations are retained: 96 catalog locations and two newly added defensive-guard locations. Changed catalog locations remain pending separate repair approval. Post-freeze MPS tests, covariate-name tests and repair #306 receive no credit in these results.
