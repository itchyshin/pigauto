# Frozen evidence for e71e5ddb

These files describe qualification source `e71e5ddbf37d5bcd1770fc02e7be7ff48290c5fe`. The campaign continues beyond this snapshot.

`cobertura.xml` and `manifest.json` retain the coverage counts and source identities. `strict-audit-summary.json` intentionally reports incomplete source identity across the runtime repairs. `reconciliation.json` and `line-accounting.csv` account for every prior location without relaxing that verifier. The compact CSV replaces each full named-test set with its SHA-256; `full-line-accounting.csv.gz` retains those complete strings. `named-test-evidence.csv.gz` retains single and collective witnesses, with hashes in `named-test-provenance.json`.

`new-current-locations.csv` accounts for added current locations. `residuals.csv` lists the exact 324 current misses, with baseline and current positions. The two condition CSVs distinguish package-check and instrumented-suite warnings/skips. Source and NEWS receipts record their narrower assertion scopes. `dependency-versions.csv` records package versions; the full runtime sessions and original library paths remain in the retained server artifacts.

Raw RDS files and output remain under Totoro's `sim-archive/pigauto/coverage100-20261109/integration-e71e5ddb-{coverage,check,source,news-source}/` directories. The raw counter RDS is approximately 148 MiB and is retained there rather than committed to Git. Decompress the CSVs before passing them to the reconciliation tool. Compare uncompressed hashes against the receipts.
