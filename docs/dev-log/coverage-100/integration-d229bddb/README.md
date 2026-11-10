# Verified CI-portability follow-up

Frozen source: `d229bddb1f280d2537e2753976fbd369af4fff4c`. This collection includes the reviewed explicit-reference-tree fixture fix. Runtime files and coverage configuration match82ce; this is fresh evidence for the updated test suite.

The full raw RDS remains on Totoro at `~/sim-archive/pigauto/coverage100-20261109/integration-d229bddb-coverage/coverage.rds`, SHA-256 `8b2998deb08b3c1862f01c6ade94c9cca7dd9009167df73772139615a609d0e5`. The provenance field `source_manifest_sha256` binds coverage `manifest.json`; `source-manifest.json` separately records the source-only gate. Raw logs and full accounting are compressed without deleting entries; dependency versions omit library paths.

The final CPU report retains97 misses. Same-source real-MPS behavior passes12 assertions separately and is not aggregated into CPU coverage. Historical100 remains unmet and is no longer a delivery requirement. The82ce oracle reviews remain relevant to unchanged runtime code; this bundle adds current-source counter and assembly reviews.
