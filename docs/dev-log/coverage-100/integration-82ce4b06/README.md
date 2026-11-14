# Final prepared-source evidence

Frozen source: `82ce4b0625a679026a6f9a9e0494aaa43da31ecf`. The full CPU raw RDS remains on Totoro at `~/sim-archive/pigauto/coverage100-20261109/integration-82ce4b06-coverage/coverage.rds`, SHA-256 `dbdbf11a3b5a45cf6c3276b9cf3cc8d4d62bcfe8f065ab18f91a6230126231c0`.

The provenance field `source_manifest_sha256` hashes coverage `manifest.json`; `source-manifest.json` separately records the source-only gate. Compressed files retain complete line accounting, named witnesses and raw logs. Machine library paths are omitted from dependency versions.

There are 97 current misses: 95 catalog locations and two new defensive guard locations. All remain in the inventory. Changed catalog locations require separate repair approval. The historical 100% gate remains unmet; Shinichi approved delivering the prepared fixes and exact remainder instead of continuing percentage-driven expansion. The focused real-MPS test passes separately and is not aggregated into CPU coverage.
