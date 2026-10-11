# Win-builder result receipt, 2026-10-09

The frozen artifact is `pigauto_0.11.0.tar.gz`, SHA-256 `f8fee9f631460a1a49fa7ef78b4256e0861ea1929f25cba5093d0ea034f5a9c4` (5,135,324 bytes; source commit `d76804e768bf77f430f43dcf591633f9cd900dab`). The two current Win-builder emails/logs report **Status: 1 ERROR** on R-release and R-devel. Each log records five failed test expectations: three torch-dependent tests report “Lantern is not loaded,” and two mocked accelerator checks expected CUDA/MPS but received `unavailable`. The emails and logs contain no source archive checksum, so these failures are diagnostic evidence and are **not proven to belong to the frozen SHA-256 above**.

- R-release log SHA-256: `575179ab3b680d9ea23e54e8c0550bd49d817dee2e652ab11aeacc5ed03af971`.
- R-devel log SHA-256: `23cc52068814203057efa2862110f9c08eaab62e4fcf4a62985390975849be54`.
- A separate candidate repair exists in another active lane; it has not been merged or validated on this frozen artifact. No Windows pass is established.

The current exact-artifact gate therefore remains open. Re-freeze after any source change and obtain Windows results with upload provenance bound to that exact archive before claiming platform passage.


Raw log and compressed-file hashes:
- `winbuilder-r-release-2026-10-09.log` raw SHA-256 `575179ab3b680d9ea23e54e8c0550bd49d817dee2e652ab11aeacc5ed03af971`; retained as `winbuilder-r-release-2026-10-09.log.gz` with compressed SHA-256 `dc69a974ea7e693a632d063ca6d205b9c4d61b2eb171f6068c4e8aff8aa19fb7`.
- `winbuilder-r-devel-2026-10-09.log` raw SHA-256 `23cc52068814203057efa2862110f9c08eaab62e4fcf4a62985390975849be54`; retained as `winbuilder-r-devel-2026-10-09.log.gz` with compressed SHA-256 `6242c79102333744dab7e14d86237614a51dea51247888cae188eda25a9fb9b2`.
