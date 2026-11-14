# Capture timing and provenance amendment

The earlier derived verdict combined capture B's supervisor duration with capture A's diagnostic duration. That combination was invalid. Two successful captures are retained:

- Capture A: supervisor 171.199674 s; diagnostic 168.258 s. Supervisor receipt SHA-256 `39c9908db1e4e8db6fee2126eecec108680d20c819fcacdb16a31d2244dc560b`; log SHA-256 `48f9938cc20133be2aeba43377b1cbc8ca51738acdce5d6e1746f335a11e4e73`.
- Capture B: supervisor 167.394193 s; diagnostic 164.807 s. Supervisor receipt SHA-256 `e2ae6e1f39ba49b86fc0fffabab19b97d26348af67d257b53d9fdf21041680af`; log SHA-256 `984e4e509a06af70427a78c044d7b8d3acb2cf1625426db1dd7b707c9003112e`.

Both records show 0 failures, 161 warnings, 86 skips and 3,195 passes. The exact local replay JSON used by the audit has SHA-256 `c1f759a30a97195184e7b52b40363bbd1a5fe67062c3a6fc70f6ffe3fe65713b` and carries capture A's diagnostic time. The local RDS input has SHA-256 `a01ed64252886a8f4f2bb0c5bb817f58b908c0d008f0c70d8f9f229e5520e852`; the surviving capture B RDS has the same bytes. Its unique process attribution cannot be resolved from the hash.

The warning audit remains PASS for the exact inputs recorded by the auditor. No replay or platform check was run for this amendment. Failed setup attempts and prior derived files are preserved. See `amendment.json` for the full input hashes and `verdict.json.before-amendment` for the original derived receipt.
