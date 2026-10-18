# Final-archive warning audit

This scratch audit reads monitor outputs without modifying them. It is pinned to archive `14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb`, 5,137,648 bytes, and direct run `39082452558`. It launches no fit or package test suite. R is used only to deserialize the actual warning objects; all process thread settings are one.

Preparation froze a comparison contract from the previously adjudicated c2 CSV and the byte-identical preserved raw RDS. That contract is not final-archive evidence. It permits the known message categories and provenance only. A changed message, test, file, source line, caller, category set or skip identity fails closed for fresh adjudication; it is not silently added to the contract.

The audit checks actual native versus replay totals, without requiring the old 3,195 passes. It compares CSV rows to actual RDS messages/provenance, checks complete IDs, retains full actual calls, and reports the consequences of each warning category. Matching counts alone is insufficient. Remote download/pre/post identities and replay installation paths must bind to the native checked archive/install. The local native command must name a rehashed final archive copy or the canonical frozen path.

Inputs are expected at:

- Local: the final archive parent's `local-macos/platform-evidence/` and `local-macos/frozen-check/pigauto.Rcheck/`.
- Remote: `/private/tmp/pigauto-final-closure-20261017/final-direct-39082452558/frozen-14ee7e5d-{windows-latest-release,windows-latest-devel,ubuntu-latest-release}/`.

Run once all four cells have finished and been downloaded:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 /private/tmp/pigauto-final-closure-20261017/final-warning-audit/audit.py selftest
PYTHONDONTWRITEBYTECODE=1 python3 /private/tmp/pigauto-final-closure-20261017/final-warning-audit/audit.py audit --output final-results
```

Every invocation requires a new output-directory name; existing results are never overwritten. `--local-root` and `--direct-root` can specify actual monitor directories if they differ. All output remains beneath this audit directory. Missing or incomplete inputs produce NOT_DONE with NOT_AVAILABLE cells, not PASS. The initial availability check found no complete cells.

Successful output includes `audit.json`, `report.md`, classified CSVs and actual object/caller JSON for each cell. Input hashes are recorded and rechecked to detect concurrent changes. Unknown-warning, changed-caller/source, missing-summary, count-mismatch and unexpected-skip controls are included in `selftest`.

The final consequence report is a warning-audit result only. It is not a fresh Grace/Rose/Pat release-panel READY vote, broad calibration/coverage validation, timing acceptance, merge permission or submission authorization. Native check timing remains separate from replay diagnostic time.

## Exact local timing NOTE diagnostic option

The opt-in `--allow-local-compare-methods-timing-note` permits only the observed local examples section `[13s/63s] NOTE`, with `compare_methods` user 2.532, system 0.992, elapsed 18.406 seconds, final `Status: 1 NOTE`, native exit exactly zero and no other NOTE/error/warning finding. It does not admit a remote NOTE. Every archive identity, environment, installation, raw CSV/RDS, warning provenance, skip and native/replay count check remains required. There is no old passing-count pin.

When admitted, the native result stays `NOTE`, adjudication `UNADJUDICATED`, with the full raw NOTE section and exact typed timings. A successful warning review has cell status `PASS_WARNING_AUDIT_NATIVE_NOTE_UNADJUDICATED`; the overall warning-only verdict explicitly retains the unadjudicated native NOTE. No release, timing or submission approval follows. Exit zero from this audit denotes warning-audit success only.

Prior script preserved verbatim as `audit-before-local-timing-note.py` with its SHA-256 file. The extension receipt binds both script versions. Selftests accept the exact local case and reject 14 changes including absent opt-in, remote NOTE, changed method, timing/body/header, additional NOTE/warning/error, nonzero or conflicting exit, incomplete summary and contradictory Status OK. Actual final input audit has not been run; wait for completed-input confirmation.

After that confirmation, use a fresh output directory:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 /private/tmp/pigauto-final-closure-20261017/final-warning-audit/audit.py audit --allow-local-compare-methods-timing-note --output actual-results-with-native-note
```
