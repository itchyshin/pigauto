# Maintainer warranty and lane-state refresh: after-task report

## 1. Goal

Record the maintainer's rights confirmation and reconcile the CRAN audit gate ledger with the current source candidate.

## 2. Implemented

Marked G0 met for the candidate source revision, documented the maintainer warranty basis and the limits of the published BirdTree terms, and refreshed the gate tally to 8 met and 3 open. Preflight confirms this is the only active pigauto lane; the 58 worktree entries are checkouts, not separate lanes.

## 3a. Decisions and Rejected Alternatives

Recorded the maintainer's confirmation as the warranty basis for bundled BirdTree-derived example trees with citation and credit. Did not describe it as a separately published BirdTree licence or written rights-holder permission. Exact shipped objects and notices remain subject to post-merge tarball inspection under G8.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/provenance/rights-and-policy.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-09-rights-refresh.md`

## 5. Checks Run

- `bash ~/shinichi-brain/tools/lane_preflight.sh <pigauto-worktree>` reported one active lane (worktree census 58).
- `python3 ~/shinichi-brain/tools/route.py pigauto` loaded the repository's current LOAD-FIRST manifest.
- The candidate-source GitHub Actions matrix at run [#37878996605](https://github.com/itchyshin/pigauto/actions/runs/37878996605) passed Ubuntu release, Ubuntu R-devel, and macOS arm64 jobs at source head `efb007e9592511fcb902112ba99625e9fd56082d`; this is not exact-tarball evidence.
- Chrome rechecked official BirdTree downloads and FAQ pages and the CRAN source-package policy. Those BirdTree pages require attribution but do not state a separate redistribution licence.
- `git diff --check` passed. Absolute-path prose assessment and after-task structure checks remain to be run for this report before commit.

## 6. Tests of the Tests

No test code changed. The candidate-source platform matrix exercised the existing package test suite; no negative control was warranted for this documentation-only update.

## 7a. Issue Ledger

Resolved: G0 is met for the current candidate source on the maintainer-warranty basis described above.

Open: G7 deployed-site verification, G8 exact post-merge tarball and checks, and G9 independent final review. The exact contents and attribution of the bundled tree objects remain to be verified in G8.

## 8. Consistency Audit

Aligned the former open-G0 language with the maintainer's subsequent confirmation. Preserved the distinction between the `megatrees` software licence and BirdTree data rights. The gate ledger records 8 met and 3 open. The prior CI run remains candidate-source evidence only and used `NOT_CRAN=true` with force-Suggests disabled.

## 9. What Did Not Go Smoothly

The first preflight invocation used the Python interpreter on a shell script and failed without changing files. Rerunning it with `bash` produced the intended one-lane census. No other operational problem arose in this slice.

## 10. Known Residuals

No merge, deployment, final tarball, Windows result, exact-artifact check, independent final verdict, or CRAN submission has occurred. A successful candidate-source matrix does not satisfy G8. Public source pages do not provide a separate BirdTree redistribution licence, so the recorded rights basis is the maintainer's warranty confirmation.

## 11. Team Learning

The current lane census is one; saved worktrees are checkouts and do not represent additional active pigauto lanes. Keep this distinction explicit when interpreting preflight output. Memory receipt: project LOAD-FIRST manifest was refreshed with `route.py`.

## 12. Cross-Product Coverage

Covers the rights-policy record and candidate-source CI evidence for this source revision. Does NOT cover the deployed website, exact frozen tarball, force-Suggests checks, Windows, final independent review, or CRAN acceptance.
