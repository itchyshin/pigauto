# GitHub R-devel source-CI receipt — 2026-10-16

- Repository: `itchyshin/pigauto`
- Workflow run: [38843648553](https://github.com/itchyshin/pigauto/actions/runs/38843648553)
- GitHub run API: `https://api.github.com/repos/itchyshin/pigauto/actions/runs/38843648553`
- Source commit: `8082ef8273177ff664d49ac4bd643f7def671f20` (`main`, push after PR #234)
- Job: `ubuntu-latest (R devel)`, job ID `116789146589`
- Runner: Ubuntu 24.04.5 LTS
- R: `R Under development (unstable) (2026-10-14 r90673)`, R 4.7.0
- Job conclusion: success; package check step conclusion: success
- Check summary: `Status: OK`; testthat `[ FAIL 0 | WARN 175 | SKIP 83 | PASS 3196 ]`
- Workflow environment: `_R_CHECK_FORCE_SUGGESTS_=false`, `NOT_CRAN=true`; check action used its workflow configuration rather than the frozen local `--as-cran` command.

This is source-CI evidence on the exact frozen source commit, not a check of the frozen tarball. It does not substitute for Windows R-devel, nor establish the exact-artifact platform gate. The full job log was fetched read-only from GitHub job ID `116789146589` on 2026-10-16; this receipt preserves the relevant log facts rather than the 279-KB raw job transcript.

## Bounded conclusion

GitHub R-devel source CI passed on the source commit used to build the frozen archive. The run's `NOT_CRAN` and force-Suggests settings differ from the exact local CRAN check. Win-builder R-devel is still pending, so platform-clean remains unproven.
