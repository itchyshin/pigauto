# Win-builder upload acknowledgement: `fe1600d7d562`

Checked 2026-10-17 00:40–00:50 UTC using Chrome to inspect the official form and the installed `devtools` form uploader to send the exact frozen file. No package rebuild occurred.

- Local input: `/private/tmp/pigauto-release-candidate/sha256-fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a/pigauto_0.11.0.tar.gz`
- Local SHA-256 before/after: `fe1600d7d5629212063802a72dd4083bbd14c453ab5786d1d0e871e0ccfeff9a`
- Local size: 5,136,083 bytes
- R-release: HTTP 200. The returned form displayed `pigauto_0.11.0.tar.gz`, 5,136,083 bytes, and `application/octet-stream`.
- R-devel: HTTP 200. The returned form displayed the same filename, byte size, and content type. The form also displayed `Currently unavailable` for R-devel.

These are upload acknowledgements, not check results. Win-builder did not provide a server-side SHA-256; byte size and filename alone cannot prove that the server checked these exact bytes. No current-hash Windows result is yet available. A separate GET after submission returned the form with both fields blank, so it is not retained as a persistent receipt. The immediate POST response bodies were captured in the task output; the acknowledgement details above transcribe those responses. The frozen archive SHA was rechecked and unchanged after both submissions.

Next: retrieve and retain any R-release result email/log, verify its filename and size against this receipt, then re-evaluate platform evidence. R-devel availability remains an external dependency. Do not promote the release rung from `tarball-clean` until valid platform results are independently reviewed and bound as far as the service permits.
