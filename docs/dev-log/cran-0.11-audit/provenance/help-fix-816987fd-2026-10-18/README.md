# Source help correction awaiting merge approval

PR #238 at `816987fd160c1d2c5c58b020b07357ff17c077e5` corrects two installed-help descriptions. `impute()` now states that default `gnn = FALSE` uses the phylogenetic baseline; GNN training/calibration requires `gnn = TRUE`. The report help likewise describes baseline-only reports and conditional GNN comparison/history. The report title now says imputation report.

Only `R/impute.R`, `R/report.R` and their generated Rd files change. All 49 executable R syntax trees match merged source `8e275ca`. The provisional installed help and novice workflow pass, including preservation of observed values/classes. Pat's bounded source-fix review is under `pat-review/`. Source CI `39087559368` passed Linux release/devel and macOS release on this exact head. PR #238 is OPEN, ready and CLEAN. The fresh local build/crawl and Rose's Chrome review passed; see `source-ci-final.json`, `pr238-final.json` and `site-review/`.

The predecessor archive `14ee7e5d` remains unchanged and is withdrawn as a submission candidate. Its original three panel verdicts and exact checks are retained in the neighboring predecessor directory. `withdrawal-and-successor.json` records the distinction. The generic validator accepts the reviewed proposed-source rung only; `release-checkpoint-refusal.log` records the strict final verifier refusing the withdrawn archive.

Approval is still needed for #238 merge and its automatic website deployment. Then verify live pages, freeze a successor archive, repeat affected artifact checks and obtain final independent verdicts. Keep #236 unmerged. CRAN upload and maintainer confirmation require separate approval after a new exact archive and submission text are ready.

`files.sha256` binds this evidence packet. Verification scripts retain their exact original scratch/worktree inputs; copying this packet does not reproduce installed libraries or native checks. This slice does NOT cover a corrected frozen artifact, deployment, submission, acceptance or new scientific evidence. Separate assessment records bind the final saved reports.
