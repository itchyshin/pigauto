# Independent pre-push checker review

**PASS for checker retargeting only.** This does not establish package, platform or release readiness.

The checker is bound to successor source `b2e6054685e6bdedb073ba46e058af2c0a619a18` through the retained freeze receipt: archive SHA-256 `a93fece7e272565a97518021ed12c809e04f41462d2e163f8b806f5060a29feb`, 5,137,837 bytes, 250 entries, inventory SHA-256 `d7b9b0e59ea051b3a22d7f347d8d925e6265b076b238ac893410e79d5b61329d`. The new public-object readback receipt records the same hash and size. No new network request was made in this review.

Each backup is byte-identical to the corresponding file at the starting checker HEAD. Independent Python syntax-tree comparison, masking only URL/SHA/SIZE/INVENTORY assignments, confirms that `artifact.py` has no other executable change. Its new constants agree with successor FINAL.json and the object-readback receipt. INVENTORY is unchanged because the package paths are unchanged; the different installed bytes require the new gzip hash and size.

`verify_results.py` changes exactly one artifact-folder label. The workflow changes exactly its concurrency-group label and retained-artifact name. Matrix cells, dependency/manual setup, direct original-gzip check, all-Suggests/CRAN settings, numerical thread caps, warning replay, pre/post identity checks, failure retention and 90-minute limit are unchanged. The older dependency-cache version label is inherited cache reuse; it does not select the checked gzip and is not a substituted artifact identity.

The retained selftest reports one positive and five negative controls. A fresh independent local readback also exercised the retargeted identity function: the real successor passes; wrong SHA, wrong size, wrong inventory, wrong entry count and same-size byte mutation are all refused. The mutation existed only in this scratch directory and was removed. The real archive hash remained unchanged afterward. No source check, native check, fit, upload, repository edit or network mutation was performed.

`readback.json` records the exact before/after file hashes and control results. `verdict.json` is scoped to those hashes. Later changes to the three reviewed files require renewed review. Successful future platform checks must retain their actual run/checker identities separately from archive source b2e6054; R-hub rebuilt-source status remains supplementary evidence.

Assessment: separate exact-draft self-assessment accompanies this report; it does not independently certify the package or release.
