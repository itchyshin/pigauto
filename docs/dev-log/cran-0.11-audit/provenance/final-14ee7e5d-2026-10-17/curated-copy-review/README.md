# Independent curated-copy mechanical review

PASS for current copy mechanics only. No final release-panel or readiness verdict.

Both source publication trees and copied repository subtrees have exactly their indexed files plus the index itself: checks-and-integration uses a `files` wrapper (40 indexed, 41 actual); external-platforms uses a direct map (162 indexed, 163 actual). All 202 indexed hashes/sizes and both index files match byte-for-byte from scratch publication to target. No unindexed files or symlinks were found.

All seven bundles were decoded independently. Their exact member sets, 228 original member hashes/sizes, deterministic regular-file metadata, and retained scratch source bytes match. Gzip decode/reencode/decode preserves every tar byte. All 108 readable views match their declared source; 51 changed views apply only CRLF-to-LF and retain exact raw bytes in .raw.gz. CSV logical records match after that declared newline normalization. The first independent reader compared literal CRLF/LF inside quoted Windows warning messages; its failed script and correction are preserved here. This was a reader correction, not altered evidence.

Exclusion checks found no release tarball, exact-source checkout, dependency library, symlink, installed package tree, credential/config file, duplicate ZIP or Python cache in these bundles. The retained platform environment capture contains the expected check/R/GitHub-path variables. The repository's `^docs$` build exclusion covers the copied target. Read-only git diff of package-bearing paths against source 8e275ca exited zero.

The final archive digest remains 14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb, 5,137,648 bytes. Artifact receipts bind package source 8e275ca738129c8d608f06ffc407b8141c1a2636. Direct identities bind run 39082452558 and executed workflow merge e94e40a2f66a7f27dcc69b515d762d69c6eadf66. Retained primary GitHub commit bytes hash to 2bdcfce8d0f76182781f617c1902fd26a67c71cfe8b36d5621227e14a1146218, with parents 8e275ca and checker head 6cd5f65. Its tree f10989373c12981c0be52c7a33a665b73e0c6fce matches local checker-head git output. Package source, checker head, and executed merge are distinct identities.

Four planted in-memory controls refused changed hashes, changed sizes, traversal and absolute paths. Original source/target/monitor evidence was read only. No R, network, source edits or ref mutations occurred; all review writes are in this scratch directory.

No mechanical publication blocker found. One reader instruction merits correction: checks-and-integration/README.md names `python3 ../bundle.py`, which is absent from the copied tree; its retained helper is `packaging-helper.py` and depends on the original scratch evidence paths. Describe that dependency and helper name accurately rather than implying the published subtree alone regenerates the bundles. Raw bundles and indices remain independently inspectable without those original paths.

This review does not adjudicate local timing NOTE, warning consequences, URL access denial, package size, or submission approval. The preserved NOTASSESSED receipts remain NOTASSESSED.
