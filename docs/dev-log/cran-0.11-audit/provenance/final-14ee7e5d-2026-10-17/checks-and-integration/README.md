# Final candidate evidence bundles

These four evidence archives preserve original proof bytes for final release archive SHA-256
14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb.
They are evidence bundles, not the R release tarball. No release READY claim is made.

`manifests/` records every original member hash and size, and each bundle hash.
`views-manifest.json` binds selected readable views to original members. Views apply
CRLF-to-LF only; changed views also have an original-byte `.raw.gz`. All original
bytes, original manifests, failed attempts and NOTASSESSED receipts remain in `raw/`.
No dependency library, symlink directory, exact-source clone or installed check tree is included.

`publication-result.json` records final archive checks before/after, deterministic
round-trip verification, negative controls, exclusions and limitations. Service proof
has 31 declared files plus its manifest; the earlier expected count of 33 was inaccurate.
The DOI access denial remains explicit. Pending native/direct/platform warning proofs
are outside this bundle. MI provenance references checker/backend sources; those
external source trees are not included in the requested manifest-and-logs bundle.

To rerun evidence-bundle packaging, copy the retained `packaging-helper.py` into
an otherwise fresh directory and run `python3 /absolute/path/to/fresh-directory/packaging-helper.py`.
The helper requires the original absolute scratch inputs under
`/private/tmp/pigauto-final-closure-20261017`, including `FINAL.json` and its named
proof directories, plus the final release archive and MI matrix at the absolute
paths recorded there. This published subtree alone is not a self-contained input set.
The helper creates `publication/` beside its copied script and refuses to overwrite
that output directory if it already exists.
Writing-style assessment: NOTASSESSED. Packaging is not a science or release review.
