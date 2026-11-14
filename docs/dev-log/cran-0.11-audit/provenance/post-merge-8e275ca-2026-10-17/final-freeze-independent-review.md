# Independent final-freeze verification

Verdict: **PASS for artifact identity, source binding, exclusion inventory and unchanged data/NOTICE bytes only.** This is not a final release-panel READY vote or evidence of completed checks on this archive.

Reviewed independently on 2026-10-18T00:10:38.122575+00:00. No producer freeze script was executed. Only this report was written.

## Exact archive

- Path: `/private/tmp/pigauto-release-candidate/sha256-14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb/pigauto_0.11.0.tar.gz`
- SHA-256, independently read twice: `14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb`
- Size: 5,137,648 bytes.
- Inventory: 250 unique paths, 236 regular files, 14 directories; no links, absolute paths or traversal members.
- Sorted inventory SHA-256: `d7b9b0e59ea051b3a22d7f347d8d925e6265b076b238ac893410e79d5b61329d`; independently derived inventory is byte-identical to inventory.txt.
- Regular uncompressed bytes: 7,150,500.
- FINAL.json and freeze-receipt.json agree. Every extracted regular file also matches its archive member.

## Tracked source binding

Source commit: `8e275ca738129c8d608f06ffc407b8141c1a2636`. The recorded checkout was at that commit and clean before and after this verification. A separate git-archive read of the commit found **215 byte-identical tracked regular package files**, only `DESCRIPTION` with expected build normalization, and exactly **20 build-generated files**. DESCRIPTION changes are whitespace normalization plus generated Author, Maintainer, NeedsCompilation and Packaged fields; all original semantic fields agree. Package/version are pigauto 0.11.0, with NeedsCompilation=no.

Generated-only files are build/partial.rdb, build/vignette.rds, and the R/Rmd/HTML outputs in inst/doc for six vignettes. No other untracked package payload was found. NEWS is shipped and source-identical.

## Exclusions and data

All archive payload paths were checked against the recorded source .Rbuildignore, including excluded parent directories. No forbidden match, docs/script/.github/BACE/dev payload or package-owned src file was present.

All six .rda streams and inst/NOTICE were independently hashed and compared byte-for-byte with the independently rehashed c2 predecessor. All seven are unchanged:

| Package path | SHA-256 |
|---|---|
| `data/avonet300.rda` | `86af337b49cbe2a81ff7acb7245b75c9a00cd48bb9131d7f4eb5d2676af131be` |
| `data/avonet_full.rda` | `29b47bccdf17233328321378b77e99e8d4aa515d9df4451e347c709b97568e81` |
| `data/ctmax_sim.rda` | `987a05458b7f3e548f22c52b75d8b7463e29a281a48895303f564effb0b9cecc` |
| `data/tree300.rda` | `6050855d2ae7de42fa0de239134aad8e8975a224610b8dffb48f3d7412bb383f` |
| `data/tree_full.rda` | `5f7b36a4b9c02fd53c98d12fe4871e0f34098e3aa7a8f34f2db38a21a8dc09f5` |
| `data/trees300.rda` | `8293e06503e89e79d9233abd1cbae98fc8d0c18e05ef5a41cca99a071fee9108` |
| `inst/NOTICE` | `9cc7fa419cacea99dbb0b822c59d6fffad514ee072e119a1ce1c7f365883ba02` |

This establishes byte identity, not an independent determination of redistribution rights. The recorded maintainer-warranty basis is unchanged.

## Limits and observation

Local exact checking and installed-MI runs were in progress outside this verification. Their results, direct platform checks, warning identities, deployment parity, final panel votes, timing acceptance and submission readiness are not established here.

The tracked shipped file `tests/testthat/_problems/test-new-features-526.R` is inherited byte-identically from c2; it is not a new archive delta and is not matched by the current exclusion rules. Its presence is recorded for final inventory review, not classified as a new failure.

Assessment: independent engineering verification report; perceived AI-like style 2/10, medium confidence. The exact-archive and limits paragraphs state measured scope without a readiness inference; no repair needed. Science gate: not applicable. Facts gate: pass within the direct checks above. References gate: not applicable; no literature claims. Sources and prior receipts were visible during review.
