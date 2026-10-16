## 1. Goal

Use current official BirdTree pages to decide whether pigauto should fetch trees automatically or direct users to obtain and import their own trees, and clarify what the pages establish about attribution and redistribution.

## 2. Implemented

Recorded the official acquisition and citation workflow and a bounded recommendation in the release ledger. No package source, data, or user documentation was changed.

## 3a. Decisions and Rejected Alternatives

Recommend that pigauto point users to BirdTree's official downloads and subset tool, then accept a locally selected tree. Keep automatic retrieval out of the package. The official citation requirements support attribution for research use but do not establish permission to redistribute the bundled tree files. Keep G0 open until redistribution permission is documented or package distribution changes under an explicit implementation decision.

## 4. Files Touched

- `docs/dev-log/cran-0.11-audit/GATES.md`
- `docs/dev-log/cran-0.11-audit/after-task-2026-10-08-birdtree-recommendation.md`

## 5. Checks Run

- Opened the official [Downloads page](https://birdtree.org/downloads/). It offers complete tree archives and user-selected subsets, requires the Jetz et al. (2012) citation for use of full or partial data, and asks users to cite BirdTree.org for web-tool use.
- Opened the official [subset tool](https://birdtree.org/subsets/). It documents user-selected species and tree distributions, a 2,500-species cap, and a downloaded sample with metadata and citations.
- Checked these pages for redistribution wording. They do not state a data redistribution licence or permission to package the tree bytes.
- Ran `git diff --check`, the after-task structure validator, and the no-slop check before staging.
- Did not run a package test or build because no package code or data changed.

## 6. Tests of the Tests

The sources were opened directly from the official BirdTree website, not inferred from a secondary description. The conclusion is limited to those pages. Absence of a redistribution grant on the checked pages is recorded as an evidence gap, not as a claim that redistribution is prohibited.

## 7a. Issue Ledger

- Confirmed: BirdTree supports direct full-tree downloads and user-selected subsets.
- Confirmed: Jetz et al. (2012) and BirdTree.org citations are specified by BirdTree.
- Open: the checked public pages do not establish redistribution rights for pigauto's bundled tree bytes.
- Recommended: direct users to official acquisition and local import; do not build automatic network retrieval into pigauto.
- Open: G0 rights evidence, final source/docs reconciliation, deployed-site review, and exact-artifact checks.

## 8. Consistency Audit

The recommendation matches the source/docs PR's stated user-directed acquisition workflow. It keeps user selection separate from the rights of pigauto to redistribute prebundled data. It does not treat the MIT licence for the `megatrees` software as a licence for the underlying BirdTree data.

## 9. What Did Not Go Smoothly

The public site provides usable download and citation instructions but no explicit redistribution term on the pages checked, so that documentation alone cannot close the CRAN rights question.

## 10. Known Residuals

No legal determination was made. The package's bundle remains in place, the maintainer's redistribution basis is not documented in the sources checked, and no source or data change was authorized in this slice.

## 11. Team Learning

User permission to use and cite a dataset does not, by itself, document permission for a package maintainer to redistribute its bytes. Keep the user acquisition path simple while tracking those rights separately.

Memory receipt: used the pigauto CRAN audit memory entry to preserve the separate rights, source, and publication gates.

Golden Set: not applicable because this slice reviewed official data-access wording and release policy.

## 12. Cross-Product Coverage

This work does NOT cover package source, bundled data bytes, legal advice, installed behavior, or any live pigauto site route. It covers only the official BirdTree acquisition and citation pages and the resulting recommendation.
