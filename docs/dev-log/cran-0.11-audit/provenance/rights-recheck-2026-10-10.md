# Rights-source recheck, 2026-10-10

This is a dated supplement to `rights-and-policy.md`; it preserves that report's earlier findings and does not independently determine legal rights.

## Sources checked

- The official [BirdTree home page](https://birdtree.org/) and [downloads page](https://birdtree.org/downloads/) describe access to full or partial phylogenetic tree data and require citation of Jetz et al. (2012); use of BirdTree's web tool also requires citing BirdTree. The pages checked do not publish a data redistribution licence.
- The official [CRAN Repository Policy](https://cran.r-project.org/web/packages/policies.html) requires clear ownership and intellectual-property rights for package components and states the maintainer's warranty for included third-party material.
- The upstream [`daijiang/megatrees` repository](https://github.com/daijiang/megatrees) currently exposes `get_tree_bird_n100()` and describes BirdTree-derived large tree distributions as release assets. Its package metadata declares MIT + file LICENSE. That software licence is provenance about the `megatrees` package; by itself it does not establish the licence or redistribution scope of the upstream BirdTree data.

## Disposition

Shinichi's direct maintainer statement, “we should cite and give credits but we can use the trees freely yes,” remains the recorded maintainer-warranty basis for including the bundled tree objects with attribution. The pages above provide citation and acquisition instructions but do not independently document a BirdTree data redistribution grant. Accordingly, this recheck adds source context and preserves the distinction between the maintainer's warranty and independent rights evidence. It does not change the existing disposition, clear a separate legal question, or close the exact-shipped-file check: the frozen tarball's tree objects and notices still require inspection under G8.

No external contact was made. The required citations and attribution remain in the package materials; this note does not authorize changing those materials or the bundled data.
