# DOI identity and access: bounded adjudication

The exact final source URL checker returned 403 Forbidden for https://doi.org/10.1126/science.1157704 in vignettes/getting-started.Rmd and the built inst/doc/getting-started.html. Preserve that result; no URL/source edit or URL-clean verdict is made.

A single independent GET of primary Crossref metadata returned HTTP 200 and the registered identifier 10.1126/science.1157704: Shannon J. Hackett as first author; A Phylogenomic Study of Birds Reveals Their Evolutionary History; Science 320(5884), 1763–1768, 27 June 2008. The fields agree with the cited reference. See https://api.crossref.org/works/10.1126/science.1157704 and NLM https://pubmed.ncbi.nlm.nih.gov/18583609/ . Raw response, digest and parsed fields are retained in crossref-raw.json and doi-metadata.json.

Conclusion: bibliographic identifier identity is supported by primary registration metadata. Automated URL accessibility remains a recorded 403; the reason for denial is not established. This is suitable for explicit URL-finding disclosure/adjudication in the final release record or cran-comments; it does not prove CRAN will accept the access result. No recheck was repeated merely to obtain a passing status.
