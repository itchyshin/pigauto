# Bounded independent spelling review

PASS for adjudication coverage of the 404 spelling flags from the final frozen archive, SHA-256 `14ee7e5df99b0474c3b9b970acc7ec2e0a82d22976827dceb85b5d00755c03bb` (5,137,648 bytes). This is not a release-readiness verdict.

All 404 distinct flagged tokens were reviewed in their frozen-source contexts. All 43 referenced source files were byte-compared directly with their archive members. The context CSV retains 1,886 reported locations and the nearest actual token lines; the spelling tool sometimes reports the paragraph start. Boundary matching excludes incidental substrings inside longer words.

| Category | Tokens |
|---|---:|
| Mathematical, code, abbreviation, DOI or reference particle | 144 |
| Proper, taxon, package, named-method or exact API token | 110 |
| Valid scientific/statistical term | 65 |
| Accepted British spelling | 52 |
| Valid computing term or productive compound | 33 |

No actual ordinary-word typo was identified among these 404 flags; none remains uncertain as a spelling token. This does not establish that unflagged prose, grammar, factual claims or references are clean. Proper-token classification is not independent bibliographic verification.

Scope is explicit per token: 154 occur only in historical NEWS, 137 only in current reader documentation, and 113 in both. Historical content has not been promoted into current evidence.

The requested odd fragments resolve as follows: `cffa` is part of commit `dc8cffa`; `decd` is part of `573decd`; `bigl`/`bigr` are LaTeX sizing commands; `mis` comes from `mis-extrapolates`/`mis-pool`; `un` comes from `un-log`/`un-restored`; `arg` is historical shorthand for argument; `passthrough` is a computing compound in historical NEWS. `writeup` is informal historical wording, optionally styled `write-up`, not an ordinary-word typo. `mondrian` is the exact quoted API value.

`token-adjudication.csv` contains one category/adjudication per input token and all source references. `source-contexts.csv` preserves exact contexts and both reported/actual lines. `receipt.json` binds the input CSV and all examined source files. `review.py` is the reproducible coverage/archive-binding checker; the category sets and contextual notes encode this review's human-style adjudication.

No package source, archive, dictionary, workflow, warning audit or monitor output was edited. No package fits or test suites were run. The highest proven statement is that this flagged-token list has no identified ordinary-word spelling error; it is not a claim that all package spelling is clean.
