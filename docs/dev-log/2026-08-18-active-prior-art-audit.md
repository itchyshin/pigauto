# Third-party prior-art audit: active acquisition and imputation

## Scope and method

This is a narrow, third-party scan supporting the Stage-A evidence programme.
NotebookLM deep research was asked to find active data-acquisition / active
feature-acquisition work in phylogenetic and mixed-type imputation, excluding
Shinichi Nakagawa's work. Its broad synthesis was treated as a source-discovery
tool only. The claims below were checked against the linked primary conference
or paper source where available.

## Verified adjacent prior art

| Work | Verified contribution | Relationship to Stage A |
|---|---|---|
| Melville et al. (2004), *Active Feature-Value Acquisition for Classifier Induction* | Selects incomplete training instances for paid feature acquisition, comparing acquisition policies against representative/random selection. | Establishes the general idea of feature-value acquisition; it is not phylogenetic trait imputation or a conditional-MVN one-cell variance rule. |
| Saar-Tsechansky, Melville & Provost (2006), *Active Feature-Value Acquisition* | Gives an expected-utility acquisition framework and reports lower acquisition cost at a target predictive accuracy. | Supports using matched random acquisition and realised post-acquisition outcomes, not any pigauto novelty claim. |
| Ma et al. (2019), EDDI | Uses a partial VAE and expected information gain to acquire costly information for target variables. | A mixed-data active-acquisition comparator in principle, but different model, target, and decision utility. |
| Gong et al. (2019), Icebreaker | Uses Bayesian latent-model uncertainty for element-wise training-data acquisition, including imputation and active prediction functions. | Confirms that element-wise acquisition can be evaluated empirically; it is not a phylogenetic comparative trait workflow. |

Primary sources checked: [Melville et al. 2004](https://www.cs.utexas.edu/~ai-lab/pub-view.php?PubID=126910),
[Saar-Tsechansky et al. 2006](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=936445),
[Ma et al. 2019](https://proceedings.mlr.press/v97/ma19c.html), and
[Gong et al. 2019](https://papers.nips.cc/paper_files/paper/2019/hash/c055dcc749c2632fd4dd806301f05ba6-Abstract.html).

## What the scan did not establish

The search did not identify a directly comparable method that selects a
single missing **phylogenetic comparative trait cell** using a conditional-MVN
variance-reduction or phylogenetic label-propagation entropy criterion and then
tests the realised improvement on a protected fixed test set. This is not proof
that no such work exists, so the package must make **no novelty claim** from
this scan. The appropriate public wording, if Stage A succeeds, is only that
pigauto provides and evaluates an active-imputation helper in its stated
single-observation, one-step scope.

## Consequence for Stage A

The adjacent literature supports the locked evaluation discipline: common
candidate set, random comparator, one acquisition, and realised downstream
loss on an untouched test set. It does not alter the DGP, policy definitions,
claim gate, or the exclusion of batch/species/mixed-type claims.
