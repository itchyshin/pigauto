# Independent review receipt

Source: implementation commit 7e542057, including the final report repair. Base 51e9a372. Reviewed on 7 November 2026.

The implementation producers used Terra medium in fresh native sessions.
A Sol high reviewer independently inspected the combined changes. A Luna
medium reviewer checked the export, manual, inventory, navigation and issue
response coverage. These are dispatch identities, not measured cost savings.

Findings corrected before publication:

1. A factor level present only among missing outcomes has an unidentified
   covariate effect. Missing-row aliases are now checked and warned about;
   equivalence is asserted only for estimable designs.
2. An unavailable lambda estimate could break prediction replay. The stored
   numerical fallback remains finite, with unavailable estimates marked in
   metadata. A public sparse fit/replay regression exercises this path.
3. The visualisation example passed preprocessed data to a plot expecting
   original-scale values. The revised guide uses held-out evaluation instead.
4. Plot labels were clipped. Adaptive margins, restored graphics state,
   larger example figures and co-missingness fractions improve readability.
5. Bare-fit reports omitted held-out metrics despite the README promise.
   They now recover stored observations and splits; a regression compares
   their metric table with the full-result report.

The Sol reviewer found no remaining P1 or P2 issue in the reviewed source.
It judged the patch suitable for a reviewable PR after final tests and rebuilt
check, not as a release or scientific calibration claim. Potential metadata
reconciliation across prediction routes was considered defensive hygiene;
no reachable wrong-status example was established.

The posted response numbers for lambda and data-only rows were checked
against smoke.log. Public GNN evidence was checked against the validation
article; private mammal results were attributed to the reporters rather than
presented as independent reproductions. All ten issue comments were read
back through the connector, and all ten issues remained open at 22:56:45 UTC.
