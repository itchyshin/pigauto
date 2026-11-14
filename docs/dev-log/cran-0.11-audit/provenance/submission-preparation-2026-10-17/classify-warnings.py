import collections
import csv
import json
from pathlib import Path

root = Path(__file__).resolve().parent
rules = [
    ('small_validation', 'Small validation set'),
    ('conformal_sample_limit', '95% split-conformal coverage is NOT achievable'),
    ('constant_matched_trait', 'At least one matched trait is constant.'),
    ('few_tree_draws', 'Few tree-sensitivity draws'),
    ('covariates_joint_baseline', "'covariates' are ignored by the joint MVN baseline"),
    ('covariates_gnn_off', 'gnn = FALSE (the default): user covariates'),
    ('imbalanced_ordinal', 'Imbalanced K>=3 ordinal trait'),
    ('signal_gate_fallback', 'phylo_signal_gate triggered'),
    ('species_dropped', '3 species in traits not found in tree -- dropped.'),
]

def classify(row):
    if row['message'] == 'the standard deviation is zero':
        assert row['last_call'].startswith('stats::cor('), row
        return 'undefined_constant_correlation'
    for category, prefix in rules:
        if row['message'].startswith(prefix):
            return category
    raise ValueError('UNCLASSIFIED WARNING: '+row['message'])

rows = list(csv.DictReader((root/'warning-details.csv').open()))
summary = json.loads((root/'warning-replay.json').read_text())
assert summary['archive_sha256'] == 'c2e768adad53baf6b772538e4c5a788161ae6d21ee029311654e9c28d3a5dfe2'
assert summary['pass'] == 3195 and summary['fail'] == 0 and summary['skips'] == 86
assert summary['warnings'] == len(rows) == 161
for row in rows:
    assert row['file'].startswith('test-') and row['test'] and int(row['line']) > 0
    row['category'] = classify(row)
try:
    classify(dict(message='planted unclassified warning',last_call=''))
except ValueError:
    negative = 'PASS: planted unknown warning refused'
else:
    raise AssertionError('unknown-warning control accepted')
with (root/'warning-classification.csv').open('w') as out:
    writer = csv.DictWriter(out, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
counts = dict(collections.Counter(row['category'] for row in rows))
report = {'scope':'one exact-archive local CRAN-layout replay; not a direct capture of Windows warning objects',
          'counts':counts,'total':len(rows),'unclassified':0,'negative_control':negative,
          'independent_consequence_review_receipt':'independent-warning-review.md'}
(root/'warning-classification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
print('ALL_161_LOCAL_WARNINGS_CLASSIFIED')
