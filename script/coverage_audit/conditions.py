#!/usr/bin/env python3
"""Inventory every warning and skip from testthat's installed-package Rout."""
import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path


def inventory(text):
    summary = re.search(r'\[ FAIL (\d+) \| WARN (\d+) \| SKIP (\d+) \| PASS (\d+) \]', text)
    if not summary:
        raise ValueError('Missing final testthat summary')
    counts = dict(zip(('fail', 'warn', 'skip', 'pass'), map(int, summary.groups())))
    skips = []
    section = re.search(r'══ Skipped tests \(\d+\).*?\n(.*?)(?=\n══|\Z)', text, re.S)
    if section:
        for block in re.split(r'(?m)^• ', section[1])[1:]:
            flat = ' '.join(block.split())
            match = re.match(r'(.*?) \((\d+)\): (.*)', flat)
            if not match:
                raise ValueError('Unrecognised skip block: ' + flat)
            reason, n, locations = match.groups()
            refs = re.findall(r"'([^']+)'", locations)
            if len(refs) != int(n):
                raise ValueError('Skip location count mismatch: ' + flat)
            if 'BIEN cache' in reason or 'fixture not found' in reason.lower():
                category = 'missing_fixture'
            elif 'HTML Tidy' in reason or reason == 'zip executable is not available':
                category = 'missing_optional_tool'
            elif 'installed' in reason and ('cannot test' in reason or 'integration test covers' in reason):
                category = 'dependency_absence_control'
            elif reason == 'NEWS.md not found' or any(s in reason for s in ('source', 'not shipped', 'script/gha/')):
                category = 'source_only_assertion'
            else:
                category = 'unclassified'
            skips.extend(dict(kind='skip', location=ref, category=category, message=reason) for ref in refs)
    warnings = []
    for match in re.finditer(r"── Warning \('([^']+)'\): (.*?) ─{2,}\n(.*?)(?=\n── Warning |\n\[ FAIL |\Z)", text, re.S):
        location, title, detail = match.groups()
        message = detail.split('\nBacktrace:')[0].strip()
        if message.startswith(('Small validation set', '95% split-conformal coverage')):
            category = 'small_calibration_fixture'
        elif any(s in message for s in ('standard deviation is zero', 'matched trait is constant')):
            category = 'constant_trait_fixture'
        elif message.startswith('fit_joint_solver:') and '(boom)' in message:
            category = 'deliberate_solver_failure_fixture'
        elif message.startswith('Few tree-sensitivity datasets'):
            category = 'small_tree_sensitivity_fixture'
        elif message.startswith(("'covariates' are omitted", 'gnn = FALSE:')):
            category = 'covariate_route_notice'
        elif message.startswith('Imbalanced K>=3 ordinal trait'):
            category = 'imbalanced_ordinal_fixture'
        elif message.startswith('phylo_signal_gate triggered'):
            category = 'signal_gate_fallback_notice'
        elif 'species in traits not found in tree' in message:
            category = 'unmatched_species_fixture'
        elif 'phytools' in message and ('not installed' in message or 'install' in message):
            category = 'optional_dependency_notice'
        elif any(s in message.lower() for s in ('singular', 'converg', 'positive definite')):
            category = 'numerical_warning_requires_review'
        else:
            category = 'unclassified'
        warnings.append(dict(kind='warning', location=location, category=category, message=message))
    if len(skips) != counts['skip'] or len(warnings) != counts['warn']:
        raise ValueError(f'Condition count mismatch: {len(skips)} skips, {len(warnings)} warnings, summary {counts}')
    return counts, skips + warnings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('rout', type=Path, nargs='?')
    parser.add_argument('out', type=Path, nargs='?')
    parser.add_argument('--selftest', action='store_true')
    args = parser.parse_args()
    if args.selftest:
        sample = ("══ Skipped tests (1) ══\n"
                  "• HTML Tidy is not installed (1): 'test-html.R:3:3'\n"
                  "══ Warnings ══\n"
                  "── Warning ('test-x.R:1:1'): constant ────────────\n"
                  "the standard deviation is zero\nBacktrace:\nfixture\n"
                  "[ FAIL 0 | WARN 1 | SKIP 1 | PASS 2 ]\n")
        counts, rows = inventory(sample)
        assert len(rows) == 2 and counts['pass'] == 2
        assert rows[0]['category'] == 'missing_optional_tool'
        assert rows[1]['message'] == 'the standard deviation is zero'
        zip_sample = sample.replace('HTML Tidy is not installed',
                                    'zip executable is not available')
        zip_counts, zip_rows = inventory(zip_sample)
        assert zip_counts == counts
        assert zip_rows[0]['category'] == 'missing_optional_tool'
        unknown_sample = sample.replace('HTML Tidy is not installed',
                                        'unrecognised prerequisite')
        assert inventory(unknown_sample)[1][0]['category'] == 'unclassified'
        for invalid in ('', sample.replace('WARN 1', 'WARN 2'),
                        sample.replace("'test-html.R:3:3'", '')):
            try:
                inventory(invalid)
            except ValueError:
                pass
            else:
                raise AssertionError('Incomplete condition inventory was accepted')
        print('CONDITIONS_SELFTEST_PASS')
        return
    if args.rout is None or args.out is None:
        parser.error('rout and out are required unless --selftest is used')
    counts, rows = inventory(args.rout.read_text())
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / 'conditions.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=('kind', 'location', 'category', 'message'))
        writer.writeheader(); writer.writerows(rows)
    result = dict(counts=counts, categories=dict(Counter(r['category'] for r in rows)),
                  source=str(args.rout.resolve()))
    (args.out / 'conditions-summary.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
