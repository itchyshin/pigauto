"""Descriptive audit of the public prediction attachment in issue #260.

Usage: python3 summarize_seed_predictions.py gnn_predictions.csv output.json
No refitting, baseline comparison, confidence intervals or ensemble validation.
"""
import csv
import hashlib
import json
import math
import statistics
import sys
from pathlib import Path


def summarize(path):
    rows = list(csv.DictReader(Path(path).open()))
    grouped = {}
    for row in rows:
        run = int(row["run"])
        seed_group = int(row["seed"]) - run
        key = (seed_group, int(row["row_n"]))
        # The issue's five folds repeat some rows in fold 5. Retain the
        # earlier fold, identically for every seed, rather than double count.
        if key in grouped:
            old = grouped[key]
            if run != 5 or old["run"] == 5:
                raise ValueError("Unexpected duplicate row/seed")
            if row["species"] != old["species"] or float(row["observed"]) != old["observed"]:
                raise ValueError("Duplicate observations disagree")
            continue
        values = {k: float(row[k]) for k in ("observed", "predicted", "lower", "upper")}
        if not all(math.isfinite(x) for x in values.values()):
            raise ValueError("Non-finite prediction or truth")
        if values["lower"] > values["upper"]:
            raise ValueError("Reversed interval")
        grouped[key] = dict(values, species=row["species"], run=run)
    seeds = sorted({k[0] for k in grouped})
    ids = sorted({k[1] for k in grouped})
    if len(seeds) != 3 or len(grouped) != len(seeds) * len(ids):
        raise ValueError("Expected a complete three-seed panel")
    for i in ids:
        first = grouped[seeds[0], i]
        for s in seeds[1:]:
            r = grouped[s, i]
            if (r["observed"], r["species"], r["run"]) != (first["observed"], first["species"], first["run"]):
                raise ValueError("Seed panels do not score identical cells")
    spreads = [max(grouped[s, i]["predicted"] for s in seeds) -
               min(grouped[s, i]["predicted"] for s in seeds) for i in ids]
    scores = {}
    for s in seeds:
        rr = [grouped[s, i] for i in ids]
        scores[str(s)] = dict(
            rmse=math.sqrt(statistics.mean((r["predicted"] - r["observed"]) ** 2 for r in rr)),
            covered=sum(r["lower"] <= r["observed"] <= r["upper"] for r in rr))
    return dict(source_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                input_rows=len(rows), scored_cells=len(ids),
                scored_species=len({grouped[seeds[0], i]["species"] for i in ids}),
                median_seed_range=statistics.median(spreads),
                cells_range_over_3=sum(x > 3 for x in spreads),
                seed_scores=scores,
                limits="Published cached predictions from old source; no current-source fit, independent validation, baseline comparison or ensemble recommendation.")


if __name__ == "__main__":
    result = summarize(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(result, indent=2) + "\n")
