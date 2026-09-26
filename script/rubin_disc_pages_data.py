#!/usr/bin/env python3
# script/rubin_disc_pages_data.py
#
# Condenses the discrete results of the Rubin study (script/rubin_study/data/discrete/, written by
# script/rubin_discrete_aggregate.R) into one JSON object for the campaign pages (report, summary, accuracy), which
# script/rubin_report_build.sh injects at __DISC__. Numbers are rounded to 5 places.
#   python3 script/rubin_disc_pages_data.py [STUDY_DATA_DIR] [OUT_JSON]
import csv, gzip, json, sys, os

data = sys.argv[1] if len(sys.argv) > 1 else "script/rubin_study/data"
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(data, "discrete", "pages.json")
disc = os.path.join(data, "discrete")


def rows(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", newline="") as f:
        return list(csv.DictReader(f))


def num(v):
    if v in (None, "", "NA", "NaN"):
        return None
    try:
        x = float(v)
    except ValueError:
        return v
    return int(x) if x.is_integer() and abs(x) < 1e6 else round(x, 5)


def pick(rs, keys):
    return [{k: num(r.get(k)) if k not in ("arm", "trait", "contrast", "measure") else r.get(k) for k in keys} for r in rs]


cells_keys = ["n", "lambda", "arm", "trait", "n_fits", "n_fits_trait_absent", "accuracy", "accuracy_se", "brier",
              "brier_se", "set_coverage", "set_coverage_se", "set_size", "set_size_se", "mae_class"]
D = {
    "cl": pick(rows(os.path.join(disc, "agg_disc_cells_l.csv")), cells_keys),
    "cr": pick(rows(os.path.join(disc, "agg_disc_cells.csv")), ["rho"] + cells_keys),
    "dn": pick(rows(os.path.join(disc, "agg_disc_down.csv")),
               ["n", "lambda", "rho", "arm", "n_fits", "n_excluded", "n_undefined", "n_missing", "truth",
                "coverage", "coverage_se", "coverage_cond", "coverage_cond_se", "complete_coverage_cond",
                "bias_cond", "bias_cond_se", "diff_cd_mean", "diff_cd_sd"]),
    "pr": pick(rows(os.path.join(disc, "agg_disc_paired.csv")),
               ["contrast", "measure", "n", "lambda", "rho", "trait", "diff", "se", "datasets"]),
}

# meta: counts, BACE failures (and overlap with the continuous run), NA draws, reproduction
fits = rows(os.path.join(disc, "fits.csv"))
fc = rows(os.path.join(disc, "fit_disc_cells.csv.gz"))
fail = sorted({r["tag"] for r in fc if r["arm"] == "bace" and r["status"] == "missing" and "bace_fit" in (r["errors"] or "")})
chain_only = sorted({r["tag"] for r in fc if r["arm"] == "bace_chain" and r["status"] == "missing"
                     and "bace_fit" not in (r["errors"] or "") and r["errors"] not in ("", "NA")})
cont = rows(os.path.join(data, "fits.csv"))
cont_fail = {r["tag"] for r in cont if r["set"] == "bace" and "bace_fit" in (r["errors"] or "") and int(float(r["n"])) <= 300}
na = {}
for arm in ("bace", "bace_chain"):
    s = [r for r in fc if r["arm"] == arm and r["trait"] == "cat3" and r["status"] == "scored"]
    draws = sum((float(r["n_cells"]) + float(r["n_cells_all_na"] or 0)) * 20 for r in s)
    na[arm] = round(sum(float(r["n_na"]) for r in s) / draws, 5) if draws else None
rp = rows(os.path.join(disc, "repro.csv"))
rb = [r for r in rp if r["set"] == "bace" and r["arm"] == "all" and r["host"] != "all"]
same = [r for r in rb if r["host"] == r["cont_host"]]
cross = [r for r in rb if r["host"] != r["cont_host"]]
rf = [r for r in rp if r["set"] == "freq" and r["arm"] == "all" and r["host"] == "all"]
truth = rows(os.path.join(disc, "truth_slope_c1_bin.csv"))
D["meta"] = {
    "bace_files": sum(1 for r in fits if r["set"] == "bace"),
    "freq_files": sum(1 for r in fits if r["set"] == "freq"),
    "bace_fit_fail": len(fail),
    "bace_fit_fail_n": {n: sum(1 for t in fail if f"_n{n}_" in t) for n in ("100", "300")},
    "bace_fit_fail_lambda1": all("_l1_" in t for t in fail),
    "fail_both_runs": len(set(fail) & cont_fail),
    "chain_only_fail": len(chain_only),
    "na_cat3": na,
    "repro": {"same_n": sum(int(float(r["n_fits"])) for r in same), "same_ident": sum(int(float(r["n_fits_identical"])) for r in same),
              "cross_n": sum(int(float(r["n_fits"])) for r in cross), "cross_ident": sum(int(float(r["n_fits_identical"])) for r in cross),
              "cross_median_min": min(float(r["median_abs_estimate"]) for r in cross) if cross else None,
              "cross_median_max": max(float(r["median_abs_estimate"]) for r in cross) if cross else None,
              "freq_n": int(float(rf[0]["n_fits"])) if rf else None, "freq_ident": int(float(rf[0]["n_fits_identical"])) if rf else None},
    "truth_R": sorted({int(float(r["R"])) for r in truth}),
}
with open(out, "w") as f:
    json.dump(D, f, separators=(",", ":"))
print(out, os.path.getsize(out), "bytes;", len(D["cl"]), "cell rows;", len(D["dn"]), "downstream rows;", json.dumps(D["meta"]))
