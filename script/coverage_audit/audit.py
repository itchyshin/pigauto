#!/usr/bin/env python3
"""Mechanical baseline inventory and fail-closed Cobertura audit for pigauto."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

BASE_SOURCE_SHA = "033649b93115c6d8191e791f06b1ddf52eb82d19"
GROUPS = {
    "presentation": {"plot.R", "plot_pigauto.R", "plot_missingness.R", "report.R"},
    "inputs": set("evaluate.R evaluate_imputation.R cross_validate.R preprocess_traits.R check_pigauto.R read_traits.R read_tree.R utils.R mask_missing.R result_accessors.R save_load.R gbif_centroids.R worldclim_covariates.R simulate_traits.R simulate_benchmark.R".split()),
    "gnn": set("utils_torch.R fit_pigauto.R fit_helpers.R zzz.R calibration_three_way.R".split()),
    "mi": set("multi_impute_analysis.R pool_mi.R mi_posterior.R draws_conditional.R with_imputations.R multi_impute_trees.R multi_impute.R".split()),
}
FIELDS = ["file", "line", "source_text", "enclosing_function", "relevant_existing_tests", "owner_group", "class", "disposition", "test_evidence"]


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def coverage(path: Path) -> tuple[dict[str, dict[int, int]], dict[str, int]]:
    root = ET.parse(path).getroot()
    files: dict[str, dict[int, int]] = {}
    for cls in root.findall(".//class"):
        name = cls.get("filename", "").replace("\\", "/").strip()
        if not name or name.startswith("/") or ".." in Path(name).parts:
            raise ValueError(f"empty or unsafe Cobertura filename: {name!r}")
        if not name.startswith("R/"):
            name = "R/" + name
        if not name.endswith(".R"):
            raise ValueError(f"non-R source path in Cobertura XML: {name}")
        dest = files.setdefault(name, {})
        for line in cls.findall("./lines/line"):
            n, hits = int(line.get("number", "0")), int(line.get("hits", "0"))
            if n < 1 or hits < 0:
                raise ValueError(f"invalid line record in {path}: {name}:{n} hits={hits}")
            # Cobertura may repeat line records in method and class summaries.
            # Merge them by source identity, taking the largest execution count.
            dest[n] = max(dest.get(n, 0), hits)
    totals = {"covered": int(root.get("lines-covered", "-1")), "valid": int(root.get("lines-valid", "-1"))}
    observed = [h for ls in files.values() for h in ls.values()]
    if totals["valid"] != len(observed) or totals["covered"] != sum(h > 0 for h in observed):
        raise ValueError(f"Cobertura summary disagrees with line records in {path}")
    return files, totals


def owner(file: str) -> str:
    base = Path(file).name
    for group, names in GROUPS.items():
        if base in names:
            return group
    return "prediction"


def function_ranges(root: Path, files: list[str]) -> dict[str, list[tuple[int, int, str]]]:
    """Return top-level function assignment ranges in one retained R parse pass."""
    rcode = r'''a <- commandArgs(TRUE)
for (path in a) {
  e <- parse(path, keep.source=TRUE); refs <- attr(e, "srcref")
  for (i in seq_along(e)) {
    x <- e[[i]]; fn <- ""
    if (is.call(x) && length(x) >= 3L && as.character(x[[1]]) %in% c("<-", "=") &&
        is.name(x[[2]]) && is.call(x[[3]]) && identical(x[[3]][[1]], as.name("function"))) fn <- as.character(x[[2]])
    r <- refs[[i]]; cat(path, "\t", r[[1]], "\t", r[[3]], "\t", if (nzchar(fn)) fn else "<toplevel>", "\n", sep="")
  }
}'''
    source_files = [root / f for f in files if (root / f).is_file()]
    if not source_files:
        return {}
    proc = subprocess.run(["Rscript", "-e", rcode, *map(str, source_files)], text=True, capture_output=True)
    if proc.returncode:
        raise ValueError(f"R parse failed: {proc.stderr.strip()}")
    out: dict[str, list[tuple[int, int, str]]] = {}
    for row in proc.stdout.splitlines():
        path, start, end, name = row.split("\t", 3)
        if name != "<toplevel>":
            rel = str(Path(path).relative_to(root))
            out.setdefault(rel, []).append((int(start), int(end), name))
    return out


def enclosing(ranges: list[tuple[int, int, str]], line: int) -> str:
    for start, end, name in ranges:
        if start <= line <= end:
            return name
    return "<toplevel>"


def test_hints(root: Path, file: str, fn: str, text: str, test_bodies: list[tuple[str, str]]) -> str:
    tests = root / "tests" / "testthat"
    if not tests.exists():
        return ""
    keys = [Path(file).stem, fn]
    keys = [k.lower() for k in keys if len(k) >= 4]
    matches = []
    for path, body in test_bodies:
        if any(k in body for k in keys):
            matches.append(path)
    return ";".join(matches[:8])


def inventory(args: argparse.Namespace) -> None:
    xml = Path(args.baseline)
    provenance_path = xml.with_name("baseline-provenance.json")
    prov = json.loads(provenance_path.read_text())
    actual_sha = sha(xml)
    if actual_sha != prov.get("baseline_xml_sha256"):
        raise ValueError("baseline XML hash does not match provenance")
    files, totals = coverage(xml)
    expected = (10845, 12749, 1904, 47)
    missing = [(f, n) for f, ls in files.items() for n, h in ls.items() if h == 0]
    missfiles = {f for f, _ in missing}
    actual = (totals["covered"], totals["valid"], len(missing), len(missfiles))
    if actual != expected or actual != tuple(prov.get(k) for k in ("expected_covered", "expected_total", "expected_missing", "expected_missing_files")):
        raise ValueError(f"baseline totals mismatch: got {actual}")
    root = Path(args.source_root or ".")
    tests_dir = root / "tests" / "testthat"
    test_bodies = [(str(p.relative_to(root)), p.read_text(errors="replace").lower()) for p in sorted(tests_dir.glob("*.R"))] if tests_dir.exists() else []
    rows = []
    range_cache = function_ranges(root, sorted({f for f, _ in missing}))
    for file, n in sorted(missing):
        src = root / file
        textlines = src.read_text(errors="replace").splitlines() if src.exists() else []
        txt = textlines[n-1].strip() if n <= len(textlines) else ""
        ranges = range_cache.get(file, [])
        fn = enclosing(ranges, n)
        # Classification remains an unreviewed routing hypothesis.
        rows.append({"file": file, "line": n, "source_text": txt, "enclosing_function": fn,
                     "relevant_existing_tests": test_hints(root, file, fn, txt, test_bodies), "owner_group": owner(file),
                     "class": "candidate reachable (provisional)", "disposition": "unresolved", "test_evidence": ""})
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS); w.writeheader(); w.writerows(rows)
    print(f"BASELINE_INVENTORY_PASS rows={len(rows)} files={len(missfiles)} covered={totals['covered']} valid={totals['valid']}")


def git_source_hashes(git_root: Path, commit: str, hashes: dict[str, str]) -> set[str]:
    if not re.fullmatch(r"[0-9a-fA-F]{40,64}", commit):
        raise ValueError("--git-root requires a full hexadecimal commit SHA")
    for rel in hashes:
        path = Path(rel)
        if not rel or path.is_absolute() or ".." in path.parts or "\\" in rel:
            raise ValueError(f"unsafe manifest source path: {rel!r}")
    listed = subprocess.run(["git", "-C", str(git_root), "ls-tree", "-r", "--name-only", commit], capture_output=True, check=True).stdout.decode().splitlines()
    tracked_r = {p for p in listed if p.startswith("R/") and p.endswith(".R")}
    required = {p for p in listed if (p.startswith("R/") and p.endswith(".R")) or (p.startswith("tests/") and p.endswith(".R")) or p in {"DESCRIPTION", "NAMESPACE", "codecov.yml", ".github/workflows/test-coverage.yaml"}}
    if not required <= set(hashes):
        missing = sorted(required - set(hashes))
        raise ValueError(f"manifest source_hashes omit files from expected git tree: {missing[:8]}")
    for rel, expected_digest in hashes.items():
        if not isinstance(expected_digest, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", expected_digest):
            raise ValueError(f"invalid SHA256 for manifest source path: {rel}")
        blob = subprocess.run(["git", "-C", str(git_root), "show", f"{commit}:{rel}"], capture_output=True)
        if blob.returncode:
            raise ValueError(f"manifest source path is absent from expected git tree: {rel}")
        if hashlib.sha256(blob.stdout).hexdigest() != expected_digest:
            raise ValueError(f"manifest source hash differs from expected git tree: {rel}")
    return tracked_r


def validate_torch_inventory(expected: object, observed: object) -> tuple[bool, dict[tuple[str, str], dict]]:
    if not isinstance(expected, list) or not expected or not isinstance(observed, list) or len(observed) != len(expected):
        return False, {}
    expected_by_key = {}
    for row in expected:
        if not isinstance(row, dict):
            return False, {}
        key = (row.get("factory"), row.get("method"))
        tracked = row.get("tracked_lines")
        first, last, file = row.get("first_line"), row.get("last_line"), row.get("file")
        if (not all(isinstance(x, str) and x for x in key) or not isinstance(file, str) or not file.startswith("R/") or
            not isinstance(first, int) or not isinstance(last, int) or first < 1 or last < first or
            not isinstance(tracked, list) or not tracked or any(not isinstance(n, int) or n < first or n > last for n in tracked)):
            return False, {}
        if key in expected_by_key:
            return False, {}
        expected_by_key[key] = row
    observed_by_key = {}
    for row in observed:
        if not isinstance(row, dict):
            return False, {}
        key = (row.get("factory"), row.get("method"))
        if key in observed_by_key:
            return False, {}
        observed_by_key[key] = row
    if set(expected_by_key) != set(observed_by_key):
        return False, {}
    for key, exp in expected_by_key.items():
        got = observed_by_key[key]
        expected_lines = set(exp["tracked_lines"])
        actual_lines = got.get("tracked_lines")
        if (not isinstance(got.get("file"), str) or got.get("file") != exp.get("file") or got.get("first_line") != exp.get("first_line") or
            got.get("last_line") != exp.get("last_line") or not isinstance(actual_lines, list) or
            any(not isinstance(n, int) or n < exp["first_line"] or n > exp["last_line"] for n in actual_lines) or not expected_lines <= set(actual_lines)):
            return False, {}
    return True, observed_by_key


def torch_method_gaps(expected: list[dict], observed_by_key: dict[tuple[str, str], dict], inventory_ok: bool,
                      candidate: dict[str, dict[int, int]], changed_sources: set[str]) -> tuple[list[str], list[dict]]:
    missing_lines = []
    missing_methods = []
    for exp in expected:
        key = (exp["factory"], exp["method"])
        obs = observed_by_key.get(key, {})
        line_numbers = sorted(set(exp.get("tracked_lines", [])) | set(obs.get("tracked_lines", [])))
        method_missing = not inventory_ok
        for line in line_numbers:
            hit = candidate.get(exp["file"], {}).get(line)
            if exp["file"] in changed_sources:
                state = "source_changed_unmapped"
                method_missing = True
                missing_lines.append({"factory": key[0], "method": key[1], "file": exp["file"], "line": line, "status": state})
            elif not inventory_ok:
                state = "inventory_unverified"
                missing_lines.append({"factory": key[0], "method": key[1], "file": exp["file"], "line": line, "status": state})
            elif hit is None or hit <= 0:
                state = "missing" if hit is None else "uncovered"
                method_missing = True
                missing_lines.append({"factory": key[0], "method": key[1], "file": exp["file"], "line": line, "status": state})
        if method_missing:
            missing_methods.append(f"{key[0]}.{key[1]}")
    return missing_methods, missing_lines


def verify(candidate: Path, manifest_path: Path, expected_sha: str, source_root: Path | None, git_root: Path | None, outdir: Path | None, require_complete: bool, allow_partial: bool) -> dict:
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("status") != "ok":
        raise ValueError("manifest status must be exactly 'ok'")
    if require_complete and not git_root:
        raise ValueError("--require-complete requires --git-root to bind source hashes to the expected commit")
    baseline = candidate.parent / "baseline.xml"
    baseline_prov_path = candidate.parent / "baseline-provenance.json"
    # Allow the canonical neighboring baseline path when the candidate is elsewhere.
    if not baseline.exists():
        baseline = Path("docs/dev-log/coverage-100/baseline.xml")
        baseline_prov_path = baseline.with_name("baseline-provenance.json")
    prov = json.loads(baseline_prov_path.read_text())
    if sha(baseline) != prov.get("baseline_xml_sha256"):
        raise ValueError("pinned baseline XML hash differs from provenance")
    if manifest.get("xml_sha256", manifest.get("xml_sha")) != sha(candidate):
        raise ValueError("candidate XML hash differs from manifest xml_sha256/xml_sha")
    if manifest.get("source_sha") != expected_sha:
        raise ValueError("manifest source_sha differs from expected SHA")
    if manifest.get("baseline_source_sha") != BASE_SOURCE_SHA or prov.get("baseline_source_sha") != BASE_SOURCE_SHA:
        raise ValueError("baseline source SHA is not the pinned 033649b baseline")
    baseline_configs = prov.get("coverage_config_sha256", {})
    if not isinstance(baseline_configs, dict) or not baseline_configs:
        raise ValueError("baseline provenance lacks coverage_config_sha256")
    approved_configs = prov.get("candidate_coverage_config_sha256")
    if approved_configs is not None:
        if (not isinstance(approved_configs, dict) or set(approved_configs) != set(baseline_configs) or
            any(not isinstance(v, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", v) for v in approved_configs.values())):
            raise ValueError("candidate_coverage_config_sha256 must whitelist all baseline config paths with SHA256 values")
        expected_configs = approved_configs
    else:
        expected_configs = baseline_configs
    hashes = manifest.get("source_hashes")
    baseline_r_hashes = prov.get("r_source_hashes")
    if not isinstance(baseline_r_hashes, dict) or not baseline_r_hashes:
        raise ValueError("baseline provenance lacks the complete r_source_hashes map")
    non_executable = set(prov.get("non_executable_sources", []))
    if not non_executable <= set(baseline_r_hashes):
        raise ValueError("baseline non_executable_sources contains unknown R paths")
    if not isinstance(hashes, dict):
        raise ValueError("manifest source_hashes must be an object")
    for rel in hashes:
        p = Path(rel)
        if not rel or p.is_absolute() or ".." in p.parts or "\\" in rel:
            raise ValueError(f"unsafe manifest source path: {rel!r}")
    base, _ = coverage(baseline)
    base_files = set(base)
    hash_files = {k for k in hashes if k.startswith("R/")} if isinstance(hashes, dict) else set()
    if not base_files <= set(baseline_r_hashes):
        raise ValueError("baseline XML contains an R source absent from r_source_hashes")
    if any(hashes.get(k) != v for k, v in expected_configs.items()):
        raise ValueError("manifest coverage config hashes differ from baseline provenance")
    if source_root:
        for rel, digest in hashes.items():
            p = source_root / rel
            if not p.is_file() or sha(p) != digest:
                raise ValueError(f"local source hash mismatch or missing file: {rel}")
    tracked_r = git_source_hashes(git_root, expected_sha, hashes) if git_root else set()
    candidate_r_hashes = {k: v for k, v in hashes.items() if k.startswith("R/")}
    same_r = {k for k in set(baseline_r_hashes) & set(candidate_r_hashes) if baseline_r_hashes[k] == candidate_r_hashes[k]}
    changed_r_files = sorted((set(baseline_r_hashes) | set(candidate_r_hashes)) - same_r)
    if git_root and tracked_r != set(baseline_r_hashes):
        changed_r_files = sorted(set(changed_r_files) | (tracked_r ^ set(baseline_r_hashes)))
    cand, totals = coverage(candidate)
    if not set(cand) <= hash_files:
        raise ValueError("candidate contains R source files absent from manifest source_hashes")
    expected_methods = prov.get("expected_torch_methods")
    torch_inventory_ok, observed_methods = validate_torch_inventory(expected_methods, manifest.get("torch_method_inventory"))
    safe_expected_methods = expected_methods if isinstance(expected_methods, list) else []
    missing_torch_methods, missing_torch_method_lines = torch_method_gaps(safe_expected_methods, observed_methods, torch_inventory_ok, cand, set(changed_r_files))
    executable_sources = set(baseline_r_hashes) - non_executable
    missing_source_files = sorted(executable_sources - set(cand))
    failures = []
    if not git_root:
        failures.append(("", 0, "git_source_binding_missing"))
    failures.extend((file, 0, "source_changed_unmapped") for file in changed_r_files)
    failures.extend((file, 0, "missing_source_file") for file in missing_source_files)
    if not torch_inventory_ok:
        failures.append(("", 0, "torch_method_inventory_unverified"))
    failures.extend((row["file"], row["line"], f"torch_method_{row['status']}") for row in missing_torch_method_lines)
    regressions = []
    for f, lines in base.items():
        for n, h in lines.items():
            if f in changed_r_files:
                continue  # No line credit without a source map after line shifts.
            if f not in cand or n not in cand[f] or cand[f][n] <= 0:
                category = "unresolved" if h == 0 else "regression"
                if category == "regression": regressions.append((f, n, "regression"))
                failures.append((f, n, category))
    # Newly executable lines are also required to be covered.
    new_uncovered = [(f, n, "new_uncovered") for f, ls in cand.items() for n, h in ls.items() if f not in base or n not in base[f] if h <= 0]
    failures.extend(new_uncovered)
    base_ids = {(f, n) for f, ls in base.items() for n in ls}
    cand_ids = {(f, n) for f, ls in cand.items() for n in ls}
    source_identity = {f for f in base_files if f not in changed_r_files}
    source_identity_ok = bool(git_root) and not changed_r_files and tracked_r == set(baseline_r_hashes)
    complete = source_identity_ok and torch_inventory_ok and not missing_torch_methods and not failures and not missing_source_files and totals["valid"] > 0 and totals["covered"] == totals["valid"] and len(base_ids - cand_ids) == 0
    summary = {"candidate_xml_sha256": sha(candidate), "source_sha": manifest["source_sha"], "covered": totals["covered"], "valid": totals["valid"],
               "baseline_lines": len(base_ids), "baseline_lines_covered": sum(1 for f,n in base_ids if f in source_identity and f in cand and cand[f].get(n, 0) > 0),
               "baseline_missed_lines": sum(1 for f,ls in base.items() for h in ls.values() if h == 0),
               "baseline_missed_lines_covered": sum(1 for f,ls in base.items() for n,h in ls.items() if h == 0 and f in source_identity and f in cand and cand[f].get(n, 0) > 0),
               "new_lines": len(cand_ids - base_ids), "failures": len(failures), "regressions": len(regressions),
               "source_changed_files": changed_r_files, "source_identity_verified": source_identity_ok,
               "missing_source_files": missing_source_files, "torch_method_inventory_verified": torch_inventory_ok,
               "missing_torch_methods": missing_torch_methods, "missing_torch_method_lines": missing_torch_method_lines,
               "complete": complete}
    if failures and allow_partial:
        summary["partial_report"] = True
    if outdir:
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        with (outdir / "disposition.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f); w.writerow(["file", "line", "disposition", "test_evidence"])
            for file, lines in base.items():
                for line, hits in lines.items():
                    if hits == 0:
                        covered = file not in changed_r_files and file in cand and cand[file].get(line, 0) > 0
                        evidence = f"integrated suite@{manifest['source_sha']}" if covered else ""
                        w.writerow([file, line, "covered" if covered else "unresolved", evidence])
            for row in regressions:
                w.writerow([*row, ""])
            for row in new_uncovered:
                w.writerow([*row, ""])
            for file in missing_source_files:
                w.writerow([file, 0, "missing_source_file", ""])
            for row in missing_torch_method_lines:
                w.writerow([row["file"], row["line"], f"torch_method_{row['status']}", ""])
    if failures and not allow_partial:
        raise ValueError(f"coverage failures={len(failures)} (see disposition.csv if --out supplied)")
    if require_complete and not summary["complete"]:
        raise ValueError("require-complete failed: candidate does not cover all executable lines")
    return summary


def selftest() -> None:
    # Small, generated XML fixtures exercise verifier semantics independently of
    # the real baseline's line count.
    with tempfile.TemporaryDirectory(prefix="pigauto-audit-") as td:
        d = Path(td)
        def xml(path: Path, entries: dict[str, dict[int,int]], declared: tuple[int,int] | None=None) -> None:
            covered=sum(h>0 for ls in entries.values() for h in ls.values()); valid=sum(len(ls) for ls in entries.values())
            if declared is not None: covered, valid=declared
            root=ET.Element("coverage", {"lines-covered":str(covered), "lines-valid":str(valid)})
            packages=ET.SubElement(root,"packages"); package=ET.SubElement(packages,"package"); classes=ET.SubElement(package,"classes")
            for fn, ls in entries.items():
                c=ET.SubElement(classes,"class",{"filename":fn}); lines=ET.SubElement(c,"lines")
                for n,h in ls.items(): ET.SubElement(lines,"line",{"number":str(n),"hits":str(h)})
            ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)
        checkout=d/"git-checkout"
        for rel, content in {"R/a.R":"a <- 1\n", "R/b.R":"b <- 2\n", "R/c.R":"c <- 3\n", "tests/testthat/test_a.R":"test_that('a', {})\n", "DESCRIPTION":"Package: fixture\n", "NAMESPACE":"", "codecov.yml":"coverage: fixture\n", ".github/workflows/test-coverage.yaml":"name: fixture\n"}.items():
            p=checkout/rel; p.parent.mkdir(parents=True,exist_ok=True); p.write_text(content)
        subprocess.run(["git","-C",str(checkout),"init","-q"],check=True)
        subprocess.run(["git","-C",str(checkout),"config","user.email","audit@example.test"],check=True)
        subprocess.run(["git","-C",str(checkout),"config","user.name","Audit Selftest"],check=True)
        subprocess.run(["git","-C",str(checkout),"add","."],check=True)
        subprocess.run(["git","-C",str(checkout),"commit","-qm","fixture"],check=True)
        git_sha=subprocess.run(["git","-C",str(checkout),"rev-parse","HEAD"],check=True,text=True,capture_output=True).stdout.strip()
        config_hashes={p:hashlib.sha256((checkout/p).read_bytes()).hexdigest() for p in ("codecov.yml",".github/workflows/test-coverage.yaml")}
        # Class-level lines are authoritative. Method summaries can repeat or
        # misattribute locations when different classes share a method name.
        noise=d/"method-summary.xml"; xml(noise,{"R/a.R":{1:1,2:0}})
        tree=ET.parse(noise); cls=tree.find(".//class")
        methods=ET.SubElement(cls,"methods"); method=ET.SubElement(methods,"method")
        lines=ET.SubElement(method,"lines")
        ET.SubElement(lines,"line",{"number":"999","hits":"10"})
        tree.write(noise)
        assert coverage(noise)[1] == {"covered":1,"valid":2}
        b=d/"baseline.xml"; xml(b,{"R/a.R":{1:1,2:0},"R/b.R":{1:0}})
        prov_path=d/"baseline-provenance.json"
        baseline_r_hashes={p:hashlib.sha256((checkout/p).read_bytes()).hexdigest() for p in ("R/a.R","R/b.R","R/c.R")}
        torch_methods=[{"factory":"Toy","method":"forward","file":"R/a.R","first_line":1,"last_line":2,"tracked_lines":[1,2]}]
        prov={"baseline_source_sha":BASE_SOURCE_SHA,"baseline_xml_sha256":sha(b),"coverage_config_sha256":config_hashes,"r_source_hashes":baseline_r_hashes,"non_executable_sources":[],"expected_torch_methods":torch_methods}
        prov_path.write_text(json.dumps(prov))
        hashes={p:hashlib.sha256((checkout/p).read_bytes()).hexdigest() for p in ("R/a.R","R/b.R","R/c.R","tests/testthat/test_a.R","DESCRIPTION","NAMESPACE","codecov.yml",".github/workflows/test-coverage.yaml")}
        def candidate(name: str, entries: dict[str,dict[int,int]], manifest_changes: dict | None=None) -> tuple[Path,Path]:
            x=d/name; xml(x,entries); m={"source_sha":"right","xml_sha256":sha(x),"baseline_source_sha":BASE_SOURCE_SHA,"source_hashes":hashes,"status":"ok","torch_method_inventory":torch_methods}; m.update(manifest_changes or {})
            mp=d/(name+".json"); mp.write_text(json.dumps(m)); return x,mp
        good=candidate("good.xml",{"R/a.R":{1:3,2:1},"R/b.R":{1:2},"R/c.R":{1:4}})
        gm=json.loads(good[1].read_text()); gm["source_sha"]=git_sha; good[1].write_text(json.dumps(gm))
        assert verify(*good,git_sha,None,checkout,None,True,False)["complete"]
        try: verify(*good,git_sha,checkout,None,None,True,False)
        except ValueError: pass
        else: raise AssertionError("--require-complete passed without --git-root")
        missing_class=candidate("missing-class.xml",{"R/a.R":{1:3,2:1},"R/b.R":{1:2}})
        missing_manifest=json.loads(missing_class[1].read_text()); missing_manifest["source_sha"]=git_sha; missing_class[1].write_text(json.dumps(missing_manifest))
        missing_result=verify(*missing_class,git_sha,None,checkout,d/"missing-class-report",False,True)
        assert not missing_result["complete"] and missing_result["missing_source_files"]==["R/c.R"]
        assert any(r["file"]=="R/c.R" and r["disposition"]=="missing_source_file" for r in csv.DictReader((d/"missing-class-report"/"disposition.csv").open()))
        partial_method=candidate("partial-method.xml",{"R/a.R":{1:3},"R/b.R":{1:2},"R/c.R":{1:4}})
        pm=json.loads(partial_method[1].read_text()); pm["source_sha"]=git_sha; partial_method[1].write_text(json.dumps(pm))
        pm_result=verify(*partial_method,git_sha,None,checkout,d/"partial-method-report",False,True)
        assert "Toy.forward" in pm_result["missing_torch_methods"] and any(x["line"]==2 and x["status"]=="missing" for x in pm_result["missing_torch_method_lines"])
        bad_method_inventory=candidate("bad-method-inventory.xml",{"R/a.R":{1:3,2:1},"R/b.R":{1:2},"R/c.R":{1:4}})
        bi=json.loads(bad_method_inventory[1].read_text()); bi["source_sha"]=git_sha; bi["torch_method_inventory"]=[]; bad_method_inventory[1].write_text(json.dumps(bi))
        bi_result=verify(*bad_method_inventory,git_sha,None,checkout,d/"bad-method-inventory-report",False,True)
        assert not bi_result["torch_method_inventory_verified"] and "Toy.forward" in bi_result["missing_torch_methods"]
        bad_tree=d/"bad-tree.json"; bad_manifest={**gm,"source_hashes":{**hashes,"R/a.R":"0"*64}}; bad_tree.write_text(json.dumps(bad_manifest))
        try: verify(good[0],bad_tree,git_sha,None,checkout,None,False,False)
        except ValueError: pass
        else: raise AssertionError("manifest content hash not bound to git tree")
        path_manifest=d/"path-traversal.json"; path_manifest.write_text(json.dumps({**gm,"source_hashes":{**hashes,"../escape":"0"*64}}))
        try: verify(good[0],path_manifest,git_sha,None,checkout,None,False,False)
        except ValueError: pass
        else: raise AssertionError("manifest traversal path unexpectedly passed")
        (checkout/"codecov.yml").write_text("coverage: reviewed-addition\n")
        subprocess.run(["git","-C",str(checkout),"add","codecov.yml"],check=True)
        subprocess.run(["git","-C",str(checkout),"commit","-qm","reviewed config fixture"],check=True)
        config_sha=subprocess.run(["git","-C",str(checkout),"rev-parse","HEAD"],check=True,text=True,capture_output=True).stdout.strip()
        reviewed_config_hash=hashlib.sha256((checkout/"codecov.yml").read_bytes()).hexdigest()
        current_hashes={**hashes,"codecov.yml":reviewed_config_hash}
        config_manifest={**gm,"source_sha":config_sha,"source_hashes":current_hashes}
        config_manifest_path=d/"config-change.json"; config_manifest_path.write_text(json.dumps(config_manifest))
        try: verify(good[0],config_manifest_path,config_sha,None,checkout,None,False,False)
        except ValueError: pass
        else: raise AssertionError("unwhitelisted config change unexpectedly passed")
        prov["candidate_coverage_config_sha256"]={**config_hashes,"codecov.yml":reviewed_config_hash}
        prov_path.write_text(json.dumps(prov))
        assert verify(good[0],config_manifest_path,config_sha,None,checkout,None,True,False)["complete"]
        (checkout/".github/workflows/test-coverage.yaml").write_text("name: unapproved\n")
        subprocess.run(["git","-C",str(checkout),"add",".github/workflows/test-coverage.yaml"],check=True)
        subprocess.run(["git","-C",str(checkout),"commit","-qm","unapproved config fixture"],check=True)
        unapproved_sha=subprocess.run(["git","-C",str(checkout),"rev-parse","HEAD"],check=True,text=True,capture_output=True).stdout.strip()
        unapproved_hashes={**current_hashes,".github/workflows/test-coverage.yaml":hashlib.sha256((checkout/".github/workflows/test-coverage.yaml").read_bytes()).hexdigest()}
        unapproved_path=d/"unapproved-config.json"; unapproved_path.write_text(json.dumps({**config_manifest,"source_sha":unapproved_sha,"source_hashes":unapproved_hashes}))
        try: verify(good[0],unapproved_path,unapproved_sha,None,checkout,None,False,False)
        except ValueError: pass
        else: raise AssertionError("config outside candidate whitelist unexpectedly passed")
        subprocess.run(["git","-C",str(checkout),"reset","--hard",config_sha],check=True,capture_output=True)
        hashes=current_hashes
        gm["source_sha"]="right"; good[1].write_text(json.dumps(gm))
        (checkout/"R/a.R").write_text("a <- 1L\n")
        subprocess.run(["git","-C",str(checkout),"add","R/a.R"],check=True)
        subprocess.run(["git","-C",str(checkout),"commit","-qm","changed source"],check=True)
        changed_sha=subprocess.run(["git","-C",str(checkout),"rev-parse","HEAD"],check=True,text=True,capture_output=True).stdout.strip()
        changed_hashes={**hashes,"R/a.R":hashlib.sha256((checkout/"R/a.R").read_bytes()).hexdigest()}
        hashes=changed_hashes
        changed_manifest={**gm,"source_sha":changed_sha,"source_hashes":changed_hashes}
        changed_manifest_path=d/"changed-source.json"; changed_manifest_path.write_text(json.dumps(changed_manifest))
        changed_result=verify(good[0],changed_manifest_path,changed_sha,None,checkout,d/"changed-source-report",False,True)
        baseline_b_lines=sum(len(v) for k,v in {"R/b.R":{1:0}}.items())
        assert changed_result["source_changed_files"]==["R/a.R"] and not changed_result["complete"]
        assert changed_result["baseline_lines_covered"]==baseline_b_lines
        changed_rows=list(csv.DictReader((d/"changed-source-report"/"disposition.csv").open()))
        assert any(r["file"]=="R/a.R" and r["line"]=="2" and r["disposition"]=="unresolved" and not r["test_evidence"] for r in changed_rows)
        assert any(x["file"]=="R/a.R" and x["status"]=="source_changed_unmapped" for x in changed_result["missing_torch_method_lines"])
        (checkout/"R/c.R").unlink()
        subprocess.run(["git","-C",str(checkout),"add","-u","R/c.R"],check=True)
        subprocess.run(["git","-C",str(checkout),"commit","-qm","remove source fixture"],check=True)
        deleted_sha=subprocess.run(["git","-C",str(checkout),"rev-parse","HEAD"],check=True,text=True,capture_output=True).stdout.strip()
        deleted_hashes={k:v for k,v in changed_hashes.items() if k!="R/c.R"}
        deleted_xml=d/"deleted-source.xml"; xml(deleted_xml,{"R/a.R":{1:3,2:1},"R/b.R":{1:2}})
        deleted_manifest={**changed_manifest,"source_sha":deleted_sha,"xml_sha256":sha(deleted_xml),"source_hashes":deleted_hashes}
        deleted_manifest_path=d/"deleted-source.json"; deleted_manifest_path.write_text(json.dumps(deleted_manifest))
        deleted_result=verify(deleted_xml,deleted_manifest_path,deleted_sha,None,checkout,d/"deleted-source-report",False,True)
        assert "R/c.R" in deleted_result["source_changed_files"] and deleted_result["missing_source_files"]==["R/c.R"]
        partial=candidate("partial.xml",{"R/a.R":{1:3,2:0},"R/b.R":{1:2}})
        partial_out=d/"partial-report"
        result=verify(*partial,"right",None,None,partial_out,False,True)
        assert result.get("partial_report") and not result["complete"]
        report=json.loads((partial_out/"summary.json").read_text())
        dispositions=list(csv.DictReader((partial_out/"disposition.csv").open()))
        assert report.get("partial_report") and any(r["disposition"]=="covered" and r["test_evidence"]=="integrated suite@right" for r in dispositions)
        missing_method=candidate("missing-method-line.xml",{"R/a.R":{1:3},"R/b.R":{1:2},"R/c.R":{1:4}})
        mr=json.loads(missing_method[1].read_text()); mr["source_sha"]=changed_sha; missing_method[1].write_text(json.dumps(mr))
        method_result=verify(*missing_method,changed_sha,None,checkout,d/"missing-method-report",False,True)
        assert "Toy.forward" in method_result["missing_torch_methods"] and any(x["line"]==2 and x["status"]=="source_changed_unmapped" for x in method_result["missing_torch_method_lines"])
        unverified_inventory=candidate("unverified-method-inventory.xml",{"R/a.R":{1:3,2:1},"R/b.R":{1:2},"R/c.R":{1:4}})
        ui=json.loads(unverified_inventory[1].read_text()); ui["source_sha"]=changed_sha; ui.pop("torch_method_inventory"); unverified_inventory[1].write_text(json.dumps(ui))
        inventory_result=verify(*unverified_inventory,changed_sha,None,checkout,d/"unverified-method-report",False,True)
        assert not inventory_result["torch_method_inventory_verified"] and "Toy.forward" in inventory_result["missing_torch_methods"]
        cases=[
            ("uncovered.xml",{"R/a.R":{1:3,2:0},"R/b.R":{1:2}},None),
            ("missingfile.xml",{"R/a.R":{1:3,2:1}},None),
            ("missingxml.xml",{"R/a.R":{1:3,2:1},"R/b.R":{}},None),
            ("badhash.xml",{"R/a.R":{1:3,2:1},"R/b.R":{1:2}}, {"xml_sha256":"bad"}),
            ("wrongsha.xml",{"R/a.R":{1:3,2:1},"R/b.R":{1:2}}, {"source_sha":"wrong"}),
            ("failedstatus.xml",{"R/a.R":{1:3,2:1},"R/b.R":{1:2}}, {"status":"failed"}),
            ("runningstatus.xml",{"R/a.R":{1:3,2:1},"R/b.R":{1:2}}, {"status":"running"}),
            ("missingstatus.xml",{"R/a.R":{1:3,2:1},"R/b.R":{1:2}}, {"status":None}),
            ("removed.xml",{"R/a.R":{1:3},"R/b.R":{1:2}},None),
            ("newuncovered.xml",{"R/a.R":{1:3,2:1,3:0},"R/b.R":{1:2}},None),
            ("forgedcounts.xml",{"R/a.R":{1:3,2:1},"R/b.R":{1:2}}, {"_declared_counts":[0,0]}),
            ("omittedhash.xml",{"R/a.R":{1:3,2:1},"R/b.R":{1:2}}, {"source_hashes":{"R/a.R":"ha","codecov.yml":"a",".github/workflows/test-coverage.yaml":"b"}}),
        ]
        for name, entries, changes in cases:
            changes=changes or {}
            declared=changes.pop("_declared_counts",None)
            args=candidate(name,entries,changes if changes else None)
            if declared is not None:
                xml(args[0],entries,tuple(declared)); m=json.loads(args[1].read_text()); m["xml_sha256"]=sha(args[0]); args[1].write_text(json.dumps(m))
            try: verify(*args,"right",None,None,None,False,False)
            except (ValueError, ET.ParseError): pass
            else: raise AssertionError(f"negative fixture unexpectedly passed: {name}")
        missing=d/"missing.xml"
        try: verify(missing,good[1],"right",None,None,None,False,False)
        except (OSError, ValueError, ET.ParseError): pass
        else: raise AssertionError("missing candidate XML unexpectedly passed")
        bad=d/"not-xml.xml"; bad.write_text("<coverage>")
        try: coverage(bad)
        except ET.ParseError: pass
        else: raise AssertionError("missing/malformed XML unexpectedly passed")
        prov_path.write_text(json.dumps({**prov,"baseline_xml_sha256":"wrong"}))
        try: verify(*good,"right",None,None,None,False,False)
        except ValueError: pass
        else: raise AssertionError("bad pinned baseline XML hash unexpectedly passed")
    print("AUDIT_SELFTEST_PASS")


def main() -> int:
    p=argparse.ArgumentParser(); sub=p.add_subparsers(dest="cmd", required=True)
    q=sub.add_parser("inventory"); q.add_argument("--baseline",required=True); q.add_argument("--out",required=True); q.add_argument("--source-root")
    q=sub.add_parser("verify"); q.add_argument("--candidate",required=True); q.add_argument("--manifest",required=True); q.add_argument("--expected-sha",required=True); q.add_argument("--source-root"); q.add_argument("--git-root"); q.add_argument("--out"); q.add_argument("--require-complete",action="store_true"); q.add_argument("--allow-partial",action="store_true")
    sub.add_parser("selftest")
    a=p.parse_args()
    try:
        if a.cmd=="inventory": inventory(a)
        elif a.cmd=="selftest": selftest()
        else:
            result=verify(Path(a.candidate),Path(a.manifest),a.expected_sha,Path(a.source_root) if a.source_root else None,Path(a.git_root) if a.git_root else None,Path(a.out) if a.out else None,a.require_complete,a.allow_partial)
            print(("AUDIT_PARTIAL " if result.get("partial_report") else "") + json.dumps(result,sort_keys=True))
        return 0
    except (OSError, ValueError, ET.ParseError, KeyError) as e:
        print(f"AUDIT_FAIL: {e}",file=sys.stderr); return 1

if __name__=="__main__": raise SystemExit(main())
