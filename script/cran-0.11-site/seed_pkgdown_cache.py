#!/usr/bin/env python3
"""Seed pkgdown's offline JS cache from a prior site after SRI checks."""

import argparse
import base64
import hashlib
import shutil
from pathlib import Path

ASSETS = {
    "headroom/0.11.0/headroom.min.js": (
        "headroom-0.11.0/headroom.min.js",
        "sha256-AsUX4SJE1+yuDu5+mAVzJbuYNPHj/WroHuZ8Ir/CkE0=",
    ),
    "headroom/0.11.0/jQuery.headroom.min.js": (
        "headroom-0.11.0/jQuery.headroom.min.js",
        "sha256-ZX/yNShbjqsohH1k95liqY9Gd8uOiE1S4vZc+9KQ1K4=",
    ),
    "bootstrap-toc/1.0.1/bootstrap-toc.min.js": (
        "bootstrap-toc-1.0.1/bootstrap-toc.min.js",
        "sha256-4veVQbu7//Lk5TSmc7YV48MxtMy98e26cf5MrgZYnwo=",
    ),
    "clipboard.js/2.0.11/clipboard.min.js": (
        "clipboard.js-2.0.11/clipboard.min.js",
        "sha512-7O5pXpc0oCRrxk8RUfDYFgn0nO1t+jLuIOQdOMRp4APB7uZ4vSjspzp5y6YDtDs4VzUSTbWzBFZ/LKJhnyFOKw==",
    ),
    "search/1.0.0/fuse.min.js": (
        "search-1.0.0/fuse.min.js",
        "sha512-KnvCNMwWBGCfxdOtUpEtYgoM59HHgjHnsVGSxxgz7QH1DYeURk+am9p3J+gsOevfE29DV0V+/Dd52ykTKxN5fA==",
    ),
    "search/1.0.0/autocomplete.jquery.min.js": (
        "search-1.0.0/autocomplete.jquery.min.js",
        "sha512-GU9ayf+66Xx2TmpxqJpliWbT5PiGYxpaG8rfnBEk1LL8l1KGkRShhngwdXK1UgqhAzWpZHSiYPc09/NwDQIGyg==",
    ),
    "search/1.0.0/mark.min.js": (
        "search-1.0.0/mark.min.js",
        "sha512-5CYOlHXGh6QpOFA/TeTylKLWfB3ftPsde7AnmhuitiTX4K5SqCLBeKro6sPS8ilsz1Q4NRx3v8Ko2IBiszzdww==",
    ),
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-site", required=True, type=Path)
    parser.add_argument("--cache-root", required=True, type=Path)
    args = parser.parse_args()

    source = args.source_site / "deps"
    for cache_rel, (source_rel, sri) in ASSETS.items():
        path = source / source_rel
        if not path.is_file():
            raise SystemExit(f"Missing cached asset: {path}")
        algorithm, expected = sri.split("-", 1)
        actual = base64.b64encode(
            hashlib.new(algorithm, path.read_bytes()).digest()
        ).decode("ascii")
        if actual != expected:
            raise SystemExit(f"Integrity mismatch: {path}")
        target = args.cache_root / "R" / "pkgdown" / cache_rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)

    print(f"PKGDOWN_CACHE_SEEDED_OK: {len(ASSETS)} verified assets")


if __name__ == "__main__":
    main()
