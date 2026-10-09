#!/usr/bin/env python3
"""Audit the cleaned local pkgdown output against the retired-route receipt."""

from __future__ import annotations

import argparse
import csv
import json
import xml.etree.ElementTree as ET
from html.parser import HTMLParser
from pathlib import Path, PurePosixPath
from urllib.parse import urljoin, urlsplit


class LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs) -> None:
        self.hrefs.extend(value for name, value in attrs if name.lower() == "href" and value)


def _route_path(value: str) -> str:
    path = urlsplit(value).path
    return "/" + str(PurePosixPath(path.lstrip("/")))


def _search_paths(value) -> list[str]:
    paths: list[str] = []
    if isinstance(value, list):
        for item in value:
            paths.extend(_search_paths(item))
    elif isinstance(value, dict):
        for key, item in value.items():
            if key == "path":
                if isinstance(item, str) and item:
                    paths.append(item)
                elif isinstance(item, list):
                    paths.extend(entry for entry in item if isinstance(entry, str) and entry)
            else:
                paths.extend(_search_paths(item))
    return paths


def _read_retired_paths(manifest: Path) -> list[str]:
    with manifest.open(encoding="utf-8", newline="") as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        if not rows.fieldnames or "path" not in rows.fieldnames:
            raise ValueError("retired-route manifest must have a path column")
        paths = [_route_path(row["path"]) for row in rows if row.get("path")]
    if not paths:
        raise ValueError("retired-route manifest contains no routes")
    return paths


def audit_candidate_site(site_dir: Path, manifest: Path) -> dict:
    site_dir = site_dir.resolve()
    required = [site_dir / name for name in ("index.html", "sitemap.xml", "search.json", "articles/index.html")]
    missing = [str(path.relative_to(site_dir)) for path in required if not path.is_file()]
    if missing:
        raise ValueError("candidate site is missing required files: " + ", ".join(missing))

    retired = _read_retired_paths(manifest)
    output_files = [route for route in retired if (site_dir / route.lstrip("/")).is_file()]
    if output_files:
        raise ValueError("retired output files: " + ", ".join(output_files))

    sitemap = ET.parse(site_dir / "sitemap.xml").getroot()
    sitemap_paths = [
        _route_path(element.text or "")
        for element in sitemap.iter()
        if element.tag.endswith("loc") and element.text
    ]
    retired_sitemap = [path for path in sitemap_paths if any(path.endswith(route) for route in retired)]
    if retired_sitemap:
        raise ValueError("retired sitemap targets: " + ", ".join(retired_sitemap))

    search = json.loads((site_dir / "search.json").read_text(encoding="utf-8"))
    search_entries = search if isinstance(search, list) else search.get("results", search.get("docs", []))
    search_paths = [_route_path(path) for path in _search_paths(search_entries)]
    retired_search = [path for path in search_paths if any(path.endswith(route) for route in retired)]
    if retired_search:
        raise ValueError("retired search targets: " + ", ".join(retired_search))

    article_html = (site_dir / "articles/index.html").read_text(encoding="utf-8")
    links = LinkParser()
    links.feed(article_html)
    article_targets = [
        _route_path(urljoin("https://candidate.invalid/articles/index.html", href))
        for href in links.hrefs
    ]
    retired_article_links = [
        path for path in article_targets if any(path.endswith(route) for route in retired)
    ]
    if retired_article_links:
        raise ValueError("retired article-index links: " + ", ".join(retired_article_links))

    homepage = (site_dir / "index.html").read_text(encoding="utf-8")
    if "[!WARNING]" in homepage:
        raise ValueError("homepage contains the literal [!WARNING] marker")
    if "Warning:" not in homepage:
        raise ValueError("homepage warning callout is missing the rendered Warning: label")

    return {
        "retired_routes": len(retired),
        "sitemap_entries": len(sitemap_paths),
        "search_entries": len(search_entries),
        "search_entries_with_paths": len(search_paths),
        "html_pages": len(list(site_dir.rglob("*.html"))),
        "message": "CANDIDATE_SITE_OK",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = audit_candidate_site(args.site_dir, args.manifest)
    except (OSError, ET.ParseError, json.JSONDecodeError, ValueError) as error:
        print(f"CANDIDATE_SITE_FAILED: {error}")
        return 1
    print(json.dumps(result, indent=2))
    print(result["message"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
