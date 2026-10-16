#!/usr/bin/env python3
"""Check deployed retired routes and pkgdown discovery files using the Python stdlib."""

from __future__ import annotations

import argparse
import json
import re
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlsplit
from html.parser import HTMLParser
from pathlib import Path


SITEMAP_RETIRED_TOKENS = (
    "/dev/",
    "simulation-study",
    "VALIDATION_LEDGER",
    "/AGENTS.html",
    "/CLAUDE.html",
    "/goodagents.html",
)


def find_retired_sitemap_targets(locations: list[str]) -> list[str]:
    return [location for location in locations if any(token in location for token in SITEMAP_RETIRED_TOKENS)]


class TitleParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.in_title = False
        self.parts: list[str] = []

    def handle_starttag(self, tag: str, attrs) -> None:
        if tag.lower() == "title":
            self.in_title = True

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() == "title":
            self.in_title = False

    def handle_data(self, data: str) -> None:
        if self.in_title:
            self.parts.append(data)


def fetch_url(url: str) -> tuple[int, bytes, str]:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "pigauto-CRAN-0.11-release-audit"},
    )
    try:
        with urllib.request.urlopen(request, timeout=20) as response:
            return response.status, response.read(), response.geturl()
    except urllib.error.HTTPError as error:
        return error.code, error.read(), error.geturl()


def fetch(base: str, path: str) -> tuple[int, bytes, str]:
    return fetch_url(base.rstrip("/") + path)


def check_retained_url(url: str) -> tuple[str, int, str, str]:
    status, body, final_url = fetch_url(url)
    title_parser = TitleParser()
    title_parser.feed(body.decode("utf-8", errors="replace"))
    title = " ".join(" ".join(title_parser.parts).split())
    return url, status, final_url, title


def is_expected_utility_route(url: str, status: int, title: str) -> bool:
    return urlsplit(url).path.endswith("/404.html") and status == 200 and "Page not found" in title


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", default="https://itchyshin.github.io/pigauto")
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[2]
    manifest = args.manifest or repo / "docs/dev-log/cran-0.11-audit/retirement-manifest.md"
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    source = manifest.read_text(encoding="utf-8")
    routes = {
        "/dev/" + name
        for name in re.findall(r"pkgdown/assets/dev/([A-Za-z0-9_.-]+\.(?:html|png))", source)
    }
    routes.update(
        {
            "/VALIDATION_LEDGER.html",
            "/AGENTS.html",
            "/CLAUDE.html",
            "/goodagents.html",
            "/articles/simulation-study.html",
        }
    )

    rows = ["path\tstatus\tfinal_url\ttitle"]
    failures = []
    for path in sorted(routes):
        status, body, final_url = fetch(args.base, path)
        title_parser = TitleParser()
        title_parser.feed(body.decode("utf-8", errors="replace"))
        title = " ".join(" ".join(title_parser.parts).split())
        rows.append(f"{path}\t{status}\t{final_url}\t{title}")
        if status != 404 or (path.lower().endswith(".html") and "Page not found" not in title):
            failures.append({"path": path, "status": status, "title": title})
    (output_dir / "retired-route-statuses.tsv").write_text(
        "\n".join(rows) + "\n", encoding="utf-8"
    )

    sitemap_status, sitemap_body, _ = fetch(args.base, "/sitemap.xml")
    locations = [
        element.text or ""
        for element in ET.fromstring(sitemap_body).iter()
        if element.tag.endswith("loc")
    ]
    sitemap_retired = find_retired_sitemap_targets(locations)

    base_parts = urlsplit(args.base)
    site_prefix = f"{base_parts.scheme}://{base_parts.netloc}{base_parts.path.rstrip('/')}/"
    if any(not location.startswith(site_prefix) for location in locations):
        raise ValueError("Sitemap contains a URL outside the requested site prefix")
    with ThreadPoolExecutor(max_workers=4) as executor:
        retained_results = list(executor.map(check_retained_url, locations))
    retained_rows = ["url\tstatus\tfinal_url\ttitle"]
    retained_failures = []
    expected_utility_routes = []
    for url, status, final_url, title in retained_results:
        retained_rows.append(f"{url}\t{status}\t{final_url}\t{title}")
        if is_expected_utility_route(url, status, title):
            expected_utility_routes.append({"url": url, "status": status, "title": title})
        elif status != 200 or "Page not found" in title:
            retained_failures.append({"url": url, "status": status, "title": title})
    (output_dir / "retained-route-statuses.tsv").write_text(
        "\n".join(retained_rows) + "\n", encoding="utf-8"
    )

    search_status, search_body, _ = fetch(args.base, "/search.json")
    search = json.loads(search_body)
    entries = search if isinstance(search, list) else search.get("results", search.get("docs", []))
    retired_tokens = ("/dev/", "simulation-study", "VALIDATION_LEDGER", "/AGENTS.html", "/CLAUDE.html", "/goodagents.html")
    search_retired = [
        entry.get("path", "")
        for entry in entries
        if any(token in entry.get("path", "") for token in retired_tokens)
    ]
    history_mentions = sorted(
        {
            entry.get("path", "")
            for entry in entries
            if "/dev/" in entry.get("text", "")
        }
    )
    summary = {
        "routes_checked": len(routes),
        "route_failures": failures,
        "retained_routes_checked": len(retained_results),
        "retained_route_failures": retained_failures,
        "expected_utility_routes": expected_utility_routes,
        "sitemap_status": sitemap_status,
        "sitemap_entries": len(locations),
        "sitemap_retired_targets": sitemap_retired,
        "search_status": search_status,
        "search_entries": len(entries),
        "search_retired_targets": search_retired,
        "search_entries_with_historical_dev_text": history_mentions,
    }
    (output_dir / "discovery-check.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    if failures or retained_failures or sitemap_status != 200 or search_status != 200 or sitemap_retired or search_retired:
        print("LIVE_SITE_RETIREMENT_CHECK_FAILED")
        return 1
    print("LIVE_SITE_RETIREMENT_CHECK_OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
