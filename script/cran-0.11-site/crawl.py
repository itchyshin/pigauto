#!/usr/bin/env python3
"""Check local site links, assets, anchors, and retired public URLs."""

import json
import re
import sys
from html import unescape
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


class Page(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.ids = set()

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if attrs.get("id"):
            self.ids.add(attrs["id"])
        if tag == "a" and attrs.get("name"):
            self.ids.add(attrs["name"])
        for key in ("href", "src"):
            if attrs.get(key):
                self.links.append(attrs[key])


def indexed_urls(path):
    if path.name == "search.json":
        for entry in json.loads(path.read_text()):
            url = entry.get("path", "")
            if isinstance(url, str):
                yield url
            elif url != []:
                raise ValueError(f"Malformed search path: {url!r}")
    else:
        for url in re.findall(r"<loc>(.*?)</loc>", path.read_text(), re.S):
            yield unescape(url)


def main(site_root):
    root = site_root.resolve()
    if not root.is_dir():
        raise ValueError(f"Site directory missing: {root}")

    pages = {}
    for path in root.rglob("*.html"):
        page = Page()
        page.feed(path.read_text(errors="replace"))
        pages[path] = page
    if not pages:
        raise ValueError("Site has no HTML pages")

    errors = []
    checked = 0
    for path, page in pages.items():
        for link in page.links:
            url = urlsplit(link)
            if url.scheme or url.netloc:
                continue
            if not url.path:
                target = path
            elif url.path.startswith("/"):
                target = (root / unquote(url.path.lstrip("/"))).resolve()
            else:
                target = (path.parent / unquote(url.path)).resolve()
            checked += 1
            if not target.is_relative_to(root):
                errors.append(f"{path.relative_to(root)}: outside site {link}")
                continue
            if target.is_dir():
                target /= "index.html"
            if not target.exists():
                errors.append(f"{path.relative_to(root)}: missing {link}")
                continue
            if (url.fragment and target.suffix == ".html" and target in pages
                    and unquote(url.fragment) not in pages[target].ids):
                errors.append(f"{path.relative_to(root)}: absent anchor {link}")

    internal = ("AGENTS", "CLAUDE", "goodagents", "VALIDATION_LEDGER")
    for stem in internal:
        for extension in ("html", "md"):
            if (root / f"{stem}.{extension}").exists():
                errors.append(f"Internal page served: {stem}.{extension}")

    repo = Path(__file__).resolve().parents[2]
    retired = list((repo / "dev/archive/cran-011-public-pages").glob("*.html"))
    for path in retired:
        if (root / "dev" / path.name).exists():
            errors.append(f"Retired page served: dev/{path.name}")

    for filename in ("search.json", "sitemap.xml"):
        path = root / filename
        if not path.exists():
            errors.append(f"Missing {filename}")
            continue
        for url in indexed_urls(path):
            indexed_path = urlsplit(url).path
            for retired_path in retired:
                if indexed_path.endswith("/dev/" + retired_path.name):
                    errors.append(f"{filename}: retired URL {retired_path.name}")
            for stem in internal:
                if indexed_path.endswith("/" + stem + ".html"):
                    errors.append(f"{filename}: internal {stem}")
            for route in (
                "articles/simulation-study.html",
                "articles/simulation-study.md",
                "articles/articles/simulation-study.html",
                "articles/articles/simulation-study.md",
            ):
                if indexed_path == route or indexed_path.endswith("/" + route):
                    errors.append(f"{filename}: retired URL {route}")

    for route in (
        "articles/simulation-study.html",
        "articles/simulation-study.md",
        "articles/articles/simulation-study.html",
        "articles/articles/simulation-study.md",
    ):
        if (root / route).exists():
            errors.append(f"Retired historical page served: {route}")

    report = {
        "html_pages": len(pages),
        "local_references_checked": checked,
        "retired_pages": len(retired),
        "errors": sorted(set(errors)),
    }
    print(json.dumps(report, indent=2))
    if errors:
        return 1
    print("SITE_CRAWL_OK")
    return 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1])))
