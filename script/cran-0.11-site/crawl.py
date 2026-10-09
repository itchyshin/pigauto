#!/usr/bin/env python3
"""Check local site links, assets, anchors, and retired public URLs."""

import json
import argparse
import re
import sys
from html import unescape
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


CSS_URL = re.compile(r"url\(\s*(?:(['\"])(.*?)\1|([^)]*?))\s*\)", re.I | re.S)
CSS_IMPORT = re.compile(r"@import\s+(['\"])(.*?)\1", re.I | re.S)
CSS_COMMENT = re.compile(r"/\*.*?\*/", re.S)


def srcset_references(value):
    """Return candidate URLs from an HTML srcset attribute."""
    references = []
    position = 0
    while position < len(value):
        while position < len(value) and (value[position].isspace() or value[position] == ","):
            position += 1
        start = position
        while position < len(value) and not value[position].isspace():
            position += 1
        raw_candidate = value[start:position]
        candidate = raw_candidate.rstrip(",")
        if candidate:
            references.append(candidate)
        if raw_candidate.endswith(","):
            continue
        while position < len(value) and value[position] != ",":
            position += 1
        if position < len(value):
            position += 1
    return references


def css_references(text):
    """Return local or remote URLs referenced by CSS url() and string @import."""
    text = CSS_COMMENT.sub("", text)
    matches = []
    for match in CSS_URL.finditer(text):
        url = (match.group(2) if match.group(1) else match.group(3)).strip().strip("'\"")
        if url:
            matches.append((match.start(), url))
    for match in CSS_IMPORT.finditer(text):
        matches.append((match.start(), match.group(2)))
    return [url for _, url in sorted(matches)]


class Page(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.ids = set()
        self.in_style = False

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if attrs.get("id"):
            self.ids.add(attrs["id"])
        if tag == "a" and attrs.get("name"):
            self.ids.add(attrs["name"])
        for key in ("href", "src"):
            if attrs.get(key):
                self.links.append(attrs[key])
        if attrs.get("srcset"):
            self.links.extend(srcset_references(attrs["srcset"]))
        if attrs.get("style"):
            self.links.extend(css_references(attrs["style"]))
        if tag == "style":
            self.in_style = True

    def handle_endtag(self, tag):
        if tag == "style":
            self.in_style = False

    def handle_data(self, data):
        if self.in_style:
            self.links.extend(css_references(data))


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


def main(site_root, base_url=None):
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
    base = urlsplit(base_url) if base_url else None
    base_path = base.path.rstrip("/") if base else ""

    errors = []
    checked = 0
    processed_css = set()
    for path, page in pages.items():
        pending = [(path, link) for link in page.links]
        while pending:
            source, link = pending.pop()
            url = urlsplit(link)
            same_origin = bool(
                base
                and url.netloc == base.netloc
                and (url.scheme or base.scheme) == base.scheme
            )
            if url.scheme or url.netloc:
                if not same_origin:
                    continue
                if base_path and url.path.rstrip("/") != base_path \
                        and not url.path.startswith(base_path + "/"):
                    continue
                relative_url_path = url.path[len(base_path):] if base_path else url.path
                if not relative_url_path.strip("/"):
                    target = root / "index.html"
                else:
                    target = (root / unquote(relative_url_path.lstrip("/"))).resolve()
            elif not url.path:
                target = source
            elif url.path.startswith("/"):
                target = (root / unquote(url.path.lstrip("/"))).resolve()
            else:
                target = (source.parent / unquote(url.path)).resolve()
            checked += 1
            if not target.is_relative_to(root):
                errors.append(f"{source.relative_to(root)}: outside site {link}")
                continue
            if target.is_dir():
                target /= "index.html"
            if not target.exists():
                errors.append(f"{source.relative_to(root)}: missing {link}")
                continue
            if (url.fragment and target.suffix == ".html" and target in pages
                    and unquote(url.fragment) not in pages[target].ids):
                errors.append(f"{source.relative_to(root)}: absent anchor {link}")
            if target.suffix.lower() == ".css" and target not in processed_css:
                processed_css.add(target)
                pending.extend(
                    (target, reference)
                    for reference in css_references(target.read_text(errors="replace"))
                )

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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("site_root", type=Path)
    parser.add_argument("--base-url", help="deployed site URL used to resolve same-origin absolute links")
    arguments = parser.parse_args()
    sys.exit(main(arguments.site_root, base_url=arguments.base_url))
