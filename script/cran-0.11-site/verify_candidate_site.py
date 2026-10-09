#!/usr/bin/env python3
"""Audit the cleaned local pkgdown output against the retired-route receipt."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from html.parser import HTMLParser
from pathlib import Path, PurePosixPath
from urllib.parse import unquote, urljoin, urlsplit


EXPECTED_RETIRED_ROUTE_COUNT = 44
EXPECTED_RETIRED_MANIFEST_SHA256 = "5414f6cb69e9d140c38d0a6308f8e0ddbe8952834eea01b6b861a78eecb807e1"
EXPECTED_WARNING_BODY = "pigauto is experimental; use at your own risk."
EXPECTED_SITE_BASE = "https://itchyshin.github.io/pigauto/"
EXPECTED_SITEMAP_ENTRIES = 62
EXPECTED_SEARCH_ENTRIES = 613
EXPECTED_SEARCH_PATHS = 568

_HIDDEN_CLASSES = {
    "collapse",
    "d-none",
    "hidden",
    "invisible",
    "opacity-0",
    "sr-only",
    "visually-hidden",
}
_HIDING_DECLARATION = re.compile(
    r"(?:display\s*:\s*none|visibility\s*:\s*(?:hidden|collapse)|"
    r"opacity\s*:\s*0(?:\.0+)?(?:e[+-]?0+)?(?:\s*!important)?(?:\s*;|\s*$)|"
    r"clip(?:-path)?\s*:|text-indent\s*:\s*-|"
    r"(?:^|[;{]\s*)(?:width|height)\s*:\s*0(?:px|rem|em|%)?\b)",
    re.IGNORECASE,
)
_CSS_RULE = re.compile(r"([^{}]+)\{([^{}]*)\}")
_CSS_IMPORT = re.compile(
    r"@import\s+(?:url\(\s*([^)]*?)\s*\)|(['\"])(.*?)\2)[^;]*;",
    re.IGNORECASE,
)


def _has_hidden_class(value: str) -> bool:
    classes = set(value.lower().split())
    return bool(
        classes.intersection(_HIDDEN_CLASSES)
        or any(re.fullmatch(r"d-(?:sm|md|lg|xl|xxl)-none", item) for item in classes)
        or ("collapse" in classes and "show" not in classes)
    )


class LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs) -> None:
        self.hrefs.extend(value for name, value in attrs if name.lower() == "href" and value)


class WarningCalloutParser(HTMLParser):
    """Find structurally valid blockquote warnings and track hidden content."""

    _VOID_TAGS = {
        "area", "base", "br", "col", "embed", "hr", "img", "input",
        "link", "meta", "param", "source", "track", "wbr",
    }
    _INERT_TAGS = {"script", "style", "template", "noscript"}

    def __init__(self) -> None:
        super().__init__()
        self.stack: list[tuple[str, bool, dict]] = []
        self.current: dict | None = None
        self.callouts: list[dict] = []
        self.in_strong = False
        self.inert_depth = 0
        self.in_style_depth = 0
        self.inline_styles: list[str] = []

    def handle_starttag(self, tag: str, attrs) -> None:
        attrs = {name.lower(): value or "" for name, value in attrs}
        style = attrs.get("style", "").lower().replace(" ", "")
        locally_hidden = (
            "hidden" in attrs
            or attrs.get("aria-hidden", "").lower() == "true"
            or bool(_HIDING_DECLARATION.search(style))
            or _has_hidden_class(attrs.get("class", ""))
        )
        inherited_hidden = any(hidden for _, hidden, _ in self.stack)
        hidden = locally_hidden or inherited_hidden
        node = {
            "tag": tag.lower(),
            "id": attrs.get("id", ""),
            "classes": set(attrs.get("class", "").lower().split()),
            "attrs": attrs,
        }

        if tag.lower() == "blockquote" and self.current is None:
            self.current = {
                "hidden": hidden,
                "text": [],
                "strong": [],
                "node_paths": [[entry[2] for entry in self.stack] + [node]],
            }
        elif self.current is not None:
            self.current["node_paths"].append(
                [entry[2] for entry in self.stack] + [node]
            )
            if hidden:
                self.current["hidden"] = True

        if tag.lower() in self._INERT_TAGS:
            self.inert_depth += 1
            if tag.lower() == "style":
                self.in_style_depth += 1
                self.inline_styles.append("")

        if tag.lower() == "strong" and self.current is not None:
            self.in_strong = True
            self.current["strong"].append([])
        if tag.lower() not in self._VOID_TAGS:
            self.stack.append((tag.lower(), hidden, node))

    def handle_endtag(self, tag: str) -> None:
        tag = tag.lower()
        if tag in self._INERT_TAGS and self.inert_depth:
            self.inert_depth -= 1
        if tag == "style" and self.in_style_depth:
            self.in_style_depth -= 1
        if tag == "strong" and self.current is not None:
            self.in_strong = False
        if tag == "blockquote" and self.current is not None:
            self.callouts.append(self.current)
            self.current = None
        for index in range(len(self.stack) - 1, -1, -1):
            if self.stack[index][0] == tag:
                del self.stack[index:]
                break

    def handle_data(self, data: str) -> None:
        if self.inert_depth:
            if self.in_style_depth and self.inline_styles:
                self.inline_styles[-1] += data
            return
        if self.current is None:
            return
        self.current["text"].append(data)
        if self.in_strong and self.current["strong"]:
            self.current["strong"][-1].append(data)


class StylesheetParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.stylesheet_hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs) -> None:
        if tag.lower() != "link":
            return
        values = {name.lower(): value or "" for name, value in attrs}
        if "stylesheet" in values.get("rel", "").lower().split() and values.get("href"):
            self.stylesheet_hrefs.append(values["href"])


def _selector_may_hide_path(selector: str, nodes: list[dict]) -> bool:
    """Match the common tag/class/id subset used by package stylesheets."""
    if ":root" in selector:
        selector = selector.replace(":root", "html")
    selector = re.sub(r"::?[\w-]+(?:\([^)]*\))?", "", selector)
    tokens = list(re.finditer(r"[^\s>+~]+|[>+~]", selector.strip()))
    if not tokens:
        return False

    compounds: list[str] = []
    combinators: list[str] = []
    for index, token in enumerate(tokens):
        if token.group() in {">", "+", "~"}:
            if not combinators and len(compounds) == 0:
                continue
            combinators.append(token.group())
            continue
        if compounds:
            previous = tokens[index - 1]
            gap = selector[previous.end():token.start()]
            if previous.group() not in {">", "+", "~"} and gap.strip() == "":
                combinators.append(" ")
        compounds.append(token.group())
    def matches(compound: str, node: dict) -> bool:
        tag_match = re.match(r"^[a-zA-Z][\w-]*", compound)
        tag = tag_match.group(0).lower() if tag_match else None
        element_id = re.search(r"#([\w-]+)", compound)
        classes = re.findall(r"\.([\w-]+)", compound)
        attribute_selectors = re.findall(r"\[([^\]]+)\]", compound)
        if tag and node["tag"] != tag:
            return False
        if element_id and node["id"] != element_id.group(1):
            return False
        if classes and not set(classes).issubset(node["classes"]):
            return False
        for selector in attribute_selectors:
            match = re.fullmatch(
                r"\s*([\w:-]+)(?:\s*(=|~=|\|=|\^=|\$=|\*=)\s*['\"]?([^'\"]*)['\"]?)?\s*",
                selector,
            )
            if not match:
                continue
            name, operator, expected = match.groups()
            if name not in node["attrs"]:
                return False
            actual = node["attrs"][name]
            if operator == "=" and actual != expected:
                return False
            if operator == "~=" and expected not in actual.split():
                return False
            if operator == "|=" and actual != expected and not actual.startswith(expected + "-"):
                return False
            if operator == "^=" and not actual.startswith(expected):
                return False
            if operator == "$=" and not actual.endswith(expected):
                return False
            if operator == "*=" and expected not in actual:
                return False
        return bool(tag or element_id or classes or attribute_selectors or "*" in compound)

    def _compound_matches_any(compound: str, candidates: list[dict]) -> bool:
        return any(matches(compound, node) for node in candidates)

    if len(combinators) != len(compounds) - 1 or any(
        relation in {"+", "~"} for relation in combinators
    ):
        # Sibling relations need sibling-tree state; fail closed when a simple
        # selector atom can target any node in the warning subtree.
        for compound in compounds:
            if _compound_matches_any(compound, nodes):
                return True
        return False

    for end in range(len(nodes)):
        if not matches(compounds[-1], nodes[end]):
            continue
        current = end
        matched = True
        for index in range(len(compounds) - 2, -1, -1):
            relation = combinators[index]
            if relation == ">":
                current -= 1
                if current < 0 or not matches(compounds[index], nodes[current]):
                    matched = False
                    break
            else:
                ancestor = current - 1
                while ancestor >= 0 and not matches(compounds[index], nodes[ancestor]):
                    ancestor -= 1
                if ancestor < 0:
                    matched = False
                    break
                current = ancestor
        if matched:
            return True
    return False


def _split_selector_list(selectors: str) -> list[str]:
    """Split a selector list without splitting commas inside functions/brackets."""
    parts: list[str] = []
    start = 0
    depth = 0
    quote: str | None = None
    escaped = False
    for index, character in enumerate(selectors):
        if escaped:
            escaped = False
            continue
        if character == "\\":
            escaped = True
            continue
        if quote:
            if character == quote:
                quote = None
            continue
        if character in {"'", '"'}:
            quote = character
        elif character in "([":
            depth += 1
        elif character in ")]":
            depth = max(0, depth - 1)
        elif character == "," and depth == 0:
            parts.append(selectors[start:index].strip())
            start = index + 1
    parts.append(selectors[start:].strip())
    return [part for part in parts if part]


def _expand_selector_functions(selector: str) -> list[str]:
    """Expand :is() and :where() argument lists into simple selector alternatives."""
    match = re.search(r":(?:is|where)\(", selector, flags=re.IGNORECASE)
    if not match:
        return [selector]
    open_paren = selector.find("(", match.start())
    depth = 1
    quote: str | None = None
    escaped = False
    close_paren = open_paren + 1
    while close_paren < len(selector) and depth:
        character = selector[close_paren]
        if escaped:
            escaped = False
        elif character == "\\":
            escaped = True
        elif quote:
            if character == quote:
                quote = None
        elif character in {"'", '"'}:
            quote = character
        elif character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
        close_paren += 1
    if depth:
        return [selector]
    arguments = selector[open_paren + 1:close_paren - 1]
    prefix = selector[:match.start()]
    suffix = selector[close_paren:]
    expanded: list[str] = []
    for argument in _split_selector_list(arguments):
        for alternative in _expand_selector_functions(prefix + argument + suffix):
            expanded.append(alternative)
    return expanded


def _css_may_hide_callout(css: str, node_paths: list[list[dict]]) -> bool:
    css = re.sub(r"/\*.*?\*/", "", css, flags=re.DOTALL)
    for selectors, declarations in _CSS_RULE.findall(css):
        if not _HIDING_DECLARATION.search(declarations):
            continue
        if any(
            _selector_may_hide_path(selector, nodes)
            for selector in _split_selector_list(selectors)
            for selector in _expand_selector_functions(selector)
            for nodes in node_paths
        ):
            return True
    return False


def _load_stylesheet(
    site_dir: Path,
    href: str,
    base_url: str = EXPECTED_SITE_BASE,
    visited: set[str] | None = None,
) -> str:
    if visited is None:
        visited = set()
    base = urlsplit(EXPECTED_SITE_BASE)
    resolved = urlsplit(urljoin(base_url, href))
    if (resolved.scheme, resolved.netloc) != (base.scheme, base.netloc):
        raise ValueError(f"candidate uses an external stylesheet that cannot be checked: {href}")
    if not resolved.path.startswith(base.path):
        raise ValueError(f"stylesheet is outside the candidate site prefix: {href}")
    relative = unquote(resolved.path[len(base.path):])
    for _ in range(8):
        decoded = unquote(relative)
        if decoded == relative:
            break
        relative = decoded
    else:
        raise ValueError(f"stylesheet path has excessive percent encoding: {href}")
    if any(part in {".", ".."} for part in relative.split("/")) or "\\" in relative:
        raise ValueError(f"unsafe candidate stylesheet path: {href}")
    css_path = (site_dir / relative).resolve()
    try:
        css_path.relative_to(site_dir)
    except ValueError as error:
        raise ValueError(f"stylesheet resolves outside candidate site: {href}") from error
    if not css_path.is_file():
        raise ValueError(f"candidate stylesheet is missing: {href}")
    resolved_url = resolved.geturl()
    if resolved_url in visited:
        raise ValueError(f"candidate stylesheet import cycle detected: {href}")
    visited.add(resolved_url)
    css = css_path.read_text(encoding="utf-8")

    def replace_import(match: re.Match) -> str:
        imported = (match.group(1) or match.group(3) or "").strip().strip("\"'")
        if not imported:
            raise ValueError(f"candidate has an empty stylesheet import: {href}")
        imported_css = _load_stylesheet(site_dir, imported, resolved_url, visited)
        return "\n" + imported_css + "\n"

    return _CSS_IMPORT.sub(replace_import, css)


def _site_stylesheets(site_dir: Path, homepage: str) -> list[str]:
    parser = StylesheetParser()
    parser.feed(homepage)
    styles = [
        _load_stylesheet(site_dir, href)
        for href in parser.stylesheet_hrefs
    ]
    return styles


def _warning_callout_matches(homepage: str) -> bool:
    parser = WarningCalloutParser()
    parser.feed(homepage)
    expected = EXPECTED_WARNING_BODY.casefold()
    for callout in parser.callouts:
        if callout["hidden"]:
            continue
        strong_labels = ["".join(parts).strip() for parts in callout["strong"]]
        body = " ".join("".join(callout["text"]).split()).casefold()
        if "Warning:" in strong_labels and expected in body:
            return True
    return False


def _route_path(value: str) -> str:
    path = unquote(urlsplit(value).path)
    for _ in range(8):
        decoded = unquote(path)
        if decoded == path:
            break
        path = decoded
    else:
        raise ValueError(f"retired route has excessive percent encoding: {value}")
    if "\\" in path or any(part in {".", ".."} for part in path.split("/")):
        raise ValueError(f"unsafe retired route: {value}")
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
        raw_paths = [row["path"] for row in rows if row.get("path")]
    paths = [_route_path(path) for path in raw_paths]
    if not paths:
        raise ValueError("retired-route manifest contains no routes")
    if len(paths) != len(set(paths)):
        raise ValueError("retired-route manifest contains duplicate retired routes")
    if len(paths) != EXPECTED_RETIRED_ROUTE_COUNT:
        raise ValueError(
            "retired-route manifest expected 44 unique retired routes; "
            f"found {len(paths)}"
        )
    return paths


def _site_target(site_dir: Path, route: str) -> Path:
    target = (site_dir / route.lstrip("/")).resolve()
    try:
        target.relative_to(site_dir)
    except ValueError as error:
        raise ValueError(f"retired route resolves outside candidate site: {route}") from error
    return target


def _site_url_path(value: str) -> str:
    expected = urlsplit(EXPECTED_SITE_BASE)
    actual = urlsplit(value)
    if (actual.scheme, actual.netloc) != (expected.scheme, expected.netloc):
        raise ValueError(f"discovery URL is outside the candidate site origin: {value}")
    prefix = expected.path.rstrip("/") + "/"
    route = _route_path(value)
    if not route.startswith(prefix):
        raise ValueError(f"discovery URL is outside the candidate site prefix: {value}")
    return route


def audit_candidate_site(
    site_dir: Path,
    manifest: Path,
    expected_manifest_sha256: str | None = None,
) -> dict:
    site_dir = site_dir.resolve()
    required = [site_dir / name for name in ("index.html", "sitemap.xml", "search.json", "articles/index.html")]
    missing = [str(path.relative_to(site_dir)) for path in required if not path.is_file()]
    if missing:
        raise ValueError("candidate site is missing required files: " + ", ".join(missing))

    if expected_manifest_sha256:
        actual_hash = hashlib.sha256(manifest.read_bytes()).hexdigest()
        if actual_hash != expected_manifest_sha256:
            raise ValueError("retired-route manifest checksum does not match the approved receipt")
    retired = _read_retired_paths(manifest)
    output_files = [route for route in retired if _site_target(site_dir, route).is_file()]
    if output_files:
        raise ValueError("retired output files: " + ", ".join(output_files))

    sitemap = ET.parse(site_dir / "sitemap.xml").getroot()
    sitemap_paths = [
        _site_url_path(element.text or "")
        for element in sitemap.iter()
        if element.tag.endswith("loc") and element.text
    ]
    if not sitemap_paths:
        raise ValueError("candidate sitemap contains no locations")
    retired_sitemap = [path for path in sitemap_paths if any(path.endswith(route) for route in retired)]
    if retired_sitemap:
        raise ValueError("retired sitemap targets: " + ", ".join(retired_sitemap))
    if len(sitemap_paths) != EXPECTED_SITEMAP_ENTRIES:
        raise ValueError(
            f"candidate sitemap expected {EXPECTED_SITEMAP_ENTRIES} locations; found {len(sitemap_paths)}"
        )
    sitemap_files = {path[len("/pigauto/"):] for path in sitemap_paths}
    html_files = {str(path.relative_to(site_dir)) for path in site_dir.rglob("*.html")}
    if sitemap_files != html_files:
        missing = sorted(html_files - sitemap_files)
        absent = sorted(sitemap_files - html_files)
        raise ValueError(
            f"candidate sitemap and HTML inventory differ; missing={missing}; absent={absent}"
        )
    search = json.loads((site_dir / "search.json").read_text(encoding="utf-8"))
    if isinstance(search, list):
        search_entries = search
    elif isinstance(search, dict):
        search_entries = search.get("results", search.get("docs"))
    else:
        search_entries = None
    if not isinstance(search_entries, list) or not search_entries:
        raise ValueError("search index has no recognized entries")
    if len(search_entries) != EXPECTED_SEARCH_ENTRIES:
        raise ValueError(
            f"candidate search index expected {EXPECTED_SEARCH_ENTRIES} entries; found {len(search_entries)}"
        )
    search_paths = [_site_url_path(path) for path in _search_paths(search_entries)]
    if not search_paths:
        raise ValueError("search index entries contain no recognized paths")
    if len(search_paths) != EXPECTED_SEARCH_PATHS:
        raise ValueError(
            f"candidate search index expected {EXPECTED_SEARCH_PATHS} paths; found {len(search_paths)}"
        )
    retired_search = [path for path in search_paths if any(path.endswith(route) for route in retired)]
    if retired_search:
        raise ValueError("retired search targets: " + ", ".join(retired_search))
    if any(path not in sitemap_paths for path in search_paths):
        raise ValueError("search index contains a path absent from the candidate sitemap")

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
    parser = WarningCalloutParser()
    parser.feed(homepage)
    stylesheets = _site_stylesheets(site_dir, homepage)
    hidden_by_stylesheet = any(
        _css_may_hide_callout(css, callout["node_paths"])
        for css in [*stylesheets, *parser.inline_styles]
        for callout in parser.callouts
        if not callout["hidden"]
    )
    if "[!WARNING]" in homepage:
        raise ValueError("homepage contains the literal [!WARNING] marker")
    if not _warning_callout_matches(homepage) or hidden_by_stylesheet:
        raise ValueError(
            "homepage is missing a rendered warning callout with the expected warning body"
        )

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
        result = audit_candidate_site(
            args.site_dir,
            args.manifest,
            expected_manifest_sha256=EXPECTED_RETIRED_MANIFEST_SHA256,
        )
    except (OSError, ET.ParseError, json.JSONDecodeError, ValueError) as error:
        print(f"CANDIDATE_SITE_FAILED: {error}")
        return 1
    print(json.dumps(result, indent=2))
    print(result["message"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
