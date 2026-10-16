"""Tests for checking a cleaned local pkgdown candidate site."""

import json
import tempfile
import unittest
from pathlib import Path

from verify_candidate_site import (
    EXPECTED_RETIRED_MANIFEST_SHA256,
    EXPECTED_SITE_BASE,
    audit_candidate_site,
)


class CandidateSiteTests(unittest.TestCase):
    def make_site(self, root: Path) -> tuple[Path, Path]:
        site = root / "site"
        site.mkdir()
        (site / "index.html").write_text(
            '<link href="extra.css" rel="stylesheet">'
            "<blockquote><p><strong>Warning:</strong> "
            "pigauto is experimental; use at your own risk.</p></blockquote>"
        )
        (site / "extra.css").write_text("")
        articles = site / "articles"
        articles.mkdir()
        (articles / "index.html").write_text("<a href='getting-started.html'>Start</a>")

        page_paths = ["index.html", "articles/index.html"]
        for index in range(60):
            path = f"page-{index}.html"
            (site / path).write_text("candidate page")
            page_paths.append(path)
        (site / "sitemap.xml").write_text(
            "<urlset>"
            + "".join(
                f"<url><loc>{EXPECTED_SITE_BASE}{path}</loc></url>"
                for path in page_paths
            )
            + "</urlset>"
        )
        search_paths = [
            f"{EXPECTED_SITE_BASE}{page_paths[index % len(page_paths)]}"
            for index in range(568)
        ]
        search_entries = [{"path": path} for path in search_paths]
        search_entries.extend({"title": f"unindexed-{index}"} for index in range(45))
        (site / "search.json").write_text(json.dumps(search_entries))

        manifest = root / "retired.tsv"
        rows = [
            f"/dev/retired-{index}.html\t404\thttps://example.org/dev/retired-{index}.html\tPage not found"
            for index in range(42)
        ]
        rows.extend(
            [
                "/dev/old.html\t404\thttps://example.org/dev/old.html\tPage not found",
                "/articles/simulation-study.html\t404\thttps://example.org/articles/simulation-study.html\tPage not found",
            ]
        )
        manifest.write_text("path\tstatus\tfinal_url\ttitle\n" + "\n".join(rows) + "\n")
        return site, manifest

    def replace_sitemap_first_location(self, site: Path, value: str) -> None:
        path = site / "sitemap.xml"
        original = path.read_text()
        start = original.index("<loc>") + len("<loc>")
        end = original.index("</loc>", start)
        path.write_text(original[:start] + value + original[end:])

    def test_accepts_clean_site_and_reports_full_inventory(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            (site / "extra.css").write_text("figure > p:empty { display: none; }")
            result = audit_candidate_site(site, manifest)
            self.assertEqual(result["retired_routes"], 44)
            self.assertEqual(result["sitemap_entries"], 62)
            self.assertEqual(result["search_entries"], 613)
            self.assertEqual(result["search_entries_with_paths"], 568)
            self.assertIn("CANDIDATE_SITE_OK", result["message"])

    def test_rejects_a_retired_route_reintroduced_into_output(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            (site / "articles" / "simulation-study.html").write_text("old page")
            with self.assertRaisesRegex(ValueError, "retired output files"):
                audit_candidate_site(site, manifest)

    def test_rejects_retired_sitemap_and_search_targets(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            self.replace_sitemap_first_location(site, f"{EXPECTED_SITE_BASE}dev/old.html")
            with self.assertRaisesRegex(ValueError, "retired sitemap targets"):
                audit_candidate_site(site, manifest)

            self.replace_sitemap_first_location(site, f"{EXPECTED_SITE_BASE}index.html")
            entries = json.loads((site / "search.json").read_text())
            entries[0]["path"] = f"{EXPECTED_SITE_BASE}dev/old.html"
            (site / "search.json").write_text(json.dumps(entries))
            with self.assertRaisesRegex(ValueError, "retired search targets"):
                audit_candidate_site(site, manifest)

    def test_rejects_stale_home_warning_marker_and_retired_article_link(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            (site / "index.html").write_text("&gt; [!WARNING] old rendering")
            with self.assertRaisesRegex(ValueError, r"\[!WARNING\]"):
                audit_candidate_site(site, manifest)

            (site / "index.html").write_text(
                '<link href="extra.css" rel="stylesheet">'
                "<blockquote><p><strong>Warning:</strong> "
                "pigauto is experimental; use at your own risk.</p></blockquote>"
            )
            (site / "articles/index.html").write_text(
                "<a href='simulation-study.html'>Old article</a>"
            )
            with self.assertRaisesRegex(ValueError, "retired article-index links"):
                audit_candidate_site(site, manifest)

    def test_rejects_warning_text_outside_a_callout(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            (site / "index.html").write_text(
                '<link href="extra.css" rel="stylesheet">'
                '<p>Warning: this appears in unrelated page text.</p>'
            )
            with self.assertRaisesRegex(ValueError, "rendered warning callout"):
                audit_candidate_site(site, manifest)

    def test_rejects_incomplete_warning_body(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            (site / "index.html").write_text(
                '<link href="extra.css" rel="stylesheet">'
                "<blockquote><p><strong>Warning:</strong> short.</p></blockquote>"
            )
            with self.assertRaisesRegex(ValueError, "warning body"):
                audit_candidate_site(site, manifest)

    def test_rejects_warning_hidden_by_attribute_or_known_class(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            for hidden_attribute in (
                "hidden",
                'class="d-none"',
                'class="d-md-none"',
                'class="visually-hidden"',
                'style="clip-path: inset(50%);"',
                'style="opacity:0.0;"',
                'style="width:0;height:0;overflow:hidden;"',
            ):
                (site / "index.html").write_text(
                    '<link href="extra.css" rel="stylesheet">'
                    f"<blockquote {hidden_attribute}><strong>Warning:</strong> "
                    "pigauto is experimental; use at your own risk.</blockquote>"
                )
                with self.assertRaisesRegex(ValueError, "rendered warning callout"):
                    audit_candidate_site(site, manifest)

    def test_rejects_attribute_selector_hiding_warning(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            (site / "extra.css").write_text("[data-hidden] { display: none; }")
            (site / "index.html").write_text(
                '<link href="extra.css" rel="stylesheet">'
                '<blockquote data-hidden="true"><strong>Warning:</strong> '
                "pigauto is experimental; use at your own risk.</blockquote>"
            )
            with self.assertRaisesRegex(ValueError, "rendered warning callout"):
                audit_candidate_site(site, manifest)

    def test_ignores_script_text_and_checks_inline_and_external_stylesheets(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            for html, css in (
                (
                    '<link href="extra.css" rel="stylesheet">'
                    "<blockquote><strong>Warning:</strong>"
                    "<script>pigauto is experimental; use at your own risk.</script>"
                    "</blockquote>",
                    "",
                ),
                (
                    "<style>blockquote { display: none !important; }</style>"
                    "<blockquote><strong>Warning:</strong> "
                    "pigauto is experimental; use at your own risk.</blockquote>",
                    "",
                ),
                (
                    '<link href="extra.css" rel="stylesheet">'
                    "<blockquote><strong>Warning:</strong> "
                    "pigauto is experimental; use at your own risk.</blockquote>",
                    "blockquote { visibility: hidden; }",
                ),
                (
                    '<link href="extra.css" rel="stylesheet">'
                    "<blockquote><strong>Warning:</strong> "
                    "pigauto is experimental; use at your own risk.</blockquote>",
                    "blockquote { opacity: 0.0; }",
                ),
                (
                    '<link href="extra.css" rel="stylesheet">'
                    "<blockquote><strong>Warning:</strong> "
                    "pigauto is experimental; use at your own risk.</blockquote>",
                    ":where(blockquote) { display: none; }",
                ),
                (
                    '<link href="extra.css" rel="stylesheet">'
                    "<blockquote><strong>Warning:</strong> "
                    "pigauto is experimental; use at your own risk.</blockquote>",
                    ":is(blockquote, pre) { display: none; }",
                ),
                (
                    '<link href="extra.css" rel="stylesheet">'
                    "<blockquote><strong>Warning:</strong> "
                    "pigauto is experimental; use at your own risk.</blockquote>",
                    '@import url("hide.css");',
                ),
            ):
                (site / "index.html").write_text(html)
                (site / "extra.css").write_text(css)
                (site / "hide.css").write_text("blockquote { display: none; }")
                with self.assertRaisesRegex(ValueError, "rendered warning callout"):
                    audit_candidate_site(site, manifest)

    def test_rejects_external_stylesheet_import(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            (site / "extra.css").write_text(
                '@import url("https://external.example/hide.css");'
            )
            with self.assertRaisesRegex(ValueError, "external stylesheet that cannot be checked"):
                audit_candidate_site(site, manifest)

    def test_rejects_empty_or_truncated_discovery_indexes(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            complete_sitemap = (site / "sitemap.xml").read_text()
            (site / "sitemap.xml").write_text("<urlset></urlset>")
            with self.assertRaisesRegex(ValueError, "sitemap contains no locations"):
                audit_candidate_site(site, manifest)

            (site / "sitemap.xml").write_text(
                f"<urlset><url><loc>{EXPECTED_SITE_BASE}index.html</loc></url></urlset>"
            )
            with self.assertRaisesRegex(ValueError, "expected 62 locations"):
                audit_candidate_site(site, manifest)

            (site / "sitemap.xml").write_text(complete_sitemap)
            (site / "search.json").write_text('{"unexpected": []}')
            with self.assertRaisesRegex(ValueError, "search index has no recognized entries"):
                audit_candidate_site(site, manifest)

            (site / "search.json").write_text(
                json.dumps([{"path": f"{EXPECTED_SITE_BASE}index.html"}])
            )
            with self.assertRaisesRegex(ValueError, "expected 613 entries"):
                audit_candidate_site(site, manifest)

    def test_rejects_external_sitemap_and_search_urls(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            self.replace_sitemap_first_location(site, "https://external.example/pigauto/index.html")
            with self.assertRaisesRegex(ValueError, "outside the candidate site origin"):
                audit_candidate_site(site, manifest)

            self.replace_sitemap_first_location(site, f"{EXPECTED_SITE_BASE}index.html")
            entries = json.loads((site / "search.json").read_text())
            entries[0]["path"] = "https://external.example/pigauto/index.html"
            (site / "search.json").write_text(json.dumps(entries))
            with self.assertRaisesRegex(ValueError, "outside the candidate site origin"):
                audit_candidate_site(site, manifest)

    def test_rejects_truncated_or_duplicate_retired_manifest(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            lines = manifest.read_text().splitlines()
            manifest.write_text("\n".join(lines[:3]) + "\n")
            with self.assertRaisesRegex(ValueError, "expected 44 unique retired routes"):
                audit_candidate_site(site, manifest)

            manifest.write_text("\n".join(lines[:2] + [lines[1]] + lines[3:]) + "\n")
            with self.assertRaisesRegex(ValueError, "duplicate retired routes"):
                audit_candidate_site(site, manifest)

    def test_decodes_encoded_retired_routes_for_sitemap_matching(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            self.replace_sitemap_first_location(site, f"{EXPECTED_SITE_BASE}dev/%6Fld.html")
            with self.assertRaisesRegex(ValueError, "retired sitemap targets"):
                audit_candidate_site(site, manifest)

    def test_rejects_traversal_and_double_encoded_paths(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            lines = manifest.read_text().splitlines()
            for unsafe in ("/../outside.html", "/%252e%252e/outside.html"):
                lines[1] = lines[1].replace("/dev/retired-0.html", unsafe)
                manifest.write_text("\n".join(lines) + "\n")
                with self.assertRaisesRegex(ValueError, "unsafe retired route"):
                    audit_candidate_site(site, manifest)
                lines[1] = lines[1].replace(unsafe, "/dev/retired-0.html")

    def test_rejects_substituted_route_manifest_with_same_count(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            lines = manifest.read_text().splitlines()
            lines[1] = lines[1].replace("/dev/retired-0.html", "/dev/substituted.html")
            manifest.write_text("\n".join(lines) + "\n")
            with self.assertRaisesRegex(ValueError, "manifest checksum"):
                audit_candidate_site(
                    site,
                    manifest,
                    expected_manifest_sha256=EXPECTED_RETIRED_MANIFEST_SHA256,
                )


if __name__ == "__main__":
    unittest.main()
