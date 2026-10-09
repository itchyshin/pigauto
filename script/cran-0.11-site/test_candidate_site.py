"""Tests for checking a cleaned local pkgdown candidate site."""

import tempfile
import unittest
from pathlib import Path

from verify_candidate_site import audit_candidate_site


class CandidateSiteTests(unittest.TestCase):
    def make_site(self, root: Path) -> tuple[Path, Path]:
        site = root / "site"
        site.mkdir()
        (site / "index.html").write_text("<blockquote><strong>Warning:</strong> note</blockquote>")
        (site / "sitemap.xml").write_text(
            "<urlset><url><loc>https://example.org/pigauto/index.html</loc></url></urlset>"
        )
        (site / "search.json").write_text('[{"path":"https://example.org/pigauto/index.html"}]')
        articles = site / "articles"
        articles.mkdir()
        (articles / "index.html").write_text("<a href='getting-started.html'>Start</a>")
        manifest = root / "retired.tsv"
        manifest.write_text(
            "path\tstatus\tfinal_url\ttitle\n"
            "/dev/old.html\t404\thttps://example.org/dev/old.html\tPage not found\n"
            "/articles/simulation-study.html\t404\thttps://example.org/articles/simulation-study.html\tPage not found\n"
        )
        return site, manifest

    def test_accepts_clean_site_and_reports_retired_route_count(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            result = audit_candidate_site(site, manifest)
            self.assertEqual(result["retired_routes"], 2)
            self.assertEqual(result["sitemap_entries"], 1)
            self.assertIn("CANDIDATE_SITE_OK", result["message"])

    def test_rejects_a_retired_route_reintroduced_into_output(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            (site / "articles" / "simulation-study.html").write_text("old page")
            with self.assertRaisesRegex(ValueError, "retired output files"):
                audit_candidate_site(site, manifest)

    def test_rejects_retired_discovery_and_stale_home_warning_marker(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            (site / "sitemap.xml").write_text(
                "<urlset><url><loc>https://example.org/pigauto/dev/old.html</loc></url></urlset>"
            )
            with self.assertRaisesRegex(ValueError, "retired sitemap targets"):
                audit_candidate_site(site, manifest)

            (site / "sitemap.xml").write_text(
                "<urlset><url><loc>https://example.org/pigauto/index.html</loc></url></urlset>"
            )
            (site / "index.html").write_text("&gt; [!WARNING] old rendering")
            with self.assertRaisesRegex(ValueError, r"\[!WARNING\]"):
                audit_candidate_site(site, manifest)

    def test_rejects_retired_search_target(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            (site / "search.json").write_text(
                '[{"path":"https://example.org/pigauto/dev/old.html"}]'
            )
            with self.assertRaisesRegex(ValueError, "retired search targets"):
                audit_candidate_site(site, manifest)

    def test_rejects_retired_article_index_link(self):
        with tempfile.TemporaryDirectory() as temp:
            site, manifest = self.make_site(Path(temp))
            (site / "articles/index.html").write_text(
                "<a href='simulation-study.html'>Old article</a>"
            )
            with self.assertRaisesRegex(ValueError, "retired article-index links"):
                audit_candidate_site(site, manifest)


if __name__ == "__main__":
    unittest.main()
