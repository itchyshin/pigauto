import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crawl


class PageReferenceTests(unittest.TestCase):
    def test_collects_srcset_and_inline_css_references(self):
        page = crawl.Page()
        page.feed(
            '<img srcset="small.png 1x, large.png 2x" '
            'style="background-image: url(\'inline.png\')">'
            '<style>.hero { background: url("hero.png") }</style>'
        )

        self.assertEqual(
            page.links,
            ["small.png", "large.png", "inline.png", "hero.png"],
        )

    def test_extracts_css_urls_and_string_imports(self):
        self.assertEqual(
            crawl.css_references(
                '@import "theme.css"; .hero { background: url(\'hero.png\') }'
            ),
            ["theme.css", "hero.png"],
        )

    def test_crawl_follows_stylesheets_and_checks_their_assets(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "index.html").write_text(
                '<link rel="stylesheet" href="assets/site.css">'
            )
            (root / "assets").mkdir()
            (root / "assets/site.css").write_text(
                '@import "theme.css"; .hero { background: url("../missing.png") }'
            )
            (root / "assets/theme.css").write_text(
                '.hero { background: url("../theme.png") }'
            )
            (root / "theme.png").write_bytes(b"image")
            (root / "search.json").write_text("[]")
            (root / "sitemap.xml").write_text("<urlset></urlset>")

            output = StringIO()
            with redirect_stdout(output):
                result = crawl.main(root)

        self.assertEqual(result, 1)
        self.assertIn("missing ../missing.png", output.getvalue())


if __name__ == "__main__":
    unittest.main()
