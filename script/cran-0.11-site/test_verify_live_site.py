"""Offline tests for the deployed-site verifier's retirement filter."""

import unittest

from verify_live_site import find_retired_sitemap_targets, is_expected_utility_route


class SitemapRetirementTests(unittest.TestCase):
    def test_all_named_retired_routes_are_filtered(self):
        locations = [
            "https://itchyshin.github.io/pigauto/AGENTS.html",
            "https://itchyshin.github.io/pigauto/CLAUDE.html",
            "https://itchyshin.github.io/pigauto/goodagents.html",
            "https://itchyshin.github.io/pigauto/VALIDATION_LEDGER.html",
            "https://itchyshin.github.io/pigauto/articles/simulation-study.html",
            "https://itchyshin.github.io/pigauto/dev/bench.html",
            "https://itchyshin.github.io/pigauto/articles/getting-started.html",
        ]
        self.assertEqual(
            find_retired_sitemap_targets(locations),
            locations[:-1],
        )

    def test_retained_page_is_not_filtered(self):
        self.assertEqual(
            find_retired_sitemap_targets([
                "https://itchyshin.github.io/pigauto/articles/getting-started.html"
            ]),
            [],
        )

    def test_custom_404_is_an_expected_utility_route(self):
        self.assertTrue(is_expected_utility_route(
            "https://itchyshin.github.io/pigauto/404.html",
            200,
            "Page not found • pigauto",
        ))
        self.assertFalse(is_expected_utility_route(
            "https://itchyshin.github.io/pigauto/articles/getting-started.html",
            200,
            "Getting started with pigauto • pigauto",
        ))
        self.assertFalse(is_expected_utility_route(
            "https://itchyshin.github.io/pigauto/404.html",
            404,
            "Page not found • pigauto",
        ))


if __name__ == "__main__":
    unittest.main()
