"""Verify exact preservation of retired pkgdown static assets."""

from hashlib import sha256
from pathlib import Path
import re

root = Path(__file__).resolve().parents[3]
manifest = root / "docs/dev-log/cran-0.11-audit/retirement-manifest.md"
rows = re.findall(
    r"^\| \x60(pkgdown/assets/dev/[^|]+)\x60 \| "
    r"\x60(dev/archive/cran-011-public-pages/[^|]+)\x60 \| "
    r"(\d+) \| \x60([0-9a-f]{64})\x60 \|$",
    manifest.read_text(),
    re.MULTILINE,
)
assert len(rows) == 39, f"expected 39 manifest rows, found {len(rows)}"
assert sum(old.endswith(".html") for old, *_ in rows) == 34
assert sum(old.endswith(".png") for old, *_ in rows) == 5
assert len({old for old, *_ in rows}) == len(rows)
assert len({new for _, new, *_ in rows}) == len(rows)
for old, new, size, digest in rows:
    assert not (root / old).exists(), old
    source = root / new
    assert source.is_file(), new
    data = source.read_bytes()
    assert len(data) == int(size), new
    assert sha256(data).hexdigest() == digest, new
assets = root / "pkgdown/assets/dev"
assert not assets.exists() or not any(assets.iterdir())
cfg = (root / "_pkgdown.yml").read_text()
assert "href: dev/" not in cfg
print("RETIREMENT_ARCHIVE_OK: 34 HTML and 5 PNG byte hashes match")
