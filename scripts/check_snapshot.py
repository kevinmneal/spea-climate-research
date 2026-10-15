#!/usr/bin/env python3
"""Check released file hashes and portable Markdown links; standard library only."""
import hashlib
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]


def main():
    entries = (ROOT / "CHECKSUMS.sha256").read_text().splitlines()
    for entry in entries:
        expected, name = entry.split("  ", 1)
        path = (ROOT / name).resolve()
        if not path.is_relative_to(ROOT) or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"Changed snapshot file: {name}")
        if path.suffix != ".md":
            continue
        for link in re.findall(r"\]\(([^)]+)\)", path.read_text()):
            parts = urlsplit(link)
            if parts.scheme or parts.netloc:
                continue
            target = (path.parent / unquote(parts.path)).resolve() if parts.path else path
            if not target.is_relative_to(ROOT) or not target.exists():
                raise ValueError(f"Missing public link in {name}: {link}")
            if parts.fragment and target.is_file() and target.suffix == ".md":
                headings = re.findall(r"^#+ (.+)$", target.read_text(), re.MULTILINE)
                slugs = {re.sub(r"[^\w\- ]", "", h.lower()).replace(" ", "-") for h in headings}
                if parts.fragment not in slugs:
                    raise ValueError(f"Missing heading in {name}: {link}")
    print(f"Verified {len(entries)} released files and their local Markdown links.")


if __name__ == "__main__":
    main()
