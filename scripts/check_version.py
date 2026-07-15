#!/usr/bin/env python3

import re
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def require_match(pattern: str, text: str, source: str) -> str:
    match = re.search(pattern, text, flags=re.MULTILINE | re.DOTALL)
    if not match:
        raise ValueError(f"Version not found in {source}")
    return match.group(1)


def read_versions() -> dict[str, str]:
    pixi_text = (REPO / "pixi.toml").read_text()
    nextflow_text = (REPO / "nextflow.config").read_text()
    changelog_text = (REPO / "CHANGELOG.md").read_text()
    return {
        "pixi.toml": require_match(
            r'^version\s*=\s*"([0-9]+\.[0-9]+\.[0-9]+)"',
            pixi_text,
            "pixi.toml",
        ),
        "nextflow.config": require_match(
            r"manifest\s*\{.*?version\s*=\s*'([0-9]+\.[0-9]+\.[0-9]+)'",
            nextflow_text,
            "nextflow.config",
        ),
        "CHANGELOG.md": require_match(
            r'^## \[([0-9]+\.[0-9]+\.[0-9]+)\]',
            changelog_text,
            "CHANGELOG.md",
        ),
    }


def main() -> None:
    versions = read_versions()
    unique_versions = set(versions.values())
    if len(unique_versions) != 1:
        details = ", ".join(f"{source}={version}" for source, version in versions.items())
        raise SystemExit(f"Version mismatch: {details}")
    print(unique_versions.pop())


if __name__ == "__main__":
    main()

