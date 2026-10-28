#!/usr/bin/env python3
"""Fail CI when a direct Python dependency has missing or changed license metadata."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, metadata


EXPECTED_LICENSE_MARKERS = {
    "numpy": ("bsd-3-clause", "bsd 3-clause", "bsd license"),
    "pandas": ("bsd-3-clause", "bsd 3-clause", "bsd license"),
    "scipy": ("bsd-3-clause", "bsd 3-clause", "bsd license"),
    "biopython": (
        "licenseref-biopython-license-agreement",
        "biopython license agreement",
        "bsd",
    ),
    "pyrodigal": (
        "gpl-3.0-or-later",
        "gpl-3.0",
        "gplv3",
        "general public license v3",
    ),
    "pyyaml": ("mit",),
}


def declared_license(distribution: str) -> str:
    """Combine standardized and legacy Core Metadata license declarations."""
    package_metadata = metadata(distribution)
    values = [
        package_metadata.get("License-Expression", ""),
        package_metadata.get("License", ""),
        *package_metadata.get_all("Classifier", []),
    ]
    return "\n".join(value for value in values if value).casefold()


def main() -> int:
    failures = []
    for distribution, markers in EXPECTED_LICENSE_MARKERS.items():
        try:
            declaration = declared_license(distribution)
        except PackageNotFoundError:
            failures.append(f"{distribution}: package is not installed")
            continue
        if not declaration:
            failures.append(f"{distribution}: no license metadata")
        elif not any(marker in declaration for marker in markers):
            failures.append(f"{distribution}: license metadata changed; review upstream")
        else:
            print(f"ok: {distribution}")

    if failures:
        print("License audit failed:")
        for failure in failures:
            print(f"- {failure}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
