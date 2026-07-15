"""Small, format-independent helpers for ranked taxonomy paths."""

from __future__ import annotations

from typing import Iterable, Sequence


def taxonomy_path(value: Sequence[str] | str) -> tuple[str, ...]:
    """Normalize a semicolon string or rank sequence into a non-empty path."""

    parts = value.split(";") if isinstance(value, str) else value
    taxonomy = tuple(str(part).strip() for part in parts)
    while taxonomy and not taxonomy[-1]:
        taxonomy = taxonomy[:-1]
    if not taxonomy or not taxonomy[0]:
        raise ValueError("Taxonomy must contain a non-empty domain")
    return taxonomy


def lowest_common_ancestor(
    taxonomies: Iterable[Sequence[str]],
) -> tuple[str, ...]:
    """Return the shared, non-empty prefix of ranked taxonomy paths."""

    paths = [tuple(taxonomy) for taxonomy in taxonomies]
    if not paths:
        return ()
    common: list[str] = []
    for values in zip(*paths):
        if not values[0] or len(set(values)) != 1:
            break
        common.append(values[0])
    return tuple(common)


def common_value(values: Iterable[str], *, conflict: str = "") -> str:
    """Return one non-empty value, an explicit conflict sentinel, or empty."""

    unique = sorted({value for value in values if value})
    if len(unique) == 1:
        return unique[0]
    return conflict if unique else ""
