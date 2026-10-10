"""Typing namespace identity shared by preparation and context comparison.

A frozen export identifies its own namespace; it does not invent a database
version. Nonempty conflicting database fingerprints are never compatible.
"""

import re
from typing import Any, Mapping

_MISSING = {"", "unknown", "none", "null", "na", "n/a", "?", "-"}


def _text(value: Any) -> str:
    return "" if value is None else str(value).strip()


def typing_scope(row: Mapping, kind: str) -> tuple[str, str] | None:
    scheme = _text(row.get(f"{kind}_scheme"))
    version = _text(row.get(f"{kind}_scheme_version"))
    if version.casefold() in _MISSING:
        digest = _text(row.get(f"{kind}_database_sha256"))
        version = (
            "sha256:" + digest.lower()
            if len(digest) == 64 and all(c in "0123456789abcdef" for c in digest.lower())
            else ""
        )
    if kind == "cglin" and version.casefold() in _MISSING:
        # A single frozen export defines its own internally consistent code
        # namespace. This is file identity, not a database version/fingerprint.
        for field in (
            "cglin_frozen_export_sha256",
            "cglin_export_sha256",
            "cglin_public_typing_export_sha256",
        ):
            digest = _text(row.get(field)).lower()
            if re.fullmatch(r"[0-9a-f]{64}", digest):
                version = "exportsha256:" + digest
                break
    if scheme.casefold() in _MISSING or version.casefold() in _MISSING:
        return None
    return scheme, version


def compatible_typing(left: Mapping, right: Mapping, kind: str) -> bool:
    if typing_scope(left, kind) != typing_scope(right, kind):
        return False
    hashes = [_text(row.get(f"{kind}_database_sha256")).lower() for row in (left, right)]
    return not (all(hashes) and hashes[0] != hashes[1])
