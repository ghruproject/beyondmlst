"""Prepared data shared by analysis routes, independent of provider retrieval.

Alleles are categories, not numbers. Binary matrices use vocabulary indices;
index zero represents null and is never an allele or an allele match.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date
from collections.abc import Mapping, Sequence
import calendar
import json
import re

import numpy as np

from chronoclade.metadata import SAFE_IDENTIFIER

SCHEMA_NAME = "chronoclade.prepared-dataset"
SCHEMA_VERSION = 1
SAMPLE_COLUMNS = (
    "sample_id",
    "label",
    "role",
    "species",
    "mlst_scheme",
    "mlst_st",
    "collection_date",
    "date_start",
    "date_end",
    "date_precision",
    "country",
    "region",
    "nuts2",
    "host",
    "isolation_source",
    "assembly_reference",
)
EVIDENCE_TABLES = ("lineages", "crosswalk", "provenance", "conflicts", "exclusions", "retrieval")
_SECRET_KEYS = {
    "api_key",
    "apikey",
    "token",
    "access_token",
    "password",
    "authorization",
    "credentials",
    "secret",
    "client_secret",
}


class DatasetError(ValueError):
    """The prepared bundle or its evidence is invalid or incomplete."""


def check_json(value, *, location="dataset"):
    """Require portable JSON evidence and keep credential material out of bundles."""
    if isinstance(value, Mapping):
        for key, item in value.items():
            if not isinstance(key, str):
                raise DatasetError(f"{location}: evidence keys must be strings")
            if key.casefold().replace("-", "_") in _SECRET_KEYS:
                raise DatasetError(f"{location}: credentials must be referenced by configuration")
            check_json(item, location=f"{location}.{key}")
    elif isinstance(value, (tuple, list)):
        for item in value:
            check_json(item, location=location)
    try:
        json.dumps(value, allow_nan=False)
    except (TypeError, ValueError) as error:
        raise DatasetError(f"{location}: evidence must be finite JSON values") from error


@dataclass(frozen=True)
class LocusCatalogue:
    """An explicitly supplied ordered locus universe and its frozen provenance.

    ``complete=False`` records an observed export honestly; callers must never
    infer completeness from a union of calls. ``source`` identifies the catalogue
    file, database, or provider evidence; it does not manufacture a database version.
    """

    scheme_id: str
    scheme_version: str
    loci: tuple[str, ...]
    source: str
    database_version: str | None = None
    complete: bool = True
    database_sha256: str | None = None

    def validate(self):
        for name in ("scheme_id", "scheme_version", "source"):
            if not isinstance(getattr(self, name), str) or not getattr(self, name).strip():
                raise DatasetError(f"Locus catalogue requires {name}")
        if type(self.complete) is not bool:
            raise DatasetError("Locus catalogue completeness must be explicit boolean")
        if self.database_version is not None and not isinstance(self.database_version, str):
            raise DatasetError("Locus catalogue database_version must be text or null")
        if self.database_sha256 is not None and (
            not isinstance(self.database_sha256, str)
            or not re.fullmatch(r"[0-9a-f]{64}", self.database_sha256)
        ):
            raise DatasetError("Locus catalogue database_sha256 requires a lowercase SHA256 digest")
        if not self.loci or any(not isinstance(x, str) or not x.strip() for x in self.loci):
            raise DatasetError("Locus catalogue requires nonempty named loci")
        if len(set(self.loci)) != len(self.loci):
            raise DatasetError("Locus catalogue repeats a locus")


def _allele(value) -> str | None:
    if value is None or type(value) is int and value == 0:
        return None
    if isinstance(value, str) and value.casefold() in {
        "",
        "0",
        "-",
        "?",
        "none",
        "null",
        "unknown",
        "na",
        "n/a",
    }:
        return None
    if isinstance(value, str) and value.strip() and value == value.strip():
        return value
    if type(value) is int and value > 0:
        return str(value)
    raise DatasetError("Allele calls must be categorical text, positive integers, or null/0")


@dataclass(frozen=True)
class AlleleMatrix:
    catalogue: LocusCatalogue
    sample_ids: tuple[str, ...]
    categories: tuple[str, ...]
    codes: np.ndarray = field(repr=False, compare=False)

    @classmethod
    def from_profiles(
        cls, catalogue: LocusCatalogue, sample_ids: Sequence[str], profiles: Mapping[str, Mapping]
    ) -> AlleleMatrix:
        """Encode calls against every supplied locus; absent cells remain null.

        Omitted samples receive an all-null row. Unexpected samples or loci are
        errors, including a called locus outside a supposedly complete catalogue.
        Novel allele hashes are preserved exactly as categorical strings.
        """
        catalogue.validate()
        ids = tuple(sample_ids)
        if set(profiles) - set(ids):
            raise DatasetError("Profiles contain samples outside the matrix sample IDs")
        positions = {locus: index for index, locus in enumerate(catalogue.loci)}
        categories, indices = [], {}
        codes = np.zeros((len(ids), len(catalogue.loci)), dtype=np.uint32)
        for index, ident in enumerate(ids):
            calls = profiles.get(ident, {})
            if not isinstance(calls, Mapping) or set(calls) - set(positions):
                raise DatasetError(f"Profile {ident!r} has unsupported loci or shape")
            for locus, value in calls.items():
                allele = _allele(value)
                if allele is not None:
                    if allele not in indices:
                        categories.append(allele)
                        indices[allele] = len(categories)
                    codes[index, positions[locus]] = indices[allele]
        matrix = cls(catalogue, ids, tuple(categories), codes)
        matrix.validate()
        codes.setflags(write=False)
        return matrix

    def validate(self):
        self.catalogue.validate()
        if not self.sample_ids or len(set(self.sample_ids)) != len(self.sample_ids):
            raise DatasetError("Allele matrix requires unique sample IDs")
        if any(not isinstance(x, str) or not SAFE_IDENTIFIER.fullmatch(x) for x in self.sample_ids):
            raise DatasetError("Allele matrix sample IDs must be safe identifiers")
        if any(_allele(x) is None or not isinstance(x, str) for x in self.categories):
            raise DatasetError("Allele vocabulary cannot contain missing categories")
        if len(set(self.categories)) != len(self.categories):
            raise DatasetError("Allele vocabulary repeats a category")
        if self.codes.dtype != np.dtype("uint32") or self.codes.shape != (
            len(self.sample_ids),
            len(self.catalogue.loci),
        ):
            raise DatasetError("Allele matrix shape or categorical code dtype is invalid")
        if self.codes.size and int(self.codes.max()) > len(self.categories):
            raise DatasetError("Allele matrix refers to an unknown category")

    def profile(self, sample_id: str) -> dict[str, str | None]:
        """Return a full-locus profile with actual nulls, never missing allele matches."""
        index = self.sample_ids.index(sample_id)
        return {
            locus: self.categories[int(code) - 1] if code else None
            for locus, code in zip(self.catalogue.loci, self.codes[index], strict=True)
        }


def _validate_sample(row):
    if not isinstance(row, dict):
        raise DatasetError("Sample rows must be dictionaries")
    ident = row.get("sample_id")
    if not isinstance(ident, str) or not SAFE_IDENTIFIER.fullmatch(ident):
        raise DatasetError("Samples require safe stable sample IDs")
    if row.get("role") not in {"input", "context"}:
        raise DatasetError(f"Sample {ident!r} requires input/context role")
    for key in SAMPLE_COLUMNS:
        value = row.get(key)
        if value is not None and not isinstance(value, str):
            raise DatasetError(f"Sample {ident!r}: {key} must be text or null")
    if not row.get("label"):
        raise DatasetError(f"Sample {ident!r} requires a readable label")
    precision = row.get("date_precision")
    if precision not in {None, "missing", "invalid", "day", "month", "year", "interval"}:
        raise DatasetError(f"Sample {ident!r}: unsupported date precision")
    bounds = (row.get("date_start"), row.get("date_end"))
    if precision in {"day", "month", "year", "interval"}:
        try:
            start, end = (date.fromisoformat(value) for value in bounds)
        except (ValueError, TypeError) as error:
            raise DatasetError(
                f"Sample {ident!r}: date interval requires valid ISO bounds"
            ) from error
        if start > end or precision == "day" and start != end:
            raise DatasetError(f"Sample {ident!r}: inconsistent date interval")
        if precision == "month" and (
            start.day != 1
            or (start.year, start.month) != (end.year, end.month)
            or end.day != calendar.monthrange(end.year, end.month)[1]
        ):
            raise DatasetError(f"Sample {ident!r}: month precision requires the full month")
        if precision == "year" and (
            start != date(start.year, 1, 1) or end != date(start.year, 12, 31)
        ):
            raise DatasetError(f"Sample {ident!r}: year precision requires the full year")
    elif any(bounds):
        raise DatasetError(f"Sample {ident!r}: missing/invalid date cannot have interval bounds")


@dataclass(frozen=True)
class PreparedDataset:
    samples: tuple[dict, ...]
    profiles: tuple[AlleleMatrix, ...] = ()
    lineages: tuple[dict, ...] = ()
    crosswalk: tuple[dict, ...] = ()
    provenance: tuple[dict, ...] = ()
    conflicts: tuple[dict, ...] = ()
    exclusions: tuple[dict, ...] = ()
    retrieval: tuple[dict, ...] = ()
    parameters: dict = field(default_factory=dict)
    dataset_id: str | None = None

    @property
    def sample_ids(self) -> tuple[str, ...]:
        return tuple(row["sample_id"] for row in self.samples)

    def validate(self):
        if not self.samples:
            raise DatasetError("Prepared dataset contains no samples")
        for row in self.samples:
            _validate_sample(row)
        if len(set(self.sample_ids)) != len(self.samples):
            raise DatasetError("Prepared dataset repeats a sample ID")
        ids, scopes = set(self.sample_ids), set()
        for matrix in self.profiles:
            matrix.validate()
            if set(matrix.sample_ids) - ids:
                raise DatasetError("Allele matrix refers to an unknown sample")
            scope = (
                matrix.catalogue.scheme_id,
                matrix.catalogue.scheme_version,
                matrix.catalogue.database_version,
                matrix.catalogue.database_sha256,
            )
            if scope in scopes:
                raise DatasetError("Prepared dataset repeats a profile scheme/database scope")
            scopes.add(scope)
        for name in EVIDENCE_TABLES:
            for row in getattr(self, name):
                if not isinstance(row, dict) or row.get("sample_id") not in ids:
                    raise DatasetError(f"{name} evidence requires an existing sample_id")
                if name == "lineages":
                    for key in ("assignment_method", "resolution_status"):
                        if not isinstance(row.get(key), str) or not row[key]:
                            raise DatasetError(f"Lineage evidence requires {key}")
                    for key in ("scheme_id", "scheme_version", "database_version"):
                        if key not in row or (
                            row[key] is not None and not isinstance(row[key], str)
                        ):
                            raise DatasetError(f"Lineage evidence requires nullable text {key}")
                    if "database_version" not in row or "assignment" not in row:
                        raise DatasetError(
                            "Lineage evidence requires database_version and full assignment"
                        )
                    if row["resolution_status"] == "resolved" and not row["assignment"]:
                        raise DatasetError("Resolved lineage evidence requires an assignment")
                check_json(row, location=name)
        check_json(self.samples, location="samples")
        if not isinstance(self.parameters, dict):
            raise DatasetError("Dataset parameters must be a dictionary")
        check_json(self.parameters, location="parameters")
