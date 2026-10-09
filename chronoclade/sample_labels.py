"""Accession labels for presentation, separate from analytical sample identity."""

from collections import Counter
from collections.abc import Iterable, Mapping
import re
import json
from pathlib import Path


def sample_label(row: Mapping) -> str:
    ident = str(row.get("sample_id") or row.get("source_genome_id") or "")
    if row.get("origin") in {"local", "query", "focal"} and not ident.startswith("PW_"):
        return ident
    values = set()
    for key in ("run_accessions", "biosample_accessions", "assembly_accessions", "aliases"):
        raw = row.get(key) or []
        values.update([raw] if isinstance(raw, str) else raw)
    for key in ("run_accession", "runAccession", "biosample", "sample_accession",
                "sampleAccession", "assembly_accession", "assemblyAccession", "accession"):
        if row.get(key):
            values.add(str(row[key]))
    for pattern in (r"[SED]RR\d+", r"SAM[NED][A-Z]?\d+", r"GC[AF]_\d+(?:\.\d+)?"):
        candidates = sorted(str(value) for value in values if re.fullmatch(pattern, str(value)))
        if candidates:
            return candidates[0]
    return ident


def sample_labels(records: Iterable[Mapping]) -> dict[str, str]:
    labels = {str(row["sample_id"]): sample_label(row) for row in records if row.get("sample_id")}
    counts = Counter(labels.values())
    return {ident: f"{label} ({ident})" if counts[label] > 1 else label
            for ident, label in labels.items()}


def label_analysis(value, labels: Mapping[str, str], field: str = ""):
    """Copy reader-facing fields without relabelling paths or scientific identifiers."""
    fields = {"sample_id", "sample_ids", "query_id", "context_id", "root", "exploratory_root",
              "tested_roots", "sample_a", "sample_b", "focal_sample", "context_sample", "nearest_focal"}
    if isinstance(value, dict):
        return {key: label_analysis(item, labels, key) for key, item in value.items()}
    if isinstance(value, list):
        return [label_analysis(item, labels, field) for item in value]
    if isinstance(value, str) and field in fields:
        return labels.get(value, value)
    return value


def read_sample_labels(directory: Path) -> dict[str, str]:
    try:
        value = json.loads((directory / "sample_labels.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return {str(key): str(label) for key, label in value.items()} if isinstance(value, dict) else {}
