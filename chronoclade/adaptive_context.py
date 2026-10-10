"""Select version-compatible Klebsiella context within ST and full cgLIN prefixes."""

from __future__ import annotations

from collections import Counter
from copy import deepcopy
import hashlib
import json
import re
from typing import Iterable, Mapping

from chronoclade.context_refinement import _lineage, _text
from chronoclade.typing_scopes import compatible_typing as _compatible, typing_scope as _scope


def _st(row: Mapping) -> str:
    value = _text(row.get("mlst_st"))
    if not value:
        match = re.fullmatch(r"ST([0-9]+)", _text(row.get("lineage")), re.IGNORECASE)
        value = match.group(1) if match else ""
    value = re.sub(r"^ST", "", value, flags=re.IGNORECASE)
    return value if value.isdigit() else ""


def _species(row: Mapping) -> str:
    return " ".join(_text(row.get("species")).replace("_", " ").casefold().split())


def _eligible(row: Mapping, depth: int = 5) -> bool:
    species = _species(row)
    return (
        (not species or "klebsiella" in species)
        and bool(_st(row))
        and _lineage(row, "cglin", depth) is not None
    )


def _cg_name(rows: list[dict]) -> str:
    names = set()
    for row in rows:
        export = row.get("cglin_export_row") or {}
        export = export if isinstance(export, Mapping) else {}
        value = _text(export.get("Clonal Group") or row.get("clonal_group"))
        if re.fullmatch(r"(?:CG)?[0-9]+", value):
            names.add("CG" + value.removeprefix("CG"))
    return next(iter(names)) if len(names) == 1 else ""


def adaptive_cglin_context(
    inputs: Iterable[Mapping], context: Iterable[Mapping], *, min_context: int = 20
) -> tuple[list[dict], list[dict], dict]:
    """Return copied inputs, selected public rows and an auditable selection summary.

    Caller supplies deduplicated public samples with inputs excluded. Supported
    inputs partition by ST and version-scoped four-component CG prefix. For each
    represented five-component subgroup, select the narrowest of levels 7, 6, 5
    having the requested minimum for *each* represented prefix. An unresolved
    prefix cannot meet that criterion. Sparse groups fall back to level 5 without
    expanding to other CGs. Unsupported/unassigned inputs retain their existing
    lineage pathway, without selecting their public rows. The caller handles
    that pathway. This performs no downloads or sampling.
    """
    if isinstance(min_context, bool) or not isinstance(min_context, int) or min_context < 1:
        raise ValueError("min_context must be a positive integer")
    queries, public = deepcopy(list(inputs)), deepcopy(list(context))
    audit = {
        "method": "adaptive-cglin-context",
        "min_context": min_context,
        "levels": [5, 6, 7],
        "datasets": [],
        "ineligible_inputs": [],
    }
    groups: dict[tuple, list[dict]] = {}
    for row in queries:
        if not _eligible(row):
            audit["ineligible_inputs"].append(
                {"sample_id": row.get("sample_id"), "reason": "unsupported_or_unassigned_cglin"}
            )
            continue
        key = (_st(row), _lineage(row, "cglin", 4), _species(row))
        groups.setdefault(key, []).append(row)
    names: dict[tuple, str] = {}
    for key, rows in groups.items():
        prefix = key[1][-1]
        cg = _cg_name(rows) or "CGprefix_" + "_".join(map(str, prefix))
        names[key] = "ST" + re.sub(r"[^A-Za-z0-9_-]", "_", key[0]) + "_" + cg
    name_counts = Counter(names.values())
    for key in names:
        if name_counts[names[key]] > 1:
            # Distinct scheme namespaces must never merge under the same CG label.
            digest = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()[:10]
            names[key] += "_" + digest
    selected: dict[tuple[str, str], dict] = {}
    for key, rows in sorted(groups.items(), key=lambda item: repr(item[0])):
        name = names[key]
        for row in rows:
            row["analysis_dataset"] = name
        candidates = [
            row
            for row in public
            if _eligible(row, 4)
            and _st(row) == key[0]
            and _species(row) == key[2]
            and _lineage(row, "cglin", 4) == key[1]
            and all(_compatible(q, row, "cglin") for q in rows)
        ]
        dataset = {
            "analysis_dataset": name,
            "mlst_st": key[0],
            "species": key[2],
            "clonal_group": _cg_name(rows),
            "cg_prefix": list(key[1][-1]),
            "cglin_scope": list(_scope(rows[0], "cglin")),
            "input_count": len(rows),
            "public_cg_pool_count": len(candidates),
            "subgroups": [],
        }
        subgroups: dict[tuple, list[dict]] = {}
        for row in rows:
            subgroups.setdefault(_lineage(row, "cglin", 5), []).append(row)
        dataset_selected = set()
        for prefix, subgroup in sorted(subgroups.items(), key=lambda item: repr(item[0])):
            matches, minima = {}, {}
            for level in (5, 6, 7):
                prefixes = {_lineage(q, "cglin", level) for q in subgroup}
                resolved = prefixes - {None}
                matches[level] = [p for p in candidates if _lineage(p, "cglin", level) in resolved]
                counts = [
                    sum(_lineage(p, "cglin", level) == group for p in candidates)
                    for group in resolved
                ]
                minima[level] = 0 if None in prefixes or not counts else min(counts)
            chosen = next((level for level in (7, 6, 5) if minima[level] >= min_context), 5)
            chosen_rows = matches[chosen]
            dataset["subgroups"].append(
                {
                    "input_prefix": list(prefix[-1]),
                    "input_count": len(subgroup),
                    "query_ids": [q.get("sample_id") for q in subgroup],
                    "context_counts": {str(level): len(matches[level]) for level in (5, 6, 7)},
                    "minimum_prefix_context_counts": {
                        str(level): minima[level] for level in (5, 6, 7)
                    },
                    "selected_level": chosen,
                    "selected_context_count": len(chosen_rows),
                    "limited_context": minima[chosen] < min_context,
                }
            )
            for row in chosen_rows:
                copied = deepcopy(row)
                copied["analysis_dataset"] = name
                copied["profile_pool_selection_reason"] = f"adaptive_cglin_level_{chosen}"
                ident = _text(row.get("sample_id") or row.get("source_genome_id"))
                selected[(name, ident)] = copied
                dataset_selected.add(ident)
        dataset["selected_context_count"] = len(dataset_selected)
        audit["datasets"].append(dataset)
    return queries, list(selected.values()), audit


def annotate_cglin_datasets(
    context: Iterable[Mapping], annotated_inputs: Iterable[Mapping]
) -> list[dict]:
    """Copy the full public catalogue with CG dataset labels where unambiguous.

    This also labels a public assignment resolved only to CG (level 4), without
    claiming it qualifies for the narrower selected context pool. Unmatched or
    incompatible rows remain unchanged.
    """
    queries = list(annotated_inputs)
    result = deepcopy(list(context))
    for row in result:
        names = {
            _text(q.get("analysis_dataset"))
            for q in queries
            if _eligible(q)
            and _eligible(row, 4)
            and _st(q) == _st(row)
            and _species(q) == _species(row)
            and _lineage(q, "cglin", 4) == _lineage(row, "cglin", 4)
            and _compatible(q, row, "cglin")
        } - {""}
        if len(names) == 1:
            row["analysis_dataset"] = next(iter(names))
    return result
