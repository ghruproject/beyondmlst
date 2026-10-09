"""Import frozen public typing against exact Pathogenwatch source identifiers.

Export hashes describe the imported file, never the typing database. Conflicting
assignments across duplicate public records remain unassigned rather than being
resolved by choosing the first record or borrowing a neighbour's code.
"""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Iterable, Mapping

from chronoclade.cglin import normalise_assignment


class PublicTypingError(ValueError):
    """A public typing export is unsafe or malformed."""


KINDS = ("cgmlst", "hiercc", "cglin")


def _text(value: Any) -> str:
    return "" if value is None else str(value).strip()


def _scope(row: Mapping, kind: str) -> None:
    for field in ("scheme", "scheme_version"):
        if not isinstance(row.get(f"{kind}_{field}"), str) or not row[f"{kind}_{field}"].strip():
            raise PublicTypingError(f"{kind} assignment requires explicit {field}")
    fingerprint = row.get(f"{kind}_database_sha256")
    if fingerprint is not None and not re.fullmatch(r"[0-9a-fA-F]{64}", str(fingerprint)):
        raise PublicTypingError(f"{kind} database fingerprint must be a SHA256")


def _normalise(row: Mapping, digest: str, index: int) -> dict:
    ident = row.get("source_genome_id")
    if not isinstance(ident, str) or not ident or ident != ident.strip():
        raise PublicTypingError("Public typing requires an exact nonempty source_genome_id")
    result = {
        "source_genome_id": ident,
        "public_typing_export_sha256": digest,
        "public_typing_export_row": index,
    }
    for kind in KINDS:
        if not any(str(key).startswith(kind + "_") for key in row):
            continue
        _scope(row, kind)
        fields = {
            key: deepcopy(value) for key, value in row.items() if str(key).startswith(kind + "_")
        }
        fields[f"{kind}_public_typing_export_sha256"] = digest
        fields[f"{kind}_source"] = fields.get(f"{kind}_source") or "imported_public_typing"
        if kind == "cgmlst":
            loci, profile = row.get("cgmlst_loci"), row.get("cgmlst_profile")
            if (
                not isinstance(loci, list)
                or not loci
                or any(not isinstance(x, str) or not x for x in loci)
                or len(set(loci)) != len(loci)
            ):
                raise PublicTypingError("cgMLST requires a complete unique locus universe")
            if not isinstance(profile, dict) or not set(profile).issubset(loci):
                raise PublicTypingError("cgMLST profile must map declared locus IDs to alleles")
            calls = {}
            for locus, allele in profile.items():
                value = _text(allele)
                if value in {"", "0", "?", "-"}:
                    continue
                if re.fullmatch(r"[0-9a-fA-F]{40}", value):
                    calls[locus] = value.lower()
                elif value.isdigit() and int(value) > 0:
                    calls[locus] = str(int(value))
                else:
                    raise PublicTypingError(
                        "cgMLST alleles must be positive integers or SHA1 hashes"
                    )
            novel = row.get("cgmlst_novel_alleles") or {}
            if (
                not isinstance(novel, dict)
                or not set(novel).issubset(loci)
                or set(novel).intersection(calls)
            ):
                raise PublicTypingError("Novel cgMLST calls must use distinct declared loci")
            if any(not re.fullmatch(r"[0-9a-fA-F]{40}", str(value)) for value in novel.values()):
                raise PublicTypingError("Novel cgMLST calls require exact SHA1 hashes")
            fields["cgmlst_novel_alleles"] = {
                key: str(value).lower() for key, value in novel.items()
            }
            fields.update(
                cgmlst_profile=calls,
                cgmlst_called_fraction=(len(calls) + len(novel)) / len(loci),
                cgmlst_status=fields.get("cgmlst_status")
                or ("resolved" if calls or novel else "unassigned"),
            )
        elif kind == "hiercc":
            codes = row.get("hiercc_codes")
            if not isinstance(codes, dict):
                raise PublicTypingError("HierCC requires a level-to-cluster mapping")
            normalised = {}
            for level, cluster in codes.items():
                match = re.fullmatch(r"(?:HC|d)(\d+)", str(level))
                value = _text(cluster)
                if not match or not value.isdigit():
                    raise PublicTypingError(
                        "HierCC levels must be HC/d integers and clusters nonnegative integers"
                    )
                key = "HC" + str(int(match[1]))
                if key in normalised:
                    raise PublicTypingError("Duplicate normalised HierCC level")
                normalised[key] = str(int(value))
            fields.update(
                hiercc_codes=normalised,
                hiercc_status=fields.get("hiercc_status")
                or ("resolved" if normalised else "unassigned"),
            )
        else:
            assignment = normalise_assignment(row)
            assignment.pop("source_genome_id", None)
            fields.update(assignment)
            if row.get("cglin_status") == "conflict":
                fields = _clear(fields, kind)
        result.update(fields)
    if len(result) == 3:
        raise PublicTypingError("Public typing row contains no cgMLST, cgLIN or HierCC assignment")
    return result


def load_public_typing(path: str | Path) -> list[dict]:
    """Load a JSON record list or {assignments: [...]} frozen export envelope."""
    try:
        data = Path(path).read_bytes()
        payload = json.loads(data)
    except (OSError, ValueError) as exc:
        raise PublicTypingError("Cannot read public typing JSON export") from exc
    rows = payload.get("assignments") if isinstance(payload, dict) else payload
    if not isinstance(rows, list) or not rows:
        raise PublicTypingError("Public typing export requires a nonempty assignment list")
    digest = hashlib.sha256(data).hexdigest()
    if any(not isinstance(row, dict) for row in rows):
        raise PublicTypingError("Public typing assignments must be JSON objects")
    return [_normalise(row, digest, index) for index, row in enumerate(rows, 1)]


def _has(row: Mapping, kind: str) -> bool:
    if kind == "cgmlst":
        return bool(row.get("cgmlst_profile") or row.get("cgmlst_novel_alleles"))
    if kind == "hiercc":
        return bool(row.get("hiercc_codes"))
    return bool(row.get("cglin_raw")) and row.get("cglin_status") not in {
        "missing",
        "unsupported",
        "malformed",
    }


def _signature(row: Mapping, kind: str) -> str:
    fields = {
        key: row.get(key)
        for key in (f"{kind}_scheme", f"{kind}_scheme_version", f"{kind}_database_sha256")
    }
    if kind == "cgmlst":
        fields.update(
            profile=row.get("cgmlst_profile"),
            novel=row.get("cgmlst_novel_alleles") or {},
            loci=sorted(row.get("cgmlst_loci") or []),
        )
    elif kind == "hiercc":
        fields["codes"] = row.get("hiercc_codes")
    else:
        fields["code"] = normalise_assignment(row).get("cglin_group_7")
        # Full code, not only the default selected prefix, detects conflicting narrow assignments.
        fields["raw"] = row.get("cglin_raw")
    fields["conflict"] = row.get(f"{kind}_status") == "conflict"
    return json.dumps(fields, sort_keys=True, separators=(",", ":"))


def _clear(row: Mapping, kind: str) -> dict:
    fields = {key: deepcopy(value) for key, value in row.items() if str(key).startswith(kind + "_")}
    fields[f"{kind}_status"] = "conflict"
    if kind == "cgmlst":
        fields.update(cgmlst_profile={}, cgmlst_novel_alleles={}, cgmlst_called_fraction=0)
    elif kind == "hiercc":
        fields["hiercc_codes"] = {}
    else:
        fields.update(cglin_raw="", cglin_code_status="conflict", cglin_resolved_depth=0)
        for depth in (5, 6, 7):
            fields[f"cglin_group_{depth}"] = ""
            fields[f"cglin_status_{depth}"] = "conflict"
    return fields


def annotate_public_typing(rows: Iterable[Mapping], assignments: Iterable[Mapping]) -> list[dict]:
    """Join only exact source IDs, checking assignments of deduplicated aliases."""
    by_id: dict[str, list[Mapping]] = {}
    for assignment in assignments:
        ident = _text(assignment.get("source_genome_id"))
        if not ident:
            raise PublicTypingError("Public assignment has no source_genome_id")
        by_id.setdefault(ident, []).append(assignment)
    output = []
    for row in rows:
        result = deepcopy(dict(row))
        ids = row.get("source_genome_ids") or [row.get("source_genome_id")]
        matches = [item for ident in ids for item in by_id.get(ident, [])]
        result["public_typing_record_count"] = len(matches)
        for kind in KINDS:
            available = [
                item for item in matches if any(str(key).startswith(kind + "_") for key in item)
            ]
            if not available:
                continue
            comparable = available + (
                [row] if _has(row, kind) or row.get(f"{kind}_status") == "conflict" else []
            )
            conflict = (
                any(item.get(f"{kind}_status") == "conflict" for item in comparable)
                or len({_signature(item, kind) for item in comparable}) > 1
            )
            if conflict:
                result.update(_clear(available[0], kind))
                result[f"{kind}_conflicting_assignments"] = deepcopy(comparable)
            else:
                result.update(
                    {
                        key: deepcopy(value)
                        for key, value in available[0].items()
                        if str(key).startswith(kind + "_")
                    }
                )
        output.append(result)
    return output


def resolve_focal_typing(
    focal_annotations: Iterable[Mapping], annotated_rows: Iterable[Mapping]
) -> list[dict]:
    """Copy public typing only through an already verified exact identity join."""
    by_id: dict[str, list[Mapping]] = {}
    for row in annotated_rows:
        for ident in row.get("source_genome_ids") or [row.get("source_genome_id")]:
            if ident:
                by_id.setdefault(ident, []).append(row)
    output = []
    for focal in focal_annotations:
        result = deepcopy(dict(focal))
        ident = _text(focal.get("cglin_matched_source_genome_id"))
        matches = by_id.get(ident, []) if ident else []
        if matches:
            joined = annotate_public_typing(
                [dict(focal, source_genome_id=ident)],
                [dict(row, source_genome_id=ident) for row in matches],
            )[0]
            for kind in ("cgmlst", "hiercc"):
                result.update(
                    {
                        key: deepcopy(value)
                        for key, value in joined.items()
                        if str(key).startswith(kind + "_")
                    }
                )
            result["public_typing_matched_source_genome_id"] = ident
        output.append(result)
    return output
