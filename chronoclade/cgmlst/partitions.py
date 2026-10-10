"""Explicit, disjoint lineage blocks; no adaptive narrowing or provider access."""

import hashlib
import json

from chronoclade.cglin import parse_code
from chronoclade.datasets import PreparedDataset
from chronoclade.typing_scopes import typing_scope


def profile_records(dataset: PreparedDataset) -> list[dict]:
    """Adapt a validated bundle to the existing categorical scientific engine."""
    records = {
        row["sample_id"]: dict(row, origin="local" if row["role"] == "input" else "context")
        for row in dataset.samples
    }
    for matrix in dataset.profiles:
        cat = matrix.catalogue
        for ident in matrix.sample_ids:
            row = records[ident]
            if "cgmlst_profile" in row:
                raise ValueError("cgmlst stage requires one profile scheme per sample")
            row.update(
                cgmlst_profile=matrix.profile(ident),
                cgmlst_loci=list(cat.loci),
                cgmlst_scheme=cat.scheme_id,
                cgmlst_scheme_version=cat.scheme_version,
                cgmlst_database_version=cat.database_version,
                cgmlst_database_sha256=cat.database_sha256,
                cgmlst_locus_universe_complete=cat.complete,
                cgmlst_status="resolved",
            )
    seen_lineages = set()
    for assignment in dataset.lineages:
        kind = assignment.get("kind")
        if kind not in {"cglin", "hiercc"}:
            raise ValueError("Lineage evidence requires kind cglin or hiercc")
        evidence = assignment.get("evidence") or {}
        if not isinstance(evidence, dict):
            raise ValueError("Lineage evidence must be a mapping")
        identity = (assignment["sample_id"], kind)
        if identity in seen_lineages:
            raise ValueError(
                f"Repeated lineage evidence for sample/kind {identity!r}; resolve conflicts before cgmlst"
            )
        seen_lineages.add(identity)
        row = records[assignment["sample_id"]]
        kind = assignment["kind"]
        row.update({key: value for key, value in evidence.items() if key.startswith(kind + "_")})
        row.update(
            {
                kind + "_scheme": assignment["scheme_id"],
                kind + "_scheme_version": assignment["scheme_version"],
                kind + "_database_version": assignment.get("database_version"),
                kind + "_status": assignment["resolution_status"],
                kind + ("_raw" if kind == "cglin" else "_codes"): assignment["assignment"],
            }
        )
    return list(records.values())


def _group_key(row, level, hiercc_level):
    kind = "hiercc" if hiercc_level else "cglin"
    if row.get(kind + "_status") in {
        "failed",
        "missing",
        "unassigned",
        "unsupported",
        "conflict",
        "malformed",
    }:
        return None
    scope = typing_scope(row, kind)
    if scope is None or not row.get("mlst_st") or not row.get("species"):
        return None
    scheme, version = scope
    if hiercc_level:
        code = (row.get("hiercc_codes") or {}).get(hiercc_level)
        if not str(code or "").isdigit():
            return None
        prefix = (str(code),)
    else:
        parts, status = parse_code(row.get("cglin_raw"))
        if status not in {"complete", "partial"} or len(parts) < level:
            return None
        prefix = tuple(parts[:level])
    species = " ".join(row["species"].replace("_", " ").casefold().split())
    # A known fingerprint and an unknown fingerprint do not share a block.
    return (
        species,
        str(row["mlst_st"]),
        kind,
        scheme,
        version,
        row.get(kind + "_database_version"),
        row.get(kind + "_database_sha256"),
        prefix,
    )


def partition_records(records, *, lin_level=5, hiercc_level=None):
    """Match every frozen context record inside an explicit input lineage block."""
    if lin_level not in {5, 6, 7}:
        raise ValueError("LIN level must be 5, 6 or 7")
    if hiercc_level is not None and (
        not hiercc_level.startswith("HC") or not hiercc_level[2:].isdigit()
    ):
        raise ValueError("HierCC level must be an explicit HC code, such as HC10")
    groups, unresolved = {}, []
    for row in records:
        if row["role"] != "input":
            continue
        key = _group_key(row, lin_level, hiercc_level)
        if key is None:
            unresolved.append(
                {"sample_id": row["sample_id"], "reason": "missing_or_unresolved_partition_typing"}
            )
        else:
            groups.setdefault(key, [])
    for row in records:
        key = _group_key(row, lin_level, hiercc_level)
        if key in groups:
            groups[key].append(row)
    blocks = []
    for key, rows in sorted(groups.items(), key=lambda item: repr(item[0])):
        digest = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()[:12]
        prefix = list(key[-1])
        blocks.append(
            {
                "block_id": f"{key[2]}_{digest}",
                "species": rows[0]["species"].replace("_", " "),
                "mlst_st": key[1],
                "kind": key[2],
                "scheme": key[3],
                "scheme_version": key[4],
                "database_version": key[5],
                "database_sha256": key[6],
                "prefix": prefix,
                "level": hiercc_level or lin_level,
                "records": sorted(rows, key=lambda r: r["sample_id"]),
            }
        )
    counts = []
    if not hiercc_level:
        for level in (5, 6, 7):
            keys = {_group_key(row, level, None) for row in records if row["role"] == "input"} - {
                None
            }
            for key in sorted(keys, key=repr):
                matches = [row for row in records if _group_key(row, level, None) == key]
                counts.append(
                    {
                        "level": level,
                        "species": key[0],
                        "mlst_st": key[1],
                        "scheme": key[3],
                        "scheme_version": key[4],
                        "prefix": list(key[-1]),
                        "input_count": sum(row["role"] == "input" for row in matches),
                        "context_count": sum(row["role"] == "context" for row in matches),
                    }
                )
    used = {row["sample_id"] for block in blocks for row in block["records"]}
    return blocks, {
        "policy": "explicit-lineage-depth",
        "lin_level": lin_level,
        "hiercc_level": hiercc_level,
        "unresolved_inputs": unresolved,
        "available_level_counts": counts,
        "unmatched_context_ids": sorted(
            row["sample_id"]
            for row in records
            if row["role"] == "context" and row["sample_id"] not in used
        ),
        "context_scope": "all matching context already present in the frozen bundle",
        "live_context_retrieval": "not_requested",
    }
