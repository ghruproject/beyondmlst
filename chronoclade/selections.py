"""Representative selection shared by profile exploration and assembly stages.

Scientific distance evidence is supplied by the caller. Stage orchestration and
report generation do not belong here.
"""

from __future__ import annotations

import csv
import math
import random
from collections import defaultdict
from pathlib import Path

from chronoclade.errors import WorkflowError


def selected_context_records(
    selection: dict, queries: list[dict], context: list[dict]
) -> list[dict]:
    """Validate a saved selection against its exact pool before acquiring assemblies."""
    by_id = {row["sample_id"]: row for row in context}
    context_ids = selection.get("selected_context_ids", [])
    sample_ids = selection.get("selected_sample_ids", [])
    query_ids = [row["sample_id"] for row in queries]
    if (
        len(by_id) != len(context)
        or selection.get("available_context_ids") != sorted(by_id)
        or selection.get("query_ids") != query_ids
        or len(context_ids) != len(set(context_ids))
        or not set(context_ids).issubset(by_id)
        or sample_ids != query_ids + context_ids
        or len(sample_ids) != len(set(sample_ids))
    ):
        raise WorkflowError("Shared context selection does not match the resolved query lineage pool")
    decisions = selection.get("decisions", [])
    if (
        [decision.get("sample_id") for decision in decisions] != context_ids
        or any(not isinstance(decision.get("reason"), str) or not decision["reason"]
               for decision in decisions)
    ):
        raise WorkflowError("Shared context selection decisions do not match selected identifiers")
    return [
        dict(by_id[decision["sample_id"]], selection_reason=decision["reason"])
        for decision in decisions
    ]


def _distances(analysis: dict) -> dict[tuple[str, str], float]:
    path = analysis.get("paths", {}).get("pairwise_distances")
    result = {}
    if path and Path(path).is_file():
        with Path(path).open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                try:
                    distance = float(row["distance"])
                except (ValueError, TypeError, KeyError):
                    continue
                if math.isfinite(distance):
                    a, b = row["sample_id_1"], row["sample_id_2"]
                    result[a, b] = result[b, a] = distance
    return result


def select_assembly_context(
    queries: list[dict],
    context: list[dict],
    analysis: dict,
    *,
    size: int = 50,
    nearest_per_query: int = 3,
    include: list[str] | None = None,
    seed: int = 42,
) -> tuple[list[dict], dict]:
    """Select fairly around queries, then across time/place and genetic diversity."""
    if size < 0 or nearest_per_query < 1:
        raise WorkflowError("Context size must be non-negative and nearest-per-query positive")
    by_id = {row["sample_id"]: row for row in context}
    if len(by_id) != len(context):
        raise WorkflowError("Duplicate contextual sample identifiers")
    include = include or []
    pinned = []
    for identifier in include:
        matches = [
            r["sample_id"]
            for r in context
            if identifier
            in {r["sample_id"], str(r.get("source_genome_id", "")), str(r.get("accession", ""))}
        ]
        if len(matches) != 1:
            raise WorkflowError(
                f"Requested context genome {identifier!r} is not uniquely in the pool"
            )
        if matches[0] not in pinned:
            pinned.append(matches[0])
    if len(pinned) > size:
        raise WorkflowError("Requested context genomes exceed --context-size; increase the budget")
    distances = _distances(analysis)
    reasons: dict[str, str] = {identifier: "user_requested" for identifier in pinned}
    selected = list(pinned)
    queues = {
        q["sample_id"]: sorted(
            [identifier for identifier in by_id if (q["sample_id"], identifier) in distances],
            key=lambda identifier: (distances[q["sample_id"], identifier], identifier),
        )[:nearest_per_query]
        for q in sorted(queries, key=lambda r: r["sample_id"])
    }
    nearest_budget = min(size, max(len(queries), math.ceil(size / 2)))
    for rank in range(nearest_per_query):
        for query_id, candidates in queues.items():
            if len(selected) >= nearest_budget:
                break
            if rank < len(candidates) and candidates[rank] not in reasons:
                identifier = candidates[rank]
                selected.append(identifier)
                reasons[identifier] = f"cgmlst_neighbour:{query_id}:rank={rank + 1}"
    groups = {
        identifier: str(group["group_id"])
        for group in analysis.get("genetic_groups", [])
        for identifier in group.get("sample_ids", [])
    }
    strata = defaultdict(list)
    for identifier, row in sorted(by_id.items()):
        key = (
            groups.get(identifier, "unassigned"),
            str(row.get("collection_date", ""))[:4] or "Unknown",
            str(row.get("nuts2") or row.get("region") or row.get("country") or "Unknown"),
        )
        strata[key].append(identifier)
    cells = sorted(strata)
    random.Random(seed).shuffle(cells)
    target = min(size, len(selected) + math.ceil(size / 4))
    for cell in cells:
        if len(selected) >= target:
            break
        candidates = [identifier for identifier in strata[cell] if identifier not in reasons]
        if candidates:
            identifier = candidates[0]
            selected.append(identifier)
            reasons[identifier] = "genetic_group_time_region_representative"
    while len(selected) < min(size, len(by_id)):
        remaining = [identifier for identifier in by_id if identifier not in reasons]
        comparable = [
            identifier
            for identifier in remaining
            if any((identifier, other) in distances for other in selected)
        ]
        if comparable:
            identifier = max(
                comparable,
                key=lambda candidate: (
                    min(
                        distances[candidate, other]
                        for other in selected
                        if (candidate, other) in distances
                    ),
                    candidate,
                ),
            )
            reason = "cgmlst_diversity_representative"
        else:
            identifier = sorted(remaining)[0]
            reason = "metadata_background_no_comparable_profile"
        selected.append(identifier)
        reasons[identifier] = reason
    records = [
        dict(by_id[identifier], selection_reason=reasons[identifier]) for identifier in selected
    ]
    return records, {
        "context_budget": size,
        "available_contexts": len(context),
        "selected_contexts": len(records),
        "nearest_per_query": nearest_per_query,
        "seed": seed,
        "decisions": [
            {"sample_id": identifier, "reason": reasons[identifier]} for identifier in selected
        ],
        "queries_without_selected_comparable_neighbour": [
            q["sample_id"]
            for q in queries
            if not any((q["sample_id"], identifier) in distances for identifier in selected)
        ],
        "scope": "Selection within the recorded profile catalogue; diversity uses cgMLST distance, not calendar time.",
    }
