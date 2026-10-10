"""Refine a frozen same-ST pool without inferring missing lineage assignments.

Lineage matches receive priority, compatible allele profiles order those matches
and may nominate additional candidates. A reserved country/year-balanced budget
retains wider ST context. Distances apply only to available, callable loci in this
pool: they are neither whole-genome distances nor global nearest neighbours.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import replace
import math
import re
from typing import Any, Iterable, Mapping

from chronoclade.typing_scopes import compatible_typing as _compatible, typing_scope as _scope
from chronoclade.cglin import DEFAULT_SCHEME, parse_code
from chronoclade.context import ContextCandidate, stratified_candidate_pool

_MISSING = {"", "unknown", "none", "null", "na", "n/a", "?", "-"}


def _text(value: Any) -> str:
    return "" if value is None else str(value).strip()


def _lineage(row: Mapping, kind: str, level: int | str) -> tuple | None:
    scope = _scope(row, kind)
    if scope is None or _text(row.get(f"{kind}_status")).casefold() in {
        "conflict",
        "malformed",
        "unsupported",
        "missing",
        "unassigned",
        "failed",
    }:
        return None
    if kind == "cglin":
        components, status = parse_code(row.get("cglin_raw"))
        if scope[0] != DEFAULT_SCHEME or status not in {"complete", "partial"}:
            return None
        return (*scope, components[:level]) if len(components) >= level else None
    codes = row.get("hiercc_codes")
    if not isinstance(codes, Mapping):
        return None
    code = _text(codes.get(level))
    # A provisional/unknown cluster is not a resolved HierCC label.
    return (*scope, code) if code.isdigit() else None


def _profile(row: Mapping) -> Mapping | None:
    if row.get("cgmlst_status") in {"conflict", "failed", "unassigned", "unsupported"}:
        return None
    profile = row.get("cgmlst_profile")
    novel = row.get("cgmlst_novel_alleles") or {}
    if not isinstance(profile, Mapping) or not isinstance(novel, Mapping):
        return None
    if set(profile).intersection(novel):
        return None
    combined = dict(profile, **novel)
    return combined or None


def _loci(row: Mapping, profile: Mapping) -> set[str] | None:
    declared = row.get("cgmlst_loci")
    if declared is None:
        # Imported dense profiles must enumerate their complete scheme.
        return set(profile)
    if not isinstance(declared, (list, tuple)) or not declared:
        return None
    loci = set(declared)
    if len(loci) != len(declared) or not set(profile).issubset(loci):
        return None
    return loci


def _allele(value: Any) -> str | None:
    value = _text(value)
    # Official allele identifiers are positive integers; zero is an absent call.
    if re.fullmatch(r"[0-9a-fA-F]{40}", value):
        return value.lower()
    return str(int(value)) if value.isdigit() and int(value) > 0 else None


def _compare(query: Mapping, public: Mapping, depth: int, level: str, overlap: float) -> dict:
    evidence: dict[str, Any] = {
        "query_id": _text(query.get("sample_id")),
        "matches": [],
        "assignment_scopes": {
            kind: {"query": _scope(query, kind), "candidate": _scope(public, kind)}
            for kind in ("cglin", "hiercc", "cgmlst")
        },
    }
    reasons = []
    for kind, selected in (("cglin", depth), ("hiercc", level)):
        left, right = _lineage(query, kind, selected), _lineage(public, kind, selected)
        if left is None or right is None:
            reasons.append(f"{kind}_unassigned_or_unversioned")
        elif not _compatible(query, public, kind):
            reasons.append(f"{kind}_incompatible_scheme_version")
        elif left == right:
            evidence["matches"].append(kind)
        else:
            reasons.append(f"{kind}_different_group")
    left, right = _profile(query), _profile(public)
    if left is None or right is None:
        reasons.append("cgmlst_profile_unavailable")
    elif _scope(query, "cgmlst") is None or _scope(public, "cgmlst") is None:
        reasons.append("cgmlst_unversioned")
    elif not _compatible(query, public, "cgmlst"):
        reasons.append("cgmlst_incompatible_scheme_version")
    elif _loci(query, left) is None or _loci(query, left) != _loci(public, right):
        reasons.append("cgmlst_incompatible_locus_set")
    else:
        pairs = [
            (_allele(left.get(key)), _allele(right.get(key))) for key in sorted(_loci(query, left))
        ]
        called = [(a, b) for a, b in pairs if a is not None and b is not None]
        shared = len(called)
        fraction = shared / len(pairs)
        evidence.update(shared_called_loci=shared, scheme_loci=len(pairs), call_overlap=fraction)
        if fraction < overlap or not shared:
            reasons.append("cgmlst_insufficient_called_overlap")
        else:
            differences = sum(a != b for a, b in called)
            evidence.update(
                allele_differences=differences, allele_difference_fraction=differences / shared
            )
    evidence["reasons"] = reasons
    evidence["priority"] = bool(evidence["matches"] or "allele_differences" in evidence)
    return evidence


def refine_candidate_pool(
    candidates: list[ContextCandidate],
    catalogue_rows: Iterable[Mapping],
    focal_assignments: Iterable[Mapping],
    *,
    limit: int,
    seed: int,
    cglin_depth: int = 7,
    hiercc_level: str = "HC1100",
    background_fraction: float = 0.25,
    min_profile_overlap: float = 0.9,
) -> tuple[list[ContextCandidate], dict]:
    """Prioritise query-relative evidence within an already eligible same-ST pool.

    Requires explicit matching scheme versions, full cgLIN prefixes and exact
    HierCC level labels. Unversioned cgLIN assignments may match only within
    the same frozen export identity (exportsha256), without asserting a database
    version or compatibility with independent exports. Profiles must enumerate the same complete locus set;
    missing allele calls never count as matches. At least ceil(limit * fraction)
    slots (subject to pool size) remain available to balanced wider-ST sampling.
    Queries take turns selecting evidence-ranked candidates, then the remaining
    slots are filled using the existing seeded country/year-balanced sampler.
    Input assignments and candidate records are not modified.
    """
    if limit < 1 or not 1 <= cglin_depth <= 10:
        raise ValueError("Positive pool limit and cgLIN depth between 1 and 10 required")
    if not math.isfinite(background_fraction) or not 0 < background_fraction <= 1:
        raise ValueError("background_fraction must be greater than zero and at most one")
    if not math.isfinite(min_profile_overlap) or not 0 < min_profile_overlap <= 1:
        raise ValueError("min_profile_overlap must be greater than zero and at most one")
    if not hiercc_level.startswith("HC") or not hiercc_level[2:].isdigit():
        raise ValueError("HierCC level must be an explicit HC label, e.g. HC1100")
    if len({c.sample_id for c in candidates}) != len(candidates):
        raise ValueError("Candidate sample IDs must be unique")
    rows: dict[str, Mapping] = {}
    for row in catalogue_rows:
        ident = _text(row.get("source_genome_id"))
        if not ident or ident in rows:
            raise ValueError("Catalogue source genome IDs must be nonempty and unique")
        rows[ident] = row
    queries = sorted(focal_assignments, key=lambda row: _text(row.get("sample_id")))
    query_ids = [_text(row.get("sample_id")) for row in queries]
    if any(not ident for ident in query_ids) or len(set(query_ids)) != len(query_ids):
        raise ValueError("Focal assignment sample IDs must be nonempty and unique")
    evidence = {}
    rankings = {ident: [] for ident in query_ids}
    for candidate in sorted(candidates, key=lambda candidate: candidate.sample_id):
        public = rows.get(candidate.source_genome_id, {})
        comparisons = [
            _compare(query, public, cglin_depth, hiercc_level, min_profile_overlap)
            for query in queries
        ]
        evidence[candidate.sample_id] = {
            "source_genome_id": candidate.source_genome_id,
            "comparisons": comparisons,
            "selected": False,
        }
        for comparison in comparisons:
            if comparison["priority"]:
                rank = (
                    0 if comparison["matches"] else 1,
                    comparison.get("allele_difference_fraction", math.inf),
                    comparison.get("allele_differences", math.inf),
                    -comparison.get("shared_called_loci", 0),
                    candidate.sample_id,
                )
                rankings[comparison["query_id"]].append((rank, candidate, comparison))
    for ranking in rankings.values():
        ranking.sort(key=lambda item: item[0])
    budget = min(limit, len(candidates))
    reserved = min(budget, math.ceil(budget * background_fraction))
    selected = []
    selected_ids = set()
    while len(selected) < budget - reserved:
        progressed = False
        for ident in query_ids:
            ranking = rankings[ident]
            while ranking and ranking[0][1].sample_id in selected_ids:
                ranking.pop(0)
            if not ranking:
                continue
            _, candidate, comparison = ranking.pop(0)
            reason = (
                "lineage_priority:" + ident if comparison["matches"] else "cgmlst_priority:" + ident
            )
            selected.append(replace(candidate, selection_reason=reason))
            selected_ids.add(candidate.sample_id)
            evidence[candidate.sample_id].update(selected=True, selection_reason=reason)
            progressed = True
            if len(selected) == budget - reserved:
                break
        if not progressed:
            break
    priority_count = len(selected)
    remaining = [c for c in candidates if c.sample_id not in selected_ids]
    if budget > len(selected):
        background = stratified_candidate_pool(remaining, limit=budget - len(selected), seed=seed)
        for candidate in background:
            reason = (
                "balanced_same_st_background" if priority_count else "balanced_same_st_fallback"
            )
            selected.append(replace(candidate, selection_reason=reason))
            evidence[candidate.sample_id].update(selected=True, selection_reason=reason)
    counts = Counter(
        reason
        for item in evidence.values()
        for comparison in item["comparisons"]
        for reason in comparison["reasons"]
    )
    audit = {
        "method": "query-lineage-priority-with-country-year-balanced-same-ST-background",
        "scope": "available eligible genomes within the frozen same-ST pool; not global nearest neighbours",
        "settings": {
            "cglin_depth": cglin_depth,
            "hiercc_level": hiercc_level,
            "background_fraction": background_fraction,
            "min_profile_overlap": min_profile_overlap,
            "seed": seed,
            "limit": limit,
        },
        "counts": {
            "eligible_same_st": len(candidates),
            "selected": len(selected),
            "priority_selected": priority_count,
            "background_selected": len(selected) - priority_count,
            "reserved_background": reserved,
            "priority_eligible": sum(
                any(c["priority"] for c in item["comparisons"]) for item in evidence.values()
            ),
            "comparison_reasons": dict(sorted(counts.items())),
        },
        "query_status": {
            ident: "has_comparable_evidence"
            if any(
                any(c["query_id"] == ident and c["priority"] for c in item["comparisons"])
                for item in evidence.values()
            )
            else "no_comparable_evidence_same_st_fallback"
            for ident in query_ids
        },
        "candidates": evidence,
    }
    return selected, audit
