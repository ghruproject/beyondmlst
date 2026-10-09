"""Exploratory allele-profile analysis; no assembly downloads or dating claims.

Allele identifiers are categorical. Every distance reports its callable denominator.
Separate complete-comparability cohorts prevent missing distances being imputed.
"""

from __future__ import annotations

import calendar
import copy
import csv
import datetime as dt
import json
from collections import Counter
from pathlib import Path

import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from chronoclade.sample_labels import sample_labels
from chronoclade.context_refinement import _allele, _compatible, _lineage, _loci, _profile, _scope


MAX_RECORDS = 1500


def draw_profile_tree(tree, records, labels, path):
    """Draw neutral branches with dataset-coloured tips and supplied metadata."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    by_id = {row["sample_id"]: row for row in records}
    colours = {"input": "#2166ac", "context": "#666666"}
    tip_labels, label_colours = {}, {}
    for tip in tree.get_terminals():
        row = by_id[tip.name]
        dataset = "input" if row.get("origin") in {"local", "query", "focal"} else "context"
        country = str(row.get("country") or "Country unknown")
        date = str(row.get("collection_date") or row.get("collection_year") or "Date unknown")
        label = f"{labels.get(tip.name, tip.name)} | {country} | {date}"
        tip_labels[tip.name] = label
        label_colours[label] = colours[dataset]

    fig, ax = plt.subplots(figsize=(13, max(4, len(records) * 0.24)))
    Phylo.draw(
        tree, axes=ax, do_show=False,
        label_func=lambda clade: tip_labels.get(clade.name),
        label_colors=label_colours,
    )
    # Draw a coloured marker at each terminal without colouring ancestral branches.
    depths = tree.depths()
    if not max(depths.values()):
        depths = tree.depths(unit_branch_lengths=True)
    terminals = tree.get_terminals()
    for index, tip in enumerate(terminals):
        ax.plot(depths[tip], index + 1, "o", markersize=4,
                color=label_colours[tip_labels[tip.name]], zorder=3)
    handles = [Line2D([], [], color=colour, marker="o", linestyle="None", label=label)
               for colour, label in ((colours["input"], "Input genomes"),
                                     (colours["context"], "Public comparisons"))]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, 1.01),
              frameon=False, ncol=2)
    ax.set_xlabel("Fraction of mismatching callable cgMLST loci")
    ax.set_ylabel("")
    ax.set_yticks([])
    left, right = ax.get_xlim()
    ax.set_xlim(left, right + (right - left) * 0.35)
    ax.set_title("Neighbour-joining tree · tip labels: genome | country | collection date", pad=38)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def _date_interval(row, *, allow_future=False):
    value = str(row.get("collection_date") or row.get("collection_year") or "").strip()
    try:
        parts = value.split("-")
        year = int(parts[0])
        if len(parts) == 1:
            start, end, precision = dt.date(year, 1, 1), dt.date(year, 12, 31), "year"
        elif len(parts) == 2:
            month = int(parts[1])
            start = dt.date(year, month, 1)
            end = dt.date(year, month, calendar.monthrange(year, month)[1])
            precision = "month"
        elif len(parts) == 3:
            start = end = dt.date.fromisoformat(value)
            precision = "day"
        else:
            return None

        if not allow_future and start > dt.date.today():
            return None

        def decimal(day):
            return day.year + (day - dt.date(day.year, 1, 1)).days / (
                366 if calendar.isleap(day.year) else 365
            )

        return {
            "start": start.isoformat(),
            "end": end.isoformat(),
            "precision": precision,
            "year_min": decimal(start),
            "year_max": decimal(end),
            "year": year,
        }
    except (ValueError, TypeError, OverflowError):
        return None


def _pair(left, right, overlap):
    if str(left.get("species", "")).casefold().replace("_", " ") != str(
        right.get("species", "")
    ).casefold().replace("_", " "):
        return None, "different_species"
    a, b = _profile(left), _profile(right)
    if a is None or b is None:
        return None, "profile_unavailable"
    if _scope(left, "cgmlst") is None or not _compatible(left, right, "cgmlst"):
        return None, "incompatible_or_unversioned_scheme"
    loci = _loci(left, a)
    if not loci or loci != _loci(right, b):
        return None, "incompatible_locus_set"
    calls = [(_allele(a.get(locus)), _allele(b.get(locus))) for locus in sorted(loci)]
    shared = [(x, y) for x, y in calls if x is not None and y is not None]
    if not shared or len(shared) / len(loci) < overlap:
        return None, "insufficient_called_overlap"
    differences = sum(x != y for x, y in shared)
    return {
        "allele_differences": differences,
        "shared_called_loci": len(shared),
        "scheme_loci": len(loci),
        "call_overlap": len(shared) / len(loci),
        "distance": differences / len(shared),
    }, None


def _clusters(matrix, threshold):
    """Deterministic complete linkage: all members satisfy the distance threshold."""
    groups = [(i,) for i in range(len(matrix))]
    linkage = matrix.copy()
    np.fill_diagonal(linkage, np.inf)
    while len(groups) > 1:
        position = int(np.argmin(linkage))
        i, j = np.unravel_index(position, linkage.shape)
        if linkage[i, j] > threshold:
            break
        if i > j:
            i, j = j, i
        groups[i] = tuple(sorted(groups[i] + groups[j]))
        groups.pop(j)
        linkage[i, :] = np.maximum(linkage[i, :], linkage[j, :])
        linkage[:, i] = linkage[i, :]
        linkage = np.delete(np.delete(linkage, j, axis=0), j, axis=1)
        np.fill_diagonal(linkage, np.inf)
    return sorted(groups)


def _csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: json.dumps(value, sort_keys=True)
                    if isinstance(value, (dict, list))
                    else value
                    for key, value in row.items()
                }
            )
    return str(path)


def _location_network(tree, rows):
    """Optimal parsimony state ambiguity across roots; no probability interpretation."""
    countries = {r["sample_id"]: str(r.get("country") or "Unknown") for r in rows}
    roots = sorted(countries)[:5]
    counts = Counter()
    ambiguous = 0
    for root in roots:
        rooted = copy.deepcopy(tree)
        rooted.root_with_outgroup(root)
        locations = sorted(set(countries.values()) - {"Unknown"})
        if not locations:
            continue
        costs = {}
        for node in rooted.find_clades(order="postorder"):
            if node.is_terminal():
                observed = countries[node.name]
                costs[node] = {
                    state: 0 if observed in {state, "Unknown"} else float("inf")
                    for state in locations
                }
            else:
                costs[node] = {
                    state: sum(
                        min(costs[child][target] + (state != target) for target in locations)
                        for child in node.clades
                    )
                    for state in locations
                }
        optimum = min(costs[rooted.root].values())
        states = {
            rooted.root: {state for state in locations if costs[rooted.root][state] == optimum}
        }
        possible = set()
        for parent in rooted.find_clades(order="preorder"):
            for child in parent.clades:
                allowed = set()
                for source in states[parent]:
                    best = min(costs[child][target] + (source != target) for target in locations)
                    for target in locations:
                        if costs[child][target] + (source != target) == best:
                            allowed.add((source, target))
                states[child] = {target for _, target in allowed}
                if len(allowed) > 1:
                    ambiguous += 1
                possible.update((source, target) for source, target in allowed if source != target)
        counts.update(possible)
    return {
        "method": "ambiguous maximum-parsimony reconstruction across alternate tip roots",
        "tested_roots": roots,
        "ambiguous_edge_reconstructions": ambiguous,
        "interpretation": "Possible location changes, not confirmed transmission. Root fractions "
        "are sensitivity summaries, not probabilities or confidence intervals.",
        "edges": [
            {
                "source": a,
                "target": b,
                "roots_with_possible_change": count,
                "roots_tested": len(roots),
                "root_fraction": count / len(roots),
                "uncertain": True,
            }
            for (a, b), count in sorted(counts.items())
        ],
    }


def _plots(output, cohort_id, coords, rows, network):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paths = {}
    fig, ax = plt.subplots(figsize=(7, 5))
    for origin in sorted({str(r.get("origin", "unknown")) for r in rows}):
        indexes = [i for i, r in enumerate(rows) if str(r.get("origin", "unknown")) == origin]
        ax.scatter(coords[indexes, 0], coords[indexes, 1], label=origin, alpha=0.75)
    ax.set(xlabel="PCoA axis 1", ylabel="PCoA axis 2", title="cgMLST distance ordination")
    ax.legend()
    fig.tight_layout()
    path = output / f"{cohort_id}_pcoa.svg"
    fig.savefig(path)
    plt.close(fig)
    paths["pcoa_figure"] = str(path)
    counts = Counter(str(r.get("country") or "Unknown") for r in rows)
    fig, ax = plt.subplots(figsize=(7, max(3, len(counts) * 0.3)))
    keys = sorted(counts, key=lambda key: (-counts[key], key))
    ax.barh(keys, [counts[key] for key in keys])
    ax.invert_yaxis()
    ax.set(xlabel="Genomes in analysed cohort", title="Available country metadata")
    fig.tight_layout()
    path = output / f"{cohort_id}_countries.svg"
    fig.savefig(path)
    plt.close(fig)
    paths["country_figure"] = str(path)
    fig, ax = plt.subplots(figsize=(7, 5))
    countries = sorted({e[key] for e in network["edges"] for key in ("source", "target")})
    xy = {
        name: (
            np.cos(2 * np.pi * i / max(1, len(countries))),
            np.sin(2 * np.pi * i / max(1, len(countries))),
        )
        for i, name in enumerate(countries)
    }
    for edge in network["edges"]:
        ax.annotate(
            "",
            xy=xy[edge["target"]],
            xytext=xy[edge["source"]],
            arrowprops={
                "arrowstyle": "->",
                "linestyle": "--",
                "alpha": 0.5,
                "connectionstyle": "arc3,rad=.15",
            },
        )
    for country, (x, y) in xy.items():
        ax.scatter([x], [y], s=100)
        ax.text(x, y + 0.08, country, ha="center")
    if not countries:
        ax.text(
            0.5,
            0.5,
            "No location changes inferred from available metadata",
            ha="center",
            transform=ax.transAxes,
        )
    ax.set_title("Exploratory possible location changes\nAll links uncertain; root-dependent")
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.3, 1.4)
    ax.axis("off")
    fig.tight_layout()
    path = output / f"{cohort_id}_location_network.svg"
    fig.savefig(path)
    plt.close(fig)
    paths["network_figure"] = str(path)
    return paths


def analyse_profiles(
    records: list[dict],
    *,
    output: Path,
    seed: int = 42,
    min_overlap: float = 0.9,
    bootstrap_replicates: int = 30,
    distance_threshold: float = 0.02,
) -> dict:
    """Write a reproducible, descriptive report using available compatible profiles.

    Distances are fractions of mismatching jointly called loci. Genetic groups use
    complete linkage and carry compatible LIN/HierCC labels as annotations. Genetic
    group persistence and time/place concentration are independent summaries.
    """
    if len(records) > MAX_RECORDS:
        raise ValueError(
            f"Profile analysis is limited to {MAX_RECORDS} records; explicitly bound "
            "the catalogue first. No records have been silently discarded."
        )
    if not 0 < min_overlap <= 1 or not 0 <= distance_threshold <= 1:
        raise ValueError("Overlap must be in (0,1] and distance threshold in [0,1]")
    if not isinstance(bootstrap_replicates, int) or not 0 <= bootstrap_replicates <= 200:
        raise ValueError("bootstrap_replicates must be an integer between 0 and 200")
    ids = [r.get("sample_id") for r in records]
    if any(not isinstance(i, str) or not i for i in ids) or len(set(ids)) != len(ids):
        raise ValueError("Sample identifiers must be nonempty and unique")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    records = [dict(row) for row in sorted(records, key=lambda r: r["sample_id"])]
    observed_universes = {}
    for row in records:
        if row.get("cgmlst_locus_universe_complete") is False:
            scope = (
                str(row.get("species", "")).casefold().replace("_", " "),
                _scope(row, "cgmlst"),
                row.get("cgmlst_database_sha256", ""),
            )
            observed_universes.setdefault(scope, set()).update(row.get("cgmlst_loci") or [])
            observed_universes[scope].update((_profile(row) or {}).keys())
    for row in records:
        if row.get("cgmlst_locus_universe_complete") is False:
            scope = (
                str(row.get("species", "")).casefold().replace("_", " "),
                _scope(row, "cgmlst"),
                row.get("cgmlst_database_sha256", ""),
            )
            row["cgmlst_loci"] = sorted(observed_universes[scope])
    pairs, incompatible, exclusions = [], Counter(), []
    distances = np.full((len(records), len(records)), np.nan)
    np.fill_diagonal(distances, 0)
    eligible = []
    for i, row in enumerate(records):
        evidence, reason = _pair(row, row, min_overlap)
        if evidence is None:
            exclusions.append({"sample_id": row["sample_id"], "reason": reason})
        else:
            eligible.append(i)
        for j in range(i):
            evidence, reason = _pair(records[j], row, min_overlap)
            if evidence is None:
                incompatible[reason] += 1
            else:
                distances[i, j] = distances[j, i] = evidence["distance"]
                pairs.append(
                    {"sample_a": records[j]["sample_id"], "sample_b": row["sample_id"], **evidence}
                )
    # Every cohort has all pairwise distances. Greedy partition is deterministic,
    # not a biological grouping, and its limitation is reported below.
    cohorts = []
    for index in eligible:
        for cohort in cohorts:
            if all(np.isfinite(distances[index, other]) for other in cohort):
                cohort.append(index)
                break
        else:
            cohorts.append([index])
    cohort_membership = {
        index: f"cohort_{number}" for number, cohort in enumerate(cohorts, 1) for index in cohort
    }
    date_exclusions = []
    for row in records:
        raw_date = _date_interval(row, allow_future=True)
        if raw_date and dt.date.fromisoformat(raw_date["start"]) > dt.date.today():
            date_exclusions.append(
                {
                    "sample_id": row["sample_id"],
                    "collection_date": row.get("collection_date") or row.get("collection_year"),
                    "reason": "future_collection_date",
                }
            )
    future_ids = {row["sample_id"] for row in date_exclusions}
    nearest = []
    for i, row in enumerate(records):
        if row.get("origin") not in {"local", "query", "focal"}:
            continue
        choices = [
            (distances[i, j], j)
            for j, other in enumerate(records)
            if other.get("origin") == "context" and np.isfinite(distances[i, j])
        ]
        if not choices:
            nearest.append({"query_id": row["sample_id"], "status": "no_comparable_context"})
            continue
        minimum = min(value for value, _ in choices)
        for value, j in choices:
            if value == minimum:
                evidence, _ = _pair(row, records[j], min_overlap)
                nearest.append(
                    {
                        "query_id": row["sample_id"],
                        "context_id": records[j]["sample_id"],
                        "status": "matched",
                        "country": records[j].get("country") or "Unknown",
                        "region": records[j].get("region") or "",
                        "nuts2": records[j].get("nuts2") or "",
                        "collection_date": records[j].get("collection_date")
                        or records[j].get("collection_year")
                        or "",
                        "year": (_date_interval(records[j]) or {}).get("year"),
                        "date_status": "future_collection_date"
                        if records[j]["sample_id"] in future_ids
                        else ("valid" if _date_interval(records[j]) else "missing_or_invalid"),
                        "cohort_id": cohort_membership.get(j),
                        "query_cohort_id": cohort_membership.get(i),
                        "tied_neighbours": int(sum(v == minimum for v, _ in choices)),
                        **evidence,
                    }
                )
    summary = {
        "schema_version": 1,
        "records_count": len(records),
        "sample_labels": sample_labels(records),
        "seed": seed,
        "min_overlap": min_overlap,
        "distance_threshold": distance_threshold,
        "bootstrap_replicates": bootstrap_replicates,
        "coverage": {
            "available_profiles": len(eligible),
            "comparable_pairs": len(pairs),
            "unavailable_pairs_by_reason": dict(incompatible),
        },
        "exclusions": exclusions,
        "date_exclusions": date_exclusions,
        "cohorts": [],
        "nearest_neighbours": nearest,
        "genetic_groups": [],
        "temporal_persistence": [],
        "time_place_concentration": [],
        "root_to_tip": [],
        "location_network": [],
        "paths": {},
        "limitations": [
            "Distances and nearest neighbours apply only to the supplied catalogue.",
            "Missing calls are excluded, never treated as matches.",
            "Complete-comparability cohorts are a deterministic partition, not genetic clusters; "
            "comparable pairs across cohorts remain in the distance table.",
            "Metadata summaries describe sampled records, not population prevalence or outbreaks.",
            "NJ date diagnostics do not assess the final recombination-corrected tree.",
        ],
    }
    rng = np.random.default_rng(seed)
    for number, indexes in enumerate(cohorts, 1):
        cohort_id = f"cohort_{number}"
        rows = [records[i] for i in indexes]
        matrix = distances[np.ix_(indexes, indexes)]
        cohort = {
            "cohort_id": cohort_id,
            "sample_ids": [r["sample_id"] for r in rows],
            "tree_units": "fraction of mismatching jointly called cgMLST loci",
            "warnings": [],
        }
        if any(row["sample_id"] in future_ids for row in rows):
            cohort["warnings"].append(
                "Impossible future collection dates excluded from temporal summaries and NJ date diagnostics; raw metadata retained."
            )
        centered = np.eye(len(rows)) - np.ones((len(rows), len(rows))) / len(rows)
        gram = -0.5 * centered @ (matrix**2) @ centered
        eigenvalues, eigenvectors = np.linalg.eigh(gram)
        order = np.argsort(eigenvalues)[::-1]
        eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]
        positive = np.maximum(eigenvalues, 0)
        coords = np.zeros((len(rows), 2))
        n = min(2, len(rows))
        coords[:, :n] = eigenvectors[:, :n] * np.sqrt(positive[:n])
        cohort["pcoa_eigenvalues"] = eigenvalues.tolist()
        cohort["pcoa_positive_axis_fraction"] = (
            (positive[:2] / positive.sum()).tolist() if positive.sum() else [0, 0]
        )
        if np.any(eigenvalues < -1e-10):
            cohort["warnings"].append(
                "Non-Euclidean distances: negative PCoA eigenvalues; scatterplot is an approximation."
            )
        groups = _clusters(matrix, distance_threshold)
        if any(r.get("cgmlst_locus_universe_complete") is False for r in rows):
            cohort["warnings"].append(
                "Coverage denominator is the union of observed export loci; full canonical scheme coverage is unknown."
            )
            cohort["coverage_denominator"] = (
                "observed export locus union; canonical completeness unknown"
            )
        else:
            cohort["coverage_denominator"] = "declared locus universe"
        profiles = [_profile(r) for r in rows]
        loci = sorted(_loci(rows[0], profiles[0]))
        calls = [[_allele(profile.get(locus)) for locus in loci] for profile in profiles]
        pair_support = np.zeros_like(matrix)
        valid_replicates = 0
        for _ in range(bootstrap_replicates):
            selected = rng.integers(0, len(loci), len(loci))
            resampled = np.zeros_like(matrix)
            shared_counts = np.zeros_like(matrix)
            differences = np.zeros_like(matrix)
            for locus in selected:
                values = np.array([profile[locus] or "" for profile in calls])
                present = values != ""
                shared = present[:, None] & present[None, :]
                shared_counts += shared
                differences += shared & (values[:, None] != values[None, :])
            valid = bool(np.all(shared_counts / len(loci) >= min_overlap))
            np.divide(differences, shared_counts, out=resampled, where=shared_counts > 0)
            if not valid:
                continue
            valid_replicates += 1
            for group in _clusters(resampled, distance_threshold):
                for i in group:
                    for j in group:
                        pair_support[i, j] += 1
        cohort["valid_bootstrap_replicates"] = valid_replicates
        for group_number, group in enumerate(groups, 1):
            group_id = f"{cohort_id}_group_{group_number}"
            members = [rows[i] for i in group]
            anchors = []
            for kind, level in (("cglin", 7), ("hiercc", "HC1100")):
                labels = [_lineage(r, kind, level) for r in members]
                if labels and labels[0] is not None and all(label == labels[0] for label in labels):
                    anchors.append({"kind": kind, "level": level, "assignment": list(labels[0])})
            support = (
                [pair_support[i, j] / valid_replicates for i in group for j in group if i < j]
                if valid_replicates
                else []
            )
            summary["genetic_groups"].append(
                {
                    "group_id": group_id,
                    "sample_ids": [r["sample_id"] for r in members],
                    "method": "complete linkage",
                    "assignment_anchors": anchors,
                    "minimum_pair_bootstrap_coassignment": min(support) if support else None,
                    "singleton": len(group) == 1,
                    "valid_bootstrap_replicates": valid_replicates,
                }
            )
            dated = [_date_interval(r) for r in members]
            years = sorted({d["year"] for d in dated if d})
            summary["temporal_persistence"].append(
                {
                    "group_id": group_id,
                    "observed_years": years,
                    "observed_year_count": len(years),
                    "year_span": years[-1] - years[0] if years else None,
                    "dated_records": sum(d is not None for d in dated),
                    "total_records": len(members),
                    "interpretation": (
                        "Recurrence observed across collection years; gaps do not establish continuous persistence."
                        if len(years) >= 2
                        else "Persistence across years cannot be assessed from zero or one collection year."
                    ),
                }
            )
            cells = Counter()
            for row, date in zip(members, dated):
                region = row.get("nuts2") or row.get("region") or row.get("country")
                if region and date:
                    cells[(str(region), date["year"])] += 1
            summary["time_place_concentration"].append(
                {
                    "group_id": group_id,
                    "cells": [
                        {"location": place, "year": year, "count": count}
                        for (place, year), count in sorted(cells.items())
                    ],
                    "records_with_time_and_place": sum(cells.values()),
                    "total_records": len(members),
                    "largest_cell_fraction_of_annotated": max(cells.values()) / sum(cells.values())
                    if cells
                    else None,
                    "interpretation": "Descriptive concentration among annotated sampled genomes; no significance or outbreak claim.",
                }
            )
        if len(rows) >= 2:
            names = [r["sample_id"] for r in rows]
            dm = DistanceMatrix(names, [matrix[i, : i + 1].tolist() for i in range(len(rows))])
            tree = DistanceTreeConstructor().nj(dm)
            tree.root_with_outgroup(names[0])
            path = output / f"{cohort_id}_nj.nwk"
            Phylo.write(tree, str(path), "newick")
            cohort["tree_path"] = str(path)
            import matplotlib.pyplot as plt

            if len(rows) <= 100:
                tree_figure = output / f"{cohort_id}_nj.svg"
                draw_profile_tree(tree, rows, summary["sample_labels"], tree_figure)
                cohort["tree_figure"] = str(tree_figure)
            else:
                cohort["warnings"].append(
                    "Static tree figure omitted above 100 tips; complete NJ Newick is available."
                )
            cohort["exploratory_root"] = names[0]
            if any((node.branch_length or 0) < 0 for node in tree.find_clades()):
                cohort["warnings"].append(
                    "NJ contains negative branches; retained and flagged, not biologically interpreted."
                )
            dates, lengths = [], []
            points = []
            for row in rows:
                date = _date_interval(row)
                if date:
                    length = tree.distance(row["sample_id"])
                    points.append({"sample_id": row["sample_id"], **date, "root_to_tip": length})
                    dates.append((date["year_min"] + date["year_max"]) / 2)
                    lengths.append(length)
            slope = correlation = None
            if len(dates) >= 3 and np.ptp(dates) > 0 and np.ptp(lengths) > 0:
                slope = float(np.polyfit(dates, lengths, 1)[0])
                correlation = float(np.corrcoef(dates, lengths)[0, 1])
            summary["root_to_tip"].append(
                {
                    "cohort_id": cohort_id,
                    "root": names[0],
                    "points": points,
                    "midpoint_date_slope": slope,
                    "midpoint_date_pearson_r": correlation,
                    "interpretation": "Exploratory, arbitrary-root NJ diagnostic using date interval midpoints; not a temporal-signal gate.",
                }
            )
            fig, ax = plt.subplots(figsize=(7, 5))
            for point in points:
                midpoint = (point["year_min"] + point["year_max"]) / 2
                ax.errorbar(
                    midpoint,
                    point["root_to_tip"],
                    xerr=[[midpoint - point["year_min"]], [point["year_max"] - midpoint]],
                    fmt="o",
                    alpha=0.65,
                )
            if not points:
                ax.text(
                    0.5,
                    0.5,
                    "No valid collection dates available",
                    transform=ax.transAxes,
                    ha="center",
                )
            ax.set(
                xlabel="Collection year (date interval shown)",
                ylabel="NJ root-to-tip allele-distance fraction",
                title="Exploratory date diagnostic; arbitrary tip root",
            )
            fig.tight_layout()
            diagnostic_path = output / f"{cohort_id}_root_to_tip.svg"
            fig.savefig(diagnostic_path)
            plt.close(fig)
            cohort["root_to_tip_figure"] = str(diagnostic_path)
            network = _location_network(tree, rows)
        else:
            network = {
                "edges": [],
                "tested_roots": [],
                "interpretation": "At least two profiles required.",
            }
            cohort["warnings"].append("Single comparable profile: no NJ tree or date diagnostic.")
        network["cohort_id"] = cohort_id
        summary["location_network"].append(network)
        cohort.update(_plots(output, cohort_id, coords, rows, network))
        cohort["pcoa_csv"] = _csv(
            output / f"{cohort_id}_pcoa.csv",
            [
                {
                    "sample_id": row["sample_id"],
                    "axis_1": float(coords[i, 0]),
                    "axis_2": float(coords[i, 1]),
                }
                for i, row in enumerate(rows)
            ],
            ["sample_id", "axis_1", "axis_2"],
        )
        summary["cohorts"].append(cohort)
    tables = {
        "distances": (
            pairs,
            [
                "sample_a",
                "sample_b",
                "allele_differences",
                "shared_called_loci",
                "scheme_loci",
                "call_overlap",
                "distance",
            ],
        ),
        "nearest_neighbours": (
            nearest,
            [
                "query_id",
                "context_id",
                "status",
                "tied_neighbours",
                "country",
                "region",
                "nuts2",
                "collection_date",
                "year",
                "date_status",
                "cohort_id",
                "query_cohort_id",
                "allele_differences",
                "shared_called_loci",
                "call_overlap",
                "distance",
            ],
        ),
        "genetic_groups": (
            summary["genetic_groups"],
            [
                "group_id",
                "sample_ids",
                "assignment_anchors",
                "minimum_pair_bootstrap_coassignment",
                "singleton",
                "valid_bootstrap_replicates",
            ],
        ),
        "temporal_persistence": (
            summary["temporal_persistence"],
            [
                "group_id",
                "observed_years",
                "observed_year_count",
                "year_span",
                "dated_records",
                "total_records",
                "interpretation",
            ],
        ),
        "time_place_concentration": (
            summary["time_place_concentration"],
            [
                "group_id",
                "cells",
                "records_with_time_and_place",
                "total_records",
                "largest_cell_fraction_of_annotated",
                "interpretation",
            ],
        ),
        "exclusions": (exclusions, ["sample_id", "reason"]),
        "date_exclusions": (date_exclusions, ["sample_id", "collection_date", "reason"]),
    }
    summary["paths"]["pairwise_distances"] = _csv(
        output / "pairwise_distances.csv",
        [
            {
                "sample_id_1": pair["sample_a"],
                "sample_id_2": pair["sample_b"],
                **{
                    key: pair[key]
                    for key in (
                        "distance",
                        "allele_differences",
                        "shared_called_loci",
                        "call_overlap",
                    )
                },
            }
            for pair in pairs
        ],
        [
            "sample_id_1",
            "sample_id_2",
            "distance",
            "allele_differences",
            "shared_called_loci",
            "call_overlap",
        ],
    )
    for name, (rows, fields) in tables.items():
        summary["paths"][name] = _csv(output / f"{name}.csv", rows, fields)
    summary["paths"]["sample_labels"] = _csv(
        output / "sample_labels.csv",
        [{"sample_id": row["sample_id"],
          "display_label": summary["sample_labels"][row["sample_id"]],
          "source_genome_id": row.get("source_genome_id", ""),
          "run_accessions": ";".join(row.get("run_accessions") or []),
          "biosample_accessions": ";".join(row.get("biosample_accessions") or [])}
         for row in records],
        ["sample_id", "display_label", "source_genome_id", "run_accessions", "biosample_accessions"],
    )
    summary["paths"]["summary"] = str(output / "profile_analysis.json")
    (output / "profile_analysis.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary
