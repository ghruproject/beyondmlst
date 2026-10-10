"""Complete large-block cgMLST analysis without silent profile sampling."""

from __future__ import annotations
from collections import Counter
import datetime as dt
import json
from pathlib import Path
import sys
import time
import numpy as np
from Bio import Phylo

from chronoclade.context_refinement import _lineage
from chronoclade.metadata_dates import date_interval
from chronoclade.sample_labels import sample_labels
from chronoclade.profile_analysis import _csv, write_profile_network_figures
from chronoclade.profile_network import build_profile_network
from chronoclade.temporal_diagnostics import date_regression, allele_unit_diagnostic
from chronoclade.report_components.ordination import write_ordination_views
from .scale_distance import prepare_records, write_distances, cohort_matrix
from .scale_tree import IndexedTree, leading_pcoa, rapidnj_tree, write_collapsible_tree


def _groups(condensed, n, cutoff):
    """Compiled complete linkage: every within-group distance meets cutoff."""
    from scipy.cluster.hierarchy import linkage, fcluster

    if n == 1:
        return [[0]], dict(method="SciPy compiled complete linkage", thresholds=[])
    linkage_matrix = linkage(condensed, method="complete")

    def cut(value):
        assignments = fcluster(linkage_matrix, value, criterion="distance")
        members = {}
        for i, group in enumerate(assignments):
            members.setdefault(int(group), []).append(i)
        return sorted(members.values())

    sensitivity = []
    for threshold in sorted(set([cutoff / 2, cutoff, min(1.0, cutoff * 2)])):
        sizes = list(map(len, cut(threshold)))
        sensitivity.append(
            dict(
                distance_fraction=threshold,
                group_count=len(sizes),
                singletons=sizes.count(1),
                largest_group=max(sizes),
                size_distribution=dict(sorted(Counter(sizes).items())),
                selected=threshold == cutoff,
            )
        )
    return cut(cutoff), dict(
        method="compiled complete linkage; every within-group pair meets cutoff",
        thresholds=sensitivity,
        calibration_status="exploratory_not_biologically_validated",
        threshold_basis="explicit descriptive mismatch fraction",
        tie_policy="SciPy deterministic input-order tie resolution; may differ from legacy merge order for exact ties",
    )


def _nearest(records, evidence, cohorts, future_ids):
    membership = {
        i: f"cohort_{number}" for number, indexes in enumerate(cohorts, 1) for i in indexes
    }
    context = np.asarray(
        [i for i, row in enumerate(records) if row.get("role", row.get("origin")) == "context"]
    )
    nearest = []
    for i, row in enumerate(records):
        if row.get("role", row.get("origin")) not in {"input", "local", "query", "focal"}:
            continue
        values = evidence.row(row["sample_id"])
        choices = context[np.isfinite(values[context])] if len(context) else []
        if not len(choices):
            nearest.append(dict(query_id=row["sample_id"], status="no_comparable_context"))
            continue
        raw = evidence.arrays["allele_differences"][i]
        minimum = int(np.min(raw[choices]))
        matches = choices[raw[choices] == minimum]
        for j in matches:
            other = records[j]
            date = date_interval(other)
            nearest.append(
                dict(
                    query_id=row["sample_id"],
                    context_id=other["sample_id"],
                    status="matched",
                    country=other.get("country") or "Unknown",
                    region=other.get("region") or "",
                    nuts2=other.get("nuts2") or "",
                    collection_date=other.get("collection_date")
                    or other.get("collection_year")
                    or "",
                    year=date["year"] if date else None,
                    date_status="future_collection_date"
                    if other["sample_id"] in future_ids
                    else ("valid" if date else "missing_or_invalid"),
                    cohort_id=membership.get(int(j)),
                    query_cohort_id=membership.get(i),
                    tied_neighbours=len(matches),
                    **{
                        k: v
                        for k, v in evidence.get(row["sample_id"], other["sample_id"]).items()
                        if k not in {"sample_a", "sample_b"}
                    },
                )
            )
    return nearest


def _date_diagnostic(
    tree, rows, labels, root, *, factor=None, fixed_callable=False, catalogue_complete=True
):
    diagnostic = date_regression(
        [
            dict(
                sample_id=row["sample_id"],
                collection_date=row.get("collection_date") or row.get("collection_year"),
                date_start=row.get("date_start"),
                date_end=row.get("date_end"),
                date_precision=row.get("date_precision"),
                country=row.get("country"),
                label=labels[row["sample_id"]],
                role="input"
                if row.get("role", row.get("origin")) in {"input", "local", "query", "focal"}
                else "context",
                root_to_tip=float(tree.distance(row["sample_id"], root)),
            )
            for row in rows
        ],
        distance_field="root_to_tip",
        distance_units="cgMLST mismatch fraction",
        label="Full profile NJ root-to-tip",
    )
    if factor is not None:
        diagnostic = allele_unit_diagnostic(
            diagnostic, factor, fixed_callable=fixed_callable, catalogue_complete=catalogue_complete
        )
    return dict(cohort_root=root, **diagnostic)


def _diagnostic_figure(diagnostic, path):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 5))
    points = diagnostic["points"]
    for point in points:
        midpoint = (point["year_min"] + point["year_max"]) / 2
        ax.errorbar(
            midpoint,
            point["root_to_tip"],
            xerr=[[midpoint - point["year_min"]], [point["year_max"] - midpoint]],
            fmt="o",
            alpha=0.4,
            markersize=3,
        )
    if diagnostic.get("slope") is not None and points:
        years = np.array([min(p["year_min"] for p in points), max(p["year_max"] for p in points)])
        ax.plot(years, diagnostic["slope"] * years + diagnostic["intercept"], color="#444")
    if not points:
        ax.text(
            0.5, 0.5, "No valid collection dates available", transform=ax.transAxes, ha="center"
        )
    ax.set(
        xlabel="Collection year (date interval shown)",
        ylabel=diagnostic.get("distance_units", "cgMLST allele differences"),
        title="Exploratory full NJ date diagnostic · declared input tip root",
    )
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return str(path)


def analyse_profiles_scale(
    records,
    *,
    output,
    seed=42,
    min_overlap=0.9,
    bootstrap_replicates=0,
    distance_threshold=0.02,
    tree_limit=None,
    batch_size=128,
    rapidnj_executable=None,
    rapidnj_memory_mb=2048,
):
    """Retain every supplied sample and every compatible pair on disk.

    Unavailable profiles remain in metadata/exclusions; complete-comparability
    cohorts contain every eligible profile and no missing/imputed distances.
    Locus-bootstrap groups are currently unavailable on this backend. Requesting
    them raises before analysis, rather than ignoring an explicit request.
    """
    start = time.perf_counter()
    if not records:
        raise ValueError("Large profile analysis needs at least one record")
    if bootstrap_replicates != 0:
        raise ValueError(
            "Large cgMLST backend currently requires explicit bootstrap_replicates=0; locus-bootstrap grouping is unavailable. No requested replicates were silently omitted."
        )
    if not 0 < min_overlap <= 1 or not 0 <= distance_threshold <= 1:
        raise ValueError("Overlap must be in (0,1] and distance threshold in [0,1]")
    if not isinstance(batch_size, int) or batch_size < 1:
        raise ValueError("Distance batch_size must be a positive integer")
    ids = [row.get("sample_id") for row in records]
    if any(not isinstance(i, str) or not i for i in ids) or len(set(ids)) != len(ids):
        raise ValueError("Sample identifiers must be nonempty and unique")
    # Biopython uses recursive Newick serialization and tree traversal. Its
    # deepest valid binary tree can have n tips along one spine.
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 4 * len(records) + 100))
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    records = prepare_records(records)
    labels = sample_labels(records)
    evidence, cohorts, exclusions = write_distances(
        records, output / "binary_distances", min_overlap=min_overlap, batch_size=batch_size
    )
    distance_seconds = time.perf_counter() - start
    date_exclusions = []
    for row in records:
        interval = date_interval(row, allow_future=True)
        if interval and dt.date.fromisoformat(interval["start"]) > dt.date.today():
            date_exclusions.append(
                dict(
                    sample_id=row["sample_id"],
                    collection_date=row.get("collection_date") or row.get("collection_year"),
                    reason="future_collection_date",
                )
            )
    future_ids = {row["sample_id"] for row in date_exclusions}
    nearest = _nearest(records, evidence, cohorts, future_ids)
    summary = dict(
        schema_version=1,
        records_count=len(records),
        sample_labels=labels,
        seed=seed,
        min_overlap=min_overlap,
        distance_threshold=distance_threshold,
        bootstrap_replicates=0,
        bootstrap_status="not_requested; large backend does not implement locus-bootstrap group support",
        coverage=dict(
            available_profiles=sum(map(len, cohorts)),
            comparable_pairs=evidence.scope["comparable_pairs"],
            unavailable_pairs_by_reason=evidence.scope["unavailable_pairs_by_reason"],
        ),
        exclusions=exclusions,
        date_exclusions=date_exclusions,
        cohorts=[],
        nearest_neighbours=nearest,
        genetic_groups=[],
        temporal_persistence=[],
        time_place_concentration=[],
        root_to_tip=[],
        location_network=[],
        paths={"distance_evidence": str(evidence.manifest_path)},
        distance_evidence=evidence.descriptor,
        scale_backend=dict(
            distance="batched compiled SciPy categorical Hamming; exact called counts",
            storage="memory-mapped square matrices and cohort condensed vectors",
            profile_count=len(records),
            silent_subsetting=False,
            batch_size=batch_size,
            distance_seconds=distance_seconds,
        ),
        limitations=[
            "Distances describe the supplied catalogue; all input IDs are retained.",
            "Missing calls are excluded, never imputed as matching loci.",
            "Complete-comparability cohorts are deterministic partitions, not genetic clusters; cross-cohort comparable pairs remain in binary evidence.",
            "Leading PCoA axes use an iterative solver; total positive inertia and the complete negative spectrum are uncomputed.",
            "Genetic thresholds are exploratory and are not validated outbreak or clone thresholds.",
            "Metadata compositions describe sampled genomes, not population prevalence.",
            "NJ temporal diagnostics and country changes do not establish transmission or assess recombination-corrected SNP temporal signal.",
        ],
    )
    for number, indexes in enumerate(cohorts, 1):
        cohort_start = time.perf_counter()
        cohort_id = f"cohort_{number}"
        rows = [records[i] for i in indexes]
        names = [row["sample_id"] for row in rows]
        matrix, condensed = cohort_matrix(evidence, indexes, output / f"{cohort_id}_matrix.npy")
        coords, ordination = leading_pcoa(matrix, seed=seed, batch_size=batch_size)
        group_lists, sensitivity = _groups(condensed, len(rows), distance_threshold)
        cohort = dict(
            cohort_id=cohort_id,
            sample_ids=names,
            tree_units="fraction of mismatching jointly called cgMLST loci",
            warnings=[],
            valid_bootstrap_replicates=0,
            group_sensitivity=sensitivity,
            **ordination,
        )
        if any(row.get("cgmlst_locus_universe_complete") is False for row in rows):
            cohort["coverage_denominator"] = (
                "observed export locus union; canonical completeness unknown"
            )
            cohort["warnings"].append(
                "Full canonical scheme coverage is unknown; denominators use the observed export locus union."
            )
        else:
            cohort["coverage_denominator"] = "declared locus universe"
        groups_by_id = {}
        for group_number, group in enumerate(group_lists, 1):
            group_id = f"{cohort_id}_group_{group_number}"
            members = [rows[i] for i in group]
            anchors = []
            for kind, level in (("cglin", 7), ("hiercc", "HC1100")):
                assignments = [_lineage(row, kind, level) for row in members]
                if (
                    assignments
                    and assignments[0] is not None
                    and all(value == assignments[0] for value in assignments)
                ):
                    anchors.append(dict(kind=kind, level=level, assignment=list(assignments[0])))
            summary["genetic_groups"].append(
                dict(
                    group_id=group_id,
                    sample_ids=[row["sample_id"] for row in members],
                    method="compiled complete linkage",
                    assignment_anchors=anchors,
                    minimum_pair_bootstrap_coassignment=None,
                    singleton=len(members) == 1,
                    valid_bootstrap_replicates=0,
                )
            )
            for row in members:
                groups_by_id[row["sample_id"]] = group_id
            dated = [date_interval(row) for row in members]
            years = sorted({date["year"] for date in dated if date})
            summary["temporal_persistence"].append(
                dict(
                    group_id=group_id,
                    observed_years=years,
                    observed_year_count=len(years),
                    year_span=years[-1] - years[0] if years else None,
                    dated_records=sum(date is not None for date in dated),
                    total_records=len(members),
                    interpretation="Recurrence across observed collection years does not establish continuous persistence."
                    if len(years) > 1
                    else "Persistence cannot be assessed from zero or one collection year.",
                )
            )
            cells = Counter(
                (str(row.get("nuts2") or row.get("region") or row.get("country")), date["year"])
                for row, date in zip(members, dated)
                if date and (row.get("nuts2") or row.get("region") or row.get("country"))
            )
            summary["time_place_concentration"].append(
                dict(
                    group_id=group_id,
                    cells=[
                        dict(location=place, year=year, count=count)
                        for (place, year), count in sorted(cells.items())
                    ],
                    records_with_time_and_place=sum(cells.values()),
                    total_records=len(members),
                    largest_cell_fraction_of_annotated=max(cells.values()) / sum(cells.values())
                    if cells
                    else None,
                    interpretation="Descriptive concentration among annotated sampled genomes; no significance or outbreak claim.",
                )
            )
        if len(rows) > 1:
            tree, nj_audit = rapidnj_tree(
                matrix,
                names,
                output / f"{cohort_id}_nj.nwk",
                executable=rapidnj_executable,
                memory_mb=rapidnj_memory_mb,
            )
            root = next(
                (
                    row["sample_id"]
                    for row in rows
                    if row.get("role", row.get("origin")) in {"input", "local", "query", "focal"}
                ),
                names[0],
            )
            tree.root_with_outgroup(root)
            Phylo.write(tree, output / f"{cohort_id}_nj.nwk", "newick")
            cohort.update(
                tree_path=str(output / f"{cohort_id}_nj.nwk"),
                nj_backend=nj_audit,
                exploratory_root=root,
                root_policy="lexicographically first input tip; first eligible tip when cohort has no input",
            )
            locus_count = int(evidence.arrays["scheme_loci"][indexes[0], indexes[0]])
            fixed_calls = all(
                np.all(evidence.arrays["shared_called_loci"][index, indexes] == locus_count)
                for index in indexes
            )
            units = (
                "cgMLST allele differences"
                if fixed_calls
                else "scheme-equivalent cgMLST allele differences"
            )
            cohort["root_to_tip_units"] = units
            cohort["root_to_tip_scale_loci"] = locus_count
            catalogue_complete = not any(
                row.get("cgmlst_locus_universe_complete") is False for row in rows
            )
            diagnostic = _date_diagnostic(
                tree,
                rows,
                labels,
                root,
                factor=locus_count,
                fixed_callable=fixed_calls,
                catalogue_complete=catalogue_complete,
            )
            summary["root_to_tip"].append(
                dict(
                    cohort_id=cohort_id,
                    root=root,
                    root_policy=cohort["root_policy"],
                    scale_loci=locus_count,
                    callable_denominators="fixed full locus count"
                    if fixed_calls
                    else "variable; mismatch fraction scaled by declared scheme locus count, not raw observed changes",
                    **diagnostic,
                    midpoint_date_slope=diagnostic["slope"],
                    midpoint_date_pearson_r=diagnostic["pearson_r"],
                )
            )
            cohort["root_to_tip_figure"] = _diagnostic_figure(
                diagnostic, output / f"{cohort_id}_root_to_tip.svg"
            )
            if any((node.branch_length or 0) < 0 for node in tree.find_clades()):
                cohort["warnings"].append("Negative NJ branch lengths are retained and flagged.")
        else:
            tree = IndexedTree(root=Phylo.BaseTree.Clade(name=names[0]))
            tree.reindex()
            cohort["warnings"].append("Singleton comparable cohort: no NJ or date regression.")
        network_start = time.perf_counter()
        network = build_profile_network(tree, rows, nearest)
        network["cohort_id"] = cohort_id
        write_profile_network_figures(output, cohort_id, network, records=())
        network["country_tree_figure"] = write_collapsible_tree(
            network, rows, labels, output / f"{cohort_id}_tree.html"
        )
        cohort["tree_figure"] = network["country_tree_figure"]
        cohort["tree_display_sample_ids"] = names
        cohort["tree_display_limit"] = None
        if len(rows) > 1:
            roots = network.get("tested_roots", [])
            summary["root_to_tip"][-1]["root_sensitivity"] = [
                dict(
                    root=alternative,
                    slope=detail["slope"],
                    pearson_r=detail["pearson_r"],
                    dated_count=detail.get("dated_count"),
                )
                for alternative in roots
                for detail in [
                    _date_diagnostic(
                        tree,
                        rows,
                        labels,
                        alternative,
                        factor=locus_count,
                        fixed_callable=fixed_calls,
                        catalogue_complete=catalogue_complete,
                    )
                ]
            ]
        network["scale_runtime_seconds"] = time.perf_counter() - network_start
        Path(network["audit_path"]).write_text(json.dumps(network, indent=2) + "\n")
        summary["location_network"].append(network)
        for label, key in (
            ("network_audit", "audit_path"),
            ("nearest_country_connections", "nearest_country_connections_path"),
            ("country_transitions", "transition_counts_path"),
            ("network_metrics", "network_metrics_path"),
        ):
            summary["paths"][f"{cohort_id}_{label}"] = network[key]
        cohort.update(
            write_ordination_views(
                output,
                cohort_id,
                coords,
                rows,
                title="Full cgMLST distance ordination",
                axis_labels=("PCoA axis 1", "PCoA axis 2"),
                groups=groups_by_id,
            )
        )
        cohort["pcoa_csv"] = _csv(
            output / f"{cohort_id}_pcoa.csv",
            [
                dict(
                    sample_id=row["sample_id"],
                    axis_1=float(coords[i, 0]),
                    axis_2=float(coords[i, 1]),
                )
                for i, row in enumerate(rows)
            ],
            ["sample_id", "axis_1", "axis_2"],
        )
        from chronoclade.report_components.geography import write_country_figure

        cohort["country_figure"] = write_country_figure(output, cohort_id, rows)
        cohort["scale_runtime_seconds"] = time.perf_counter() - cohort_start
        summary["cohorts"].append(cohort)
    for key in (
        "nearest_neighbours",
        "genetic_groups",
        "temporal_persistence",
        "time_place_concentration",
        "exclusions",
        "date_exclusions",
    ):
        fields = list(dict.fromkeys(field for item in summary[key] for field in item))
        summary["paths"][key] = _csv(output / f"{key}.csv", summary[key], fields or ["sample_id"])
    summary["paths"]["sample_labels"] = _csv(
        output / "sample_labels.csv",
        [
            dict(sample_id=row["sample_id"], display_label=labels[row["sample_id"]])
            for row in records
        ],
        ["sample_id", "display_label"],
    )
    summary["paths"]["summary"] = str(output / "profile_analysis.json")
    summary["scale_backend"]["total_seconds"] = time.perf_counter() - start
    (output / "profile_analysis.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary
