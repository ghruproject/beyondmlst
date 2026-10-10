"""Compare immutable independent selections on explicitly shared anchors."""

import csv
import json
from itertools import combinations
from pathlib import Path
from html import escape


from chronoclade.artifacts import write_json
from chronoclade.errors import WorkflowError
from chronoclade.selection_manifest import load_selection_ensemble, _resolve, _read
from chronoclade.stage_artifacts import (
    load_result,
    load_tree_result,
    reference,
    header,
    read_dated_tree,
)
from chronoclade.report_components.stage import write_report


def _spread(values):
    values = [float(value) for value in values if value is not None]
    return {
        "assessable_runs": len(values),
        "minimum": min(values) if values else None,
        "maximum": max(values) if values else None,
        "range": max(values) - min(values) if values else None,
    }


def _anchor_dates(path, result, anchors):
    if result["dated_tree"] is None:
        return None
    dated_path = _resolve(path.parent, result["dated_tree"], "dated_tree")
    biological = read_dated_tree(dated_path)
    mrca = biological.common_ancestor(anchors)
    date_files = [spec for name, spec in result["artifacts"].items() if name == "node_dates.csv"]
    if not date_files:
        return None
    with _resolve(path.parent, date_files[0], "node_dates").open() as stream:
        by_node = {row["node"]: row for row in csv.DictReader(stream)}
    row = by_node.get(str(mrca.name))
    return float(row["numeric_date"]) if row else None


def _clusters(tree_path, tree, ids, cutoff):
    spec = tree["artifacts"].get("snp_distances.json")
    if not spec:
        raise WorkflowError("Tree has no corrected SNP-distance evidence")
    distances = json.loads(_resolve(tree_path.parent, spec, "snp_distances").read_text())
    parent = {ident: ident for ident in tree["selected_sample_ids"]}

    def find(ident):
        while parent[ident] != ident:
            ident = parent[ident]
        return ident

    for row in distances:
        if row["callable_sites"] > 0 and row["clonal_snps"] <= cutoff:
            parent[find(row["sample_2"])] = find(row["sample_1"])
    callability = {
        frozenset((row["sample_1"], row["sample_2"])): row["callable_sites"] for row in distances
    }
    return {
        "|".join((left, right)): (find(left) == find(right))
        if callability.get(frozenset((left, right)), 0) > 0
        else None
        for left, right in combinations(ids, 2)
    }


def compare_time_runs(ensemble_manifest, run_manifests, *, output, cluster_snp_cutoff=10):
    if cluster_snp_cutoff < 0:
        raise WorkflowError("Cluster SNP cutoff must be nonnegative")
    ensemble_path = Path(ensemble_manifest).resolve()
    ensemble = load_selection_ensemble(ensemble_path)
    output = Path(output).resolve()
    if (output / "comparison.json").exists():
        raise WorkflowError("Comparison output already exists; use an independent directory")
    output.mkdir(parents=True, exist_ok=True)
    available = {}
    for raw_path in run_manifests:
        path = Path(raw_path).resolve()
        raw = _read(path)
        if raw.get("stage") == "time" and raw.get("status") in {"failed", "pending", "running"}:
            result = raw
            if not isinstance(result.get("selection_id"), str) or not isinstance(
                result.get("source"), dict
            ):
                raise WorkflowError(
                    "Failed/pending time job requires source tree and selection identity"
                )
        else:
            result = load_result(path, "time")
        ident = result["selection_id"]
        if ident in available:
            raise WorkflowError("Duplicate time runs for one selection")
        available[ident] = (path, result)
    expected = {entry["selection_id"] for entry in ensemble["selections"]}
    if set(available) - expected:
        raise WorkflowError("Time run does not belong to declared selection ensemble")
    runs = []
    anchors, targets = ensemble["common_anchor_ids"], ensemble["target_ids"]
    for selection_entry in ensemble["selections"]:
        ident = selection_entry["selection_id"]
        if ident not in available:
            runs.append(
                {
                    "selection_id": ident,
                    "status": "unavailable",
                    "reason": "No complete time manifest supplied; run may be pending or failed",
                }
            )
            continue
        path, result = available[ident]
        tree_path = _resolve(path.parent, result["source"]["tree"], "source.tree", external=True)
        tree = load_tree_result(tree_path)
        if tree["source"]["selection"]["sha256"] != selection_entry["sha256"]:
            raise WorkflowError("Time tree references a different selection checksum")
        if result.get("stage") == "time":
            runs.append(
                {
                    "selection_id": ident,
                    "status": result["status"],
                    "time_job": reference(output, path),
                    "reason": result.get("reason") or "Time job has not completed",
                }
            )
            continue
        confidence = result.get("clock_confidence") or {}
        runs.append(
            {
                "selection_id": ident,
                "status": result["dating_status"],
                "time": reference(output, path),
                "reason": result["temporal_assessment"]["reason"],
                "rate": confidence.get("rate"),
                "root_date": confidence.get("root_numeric_date"),
                "common_anchor_mrca_date": _anchor_dates(path, result, anchors),
                "target_clusters": _clusters(tree_path, tree, targets, cluster_snp_cutoff),
                "dating_override": result["dating_override"],
            }
        )
    assessable = [run for run in runs if run.get("status") == "dated"]
    pair_labels = sorted({pair for run in runs for pair in run.get("target_clusters", {})})
    agreements = [
        {
            "target_pair": pair,
            "assessable_runs": sum(
                run.get("target_clusters", {}).get(pair) is not None for run in runs
            ),
            "co_clustered_runs": sum(
                run.get("target_clusters", {}).get(pair) is True for run in runs
            ),
        }
        for pair in pair_labels
    ]
    result = {
        **header("time-comparison"),
        "source": {"ensemble": reference(output, ensemble_path)},
        "common_anchor_ids": anchors,
        "target_ids": targets,
        "requested_runs": len(runs),
        "dated_runs": len(assessable),
        "runs": runs,
        "sensitivity": {
            field: _spread([run.get(field) for run in assessable])
            for field in ("rate", "root_date", "common_anchor_mrca_date")
        },
        "cluster_policy": {
            "method": "single-linkage corrected SNP distances",
            "snp_cutoff": cluster_snp_cutoff,
            "calibration": "user-declared descriptive sensitivity threshold; no universal biological cutoff",
        },
        "target_pair_agreement": agreements,
        "interpretation": "Between-selection spread measures selection sensitivity, not confidence intervals. Shared-anchor MRCA dates are comparable; whole-selection roots may represent different ancestors. Cluster comparison uses exact target pairs, never run-local cluster numbers.",
    }
    write_json(output / "comparison.json", result)

    def table(columns, items):
        heads = "".join("<th>" + escape(label) + "</th>" for _, label in columns)
        body = "".join(
            "<tr>"
            + "".join(
                "<td>"
                + escape(str(item.get(key) if item.get(key) is not None else "Unavailable"))
                + "</td>"
                for key, _ in columns
            )
            + "</tr>"
            for item in items
        )
        return (
            '<div class="table-scroll"><table><thead><tr>'
            + heads
            + "</tr></thead><tbody>"
            + body
            + "</tbody></table></div>"
        )

    extra = "<h3>Independent runs</h3>" + table(
        [
            ("selection_id", "Selection"),
            ("status", "Dating status"),
            ("rate", "Rate (substitutions/site/year)"),
            ("common_anchor_mrca_date", "Shared-anchor MRCA year"),
            ("root_date", "Whole-selection root year"),
            ("dating_override", "Override"),
            ("reason", "Assessment / failure"),
        ],
        runs,
    )
    extra += "<h3>Selection sensitivity</h3>" + table(
        [
            ("metric", "Metric"),
            ("assessable_runs", "Assessable runs"),
            ("minimum", "Minimum"),
            ("maximum", "Maximum"),
            ("range", "Range"),
        ],
        [{"metric": field, **spread} for field, spread in result["sensitivity"].items()],
    )
    extra += "<h3>Shared-target clustering</h3>" + table(
        [
            ("target_pair", "Exact target pair"),
            ("assessable_runs", "Assessable runs"),
            ("co_clustered_runs", "Co-clustered runs"),
        ],
        agreements,
    )
    extra += "<p>" + escape(result["interpretation"]) + "</p>"
    write_report(
        output / "report.html",
        title="Selection sensitivity",
        overview=f"{len(anchors)} shared anchor genomes · {len(runs)} requested selections",
        result={
            "selected_sample_ids": anchors,
            "temporal_assessment": {
                "code": "sensitivity",
                "supported": False,
                "reason": f"{len(assessable)} of {len(runs)} requested runs produced dates",
            },
        },
        extra=extra,
    )
    return result
