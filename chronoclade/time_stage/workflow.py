"""Date a saved corrected tree without constructing or changing phylogenies."""

from pathlib import Path
import shutil

from Bio import Phylo

from chronoclade.artifacts import write_json
from chronoclade.errors import WorkflowError
from chronoclade.selection_manifest import _resolve
from chronoclade.stage_artifacts import (
    load_tree_result,
    load_result,
    reference,
    fingerprint,
    tools,
    require_tools,
    prepare_job,
    publish,
    fail_job,
    read_dated_tree,
)
from chronoclade.report_components.stage import write_report


def run_time(tree_manifest, *, output, force=False, allow_unsupported=False):
    tree_path = Path(tree_manifest).expanduser().resolve()
    tree = load_tree_result(tree_path)
    output = Path(output).expanduser().resolve()
    tool_evidence = tools(["treetime"])
    digest = fingerprint(
        {
            "tree": reference(output, tree_path)["sha256"],
            "tools": tool_evidence,
            "allow_unsupported": allow_unsupported,
        }
    )
    if (output / "time.json").is_file() and not force:
        import json

        previous = json.loads((output / "time.json").read_text())
        if previous.get("fingerprint") == digest:
            return load_result(output / "time.json", "time")
    output, job = prepare_job(output, digest, "time", force)
    write_json(
        job / "job.json",
        {
            "stage": "time",
            "status": "running",
            "selection_id": tree["selection_id"],
            "source": {"tree": reference(job, tree_path)},
            "fingerprint": digest,
        },
    )
    try:
        frozen_tree = tree_path.with_name(
            "tree-" + reference(output, tree_path)["sha256"] + ".json"
        )
        if not frozen_tree.exists():
            write_json(frozen_tree, tree)
        load_tree_result(frozen_tree)
        rows = _read_rows(_resolve(tree_path.parent, tree["metadata"], "metadata"))
        assessment = tree["temporal_assessment"]
        result = {
            "fingerprint": digest,
            "source": {"tree": reference(output, frozen_tree)},
            "input_tree_sha256": reference(output, tree_path)["sha256"],
            "selection_id": tree["selection_id"],
            "selected_sample_ids": tree["selected_sample_ids"],
            "temporal_assessment": assessment,
            "parameters": {"allow_unsupported": allow_unsupported},
            "tools": tool_evidence,
            "dating_status": "unsupported",
            "dated_tree": None,
            "dating_override": bool(allow_unsupported and not assessment["supported"]),
            "branch_units": None,
            "clock_confidence": None,
        }
        figures = []
        network = None
        if assessment["supported"] or allow_unsupported:
            if len(tree["dated_sample_ids"]) < 3:
                raise WorkflowError(
                    "Dating requires at least three saved eligible dated samples, even with an override"
                )
            require_tools(["treetime"])
            from chronoclade.lineage import LineageFiles, _write_csv, _run_dated_tree
            from chronoclade.temporal_report import (
                write_timetree_confidence,
                write_timetree_figures,
            )

            files = LineageFiles.in_directory(job)
            files.clock_dir.mkdir(exist_ok=True)
            rooted = _resolve(tree_path.parent, tree["rooted_tree"], "rooted_tree")
            shutil.copyfile(rooted, files.clock_dir / "rerooted.newick")
            dated = set(tree["dated_sample_ids"])
            _write_csv(
                job / "clock_dates.csv",
                ["sample_id", "collection_date"],
                (
                    {
                        "sample_id": row["sample_id"],
                        "collection_date": row.get("collection_date") or "",
                    }
                    for row in rows
                    if row["sample_id"] in dated
                ),
            )
            _run_dated_tree(files, tree["complete_alignment_sites"], force)
            tips = [tip.name for tip in Phylo.read(files.time_tree, "nexus").get_terminals()]
            if len(tips) != len(tree["selected_sample_ids"]) or set(tips) != set(
                tree["selected_sample_ids"]
            ):
                raise WorkflowError("Dated tree does not retain exact selected sample identifiers")
            from chronoclade.location_network.reconstruction import build_location_network
            from chronoclade.location_network.tree import draw_country_tree

            biological = read_dated_tree(files.time_tree)
            network = build_location_network(biological, rows, tree_basis="sequence")
            network.update(branch_units="years", root_policy="saved corrected tree root retained")
            write_json(job / "location_network.json", network)
            draw_country_tree(
                network,
                rows,
                job / "dated_country_tree.svg",
                title="Dated biological tree location history",
                distance_label="Years",
                dating_note="TreeTime time-scaled biological tree",
            )
            confidence = write_timetree_confidence(job)
            write_timetree_figures(job, confidence)
            result.update(
                dating_status="dated",
                branch_units="years",
                dated_tree=reference(output, files.time_tree),
                clock_confidence=confidence,
                root_policy="saved tree root retained",
            )
            figures = [
                ("dated_country_tree.svg", "Dated-tree location history; years"),
                (
                    "timetree_with_confidence.svg",
                    "TreeTime node-date uncertainty; 90% max-posterior regions",
                ),
            ]
        write_json(job / "assessment.json", assessment)
        write_json(job / "selected_metadata.json", rows)
        note = "<p>Dating override explicitly requested.</p>" if result["dating_override"] else ""
        result["report"] = reference(
            output,
            write_report(
                job / "report.html",
                title="Time-scaled phylogeny" if result["dated_tree"] else "Dating assessment",
                result=result,
                rows=rows,
                network=network,
                figures=figures,
                extra=note,
            ),
        )
        return publish(output, job, "time", result)
    except Exception as error:
        fail_job(job, "time", error)
        raise


def _read_rows(path):
    import json

    rows = json.loads(Path(path).read_text())
    if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
        raise WorkflowError("Selected metadata must be a table of objects")
    return rows
