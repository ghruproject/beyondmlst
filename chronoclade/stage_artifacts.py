"""Validated, checksum-bound tree/time result contracts and job publication."""

import hashlib
import json
import os
from pathlib import Path
import shutil

from chronoclade import __version__
from chronoclade.artifacts import file_sha256, write_json
from chronoclade.errors import WorkflowError
from chronoclade.selection_manifest import _resolve, _read


def reference(base, path):
    path = Path(path).resolve()
    return {"path": Path(os.path.relpath(path, base)).as_posix(), "sha256": file_sha256(path)}


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def tools(names):
    result = {}
    for name in names:
        executable = shutil.which(name)
        evidence = {
            "path": executable,
            "sha256": file_sha256(Path(executable)) if executable else None,
            "package_version": None,
            "version_status": "unavailable",
        }
        if executable:
            records = Path(executable).parent.parent / "conda-meta"
            for record_path in records.glob("*.json"):
                record = _read(record_path)
                if "bin/" + name in record.get("files", []):
                    evidence.update(
                        package_version=record.get("version"),
                        package_name=record.get("name"),
                        package_build=record.get("build"),
                        version_status="conda_manifest",
                    )
                    break
        result[name] = evidence
    return result


def require_tools(names):
    missing = [name for name in names if not shutil.which(name)]
    if missing:
        raise WorkflowError(
            "Missing stage tools: "
            + ", ".join(missing)
            + ". Run through the configured Pixi environment."
        )


def header(stage):
    return {
        "schema": f"chronoclade.{stage}",
        "schema_version": 1,
        "software": {"name": "chronoclade", "version": __version__},
        "status": "complete",
    }


def artifacts(root, job):
    return {
        p.relative_to(job).as_posix(): reference(root, p)
        for p in sorted(job.rglob("*"))
        if p.is_file() and p.name != "job.json"
    }


def load_result(path, stage):
    path = Path(path).expanduser().resolve()
    result = _read(path)
    if result.get("schema") != f"chronoclade.{stage}" or result.get("schema_version") != 1:
        raise WorkflowError(f"Expected versioned chronoclade.{stage} result")
    if result.get("status") != "complete" or not result.get("software", {}).get("version"):
        raise WorkflowError("Stage result must be complete with a software version")
    ids = result.get("selected_sample_ids")
    if (
        not isinstance(ids, list)
        or not ids
        or any(not isinstance(i, str) or not i for i in ids)
        or len(set(ids)) != len(ids)
    ):
        raise WorkflowError("Stage result requires unique exact selected_sample_ids")
    source = result.get("source", {})
    if not source:
        raise WorkflowError("Stage result requires source manifest references")
    for name, spec in source.items():
        _resolve(path.parent, spec, f"source.{name}", external=True)
    entries = result.get("artifacts")
    if not isinstance(entries, dict) or not entries:
        raise WorkflowError("Stage result requires hashed artifacts")
    for name, spec in entries.items():
        _resolve(path.parent, spec, f"artifacts.{name}")
    assessment = result.get("temporal_assessment")
    if not isinstance(assessment, dict) or type(assessment.get("supported")) is not bool:
        raise WorkflowError("Stage result requires an explicit temporal assessment")
    if stage == "time":
        from Bio import Phylo

        tree_path = _resolve(
            path.parent, result["source"].get("tree"), "source.tree", external=True
        )
        tree = load_tree_result(tree_path)
        if ids != tree["selected_sample_ids"] or result.get("selection_id") != tree["selection_id"]:
            raise WorkflowError("Time identifiers disagree with saved tree")
        if result["temporal_assessment"] != tree["temporal_assessment"]:
            raise WorkflowError("Time assessment differs from the saved tree assessment")
        override = result.get("parameters", {}).get("allow_unsupported")
        if type(override) is not bool or result.get("dating_override") is not (
            override and not assessment["supported"]
        ):
            raise WorkflowError("Dating override does not match declared parameters and gate")
        if result.get("dating_status") == "dated":
            if not assessment["supported"] and not override:
                raise WorkflowError("Dated tree bypasses unsupported temporal gate")
            dated_path = _resolve(path.parent, result.get("dated_tree"), "dated_tree")
            tips = [tip.name for tip in Phylo.read(dated_path, "nexus").get_terminals()]
            if len(tips) != len(ids) or set(tips) != set(ids):
                raise WorkflowError("Dated tree does not preserve exact selected identifiers")
        elif result.get("dating_status") != "unsupported" or result.get("dated_tree") is not None:
            raise WorkflowError("Time must be dated or unsupported without a dated tree")
    return result


def load_tree_result(path):
    from Bio import Phylo
    from chronoclade.evidence import read_alignment
    from chronoclade.selection_manifest import load_selection_manifest

    result = load_result(path, "tree")
    base = Path(path).resolve().parent
    selection = load_selection_manifest(
        _resolve(base, result["source"]["selection"], "source.selection", external=True)
    )
    if result["selected_sample_ids"] != selection["selected_sample_ids"]:
        raise WorkflowError("Tree selected identifiers disagree with the exact selection")
    ids = result["selected_sample_ids"]
    for name in ("corrected_tree", "raw_tree"):
        tips = [
            tip.name
            for tip in Phylo.read(_resolve(base, result[name], name), "newick").get_terminals()
        ]
        if len(tips) != len(ids) or set(tips) != set(ids):
            raise WorkflowError(f"{name} does not contain exact selected identifiers")
    sequences = read_alignment(_resolve(base, result["alignment"], "alignment"))
    if set(sequences) != set(ids):
        raise WorkflowError("Corrected alignment does not contain exact selected identifiers")
    import json

    metadata = json.loads(_resolve(base, result["metadata"], "metadata").read_text())
    if (
        not isinstance(metadata, list)
        or [row.get("sample_id") for row in metadata if isinstance(row, dict)] != ids
    ):
        raise WorkflowError(
            "Selected metadata does not preserve exact ordered selected identifiers"
        )
    rooted = Phylo.read(_resolve(base, result["rooted_tree"], "rooted_tree"), "newick")
    tips = [tip.name for tip in rooted.get_terminals()]
    if len(tips) != len(ids) or set(tips) != set(ids):
        raise WorkflowError("Rooted tree does not preserve exact selected identifiers")
    sites = sum(
        all(base in "ACGT" for base in column) for column in zip(*sequences.values(), strict=True)
    )
    if result.get("complete_alignment_sites") != sites:
        raise WorkflowError("Saved callable-site count disagrees with alignment")
    assessment_path = _resolve(base, result["assessment"], "assessment")
    if _read(assessment_path) != result["temporal_assessment"]:
        raise WorkflowError("Temporal assessment artifact disagrees with manifest")
    return result


def prepare_job(output, digest, stage, force):
    output = Path(output).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    # A forced replacement uses a new namespace, leaving the last valid result intact.
    import uuid

    suffix = "-" + uuid.uuid4().hex[:8] if force else ""
    job = output / "jobs" / (digest + suffix)
    job.mkdir(parents=True, exist_ok=True)
    write_json(job / "job.json", {"stage": stage, "status": "running", "fingerprint": digest})
    return output, job


def publish(output, job, stage, result):
    result["artifacts"] = artifacts(output, job)
    result.update(header(stage))
    write_json(
        job / "job.json",
        {"stage": stage, "status": "complete", "fingerprint": result["fingerprint"]},
    )
    write_json(output / f"{stage}.json", result)
    return result


def fail_job(job, stage, error):
    previous = _read(job / "job.json") if (job / "job.json").is_file() else {}
    write_json(
        job / "job.json", {**previous, "stage": stage, "status": "failed", "reason": str(error)}
    )


def read_dated_tree(path):
    """Decode TreeTime NEXUS node labels retained by Bio.Nexus as confidence."""
    from Bio import Phylo

    tree = Phylo.read(path, "nexus")
    for node in tree.find_clades():
        if isinstance(node.confidence, str):
            node.name = node.confidence if node.name is None else node.name
            node.confidence = None
    return tree
