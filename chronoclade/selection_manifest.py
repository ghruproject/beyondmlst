"""Portable, immutable selection bundles consumed by downstream tree stages.

A selection bundle contains its distance evidence and individual selections. Its
prepared dataset (and optional partition) remain independent immutable bundles:
move their common parent together to preserve relative source references.
"""

from __future__ import annotations

import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import tempfile

from chronoclade import __version__
from chronoclade.artifacts import file_sha256, write_json
from chronoclade.datasets import DatasetError, load_dataset
from chronoclade.errors import WorkflowError
from chronoclade.selections import read_distance_evidence, select_context_ensemble

SCHEMA_VERSION = 1
SELECTION_SCHEMA = "chronoclade.selection"
ENSEMBLE_SCHEMA = "chronoclade.selection-ensemble"


def _read(path):
    try:
        result = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise WorkflowError(f"Cannot read manifest {path}: {error}") from error
    if not isinstance(result, dict):
        raise WorkflowError(f"Manifest {path} must be an object")
    return result


def _reference(base, path, actual=None):
    actual = actual or path
    return {"path": Path(os.path.relpath(path, base)).as_posix(), "sha256": file_sha256(actual)}


def _resolve(base, spec, field, *, external=False):
    if not isinstance(spec, dict):
        raise WorkflowError(f"{field} requires a path/checksum reference")
    raw, digest = spec.get("path"), spec.get("sha256")
    if (
        not isinstance(raw, str)
        or not raw
        or "\\" in raw
        or PurePosixPath(raw).is_absolute()
        or not external
        and ".." in PurePosixPath(raw).parts
    ):
        raise WorkflowError(f"{field}.path must be a relative POSIX path")
    if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
        raise WorkflowError(f"{field}.sha256 requires a SHA256 checksum")
    path = (base / raw).resolve()
    if not external and not path.is_relative_to(base.resolve()):
        raise WorkflowError(f"{field}.path escapes its bundle")
    if not path.is_file():
        raise WorkflowError(f"{field}.path does not exist: {raw}")
    if file_sha256(path) != digest:
        raise WorkflowError(f"{field}.sha256 does not match {raw}")
    return path


def _header(manifest, schema):
    if manifest.get("schema") != schema:
        raise WorkflowError(f"Expected {schema} manifest")
    if (
        type(manifest.get("schema_version")) is not int
        or manifest["schema_version"] != SCHEMA_VERSION
    ):
        raise WorkflowError("Unsupported selection schema_version")
    if manifest.get("status") != "complete":
        raise WorkflowError("Selection manifest status must be complete")
    software = manifest.get("software")
    if not isinstance(software, dict) or not software.get("version"):
        raise WorkflowError("Selection manifest software.version is required")


def _ids(manifest, field):
    values = manifest.get(field)
    if not isinstance(values, list) or any(not isinstance(i, str) or not i for i in values):
        raise WorkflowError(f"{field} requires a list of nonempty sample identifiers")
    if len(values) != len(set(values)):
        raise WorkflowError(f"{field} repeats sample identifiers")
    return values


def _partition_block(path, dataset_id, dataset_hash, block_id, pool=None):
    partition = _read(path)
    if (
        partition.get("schema") not in {"chronoclade.cgmlst.partitions", "chronoclade.esm2.partitions"}
        or partition.get("schema_version") != 1
        or partition.get("dataset_id") != dataset_id
        or partition.get("dataset_sha256") != dataset_hash
    ):
        raise WorkflowError("source.partition must reference the same versioned prepared dataset")
    blocks = partition.get("blocks")
    if not isinstance(blocks, list) or any(not isinstance(block, dict) for block in blocks):
        raise WorkflowError("source.partition.blocks must be a list of objects")
    matches = (
        [block for block in blocks if block.get("block_id") == block_id]
        if block_id
        else [block for block in blocks if set(_ids(block, "sample_ids")) == pool]
    )
    if len(matches) != 1 or not isinstance(matches[0].get("block_id"), str):
        raise WorkflowError("source.partition_block_id must identify exactly one partition block")
    members = set(_ids(matches[0], "sample_ids"))
    if pool is not None and pool != members:
        raise WorkflowError("Selection pool must match exact source.partition block sample_ids")
    return matches[0]["block_id"]


def _sources(manifest, base):
    source = manifest.get("source")
    if not isinstance(source, dict):
        raise WorkflowError("source requires dataset and distance evidence references")
    dataset_path = _resolve(base, source.get("dataset"), "source.dataset", external=True)
    try:
        dataset = load_dataset(dataset_path)
    except DatasetError as error:
        raise WorkflowError(f"source.dataset is invalid: {error}") from error
    if source.get("dataset_id") != dataset.dataset_id:
        raise WorkflowError("source.dataset_id does not match the prepared dataset")
    if source.get("partition") is not None:
        if (
            not isinstance(source.get("partition_block_id"), str)
            or not source["partition_block_id"]
        ):
            raise WorkflowError("source.partition_block_id is required for a partition reference")
        partition_path = _resolve(base, source["partition"], "source.partition", external=True)
        pool = (
            set(_ids(manifest, "query_ids") + _ids(manifest, "available_context_ids"))
            if manifest.get("schema") == SELECTION_SCHEMA
            else None
        )
        _partition_block(
            partition_path,
            dataset.dataset_id,
            source["dataset"]["sha256"],
            source.get("partition_block_id"),
            pool,
        )
    evidence = _resolve(
        base, source.get("distance_evidence"), "source.distance_evidence", external=True
    )
    if not isinstance(source.get("method"), str) or not source["method"]:
        raise WorkflowError("source.method is required")
    if not isinstance(source.get("distance_definition"), str) or not source["distance_definition"]:
        raise WorkflowError("source.distance_definition is required")
    distance_data = _read(evidence)
    if (
        distance_data.get("method") != source["method"]
        or distance_data.get("distance_definition") != source["distance_definition"]
        or not isinstance(distance_data.get("records"), list)
    ):
        raise WorkflowError(
            "source.distance_evidence must preserve its method, definition and records"
        )
    dataset_ids = set(dataset.sample_ids)
    for row in distance_data["records"]:
        if not isinstance(row, dict) or any(
            not isinstance(row.get(field), str) or row[field] not in dataset_ids
            for field in ("sample_id_1", "sample_id_2")
        ):
            raise WorkflowError(
                "source.distance_evidence contains identifiers outside source.dataset"
            )
    return dataset, evidence


def load_selection_manifest(path: str | Path) -> dict:
    """Validate saved identities, mandatory membership and every source checksum.

    Return the persisted dictionary unchanged. Tree stages must consume
    selected_sample_ids in their recorded order, without choosing replacements.
    """
    path = Path(path).expanduser().resolve()
    manifest = _read(path)
    _header(manifest, SELECTION_SCHEMA)
    dataset, _ = _sources(manifest, path.parent)
    query = _ids(manifest, "query_ids")
    contexts = _ids(manifest, "selected_context_ids")
    selected = _ids(manifest, "selected_sample_ids")
    pool = _ids(manifest, "available_context_ids")
    mandatory = _ids(manifest, "mandatory_sample_ids")
    pinned = _ids(manifest, "pinned_context_ids")
    nearest = _ids(manifest, "required_nearest_context_ids")
    if selected != query + contexts or set(query) & set(pool) or not set(contexts).issubset(pool):
        raise WorkflowError(
            "selected_sample_ids must equal inputs plus selected contexts from the pool"
        )
    if not set(query + pool).issubset(dataset.sample_ids):
        raise WorkflowError("Selection pool contains identifiers outside source.dataset")
    input_ids = {row["sample_id"] for row in dataset.samples if row["role"] == "input"}
    if set(query) != input_ids.intersection(query + pool):
        raise WorkflowError(
            "query_ids must retain source.dataset input roles within the selection pool"
        )
    if set(mandatory) != set(query + pinned + nearest) or not set(mandatory).issubset(selected):
        raise WorkflowError(
            "mandatory_sample_ids must retain all inputs, pins and required nearest IDs"
        )
    if not set(pinned + nearest).issubset(contexts):
        raise WorkflowError("Pinned and required nearest identifiers must be selected contexts")
    decisions = manifest.get("decisions")
    inputs = manifest.get("input_decisions")
    if not isinstance(decisions, list) or not isinstance(inputs, list):
        raise WorkflowError("decisions and input_decisions are required")
    if any(not isinstance(d, dict) for d in decisions + inputs):
        raise WorkflowError("decisions must contain objects")
    if [d.get("sample_id") for d in decisions] != contexts or [
        d.get("sample_id") for d in inputs
    ] != query:
        raise WorkflowError("decisions must cover exact ordered selected identifiers")
    for decision in inputs + decisions:
        reasons = decision.get("reasons")
        if (
            not isinstance(reasons, list)
            or not reasons
            or any(not isinstance(r, str) or not r for r in reasons)
        ):
            raise WorkflowError("decisions.reasons must be nonempty strings")
        if decision.get("mandatory") is not (decision["sample_id"] in mandatory):
            raise WorkflowError("decisions.mandatory disagrees with mandatory_sample_ids")
    for decision in decisions:
        if decision.get("reason") != decision["reasons"][0]:
            raise WorkflowError("decisions.reason must match the first recorded reason")
        if decision["sample_id"] in pinned and "user_requested" not in decision["reasons"]:
            raise WorkflowError("Pinned decisions must record user_requested")
    if any(d["reasons"] != ["input"] for d in inputs):
        raise WorkflowError("input_decisions must record the input reason")
    neighbour_records = manifest.get("nearest_neighbours")
    if (
        not isinstance(neighbour_records, list)
        or any(not isinstance(n, dict) for n in neighbour_records)
        or any(not isinstance(n.get("query_id"), str) for n in neighbour_records)
        or sorted(n.get("query_id", "") for n in neighbour_records) != sorted(query)
    ):
        raise WorkflowError("nearest_neighbours must cover every query exactly once")
    required_ids = set()
    for record in neighbour_records:
        required_ids.update(_ids(record, "required_ids"))
        if not set(_ids(record, "boundary_tie_ids")).issubset(record["required_ids"]):
            raise WorkflowError("nearest_neighbours.boundary_tie_ids must be retained")
    if required_ids != set(nearest):
        raise WorkflowError("required_nearest_context_ids must retain every recorded neighbour/tie")
    budget = manifest.get("context_budget")
    if type(budget) is not int or budget < 0:
        raise WorkflowError("context_budget must be a non-negative integer")
    expected_overrun = max(0, len(set(mandatory) - set(query)) - budget)
    if manifest.get("budget_overrun") != expected_overrun:
        raise WorkflowError("budget_overrun does not match mandatory contexts and budget")
    if len(contexts) > budget + expected_overrun:
        raise WorkflowError("Selected optional contexts exceed context_budget")
    if manifest.get("selected_contexts") != len(contexts) or manifest.get(
        "available_contexts"
    ) != len(pool):
        raise WorkflowError("Selection context counts disagree with exact identifiers")
    if type(manifest.get("seed")) is not int:
        raise WorkflowError("Selection seed must be an integer")
    return manifest


def load_selection_ensemble(path: str | Path) -> dict:
    """Validate every independent run, including duplicate alternatives."""
    path = Path(path).expanduser().resolve()
    manifest = _read(path)
    _header(manifest, ENSEMBLE_SCHEMA)
    _sources(manifest, path.parent)
    anchors = _ids(manifest, "common_anchor_ids")
    targets = _ids(manifest, "target_ids")
    runs = manifest.get("selections")
    if not isinstance(runs, list) or not runs or manifest.get("requested_replicates") != len(runs):
        raise WorkflowError("selections must include every requested replicate")
    seen, selections = set(), []
    for index, run in enumerate(runs):
        selection_path = _resolve(path.parent, run, f"selections[{index}]")
        selection = load_selection_manifest(selection_path)
        if run.get("selection_id") != selection.get("selection_id") or run["selection_id"] in seen:
            raise WorkflowError("selections contain missing or duplicate selection_id")
        seen.add(run["selection_id"])
        if not set(anchors + targets).issubset(selection["selected_sample_ids"]):
            raise WorkflowError("Selection does not retain ensemble anchors/targets")
        if anchors != selection["mandatory_sample_ids"] or targets != selection["query_ids"]:
            raise WorkflowError("Selection mandatory inputs must match ensemble anchors/targets")
        for field in ("dataset", "partition", "distance_evidence"):
            left, right = manifest["source"].get(field), selection["source"].get(field)
            if (left is None) != (right is None) or left and left["sha256"] != right["sha256"]:
                raise WorkflowError(f"Selection source.{field} must match ensemble source")
        selections.append(selection)
    overlaps = []
    for index, left in enumerate(selections):
        a = set(left["selected_sample_ids"])
        for right in selections[index + 1 :]:
            b = set(right["selected_sample_ids"])
            overlaps.append(
                {
                    "left": left["selection_id"],
                    "right": right["selection_id"],
                    "shared_ids": sorted(a & b),
                    "intersection": len(a & b),
                    "union": len(a | b),
                    "jaccard": len(a & b) / len(a | b) if a | b else 1.0,
                }
            )
    if manifest.get("pairwise_overlap") != overlaps:
        raise WorkflowError("pairwise_overlap must match exact selected identifiers")
    return manifest


def write_selection_ensemble(
    queries: list[dict],
    context: list[dict],
    analysis: dict,
    *,
    dataset_manifest: Path,
    output: Path,
    replicates: int = 1,
    size: int = 50,
    nearest_per_query: int = 3,
    include: list[str] | None = None,
    seed: int = 42,
    partition_manifest: Path | None = None,
) -> Path:
    """Atomically publish a fresh bundle; never overwrite a saved selection."""
    dataset_manifest, output = Path(dataset_manifest).resolve(), Path(output).resolve()
    dataset = load_dataset(dataset_manifest)
    evidence = read_distance_evidence(analysis)
    dataset_ids = set(dataset.sample_ids)
    if any(
        row.get(field) not in dataset_ids
        for row in evidence.records
        for field in ("sample_id_1", "sample_id_2")
    ):
        raise WorkflowError("Distance evidence contains identifiers outside source.dataset")
    selections, ensemble = select_context_ensemble(
        queries,
        context,
        dict(analysis, distance_evidence=evidence),
        replicates=replicates,
        size=size,
        nearest_per_query=nearest_per_query,
        include=include,
        seed=seed,
    )
    pool = set(selections[0]["query_ids"] + selections[0]["available_context_ids"])
    if not pool.issubset(dataset.sample_ids):
        raise WorkflowError("Selection pool contains identifiers outside source.dataset")
    input_ids = {row["sample_id"] for row in dataset.samples if row["role"] == "input"}
    if set(selections[0]["query_ids"]) != input_ids & pool:
        raise WorkflowError(
            "query_ids must retain source.dataset input roles within the selection pool"
        )
    block_id = None
    if partition_manifest:
        block_id = _partition_block(
            Path(partition_manifest),
            dataset.dataset_id,
            file_sha256(dataset_manifest),
            analysis.get("partition", {}).get("block_id"),
            pool,
        )
    if output.exists():
        raise WorkflowError(f"Selection output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.pending-", dir=output.parent))
    try:
        evidence_path = staging / "distance_evidence.json"
        write_json(
            evidence_path,
            {
                "method": evidence.method,
                "distance_definition": evidence.distance_definition,
                "records": list(evidence.records),
            },
        )
        header = {
            "schema_version": SCHEMA_VERSION,
            "software": {"name": "chronoclade", "version": __version__},
            "status": "complete",
        }
        parameters = {
            "size": size,
            "nearest_per_query": nearest_per_query,
            "include": include or [],
            "seed": seed,
            "replicates": replicates,
        }

        def source(base):
            return {
                "dataset": _reference(base, dataset_manifest),
                "dataset_id": dataset.dataset_id,
                "partition": _reference(base, Path(partition_manifest).resolve())
                if partition_manifest
                else None,
                "partition_block_id": block_id,
                "distance_evidence": _reference(base, output / evidence_path.name, evidence_path),
                "method": evidence.method,
                "distance_definition": evidence.distance_definition,
            }

        runs = []
        for selection in selections:
            relative = Path("selections") / f"{selection['selection_id']}.json"
            selection.update(
                header,
                schema=SELECTION_SCHEMA,
                parameters=parameters,
                source=source((output / relative).parent),
            )
            write_json(staging / relative, selection)
            runs.append(
                {
                    **_reference(output, output / relative, staging / relative),
                    "selection_id": selection["selection_id"],
                    "status": "complete",
                    "alternative_status": selection["alternative_status"],
                    "duplicate_of": selection["duplicate_of"],
                }
            )
        ensemble.update(
            header,
            schema=ENSEMBLE_SCHEMA,
            parameters=parameters,
            source=source(output),
            selections=runs,
        )
        write_json(staging / "ensemble.json", ensemble)
        for artifact in staging.rglob("*"):
            if artifact.is_file():
                with artifact.open("rb") as handle:
                    os.fsync(handle.fileno())
        staging.rename(output)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return output / "ensemble.json"
