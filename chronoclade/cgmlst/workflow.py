"""Own a complete profile-stage output without acquiring assemblies or dating."""

from collections import Counter
from dataclasses import dataclass
import os
from pathlib import Path

from chronoclade import __version__
from chronoclade.artifacts import file_sha256, write_json
from chronoclade.datasets import load_dataset
from chronoclade.profile_analysis import analyse_profiles
from chronoclade.profile_report import write_profile_report
from .report import write_partition_index
from chronoclade.selection_manifest import write_selection_ensemble
from chronoclade.context_partitions import partition_records, profile_records


@dataclass(frozen=True)
class CGMLSTResult:
    manifest_path: Path
    report_path: Path


def _relative_artifacts(value, root):
    """Keep legacy stage evidence paths portable in the saved analysis adapter."""
    if isinstance(value, dict):
        return {key: _relative_artifacts(item, root) for key, item in value.items()}
    if isinstance(value, list):
        return [_relative_artifacts(item, root) for item in value]
    if isinstance(value, str) and value.startswith(str(root) + os.sep):
        return Path(value).relative_to(root).as_posix()
    return value


def run_cgmlst(
    dataset_manifest,
    output,
    *,
    lin_level=5,
    hiercc_level=None,
    replicates=1,
    context_size=50,
    nearest_per_query=3,
    include=None,
    seed=42,
    bootstrap_replicates=30,
    min_overlap=0.9,
    distance_threshold=0.02,
):
    dataset_manifest, output = Path(dataset_manifest).resolve(), Path(output).resolve()
    dataset = load_dataset(dataset_manifest)
    records = profile_records(dataset)
    blocks, audit = partition_records(records, lin_level=lin_level, hiercc_level=hiercc_level)
    if not any(row["role"] == "input" for row in records):
        raise ValueError("cgmlst requires at least one input sample")
    # Immutable namespaces avoid overwriting another stage or earlier completed result.
    if output.exists() and any(output.iterdir()):
        raise ValueError("cgmlst output must be empty; choose a separate output directory")
    if replicates < 1 or context_size < 0 or nearest_per_query < 1:
        raise ValueError(
            "Selections require positive replicates/neighbours and nonnegative context size"
        )
    requested = {}
    selectable = [row for block in blocks for row in block["records"] if row["role"] == "context"]
    for identifier in include or []:
        matches = [
            row
            for row in selectable
            if identifier
            in {
                row["sample_id"],
                str(row.get("source_genome_id") or ""),
                str(row.get("accession") or ""),
            }
        ]
        if len(matches) != 1:
            raise ValueError(f"Pinned context {identifier!r} is not uniquely in a selected block")
        requested[identifier] = matches[0]["sample_id"]
    output.mkdir(parents=True, exist_ok=True)
    partition_path = output / "partitions.json"
    partition_evidence = {
        "schema": "chronoclade.cgmlst.partitions",
        "schema_version": 1,
        "dataset_id": dataset.dataset_id,
        "dataset_sha256": file_sha256(dataset_manifest),
        **audit,
        "blocks": [
            {key: value for key, value in block.items() if key != "records"}
            | {"sample_ids": [row["sample_id"] for row in block["records"]]}
            for block in blocks
        ],
    }
    write_json(partition_path, partition_evidence)
    products = []
    for block in blocks:
        directory = output / block["block_id"]
        rows = block["records"]
        queries = [row for row in rows if row["role"] == "input"]
        contexts = [row for row in rows if row["role"] == "context"]
        analysis = analyse_profiles(
            rows,
            output=directory,
            seed=seed,
            min_overlap=min_overlap,
            bootstrap_replicates=bootstrap_replicates,
            distance_threshold=distance_threshold,
            tree_limit=max(2, len(rows)),
        )
        analysis.update(
            species=block["species"],
            lineage="ST"
            + block["mlst_st"]
            + " · "
            + block["kind"]
            + " "
            + str(block["level"])
            + " · "
            + ".".join(map(str, block["prefix"])),
            partition={key: value for key, value in block.items() if key != "records"},
        )
        composition = Counter(
            (
                row["origin"],
                row.get("country") or "Unknown",
                row.get("region") or "Unknown",
                row.get("nuts2") or "Unknown",
            )
            for row in rows
        )
        analysis["metadata_geography"] = [
            dict(origin=key[0], country=key[1], region=key[2], nuts2=key[3], count=count)
            for key, count in sorted(composition.items())
        ]
        ensemble_path = write_selection_ensemble(
            queries,
            contexts,
            analysis,
            dataset_manifest=dataset_manifest,
            output=directory / "selection_bundle",
            replicates=replicates,
            size=context_size,
            nearest_per_query=nearest_per_query,
            include=[
                alias
                for alias, ident in requested.items()
                if ident in {row["sample_id"] for row in contexts}
            ],
            seed=seed,
            partition_manifest=partition_path,
        )
        provenance = {
            "source_dataset": dataset.dataset_id,
            "context_scope": audit["context_scope"],
            "coverage": {
                label: {
                    "total": len(group),
                    "profiles_available": sum(
                        any(value is not None for value in row.get("cgmlst_profile", {}).values())
                        for row in group
                    ),
                }
                for label, group in (("queries", queries), ("context", contexts))
            },
        }
        report = write_profile_report(
            analysis,
            directory=directory,
            provenance=provenance,
            report_label="cgMLST profile report",
        )
        write_json(directory / "profile_analysis.json", _relative_artifacts(analysis, directory))
        products.append(
            {
                "block_id": block["block_id"],
                "label": "ST"
                + block["mlst_st"]
                + " · "
                + block["kind"]
                + " "
                + str(block["level"])
                + " · "
                + ".".join(map(str, block["prefix"])),
                "input_count": len(queries),
                "context_count": len(contexts),
                "report": report.relative_to(output).as_posix(),
                "ensemble": ensemble_path.relative_to(output).as_posix(),
                "report_sha256": file_sha256(report),
                "ensemble_sha256": file_sha256(ensemble_path),
            }
        )
    status = (
        "complete"
        if blocks and not audit["unresolved_inputs"]
        else "partial"
        if blocks
        else "unavailable"
    )
    report_path = write_partition_index(
        products,
        audit,
        records_count=len(records),
        lin_level=lin_level,
        hiercc_level=hiercc_level,
        status=status,
        output=output,
    )
    manifest = {
        "schema": "chronoclade.cgmlst.analysis",
        "schema_version": 1,
        "software_version": __version__,
        "status": status,
        "dataset_id": dataset.dataset_id,
        "dataset": {
            "path": os.path.relpath(dataset_manifest, output),
            "sha256": file_sha256(dataset_manifest),
        },
        "parameters": {
            "lin_level": lin_level,
            "hiercc_level": hiercc_level,
            "replicates": replicates,
            "context_size": context_size,
            "nearest_per_query": nearest_per_query,
            "include": include or [],
            "seed": seed,
            "bootstrap_replicates": bootstrap_replicates,
            "min_overlap": min_overlap,
            "distance_threshold": distance_threshold,
        },
        "partitions": {"path": partition_path.name, "sha256": file_sha256(partition_path)},
        "blocks": products,
        "report": {"path": "index.html", "sha256": file_sha256(report_path)},
    }
    manifest_path = output / "cgmlst.json"
    write_json(manifest_path, manifest)
    return CGMLSTResult(manifest_path, report_path)
