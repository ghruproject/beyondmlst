"""Prepared-data ESM2 stage: real sequences, frozen vectors and shared selections."""
from __future__ import annotations

from collections import Counter
import csv
import hashlib
import json
from dataclasses import dataclass
import os
from pathlib import Path
import shutil
import tempfile

from chronoclade import __version__
from chronoclade.artifacts import file_sha256, write_json
from chronoclade.context_partitions import partition_records, profile_records
from chronoclade.datasets import AlleleMatrix, PreparedDataset, load_dataset
from chronoclade.selection_manifest import write_selection_ensemble, load_selection_ensemble
from chronoclade.selections import DistanceEvidence
from .allele_sequences import prepare_allele_sequences
from .engine import run_embeddings
from .genome_analysis import analyse_genomes, add_dna_comparison, DISTANCE_DEFINITION
from .genome_report import write_genome_report, write_genome_index
from .proteins import EmbeddingError
from .storage import write_npz


@dataclass(frozen=True)
class ESM2Result:
    manifest_path: Path
    report_path: Path


def _portable(value, base, root):
    """Convert owned absolute evidence paths before immutable publication."""
    if isinstance(value, dict):
        return {key: _portable(item, base, root) for key, item in value.items()}
    if isinstance(value, list):
        return [_portable(item, base, root) for item in value]
    if isinstance(value, str) and value.startswith(str(root) + os.sep):
        return Path(os.path.relpath(value, base)).as_posix()
    return value


def _subset(dataset, ids):
    wanted = set(ids)
    matrices = [matrix for matrix in dataset.profiles if wanted & set(matrix.sample_ids)]
    if len(matrices) != 1:
        raise EmbeddingError("Each ESM2 lineage block requires exactly one compatible profile matrix")
    matrix = matrices[0]
    selected = [ident for ident in matrix.sample_ids if ident in wanted]
    return PreparedDataset(tuple(row for row in dataset.samples if row["sample_id"] in wanted),
                           (AlleleMatrix.from_profiles(matrix.catalogue, selected,
                                                       {ident: matrix.profile(ident) for ident in selected}),),
                           dataset_id=dataset.dataset_id)


def _mapping_subset(source, out, ids, dataset_ids):
    # Validate column schema here; full cell/allele consistency is checked later.
    with Path(source).open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ["sample_id", "locus", "record_id"]:
            raise EmbeddingError("Mapping CSV requires sample_id,locus,record_id in order")
        all_rows = list(reader)
        if any(None in row or any(not row.get(key) for key in reader.fieldnames)
               or row["sample_id"] not in dataset_ids for row in all_rows):
            raise EmbeddingError("Mapping has malformed rows or unknown prepared sample IDs")
        rows = [row for row in all_rows if row["sample_id"] in set(ids)]
    with Path(out).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=reader.fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return out


def _copy_embeddings(source, output):
    from .sample_distances import _load_embeddings
    manifest, *_ = _load_embeddings(source)
    output.mkdir(parents=True, exist_ok=True)
    descriptor = manifest["artifacts"]["embeddings"]
    shutil.copyfile(Path(source).parent / descriptor["path"], output / "vectors.npz")
    manifest["artifacts"]["embeddings"]["path"] = "vectors.npz"
    manifest["reused_source"] = {"sha256": file_sha256(source)}
    write_json(output / "embeddings.json", manifest)
    return output / "embeddings.json"


def run_esm2(dataset_manifest, output, *, catalogue_manifest=None, embeddings_manifest=None,
             mapping_csv=None, model="8M", checkpoint=None, allow_download=False, device="auto",
             cache_dir=None, token_budget=4096, max_length=1022, expected_checkpoint_sha256=None,
             panel_loci=None, reference_sample_id=None, lin_level=5, hiercc_level=None,
             selection_size=50, selection_runs=5, pinned_ids=(), nearest_per_input=3,
             selection_seed=42, neighbour_count=3, group_threshold=0.0,
             baseline=True, terminal_stop="require"):
    """Publish an independent immutable analysis, completion manifest last.

    Provide a matching frozen DNA catalogue or both saved embeddings and a mapping
    CSV. No network acquisition or allele-ID-to-sequence guessing is performed.
    Missing contexts remain audited; missing inputs or pins fail the block's
    selection instead of silently dropping mandatory genomes.
    """
    if (catalogue_manifest is None) == (embeddings_manifest is None):
        raise EmbeddingError("Supply an allele catalogue or saved embeddings, choosing one")
    if (embeddings_manifest is not None) != (mapping_csv is not None):
        raise EmbeddingError("Saved embeddings require an explicit sample/locus mapping CSV")
    if (type(selection_size) is not int or selection_size < 0
        or type(selection_runs) is not int or selection_runs < 1
        or type(nearest_per_input) is not int or nearest_per_input < 1
        or type(selection_seed) is not int):
        raise EmbeddingError("Selection size/runs/neighbours/seed are invalid")
    dataset_manifest, output = Path(dataset_manifest).resolve(), Path(output).resolve()
    dataset = load_dataset(dataset_manifest)
    records = profile_records(dataset)
    blocks, audit = partition_records(records, lin_level=lin_level, hiercc_level=hiercc_level)
    if not any(row["role"] == "input" for row in records):
        raise EmbeddingError("ESM2 requires at least one prepared input sample")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise EmbeddingError("ESM2 output must be new or empty")
    all_ids = set(dataset.sample_ids)
    if reference_sample_id is not None and reference_sample_id not in all_ids:
        raise EmbeddingError("Reference sample must be an exact prepared sample ID")
    if set(pinned_ids) - all_ids:
        raise EmbeddingError("Pinned IDs must be exact prepared sample IDs")
    block_ids = {r["sample_id"] for block in blocks for r in block["records"]}
    if set(pinned_ids) - block_ids:
        raise EmbeddingError("Pinned IDs must belong to a resolved input lineage block")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.pending-", dir=output.parent))
    products, partition_blocks = [], []
    try:
        # Partition source uses eligible selection pools while retaining original
        # membership and every excluded genome in separately named fields.
        prepared_blocks = []
        for block in blocks:
            directory = staging / block["block_id"]
            directory.mkdir()
            subset = _subset(dataset, [row["sample_id"] for row in block["records"]])
            if catalogue_manifest:
                mapped = prepare_allele_sequences(subset, Path(catalogue_manifest), directory / "sequences",
                                                  terminal_stop=terminal_stop, max_length=max_length)
                embeddings = run_embeddings(mapped.fasta_path, directory / "embeddings", model=model,
                                            checkpoint=checkpoint, allow_download=allow_download, device=device,
                                            cache_dir=cache_dir, token_budget=token_budget, max_length=max_length,
                                            expected_checkpoint_sha256=expected_checkpoint_sha256).manifest_path
                mapping = mapped.mapping_path
            else:
                embeddings = _copy_embeddings(Path(embeddings_manifest), directory / "embeddings")
                mapping = _mapping_subset(mapping_csv, directory / "sample_loci.csv", subset.sample_ids, all_ids)
            embedding_record = json.loads(embeddings.read_text(encoding="utf-8"))
            write_json(embeddings, _portable(embedding_record, embeddings.parent, staging))
            requested_reference = reference_sample_id if reference_sample_id in subset.sample_ids else None
            summary, arrays, pairs, date_diagnostic = analyse_genomes(
                subset, embeddings, mapping, panel_loci=panel_loci,
                reference_sample_id=requested_reference, neighbour_count=neighbour_count,
                group_threshold=group_threshold)
            if catalogue_manifest:
                from .allele_sequences import _read_catalogue
                sequences, _, _ = _read_catalogue(catalogue_manifest, subset, subset.profiles[0])
                add_dna_comparison(summary, pairs, subset, sequences, neighbour_count)
            eligible = set(summary["sample_ids"])
            mandatory = {row["sample_id"] for row in subset.samples if row["role"] == "input"} | (set(pinned_ids) & set(subset.sample_ids))
            missing_mandatory = sorted(mandatory - eligible)
            summary["selection_status"] = "unavailable" if missing_mandatory else "complete"
            summary["missing_mandatory_ids"] = missing_mandatory
            partition_blocks.append({key: value for key, value in block.items() if key != "records"} |
                                    {"sample_ids": sorted(eligible), "source_sample_ids": list(subset.sample_ids),
                                     "embedding_exclusions": summary["exclusions"]})
            prepared_blocks.append((block, subset, summary, arrays, pairs, date_diagnostic, directory))
        partition_path = staging / "partitions.json"
        write_json(partition_path, {"schema": "chronoclade.esm2.partitions", "schema_version": 1,
                                   "dataset_id": dataset.dataset_id, "dataset_sha256": file_sha256(dataset_manifest),
                                   **audit, "blocks": partition_blocks})
        for block, subset, summary, arrays, pairs, date_diagnostic, directory in prepared_blocks:
            from chronoclade.report_components.ordination import write_ordination_views, ordination_viewer
            eligible = set(summary["sample_ids"])
            rows = [row for row in block["records"] if row["sample_id"] in eligible]
            by_id = {r["sample_id"]: r for r in rows}
            ordered = [by_id[ident] for ident in summary["sample_ids"]]
            plots = write_ordination_views(directory, "protein", arrays["coordinates"], ordered,
                                            title="ESM2 genome protein ordination",
                                            axis_labels=("Protein chord axis 1", "Protein chord axis 2"),
                                            groups=summary["groups"]["by_sample"])
            html = ordination_viewer([plots], directory)
            write_npz(directory / "genome_vectors.npz", **arrays)
            with (directory / "pairwise_distances.csv").open("w", encoding="utf-8", newline="") as stream:
                fields = list(pairs[0]) if pairs else ["sample_id_1", "sample_id_2", "distance", "panel_loci"]
                writer = csv.DictWriter(stream, fieldnames=fields)
                writer.writeheader()
                writer.writerows(pairs)
            ensemble = None
            if not summary["missing_mandatory_ids"]:
                ensemble = write_selection_ensemble(
                    [row for row in rows if row["role"] == "input"],
                    [row for row in rows if row["role"] == "context"],
                    {"distance_evidence": DistanceEvidence(tuple(pairs), "esm2",
                        DISTANCE_DEFINITION + "; checkpoint_sha256=" + summary["model"]["model_sha256"]
                        + "; panel_sha256=" + hashlib.sha256(json.dumps(summary["panel"]["loci"], separators=(",", ":")).encode()).hexdigest()
                        + "; vector_sha256=" + summary["provenance"]["vector_artifact"]["sha256"]),
                     "genetic_groups": [{"group_id": group, "sample_ids": [ident for ident, value in summary["groups"]["by_sample"].items() if value == group]}
                                        for group in sorted(set(summary["groups"]["by_sample"].values()))],
                     "partition": {"block_id": block["block_id"]}},
                    dataset_manifest=dataset_manifest, output=directory / "selection_bundle",
                    replicates=selection_runs, size=selection_size, nearest_per_query=nearest_per_input,
                    include=[row["sample_id"] for row in rows if row["role"] == "context" and row["sample_id"] in pinned_ids], seed=selection_seed,
                    partition_manifest=partition_path)
                load_selection_ensemble(ensemble)
            date_diagnostic["provenance"]["embeddings_manifest"] = "../embeddings/embeddings.json"
            date_diagnostic["provenance"]["mapping_csv"] = "../sequences/sample_loci.csv" if catalogue_manifest else "../sample_loci.csv"
            guide_report, temporal = None, None
            if baseline:
                from chronoclade.profile_analysis import analyse_profiles
                from chronoclade.profile_report import write_profile_report
                from chronoclade.temporal_diagnostics import cgmlst_root_to_tip
                from .temporal_report import write_temporal_report
                guide = analyse_profiles(block["records"], output=directory / "cgmlst_baseline",
                                         tree_limit=max(2, len(block["records"])), bootstrap_replicates=0)
                composition = Counter((row["origin"], row.get("country") or "Unknown",
                                       row.get("region") or "Unknown", row.get("nuts2") or "Unknown")
                                      for row in block["records"])
                guide["metadata_geography"] = [dict(origin=key[0], country=key[1], region=key[2], nuts2=key[3], count=count)
                                               for key, count in sorted(composition.items())]
                coverage = {}
                for role, label in (("input", "queries"), ("context", "context")):
                    group = [row for row in block["records"] if row["role"] == role]
                    called = sum(any(value is not None for value in row.get("cgmlst_profile", {}).values()) for row in group)
                    coverage[label] = {"total": len(group), "profiles_available": called, "profiles_missing": len(group) - called}
                guide_report = write_profile_report(guide, directory=directory / "cgmlst_baseline",
                    provenance={"source_dataset": dataset.dataset_id, "network_source": "conventional cgMLST NJ", "coverage": coverage},
                    report_label="Conventional cgMLST baseline for ESM2 exploration")
                write_json(directory / "cgmlst_baseline" / "profile_analysis.json",
                           _portable(guide, directory / "cgmlst_baseline", staging))
                for saved in (directory / "cgmlst_baseline").rglob("*.json"):
                    write_json(saved, _portable(json.loads(saved.read_text(encoding="utf-8")), saved.parent, staging))
                temporal = write_temporal_report(date_diagnostic,
                    cgmlst_root_to_tip(subset, reference_sample_id=date_diagnostic["reference"]["sample_id"]),
                    directory / "date_comparison")
            # External source paths belong in the stage's input references; saved
            # block diagnostics contain checksums and relative owned paths.
            summary["provenance"]["embeddings_manifest"] = "embeddings/embeddings.json"
            summary["provenance"]["mapping_csv"] = "sample_loci.csv" if mapping_csv else "sequences/sample_loci.csv"
            write_json(directory / "genome_analysis.json", summary)
            report = write_genome_report(summary, directory, ordination_html=html,
                                         ensemble_path=ensemble, baseline_report=guide_report,
                                         temporal_report=temporal)
            products.append({"block_id": block["block_id"], "status": summary["selection_status"],
                             "eligible_genomes": len(eligible), "report": report.relative_to(staging).as_posix(),
                             "analysis": (directory / "genome_analysis.json").relative_to(staging).as_posix(),
                             "ensemble": ensemble.relative_to(staging).as_posix() if ensemble else None,
                             "missing_mandatory_ids": summary["missing_mandatory_ids"]})
        report = write_genome_index(products, audit, staging)
        status = "complete" if products and not audit["unresolved_inputs"] and all(p["status"] == "complete" for p in products) else "partial" if products else "unavailable"
        inputs = {"dataset": dataset_manifest}
        if catalogue_manifest:
            inputs["allele_catalogue"] = Path(catalogue_manifest)
        else:
            inputs.update(embeddings=Path(embeddings_manifest), mapping=Path(mapping_csv))
        artifacts = {p.relative_to(staging).as_posix(): {"path": p.relative_to(staging).as_posix(), "sha256": file_sha256(p)}
                     for p in staging.rglob("*") if p.is_file()}
        write_json(staging / "esm2.json", {"schema": "chronoclade.esm2.analysis", "schema_version": 1,
                   "software_version": __version__, "status": status, "experimental": True,
                   "dataset_id": dataset.dataset_id,
                   "inputs": {name: {"path": os.path.relpath(path.resolve(), output), "sha256": file_sha256(path)} for name, path in inputs.items()},
                   "parameters": {"model": model, "device": device, "panel_loci": panel_loci,
                                  "lin_level": lin_level, "hiercc_level": hiercc_level,
                                  "selection_size": selection_size, "selection_runs": selection_runs,
                                  "selection_seed": selection_seed, "nearest_per_input": nearest_per_input,
                                  "pinned_ids": list(pinned_ids), "group_threshold": group_threshold,
                                  "baseline": baseline},
                   "partition_audit": audit, "blocks": products, "artifacts": artifacts})
        if output.exists():
            output.rmdir()
        staging.rename(output)
        return ESM2Result(output / "esm2.json", output / "index.html")
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
