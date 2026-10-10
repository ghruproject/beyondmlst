"""Complete accessible public profile context, with explicit lineage scope and policy."""

from __future__ import annotations

import json
import shutil
from dataclasses import asdict
from pathlib import Path

from chronoclade.artifacts import file_sha256
from chronoclade.cglin import parse_code
from chronoclade.datasets import DatasetError, from_profile_records, load_dataset
from chronoclade.pathogenwatch import PathogenwatchClient
from chronoclade.prepare_provider import _has_profile, _save, catalogues_for_records
from chronoclade.profile_inputs import _canonical, _frozen_cglin, _grouped_analysis_exports
from chronoclade.typing_scopes import compatible_typing


def dataset_records(dataset):
    """Adapt authoritative typing without allowing evidence to overwrite identity/metadata."""
    from chronoclade.context_partitions import profile_records

    rows = profile_records(dataset)
    by_id = {row["sample_id"]: row for row in rows}
    allowed = {
        "source_genome_id",
        "aliases",
        "accession",
        "run_accession",
        "runAccession",
        "sample_accession",
        "sampleAccession",
        "assembly_accession",
        "assemblyAccession",
        "biosample",
        "run_accessions",
        "biosample_accessions",
        "assembly_accessions",
        "identity_candidates",
        "identity_resolved",
        "identity_conflict",
    }
    for evidence in dataset.crosswalk:
        by_id[evidence["sample_id"]].update({k: v for k, v in evidence.items() if k in allowed})
    return rows


def _matching(query, candidate, kind, level):
    if (
        not compatible_typing(query, candidate, kind)
        or query.get(kind + "_database_version") != candidate.get(kind + "_database_version")
        or query.get(kind + "_database_sha256") != candidate.get(kind + "_database_sha256")
    ) or any(
        row.get(kind + "_status")
        in {"unassigned", "unsupported", "conflict", "failed", "unverified", "unavailable"}
        for row in (query, candidate)
    ):
        return False
    if kind == "cglin":
        left, right = (
            parse_code(query.get("cglin_raw"))[0],
            parse_code(candidate.get("cglin_raw"))[0],
        )
        return (
            len(left) >= level
            and len(right) >= level
            and left[:level] == right[:level]
            and not query.get("cglin_provisional")
            and not candidate.get("cglin_provisional")
        )
    left, right = query.get("hiercc_codes") or {}, candidate.get("hiercc_codes") or {}
    return bool(left.get(level) and left.get(level) == right.get(level))


def discover_profile_context(
    dataset_manifest,
    output,
    *,
    lin_level=5,
    hiercc_level=None,
    client=None,
    public_typing=None,
    cglin_export=None,
    catalogues=None,
):
    """Publish all accessible compatible block profiles; never fetch context assemblies."""
    from tempfile import TemporaryDirectory

    from chronoclade.prepare_stage import _catalogues, run_prepare

    dataset = load_dataset(dataset_manifest, mmap_mode=None)
    if lin_level not in {5, 6, 7}:
        raise DatasetError("Context LIN level must be explicitly 5, 6 or 7")
    rows = dataset_records(dataset)
    for row in rows:
        reference = row.get("assembly_reference")
        if reference and "://" not in reference:
            source_root = Path(dataset_manifest).resolve()
            if source_root.is_file():
                source_root = source_root.parent
            row["assembly_reference"] = str((source_root / reference).resolve())
    queries = [r for r in rows if r["role"] == "input"]
    if any(
        r.get("species", "").casefold() not in {"klebsiella pneumoniae", "escherichia coli"}
        for r in queries
    ):
        raise DatasetError(
            "Public lineage context supports Klebsiella pneumoniae and Escherichia coli"
        )
    if any("coli" in r.get("species", "").casefold() for r in queries) and not hiercc_level:
        raise DatasetError("E. coli public context requires an explicit --hiercc-level")
    client = client or PathogenwatchClient()
    target = Path(output).resolve()
    if target.exists():
        raise DatasetError(f"Preparation output already exists: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".context-provider-", dir=target.parent) as temporary:
        work = Path(temporary)
        found, pools = {}, []
        scopes = {(r["species"], str(r.get("mlst_st") or "")) for r in queries}
        organisms = client.supported_organisms()
        for index, (species, st) in enumerate(sorted(scopes)):
            organism = next(
                (r for r in organisms if r.get("fullName", "").casefold() == species.casefold()),
                None,
            )
            if not st.isdigit() or organism is None:
                pools.append(
                    {
                        "species": species,
                        "st": st,
                        "status": "unavailable",
                        "reason": "missing_supported_species_or_exact_mlst_st",
                    }
                )
                continue
            envelope = client.freeze_catalogue(
                work / f"public_{index}.json", organism_id=str(organism["organismId"]), st=st
            )
            pool = [
                _canonical(
                    dict(r, sample_id="PW_" + r["source_genome_id"]),
                    origin="context",
                    species=species,
                )
                for r in envelope["rows"]
            ]
            export = {}
            pool = _grouped_analysis_exports(
                pool,
                client,
                work / f"lineages_{index}",
                export,
                download_names={"klebsiella-lincodes"},
            )
            if public_typing:
                from chronoclade.public_typing import annotate_public_typing, load_public_typing

                pool = annotate_public_typing(pool, load_public_typing(Path(public_typing)))
            if cglin_export:
                from chronoclade.cglin import load_cglin_export

                pool = _frozen_cglin(pool, load_cglin_export(Path(cglin_export)))
            relevant = [
                r for r in queries if r["species"] == species and str(r.get("mlst_st") or "") == st
            ]
            kind = "cglin" if "klebsiella" in species.casefold() else "hiercc"
            counts = {
                str(level): sum(any(_matching(q, r, kind, level) for q in relevant) for r in pool)
                for level in ([5, 6, 7] if kind == "cglin" else [hiercc_level])
            }
            matched = [
                r
                for r in pool
                if any(
                    _matching(q, r, kind, lin_level if kind == "cglin" else hiercc_level)
                    for q in relevant
                )
            ]
            requested = len(matched)
            matched = _grouped_analysis_exports(matched, client, work / f"profiles_{index}", export)
            for row in matched:
                row.pop("assembly", None)
                row.pop("assembly_reference", None)
                found[row["source_genome_id"]] = row
            pools.append(
                {
                    "species": species,
                    "st": st,
                    "public_search": envelope["provenance"],
                    "level_counts": counts,
                    "discovered": len(pool),
                    "requested": requested,
                    "retrieved": sum(_has_profile(r) for r in matched),
                    "exports": export,
                }
            )
        aliases = set()
        for row in rows:
            aliases.update(row.get("aliases") or [])
            for key in (
                "source_genome_id",
                "sample_id",
                "biosample",
                "sample_accession",
                "assembly_accession",
            ):
                if row.get(key):
                    aliases.add(row[key])
        context, exclusions = [], []
        for source_id, row in sorted(found.items()):
            if source_id in aliases or set(row.get("aliases") or []) & aliases:
                exclusions.append(
                    {"source_genome_id": source_id, "reason": "already_present_exact_identity"}
                )
            elif not _has_profile(row):
                exclusions.append({"source_genome_id": source_id, "reason": "profile_unavailable"})
            elif not any(
                compatible_typing(q, row, "cgmlst")
                and q.get("cgmlst_database_version") == row.get("cgmlst_database_version")
                and (q.get("cgmlst_database_sha256") or None)
                == (row.get("cgmlst_database_sha256") or None)
                for q in queries
            ):
                exclusions.append(
                    {"source_genome_id": source_id, "reason": "incompatible_profile_scope"}
                )
            else:
                row.update(role="context", identity_resolved=True)
                context.append(row)
        all_rows = rows + context
        explicit = [m.catalogue for m in dataset.profiles] + list(
            _catalogues(Path(catalogues)) if catalogues else []
        )
        catalogue_values = catalogues_for_records(all_rows, explicit)
        audit = {
            "policy": {"lin_level": lin_level, "hiercc_level": hiercc_level},
            "pools": pools,
            "accessible_matches": len(found),
            "usable": len(context),
            "exclusions": exclusions,
            "source_dataset_sha256": file_sha256(
                Path(dataset_manifest) / "dataset.json"
                if Path(dataset_manifest).is_dir()
                else Path(dataset_manifest)
            ),
            "context_assemblies": "not_requested",
            "global_completeness": "not_claimed",
        }
        _save(work / "catalogues.json", [asdict(c) for c in catalogue_values])
        _save(
            work / "profiles.json",
            {
                "records": all_rows,
                "lineages": list(dataset.lineages)
                + list(from_profile_records(context, catalogue_values).lineages)
                if context
                else list(dataset.lineages),
                "provenance": list(dataset.provenance)
                + [
                    {
                        "sample_id": r["sample_id"],
                        "field": "source_genome_id",
                        "source": "public_pathogenwatch",
                        "value": r["source_genome_id"],
                    }
                    for r in context
                ],
                "crosswalk": list(dataset.crosswalk)
                + list(from_profile_records(context, catalogue_values).crosswalk)
                if context
                else list(dataset.crosswalk),
                "retrieval": list(dataset.retrieval)
                + list(from_profile_records(context, catalogue_values).retrieval)
                if context
                else list(dataset.retrieval),
                "exclusions": list(dataset.exclusions),
                "conflicts": list(dataset.conflicts),
                "parameters": dict(dataset.parameters, profile_context=audit),
            },
        )
        source_root = Path(dataset_manifest).resolve()
        if source_root.is_file():
            source_root = source_root.parent
        source_copy = work / "source_dataset"
        source_copy.mkdir()
        raw_manifest = json.loads((source_root / "dataset.json").read_text())
        source_paths = {"dataset.json"}

        def collect_paths(value):
            if isinstance(value, dict):
                if "path" in value and "sha256" in value:
                    source_paths.add(value["path"])
                for child in value.values():
                    collect_paths(child)
            elif isinstance(value, list):
                for child in value:
                    collect_paths(child)

        # Copy only validated authoritative dataset descriptors, excluding legacy
        # parameter paths and deferred assembly references.
        collect_paths(raw_manifest["tables"])
        collect_paths(raw_manifest["profiles"])
        for name in sorted(source_paths):
            destination = source_copy / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source_root / name, destination)
        if not catalogue_values:
            from chronoclade.datasets import write_dataset
            from chronoclade.prepare_stage import _frozen_dataset

            sample_only = _frozen_dataset(
                json.loads((work / "profiles.json").read_text()), (), None
            )
            manifest = write_dataset(sample_only, work / "sample_only")
            return run_prepare(manifest, target, _provider_sources=work)
        result = run_prepare(
            work / "profiles.json",
            target,
            catalogues=work / "catalogues.json",
            _provider_sources=work,
        )
    return result
