"""Profile-first stages; assemblies are acquired only after context selection."""

from __future__ import annotations

import csv
import json
import math
import random
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Any

from chronoclade.errors import WorkflowError
from chronoclade.pathogenwatch import content_hash


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _distances(analysis: dict) -> dict[tuple[str, str], float]:
    path = analysis.get("paths", {}).get("pairwise_distances")
    result = {}
    if path and Path(path).is_file():
        with Path(path).open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                try:
                    distance = float(row["distance"])
                except (ValueError, TypeError, KeyError):
                    continue
                if math.isfinite(distance):
                    a, b = row["sample_id_1"], row["sample_id_2"]
                    result[a, b] = result[b, a] = distance
    return result


def select_assembly_context(
    queries: list[dict],
    context: list[dict],
    analysis: dict,
    *,
    size: int = 50,
    nearest_per_query: int = 3,
    include: list[str] | None = None,
    seed: int = 42,
) -> tuple[list[dict], dict]:
    """Select fairly around queries, then across time/place and genetic diversity."""
    if size < 0 or nearest_per_query < 1:
        raise WorkflowError("Context size must be non-negative and nearest-per-query positive")
    by_id = {row["sample_id"]: row for row in context}
    if len(by_id) != len(context):
        raise WorkflowError("Duplicate contextual sample identifiers")
    include = include or []
    pinned = []
    for identifier in include:
        matches = [
            r["sample_id"]
            for r in context
            if identifier
            in {r["sample_id"], str(r.get("source_genome_id", "")), str(r.get("accession", ""))}
        ]
        if len(matches) != 1:
            raise WorkflowError(
                f"Requested context genome {identifier!r} is not uniquely in the pool"
            )
        if matches[0] not in pinned:
            pinned.append(matches[0])
    if len(pinned) > size:
        raise WorkflowError("Requested context genomes exceed --context-size; increase the budget")
    distances = _distances(analysis)
    reasons: dict[str, str] = {identifier: "user_requested" for identifier in pinned}
    selected = list(pinned)
    queues = {
        q["sample_id"]: sorted(
            [identifier for identifier in by_id if (q["sample_id"], identifier) in distances],
            key=lambda identifier: (distances[q["sample_id"], identifier], identifier),
        )[:nearest_per_query]
        for q in sorted(queries, key=lambda r: r["sample_id"])
    }
    nearest_budget = min(size, max(len(queries), math.ceil(size / 2)))
    for rank in range(nearest_per_query):
        for query_id, candidates in queues.items():
            if len(selected) >= nearest_budget:
                break
            if rank < len(candidates) and candidates[rank] not in reasons:
                identifier = candidates[rank]
                selected.append(identifier)
                reasons[identifier] = f"cgmlst_neighbour:{query_id}:rank={rank + 1}"
    groups = {
        identifier: str(group["group_id"])
        for group in analysis.get("genetic_groups", [])
        for identifier in group.get("sample_ids", [])
    }
    strata = defaultdict(list)
    for identifier, row in sorted(by_id.items()):
        key = (
            groups.get(identifier, "unassigned"),
            str(row.get("collection_date", ""))[:4] or "Unknown",
            str(row.get("nuts2") or row.get("region") or row.get("country") or "Unknown"),
        )
        strata[key].append(identifier)
    cells = sorted(strata)
    random.Random(seed).shuffle(cells)
    target = min(size, len(selected) + math.ceil(size / 4))
    for cell in cells:
        if len(selected) >= target:
            break
        candidates = [identifier for identifier in strata[cell] if identifier not in reasons]
        if candidates:
            identifier = candidates[0]
            selected.append(identifier)
            reasons[identifier] = "genetic_group_time_region_representative"
    while len(selected) < min(size, len(by_id)):
        remaining = [identifier for identifier in by_id if identifier not in reasons]
        comparable = [
            identifier
            for identifier in remaining
            if any((identifier, other) in distances for other in selected)
        ]
        if comparable:
            identifier = max(
                comparable,
                key=lambda candidate: (
                    min(
                        distances[candidate, other]
                        for other in selected
                        if (candidate, other) in distances
                    ),
                    candidate,
                ),
            )
            reason = "cgmlst_diversity_representative"
        else:
            identifier = sorted(remaining)[0]
            reason = "metadata_background_no_comparable_profile"
        selected.append(identifier)
        reasons[identifier] = reason
    records = [
        dict(by_id[identifier], selection_reason=reasons[identifier]) for identifier in selected
    ]
    return records, {
        "context_budget": size,
        "available_contexts": len(context),
        "selected_contexts": len(records),
        "nearest_per_query": nearest_per_query,
        "seed": seed,
        "decisions": [
            {"sample_id": identifier, "reason": reasons[identifier]} for identifier in selected
        ],
        "queries_without_selected_comparable_neighbour": [
            q["sample_id"]
            for q in queries
            if not any((q["sample_id"], identifier) in distances for identifier in selected)
        ],
        "scope": "Selection within the recorded profile catalogue; diversity uses cgMLST distance, not calendar time.",
    }


def _stage_landing(output: Path, name: str, records: list[dict]) -> Path:
    from html import escape

    path = output / f"{name}.html"
    links = []
    for row in records:
        report = row.get("outputs", {}).get("html_report")
        if report and Path(report).is_file():
            relative = Path(report).resolve().relative_to(output.resolve())
            links.append(
                f'<li><a href="{escape(relative.as_posix(), quote=True)}">'
                f"{escape(str(row['species']))} {escape(str(row['lineage']))}</a></li>"
            )
        else:
            links.append(
                f"<li>{escape(str(row.get('species', '')))} "
                f"{escape(str(row.get('lineage', '')))}: "
                f"{escape(str(row.get('status', 'not analysed')))}</li>"
            )
    path.write_text(
        '<!doctype html><html lang="en"><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        f"<title>ChronoClade {name}</title><body><h1>{escape(name.title())} results</h1>"
        '<a href="index.html">All completed stages</a><ul>'
        + "".join(links)
        + "</ul></body></html>",
        encoding="utf-8",
    )
    return path


def _snapshot_stage(record: dict, stage: str) -> None:
    """Freeze reader-facing assets so later analysis cannot change an earlier report."""
    directory = Path(record["outputs"]["html_report"]).parent
    destination = directory / "stages" / stage
    shutil.rmtree(destination, ignore_errors=True)
    destination.mkdir(parents=True)
    allowed = {".html", ".csv", ".tsv", ".json", ".txt", ".svg", ".png",
               ".nexus", ".newick", ".nwk", ".fasta", ".zip"}
    excluded = {"logs", "randomisations", "stages"}
    if stage == "full":
        excluded.add("clock")
    dated_tree_supported = bool(record.get("temporal_signal_supported"))
    if stage == "full" or not dated_tree_supported:
        excluded.add("timetree")
    for path in sorted(directory.rglob("*")):
        relative = path.relative_to(directory)
        if (not path.is_file() or excluded.intersection(relative.parts)
                or path.suffix.lower() not in allowed
                or path.name.startswith(("report.full", "report.finish", "supporting_results.full",
                                         "supporting_results.finish"))
                or (stage == "full" and path.name.startswith(("temporal_signal", "date_randomisation")))
                or ((stage == "full" or not dated_tree_supported)
                    and path.name.startswith(("timetree", "node_dates")))):
            continue
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        if target.suffix == ".html":
            target.write_text(target.read_text().replace(
                "../../fast/", "../../../../fast/"))
    def relocated(value: Any) -> Any:
        if isinstance(value, dict):
            return {key: relocated(item) for key, item in value.items()}
        if isinstance(value, list):
            return [relocated(item) for item in value]
        if isinstance(value, str) and value.startswith(str(directory) + "/"):
            candidate = destination / Path(value).relative_to(directory)
            if candidate.is_file():
                return str(candidate)
        return value

    snapshot_record = relocated(record)
    _write_json(destination / "report.json", snapshot_record)
    from chronoclade.report import write_supporting_bundle
    write_supporting_bundle(destination)
    # Keep previous file names as convenient entry points, with immutable asset URLs.
    html = (destination / "report.html").read_text()
    import re
    html = re.sub(r'(href|src)="([^"#:]+)"',
                  lambda m: m[0] if "://" in m[2] or m[2].startswith("/")
                  else f'{m[1]}="stages/{stage}/{m[2]}"', html)
    (directory / f"report.{stage}.html").write_text(html)
    shutil.copy2(destination / "report.json", directory / f"report.{stage}.json")
    shutil.copy2(destination / "supporting_results.zip", directory / f"supporting_results.{stage}.zip")
    record["outputs"]["html_report"] = str(destination / "report.html")
    record["outputs"]["supporting_results"] = str(destination / "supporting_results.zip")


def run_staged_workflow(
    metadata: Path | None,
    *,
    output: Path,
    mode: str = "fast",
    collection: str | None = None,
    accessions: Path | None = None,
    species: str | None = None,
    catalogue: Path | None = None,
    public_typing: Path | None = None,
    query_typing: Path | None = None,
    typing_config: Path | None = None,
    cglin_export: Path | None = None,
    profile_limit: int = 500,
    context_size: int = 50,
    nearest_per_query: int = 3,
    include_genomes: list[str] | None = None,
    bootstrap_replicates: int = 10,
    distance_threshold: float = 0.02,
    threads: int = 4,
    lineage_jobs: int = 2,
    randomisation_jobs: int = 4,
    randomisations: int = 100,
    temporal_p_value: float = 0.05,
    min_samples: int = 4,
    seed: int = 42,
    force: bool = False,
    date_randomisation_method: str = "root_to_tip",
    context_manifest: Path | None = None,
    client: Any = None,
    dry_run: bool = False,
) -> dict:
    """Execute cumulative stages while preserving reports and selected input hashes."""
    from chronoclade.profile_analysis import analyse_profiles
    from chronoclade.profile_inputs import materialise_assemblies, resolve_profile_inputs
    from chronoclade.profile_report import write_profile_report, write_stage_index
    from chronoclade.report import write_supporting_bundle
    from chronoclade.workflow import run_workflow

    if mode not in {"fast", "full", "finish"}:
        raise WorkflowError("--mode must be fast, full or finish")
    if metadata is None and collection is None and accessions is None:
        raise WorkflowError("Supply a metadata CSV, --accessions or --collection")
    output = output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    if dry_run:
        # Planning must not fetch exports, type queries or download any assembly.
        result = {
            "analysis_mode": mode,
            "dry_run": True,
            "collection": collection,
            "metadata": str(metadata) if metadata else None,
            "accessions": str(accessions) if accessions else None,
            "context_budget_per_lineage": context_size,
            "stages": ["fast"]
            + (["full"] if mode != "fast" else [])
            + (["finish"] if mode == "finish" else []),
        }
        _write_json(output / "plan.json", result)
        return result
    inputs = resolve_profile_inputs(
        metadata,
        collection=collection,
        accessions=accessions,
        species=species,
        output=output / "inputs",
        client=client,
        catalogue=catalogue,
        public_typing=public_typing,
        query_typing=query_typing,
        typing_config=typing_config,
        cglin_export=cglin_export,
        profile_limit=profile_limit,
        seed=seed,
    )
    queries, context = inputs["queries"], inputs["context"]
    analysis = analyse_profiles(
        queries + context,
        output=output / "fast",
        seed=seed,
        bootstrap_replicates=bootstrap_replicates,
        distance_threshold=distance_threshold,
    )
    def geography(rows: list[dict]) -> list[dict]:
        counts: dict[tuple, int] = defaultdict(int)
        for row in rows:
            counts[(row.get("origin", "context"), row.get("country") or "Unknown",
                    row.get("region") or "Unknown", row.get("nuts2") or "Unknown")] += 1
        return [dict(origin=key[0], country=key[1], region=key[2], nuts2=key[3], count=count)
                for key, count in sorted(counts.items())]

    analysis["metadata_geography"] = geography(queries + context)
    analysis["public_catalogue_geography"] = geography(inputs.get("catalogue_rows", []))
    _write_json(output / "fast" / "metadata_geography.json", {
        "analysis_records": analysis["metadata_geography"],
        "public_catalogue": analysis["public_catalogue_geography"],
    })
    _write_json(output / "fast" / "profile_analysis.json", analysis)
    profile_report = write_profile_report(
        analysis, directory=output / "fast", provenance=inputs["provenance"]
    )
    write_supporting_bundle(output / "fast")
    fingerprint = content_hash(
        {
            "resolved_records_sha256": inputs["provenance"].get(
                "resolved_records_sha256", content_hash({"queries": queries, "context": context})
            ),
            "context_size": context_size,
            "seed": seed,
            "nearest_per_query": nearest_per_query,
            "include_genomes": include_genomes or [],
            "distance_threshold": distance_threshold,
            "bootstrap_replicates": bootstrap_replicates,
            "randomisations": randomisations,
            "temporal_p_value": temporal_p_value,
            "date_randomisation_method": date_randomisation_method,
            "min_samples": min_samples,
        }
    )
    previous_path = output / "stages.json"
    previous = json.loads(previous_path.read_text()) if previous_path.is_file() else {}
    stages = previous.get("stages", {}) if previous.get("fingerprint") == fingerprint else {}
    if previous and previous.get("fingerprint") != fingerprint:
        # Invalidate reader-facing dated products. Native stage fingerprints
        # independently validate the exact alignment/tree inputs before reuse.
        for path in (output / "assembly").glob("*/report.finish.*"):
            path.unlink()
        for path in (output / "assembly").glob("*/supporting_results.finish.zip"):
            path.unlink()
        for lineage_directory in (output / "assembly").glob("*"):
            if not lineage_directory.is_dir():
                continue
            for name in ("timetree", "clock", "randomisations", "stages"):
                shutil.rmtree(lineage_directory / name, ignore_errors=True)
            for name in ("temporal_signal.json", "report.full.html", "report.full.json",
                         "supporting_results.full.zip"):
                (lineage_directory / name).unlink(missing_ok=True)
        for name in ("finish.html", "finish-summary.json", "full.html", "full-summary.json"):
            (output / name).unlink(missing_ok=True)
    stages["fast"] = {"status": "completed", "report": str(profile_report)}
    result = {
        "workflow": "chronoclade",
        "analysis_mode": mode,
        "fingerprint": fingerprint,
        "stages": stages,
        "provenance": inputs["provenance"],
        "context_selections": [],
    }
    _write_json(previous_path, result)
    write_stage_index(output, stages)
    if mode == "fast":
        return result
    query_groups = defaultdict(list)
    for row in queries:
        query_groups[(row["species"], row.get("lineage") or "Unassigned")].append(row)
    requested = set(include_genomes or [])
    selectable = [r for r in context
                  if (r["species"], r.get("lineage") or "Unassigned") in query_groups]
    for identifier in requested:
        matches = [r for r in selectable if identifier in {
            r["sample_id"], str(r.get("source_genome_id", "")), str(r.get("accession", ""))}]
        if len(matches) != 1:
            raise WorkflowError(
                f"Requested context genome {identifier!r} is not uniquely in a query lineage pool"
            )
    samples = []
    selected_records = []
    used = set()
    for key, group in sorted(query_groups.items()):
        candidates = [r for r in context if (r["species"], r.get("lineage") or "Unassigned") == key]
        selected, selection = select_assembly_context(
            group,
            candidates,
            analysis,
            size=context_size,
            nearest_per_query=nearest_per_query,
            include=[
                identifier
                for identifier in include_genomes or []
                if any(
                    identifier
                    in {
                        r["sample_id"],
                        str(r.get("source_genome_id", "")),
                        str(r.get("accession", "")),
                    }
                    for r in candidates
                )
            ],
            seed=seed,
        )
        result["context_selections"].append({"species": key[0], "lineage": key[1], **selection})
        members = group + selected
        selected_records.extend(selected)
        for row in members:
            if row["sample_id"] in used:
                raise WorkflowError("Selected sample occurs in more than one analysis lineage")
            used.add(row["sample_id"])
        samples.extend(materialise_assemblies(members, output=output / "assemblies", client=client))
    _write_json(output / "context_selection.json", result["context_selections"])
    # Preserve selection reasons even when no legacy context manifest was supplied.
    if context_manifest is None and selected_records:
        context_manifest = output / "selected_context.tsv"
        manifest_rows = []
        for key, group in sorted(query_groups.items()):
            pool = [r for r in inputs.get("catalogue_rows", [])
                    if (r["species"], r.get("lineage") or "Unassigned") == key]
            selected = [r for r in selected_records
                        if (r["species"], r.get("lineage") or "Unassigned") == key]
            catalogue_path = output / "inputs" / (content_hash(key)[:16] + "_catalogue.json")
            payload = {"rows": pool, "focal_rows": group,
                       "provenance": {"scope": "Frozen public same-ST metadata catalogue",
                                      "input_provenance": inputs["provenance"]}}
            digest = content_hash(payload)
            if pool:
                _write_json(catalogue_path, {**payload, "snapshot_sha256": digest})
            public_ids = {r.get("source_genome_id") for r in pool}
            for row in selected:
                public = bool(pool and row.get("source_genome_id") in public_ids)
                manifest_rows.append({**row, "source": "pathogenwatch" if public else "provided",
                                      "catalogue_path": str(catalogue_path) if public else "",
                                      "catalogue_sha256": digest if public else ""})
        fields = ["sample_id", "species", "lineage", "country", "collection_date",
                  "source_genome_id", "source", "selection_reason", "catalogue_path",
                  "catalogue_sha256"]
        with context_manifest.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(manifest_rows)
    common = dict(
        output=output / "assembly",
        threads=threads,
        lineage_jobs=lineage_jobs,
        randomisation_jobs=randomisation_jobs,
        randomisations=randomisations,
        temporal_p_value=temporal_p_value,
        min_samples=min_samples,
        seed=seed,
        force=force,
        context_manifest=context_manifest,
        date_randomisation_method=date_randomisation_method,
    )
    corrected = run_workflow(samples, mode="corrected", **common)
    for record in corrected["lineages"]:
        if "outputs" not in record:
            continue
        _snapshot_stage(record, "full")
    stages["full"] = {
        "status": "completed",
        "report": str(_stage_landing(output, "full", corrected["lineages"])),
    }
    _write_json(output / "full-summary.json", corrected)
    result["lineages"] = corrected["lineages"]
    _write_json(previous_path, result)
    write_stage_index(output, stages)
    if mode == "finish":
        finished = run_workflow(samples, mode="full", **common)
        for record in finished["lineages"]:
            if "outputs" not in record:
                continue
            _snapshot_stage(record, "finish")
        stages["finish"] = {
            "status": "completed",
            "report": str(_stage_landing(output, "finish", finished["lineages"])),
        }
        _write_json(output / "finish-summary.json", finished)
        result["lineages"] = finished["lineages"]
    _write_json(previous_path, result)
    write_stage_index(output, stages)
    return result
