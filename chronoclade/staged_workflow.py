"""Profile-first stages; assemblies are acquired only after context selection."""

from __future__ import annotations

import csv
import json
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Any

from chronoclade.errors import WorkflowError
from chronoclade.metadata import slugify_lineage
from chronoclade.pathogenwatch import content_hash
from chronoclade.sample_labels import sample_labels
from chronoclade.selections import selected_context_records


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _stage_landing(output: Path, name: str, records: list[dict]) -> Path:
    from html import escape
    from chronoclade.report_components.styles import report_styles

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
        f"<title>ChronoClade {name}</title><style>{report_styles()}</style><body><main class=\"shell\">"
        '<header class="identity"><div class="identity-mark">ChronoClade · ANALYSIS REPORTS</div><div class="identity-copy">'
        f"<h1>{escape(name.title())} results</h1></div></header>"
        '<section class="stage"><div class="stage-body"><a href="index.html">All completed stages</a><ul>'
        + "".join(links)
        + "</ul></div></section></main></body></html>",
        encoding="utf-8",
    )
    return path


def _snapshot_stage(record: dict, stage: str) -> None:
    """Freeze reader-facing assets so later analysis cannot change an earlier report."""
    source_report = Path(record["outputs"]["html_report"])
    directory = source_report.parent
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

    # Writers may use a stage-specific filename. The advertised source, rather
    # than a stale report.html left by another mode, supplies the snapshot entry.
    entry = destination / "report.html"
    copied_report = destination / source_report.name
    if copied_report != entry:
        shutil.copy2(copied_report, entry)
    snapshot_record = relocated(record)
    snapshot_record["outputs"]["html_report"] = str(entry)
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
    tree_limit: int = 80,
    lin_min_context: int = 20,
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
    from chronoclade.profile_inputs import materialise_assemblies, resolve_profile_inputs
    from chronoclade.profile_report import write_stage_index
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
        lin_min_context=lin_min_context,
        seed=seed,
    )
    from chronoclade.fast_workflow import route_datasets, run_fast_datasets

    inputs = route_datasets(inputs)
    queries, context = inputs["queries"], inputs["context"]
    analysis, profile_report, fast_datasets = run_fast_datasets(
        inputs, output / "fast", seed=seed,
        bootstrap_replicates=bootstrap_replicates,
        distance_threshold=distance_threshold, tree_limit=tree_limit,
        context_size=context_size, nearest_per_query=nearest_per_query,
        include_genomes=include_genomes,
    )
    query_groups = defaultdict(list)
    for row in queries:
        query_groups[(row["species"], row.get("lineage") or "Unassigned")].append(row)
    selections = analysis.get("context_selections", [])
    selections_by_dataset = {}
    selected_by_dataset = {}
    for selection in selections:
        key = (selection["species"], selection["lineage"])
        if key in selections_by_dataset or key not in query_groups:
            raise WorkflowError("Shared context selection has duplicate or unknown datasets")
        group = query_groups[key]
        candidates = [r for r in context if (r["species"], r.get("lineage") or "Unassigned") == key]
        selected_by_dataset[key] = selected_context_records(selection, group, candidates)
        selections_by_dataset[key] = selection
    if set(selections_by_dataset) != set(query_groups):
        raise WorkflowError("Shared context selection is missing a query dataset")
    fingerprint = content_hash(
        {
            "resolved_records_sha256": inputs["provenance"].get(
                "resolved_records_sha256", content_hash({"queries": queries, "context": context})
            ),
            "context_size": context_size,
            "shared_selection": [{"species": row["species"], "lineage": row["lineage"],
                                  "selected_sample_ids": row["selected_sample_ids"]}
                                 for row in selections],
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
        "context_selections": selections,
        "fast_datasets": fast_datasets,
    }
    _write_json(output / "context_selection.json", selections)
    _write_json(previous_path, result)
    write_stage_index(output, stages)
    if mode == "fast":
        return result
    samples = []
    selected_records = []
    used = set()
    for key, group in sorted(query_groups.items()):
        selected = selected_by_dataset[key]
        members = group + selected
        labels_path = output / "assembly" / slugify_lineage(*key) / "sample_labels.json"
        _write_json(labels_path, sample_labels(members))
        selected_records.extend(selected)
        for row in members:
            if row["sample_id"] in used:
                raise WorkflowError("Selected sample occurs in more than one analysis lineage")
            used.add(row["sample_id"])
        samples.extend(materialise_assemblies(members, output=output / "assemblies", client=client))
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
