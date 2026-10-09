"""Pathogenwatch provider boundary for frozen context preparation."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import fields
from datetime import date
from pathlib import Path

from chronoclade.context import (
    ContextCandidate,
    ContextError,
    annotate_screening_distances,
    filter_candidates,
    parse_ska_distances,
    run_ska_screen,
    select_context,
    stratified_candidate_pool,
    write_audit,
    write_candidate_table,
    write_combined_metadata,
    write_ska_inputs,
)
from chronoclade.pathogenwatch import (
    PathogenwatchClient,
    PathogenwatchError,
    catalogue_summary,
    content_hash,
    deduplicate_catalogue,
    load_api_key,
    load_catalogue,
)
from chronoclade.pathogenwatch_download import download_assemblies


def _decimal(value: str) -> float:
    day = date.fromisoformat(value)
    start = date(day.year, 1, 1)
    return day.year + (day - start).days / (date(day.year + 1, 1, 1) - start).days


def cohort_date(row: dict) -> str:
    """Preserve partial dates; represent interval endpoints as TreeTime date ranges."""
    if not row.get("dated_cohort_eligible"):
        return ""
    if row.get("date_precision") == "interval":
        return f"[{_decimal(row['date_start']):.8f}:{_decimal(row['date_end']):.8f}]"
    return str(row.get("collection_date", ""))


def _candidate(
    row: dict,
    *,
    species: str,
    lineage: str,
    scheme: str,
    st: str,
    snapshot: str,
    catalogue_path: Path,
) -> ContextCandidate:
    names = {field.name for field in fields(ContextCandidate)}
    values = {
        key: str(value) if value is not None else ""
        for key, value in row.items()
        if key in names and not isinstance(value, (dict, list))
    }
    values.update(
        sample_id="PW_" + row["source_genome_id"],
        species=species,
        lineage=lineage,
        mlst_scheme=scheme,
        mlst_st=st,
        collection_date=cohort_date(row),
        source="pathogenwatch",
        source_numeric_id=str(row.get("numeric_source_id", "")),
        sample_accession=row.get("biosample", ""),
        run_accession="|".join(row.get("run_accessions", [])),
        assembly_accession="|".join(row.get("assembly_accessions", [])),
        study_accession="|".join(row.get("study_accessions", [])),
        country_raw=json.dumps(row.get("country_raw", []), sort_keys=True),
        counting_unit=row.get("sample_unit_id", ""),
        pathogenwatch_qc=str(row.get("qc_pass", "")),
        genome_size=str(row.get("source_length") or ""),
        contig_n50=str(row.get("source_n50") or ""),
        catalogue_sha256=snapshot,
        catalogue_path=str(catalogue_path),
    )
    return ContextCandidate(**values)


def _read_crosswalk(path: Path | None) -> list[dict]:
    if path is None:
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or not all(row.get("sample_id") for row in rows):
        raise ContextError("Focal crosswalk requires a sample_id column and nonempty rows")
    return rows


def prepare_pathogenwatch_context(
    focal,
    *,
    species,
    lineage,
    scheme,
    st,
    output,
    cache_dir,
    countries,
    year_from,
    year_to,
    host,
    isolation_source,
    candidate_pool,
    max_context,
    nearest_per_focal,
    seed,
    threads,
    dry_run,
    ska_executable="ska",
    catalogue=None,
    cglin_export=None,
    focal_crosswalk=None,
    refresh_catalogue=False,
):
    """Freeze, annotate, filter, fetch a bounded pool and screen with existing SKA gates."""
    from chronoclade.cglin import (
        CGLINError,
        annotate_catalogue,
        download_cglin_export,
        load_cglin_export,
        resolve_focal_assignments,
    )
    from chronoclade.context_geography import generate_context_geography

    if not focal:
        raise ContextError("No focal samples were supplied for context preparation")
    # Initial supported route is explicit; genome availability is separate from cgLIN support.
    normalized = " ".join(species.replace("_", " ").casefold().split())
    if normalized != "klebsiella pneumoniae" or scheme not in {"klebsiella", "mlst"}:
        raise ContextError(
            "Pathogenwatch context currently supports Klebsiella pneumoniae "
            "with --scheme klebsiella (cgLIN is organism-specific); "
            "use --context-source atb for the legacy route"
        )
    if min(candidate_pool, max_context, nearest_per_focal, threads) < 1:
        raise ContextError("Context pool, selection and thread settings must be positive")
    output = output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    # Never leave a previous selected cohort masquerading as a new successful run.
    for name in ("context_manifest.tsv", "combined_metadata.csv", "screened_candidates.tsv"):
        (output / name).unlink(missing_ok=True)
    frozen = (
        catalogue.expanduser().resolve() if catalogue else output / "pathogenwatch_catalogue.json"
    )
    try:
        if frozen.is_file() and not refresh_catalogue:
            envelope = load_catalogue(frozen)
        elif catalogue:
            raise ContextError("Supplied frozen catalogue does not exist; do not refresh imports")
        else:
            envelope = PathogenwatchClient().freeze_catalogue(frozen, organism_id="573", st=st)
    except PathogenwatchError as error:
        raise ContextError(str(error)) from None
    provenance = envelope["provenance"]
    query = provenance["query"]
    if str(query.get("organismId")) != "573" or query.get("mlst") != [str(st)]:
        raise ContextError("Frozen catalogue query does not match the requested organism/ST")
    rows = envelope["rows"]
    assignments = []
    annotation_capability = "unavailable_without_credentials"
    try:
        if cglin_export:
            assignments = load_cglin_export(cglin_export)
            annotation_capability = "imported_validated_export"
        else:
            key = load_api_key()
            export_dir = output / "cglin" / envelope["snapshot_sha256"]
            export_manifest = export_dir / "cglin_export.json"
            if export_manifest.is_file():
                assignments = load_cglin_export(export_manifest)
                annotation_capability = "frozen_export"
            elif key:
                exported = download_cglin_export(rows, export_dir, api_key=key)
                assignments = exported["assignments"]
                annotation_capability = "downloaded_export"
    except CGLINError as error:
        raise ContextError(str(error)) from None
    rows = annotate_catalogue(rows, assignments)
    # Public composition is deduplicated before dated/selection/focal exclusions.
    public, duplicate_audit = deduplicate_catalogue(
        [row for row in rows if row.get("qc_pass") is True]
    )
    crosswalk = _read_crosswalk(focal_crosswalk)
    focal_rows = [
        {
            "sample_id": sample.sample_id,
            "biosample": sample.sample_id if sample.sample_id.upper().startswith("SAM") else "",
            "country": sample.location,
            "species": species,
            "lineage": lineage,
            **next((r for r in crosswalk if r["sample_id"] == sample.sample_id), {}),
        }
        for sample in focal
    ]
    focal_annotations, focal_audit = resolve_focal_assignments(
        focal_rows, rows, crosswalk=crosswalk or None
    )
    focal_aliases = {sample.sample_id for sample in focal}
    for row in focal_annotations:
        for key in (
            "source_genome_id",
            "cglin_matched_source_genome_id",
            "biosample",
            "sample_accession",
            "assembly_accession",
        ):
            if row.get(key):
                focal_aliases.add(row[key])
    eligible, focal_duplicate_audit = deduplicate_catalogue(rows, focal_aliases=focal_aliases)
    qc = [row for row in eligible if row.get("qc_pass") is True]
    dated = [row for row in qc if row.get("dated_cohort_eligible")]
    annotated_path = output / "context_catalogue.json"
    annotated_payload = {"rows": public, "focal_rows": focal_annotations, "provenance": provenance}
    snapshot = content_hash(annotated_payload)
    annotated_payload["snapshot_sha256"] = snapshot
    write_audit(annotated_path, annotated_payload)
    candidates = [
        _candidate(
            row,
            species=species,
            lineage=lineage,
            scheme=scheme,
            st=st,
            snapshot=snapshot,
            catalogue_path=annotated_path,
        )
        for row in dated
    ]
    filtered = filter_candidates(
        candidates,
        countries=countries,
        year_from=year_from,
        year_to=year_to,
        host=host,
        isolation_source=isolation_source,
    )
    pool = stratified_candidate_pool(filtered, limit=candidate_pool, seed=seed)
    write_candidate_table(output / "candidate_pool.tsv", pool)
    settings = dict(
        countries=countries or [],
        year_from=year_from,
        year_to=year_to,
        host=host,
        isolation_source=isolation_source,
        seed=seed,
        candidate_pool=candidate_pool,
        max_context=max_context,
        nearest_per_focal=nearest_per_focal,
    )
    fingerprint = content_hash(
        {
            "snapshot": snapshot,
            "settings": settings,
            "focal_assemblies": {
                s.sample_id: hashlib.sha256(s.assembly.read_bytes()).hexdigest() for s in focal
            },
        }
    )
    audit = dict(
        context_source="pathogenwatch",
        species=species,
        lineage=lineage,
        mlst_scheme=scheme,
        mlst_st=st,
        context_metadata_snapshot=provenance,
        catalogue_sha256=snapshot,
        selection_fingerprint=fingerprint,
        filters=settings,
        seed=seed,
        same_st_accessions=len(rows),
        metadata_completeness=catalogue_summary(rows),
        deduplication=duplicate_audit,
        focal_exclusions=focal_duplicate_audit,
        focal_cglin=focal_audit,
        annotation_capability=annotation_capability,
        qc_pass_candidates=len(qc),
        dated_hq_metadata_candidates=len(dated),
        metadata_filtered_candidates=len(filtered),
        candidate_pool=len(pool),
        candidate_pool_limit=candidate_pool,
        dry_run=dry_run,
        stage_losses={
            "qc_or_unknown": len(eligible) - len(qc),
            "undated_or_invalid": len(qc) - len(dated),
            "explicit_filters": len(dated) - len(filtered),
            "bounded_pool": len(filtered) - len(pool),
        },
    )
    scope = dict(
        description=f"Public Pathogenwatch {species} ST{st} same-ST catalogue",
        snapshot=provenance["retrieved_at"],
        filters="QC pass; accession-deduplicated; "
        "includes undated records before country/year balancing and SKA selection",
    )

    def figures(selected=()):
        result = generate_context_geography(
            public,
            output / "context_geography",
            selected_source_ids=selected,
            focal_rows=focal_annotations,
            scope=scope,
        )
        audit["geography"] = {
            key: value for key, value in result.items() if key in ("summaries", "outputs")
        }

    figures()
    write_audit(output / "context_selection.json", audit)
    if dry_run:
        return audit
    if not pool:
        raise ContextError(
            "No QC-passing dated same-ST candidates remain; catalogue geography "
            "and exclusion audit are available"
        )
    key = load_api_key()
    if not key:
        raise ContextError(
            "Pathogenwatch FASTA download requires PATHOGENWATCH_API_KEY "
            "or protected ~/.config/chronoclade/pathogenwatch.json"
        )
    by_source = {row["source_genome_id"]: row for row in qc}
    assemblies, ledger = download_assemblies(
        [by_source[c.source_genome_id] for c in pool],
        output=output / "assemblies",
        cache_dir=cache_dir.expanduser().resolve(),
        api_key=key,
        workers=min(threads, 4),
    )
    write_audit(output / "download_ledger.json", {"records": ledger})
    audit["download_ledger"] = ledger
    downloaded = []
    for candidate in pool:
        if candidate.source_genome_id in assemblies:
            candidate.assembly = str(assemblies[candidate.source_genome_id])
            candidate.assembly_sha256 = hashlib.sha256(
                Path(candidate.assembly).read_bytes()
            ).hexdigest()
            downloaded.append(candidate)
    audit["downloaded_candidates"] = len(downloaded)
    audit["missing_downloads"] = [r["source_genome_id"] for r in ledger if r["status"] == "failed"]
    audit["stage_losses"]["download_failure"] = len(pool) - len(downloaded)
    write_audit(output / "context_selection.json", audit)
    if not downloaded:
        raise ContextError("No validated Pathogenwatch FASTAs downloaded; see failure ledger")
    inputs = write_ska_inputs(output / "ska_inputs.tsv", focal, downloaded)
    distance_path = run_ska_screen(
        inputs=inputs, output_dir=output, threads=threads, executable=ska_executable
    )
    distances = parse_ska_distances(distance_path)
    selected, selection = select_context(
        downloaded,
        [s.sample_id for s in focal],
        distances,
        max_context=max_context,
        nearest_per_focal=nearest_per_focal,
        seed=seed,
    )
    audit.update(selection)
    audit["stage_losses"]["missing_ska_comparisons"] = (
        len(downloaded) - selection["screened_candidates"]
    )
    audit["stage_losses"]["selection_limit"] = selection["screened_candidates"] - len(selected)
    if not selected:
        write_audit(output / "context_selection.json", audit)
        raise ContextError("SKA screening yielded no contextual genomes; see selection audit")
    write_candidate_table(
        output / "screened_candidates.tsv",
        annotate_screening_distances(downloaded, [s.sample_id for s in focal], distances),
    )
    manifest = write_candidate_table(output / "context_manifest.tsv", selected)
    combined = write_combined_metadata(output / "combined_metadata.csv", focal, selected)
    figures([c.source_genome_id for c in selected])
    audit["outputs"] = dict(
        manifest=str(manifest),
        combined_metadata=str(combined),
        ska_distances=str(distance_path),
        catalogue=str(annotated_path),
    )
    write_audit(output / "context_selection.json", audit)
    return audit
