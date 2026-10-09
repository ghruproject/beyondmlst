"""Reproducible acquisition and down-selection of public context genomes."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import random
import re
import subprocess
from collections import defaultdict
from dataclasses import asdict, dataclass, fields
from datetime import datetime
from pathlib import Path

import pyarrow.parquet as pq

from chronoclade.metadata import Sample


DEFAULT_CONTEXT_METADATA = Path(__file__).parent / "data" / "atb_context_202505.parquet"


class ContextError(RuntimeError):
    """Raised when contextual-genome preparation cannot be completed safely."""


@dataclass
class ContextCandidate:
    """One same-ST public genome considered for contextual analysis."""

    sample_id: str
    species: str
    lineage: str
    mlst_scheme: str
    mlst_st: str
    collection_date: str = ""
    country: str = ""
    host: str = ""
    isolation_source: str = ""
    hq_filter: str = ""
    aws_url: str = ""
    completeness: str = ""
    contamination: str = ""
    genome_size: str = ""
    contig_n50: str = ""
    assembly: str = ""
    nearest_focal: str = ""
    min_ska_distance: float | None = None
    mismatch_proportion: float | None = None
    selection_reason: str = ""

    source: str = "atb"
    source_genome_id: str = ""
    source_numeric_id: str = ""
    sample_accession: str = ""
    run_accession: str = ""
    assembly_accession: str = ""
    study_accession: str = ""
    counting_unit: str = ""
    country_raw: str = ""
    date_start: str = ""
    date_end: str = ""
    date_precision: str = ""
    pathogenwatch_qc: str = ""
    assembly_sha256: str = ""
    catalogue_sha256: str = ""
    catalogue_path: str = ""
    cglin_raw: str = ""
    cglin_status: str = ""
    cglin_resolved_depth: str = ""
    cgst: str = ""
    cglin_provisional: str = ""
    cglin_source: str = ""
    cglin_scheme: str = ""
    cglin_scheme_version: str = ""
    cglin_export_sha256: str = ""
    cglin_retrieved_at: str = ""
    cglin_group_5: str = ""
    cglin_group_6: str = ""
    cglin_group_7: str = ""
    cglin_status_5: str = ""
    cglin_status_6: str = ""
    cglin_status_7: str = ""

    @property
    def year(self) -> str:
        if self.date_start:
            return self.date_start[:4]
        match = re.match(r"^\[?(\d{4})", self.collection_date)
        return match.group(1) if match else "Unknown"

    @property
    def stratum(self) -> tuple[str, str]:
        return (self.country or "Unknown", self.year)


MANIFEST_FIELDS = [
    "sample_id",
    "species",
    "lineage",
    "mlst_scheme",
    "mlst_st",
    "collection_date",
    "country",
    "host",
    "isolation_source",
    "hq_filter",
    "completeness",
    "contamination",
    "genome_size",
    "contig_n50",
    "assembly",
    "nearest_focal",
    "min_ska_distance",
    "mismatch_proportion",
    "selection_reason",
    "aws_url",
]

MANIFEST_FIELDS.extend(field.name for field in fields(ContextCandidate)
                       if field.name not in MANIFEST_FIELDS)


def _run_capture(
    command: list[str], *, log: Path | None = None
) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if log is not None:
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(
            "COMMAND\n"
            + " ".join(command)
            + "\n\nSTDOUT\n"
            + completed.stdout
            + "\nSTDERR\n"
            + completed.stderr,
            encoding="utf-8",
        )
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip() or "no diagnostic output"
        raise ContextError(f"Command failed ({command[0]}): {detail}")
    return completed


def atbfetcher_version(executable: str = "atbfetcher") -> str:
    """Return the installed atbfetcher version string."""

    completed = _run_capture([executable, "--version"])
    return completed.stdout.strip() or completed.stderr.strip() or "unknown"


def resolve_context_metadata(requested: Path | None = None) -> Path:
    """Resolve the bundled or explicitly supplied compact ATB metadata table."""

    result = (requested or DEFAULT_CONTEXT_METADATA).expanduser().resolve()
    if not result.is_file():
        raise ContextError(f"Context metadata table not found: {result}")
    return result


def context_metadata_provenance(path: Path) -> dict[str, object]:
    """Return an auditable description of one compact context snapshot."""

    manifest_path = path.with_suffix(".json")
    manifest: dict[str, object] = {}
    if manifest_path.is_file():
        try:
            loaded = json.loads(manifest_path.read_text(encoding="utf-8"))
            if isinstance(loaded, dict):
                manifest = loaded
        except (OSError, json.JSONDecodeError):
            manifest = {}
    with path.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    try:
        rows = pq.ParquetFile(path).metadata.num_rows
    except (OSError, ValueError) as error:
        raise ContextError(f"Could not read context metadata table {path}: {error}") from error
    return {
        "path": str(path),
        "name": path.name,
        "size_bytes": path.stat().st_size,
        "sha256": digest,
        "rows": rows,
        "release": manifest.get("atb_release", "unknown"),
        "species": manifest.get("species", "unknown"),
        "mlst_scheme_count": manifest.get("mlst_scheme_count", "unknown"),
        "usable_collection_date_rows": manifest.get("usable_collection_date_rows", "unknown"),
        "source_manifest": str(manifest_path) if manifest_path.is_file() else None,
    }


def normalise_collection_date(value: str) -> str:
    """Return a ChronoClade-compatible public collection date or an empty string."""

    value = value.strip()
    if not value:
        return ""
    value = value.split("T", 1)[0].split(" ", 1)[0]
    patterns = (
        (r"^(\d{4})$", "%Y"),
        (r"^(\d{4})-(\d{2})$", "%Y-%m"),
        (r"^(\d{4})-(\d{2})-(\d{2})$", "%Y-%m-%d"),
    )
    for pattern, date_format in patterns:
        if not re.fullmatch(pattern, value):
            continue
        try:
            parsed = datetime.strptime(value, date_format)
        except ValueError:
            return ""
        if 1800 <= parsed.year <= 2200:
            return value
    return ""


def _normalise_species(value: str) -> str:
    """Normalise GTDB suffixes and underscore-separated user metadata."""

    without_gtdb_suffix = re.sub(r"_[A-Z]\b", "", value)
    return " ".join(without_gtdb_suffix.replace("_", " ").casefold().split())


def _text(value: object) -> str:
    return "" if value is None else str(value)


def load_same_st_candidates(
    metadata_table: Path,
    *,
    species: str,
    lineage: str,
    scheme: str,
    st: str,
) -> tuple[list[str], list[ContextCandidate]]:
    """Load all and dated same-ST candidates from the compact ATB snapshot."""

    try:
        table = pq.read_table(
            metadata_table,
            filters=[
                ("mlst_scheme", "=", scheme),
                ("mlst_st", "=", str(st)),
            ],
        )
    except (OSError, ValueError) as error:
        raise ContextError(
            f"Could not query context metadata table {metadata_table}: {error}"
        ) from error

    required = {
        "sample_id",
        "species",
        "mlst_scheme",
        "mlst_st",
        "collection_date",
        "country",
        "host",
        "isolation_source",
        "hq_filter",
        "completeness",
        "contamination",
        "genome_size",
        "contig_n50",
        "aws_url",
    }
    missing = required - set(table.column_names)
    if missing:
        raise ContextError(
            "Context metadata table is missing required columns: " + ", ".join(sorted(missing))
        )

    rows = [
        row
        for row in table.to_pylist()
        if _normalise_species(_text(row.get("species"))) == _normalise_species(species)
        and _text(row.get("hq_filter")) == "PASS"
    ]
    accessions = sorted({_text(row.get("sample_id")) for row in rows if row.get("sample_id")})
    if not accessions:
        raise ContextError(
            f"No high-quality {species} ST{st} records were found in {metadata_table.name}"
        )
    candidates: list[ContextCandidate] = []
    for row in rows:
        collection_date = normalise_collection_date(_text(row.get("collection_date")))
        if not collection_date:
            continue
        candidates.append(
            ContextCandidate(
                sample_id=_text(row.get("sample_id")),
                species=species,
                lineage=lineage,
                mlst_scheme=scheme,
                mlst_st=str(st),
                collection_date=collection_date,
                country=_text(row.get("country")),
                host=_text(row.get("host")),
                isolation_source=_text(row.get("isolation_source")),
                hq_filter=_text(row.get("hq_filter")),
                aws_url=_text(row.get("aws_url")),
                completeness=_text(row.get("completeness")),
                contamination=_text(row.get("contamination")),
                genome_size=_text(row.get("genome_size")),
                contig_n50=_text(row.get("contig_n50")),
            )
        )
    candidates.sort(key=lambda candidate: candidate.sample_id)
    return accessions, candidates


def filter_candidates(
    candidates: list[ContextCandidate],
    *,
    countries: list[str] | None = None,
    year_from: int | None = None,
    year_to: int | None = None,
    host: str | None = None,
    isolation_source: str | None = None,
) -> list[ContextCandidate]:
    """Apply explicit epidemiological filters to same-ST candidates."""

    country_values = [value.casefold() for value in countries or []]
    host_value = host.casefold() if host else ""
    source_value = isolation_source.casefold() if isolation_source else ""
    selected: list[ContextCandidate] = []
    for candidate in candidates:
        year = int(candidate.year)
        if country_values and not any(
            candidate.country.casefold().startswith(value) for value in country_values
        ):
            continue
        if year_from is not None and year < year_from:
            continue
        if year_to is not None and year > year_to:
            continue
        if host_value and host_value not in candidate.host.casefold():
            continue
        if source_value and source_value not in candidate.isolation_source.casefold():
            continue
        selected.append(candidate)
    return selected


def stratified_candidate_pool(
    candidates: list[ContextCandidate], *, limit: int, seed: int
) -> list[ContextCandidate]:
    """Select a reproducible country/year-balanced candidate pool before download."""

    if limit < 1:
        raise ValueError("candidate pool limit must be at least 1")
    if len(candidates) <= limit:
        return sorted(candidates, key=lambda candidate: candidate.sample_id)

    rng = random.Random(seed)
    strata: dict[tuple[str, str], list[ContextCandidate]] = defaultdict(list)
    for candidate in candidates:
        strata[candidate.stratum].append(candidate)
    for members in strata.values():
        members.sort(key=lambda candidate: candidate.sample_id)
        rng.shuffle(members)
    keys = sorted(strata)
    rng.shuffle(keys)

    result: list[ContextCandidate] = []
    while len(result) < limit:
        progressed = False
        for key in keys:
            if strata[key]:
                result.append(strata[key].pop())
                progressed = True
                if len(result) == limit:
                    break
        if not progressed:
            break
    return result


def write_accessions(path: Path, candidates: list[ContextCandidate]) -> Path:
    path.write_text(
        "".join(f"{candidate.sample_id}\n" for candidate in candidates), encoding="utf-8"
    )
    return path


def download_candidates(
    candidates: list[ContextCandidate],
    *,
    accessions_path: Path,
    output: Path,
    cache_dir: Path,
    source: str,
    threads: int,
    executable: str = "atbfetcher",
    log: Path | None = None,
) -> None:
    """Fetch a frozen candidate accession list with atbfetcher."""

    if not candidates:
        raise ContextError("No context candidates remain to download")
    output.mkdir(parents=True, exist_ok=True)
    command = [
        executable,
        "accessions",
        str(accessions_path),
        "--output",
        str(output),
        "--cache-dir",
        str(cache_dir),
        "--source",
        source,
        "--threads",
        str(threads),
    ]
    _run_capture(command, log=log)


def _assembly_identifier(path: Path) -> str:
    name = path.name
    if name.endswith(".gz"):
        name = name[:-3]
    for extension in (".fasta", ".fna", ".fa"):
        if name.endswith(extension):
            return name[: -len(extension)]
    return ""


def attach_downloaded_assemblies(
    candidates: list[ContextCandidate], output: Path
) -> tuple[list[ContextCandidate], list[str]]:
    """Attach downloaded FASTA paths and report missing accessions."""

    paths: dict[str, Path] = {}
    for path in output.rglob("*"):
        if not path.is_file():
            continue
        identifier = _assembly_identifier(path)
        if identifier:
            paths.setdefault(identifier, path.resolve())
    found: list[ContextCandidate] = []
    missing: list[str] = []
    for candidate in candidates:
        assembly = paths.get(candidate.sample_id)
        if assembly is None:
            missing.append(candidate.sample_id)
            continue
        candidate.assembly = str(assembly)
        found.append(candidate)
    return found, missing


def write_ska_inputs(path: Path, focal: list[Sample], candidates: list[ContextCandidate]) -> Path:
    rows = [f"{sample.sample_id}\t{sample.assembly}\n" for sample in focal]
    rows.extend(f"{candidate.sample_id}\t{candidate.assembly}\n" for candidate in candidates)
    path.write_text("".join(rows), encoding="utf-8")
    return path


def run_ska_screen(
    *, inputs: Path, output_dir: Path, threads: int, executable: str = "ska"
) -> Path:
    """Build one SKA index and calculate pairwise screening distances."""

    prefix = output_dir / "context_screen"
    skf = output_dir / "context_screen.skf"
    distances = output_dir / "context_distances.tsv"
    logs = output_dir / "logs"
    _run_capture(
        [executable, "build", "-f", str(inputs), "-o", str(prefix), "--threads", str(threads)],
        log=logs / "ska_build.log",
    )
    if not skf.is_file():
        raise ContextError(f"SKA did not create its expected index: {skf}")
    _run_capture(
        [
            executable,
            "distance",
            "--threads",
            str(threads),
            "-o",
            str(distances),
            str(skf),
        ],
        log=logs / "ska_distance.log",
    )
    if not distances.is_file():
        raise ContextError(f"SKA did not create its distance table: {distances}")
    return distances


def parse_ska_distances(path: Path) -> dict[tuple[str, str], tuple[float, float]]:
    """Parse SKA2 pairwise SNP and mismatch-proportion output."""

    result: dict[tuple[str, str], tuple[float, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Sample1", "Sample2", "Distance", "Mismatches (proportion)"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ContextError(f"Unrecognised SKA distance format: {path}")
        for row in reader:
            first = str(row["Sample1"])
            second = str(row["Sample2"])
            value = (float(row["Distance"]), float(row["Mismatches (proportion)"]))
            result[(first, second)] = value
            result[(second, first)] = value
    return result


def annotate_screening_distances(
    candidates: list[ContextCandidate],
    focal_ids: list[str],
    distances: dict[tuple[str, str], tuple[float, float]],
) -> list[ContextCandidate]:
    """Annotate every candidate with its nearest focal SKA screening distance."""

    annotated: list[ContextCandidate] = []
    for candidate in candidates:
        comparisons = [
            (distances[(candidate.sample_id, focal_id)], focal_id)
            for focal_id in focal_ids
            if (candidate.sample_id, focal_id) in distances
        ]
        if not comparisons:
            continue
        (distance, mismatch), focal_id = min(
            comparisons, key=lambda value: (value[0][0], value[0][1], value[1])
        )
        candidate.min_ska_distance = distance
        candidate.mismatch_proportion = mismatch
        candidate.nearest_focal = focal_id
        annotated.append(candidate)
    return annotated


def select_context(
    candidates: list[ContextCandidate],
    focal_ids: list[str],
    distances: dict[tuple[str, str], tuple[float, float]],
    *,
    max_context: int,
    nearest_per_focal: int,
    seed: int,
) -> tuple[list[ContextCandidate], dict[str, object]]:
    """Retain nearest neighbours, then fill with a balanced background."""

    if max_context < 1 or nearest_per_focal < 1:
        raise ValueError("max_context and nearest_per_focal must be at least 1")
    candidates = annotate_screening_distances(candidates, focal_ids, distances)
    by_id = {candidate.sample_id: candidate for candidate in candidates}
    nearest_for: dict[str, list[str]] = defaultdict(list)
    support: dict[str, int] = defaultdict(int)
    for focal_id in focal_ids:
        ranked = sorted(
            (
                (distances[(candidate.sample_id, focal_id)][0], candidate.sample_id)
                for candidate in candidates
                if (candidate.sample_id, focal_id) in distances
            ),
            key=lambda value: (value[0], value[1]),
        )[:nearest_per_focal]
        for _, candidate_id in ranked:
            nearest_for[candidate_id].append(focal_id)
            support[candidate_id] += 1

    selected_ids = sorted(
        nearest_for,
        key=lambda candidate_id: (
            -support[candidate_id],
            (
                by_id[candidate_id].min_ska_distance
                if by_id[candidate_id].min_ska_distance is not None
                else float("inf")
            ),
            candidate_id,
        ),
    )[:max_context]
    selected_set = set(selected_ids)
    for candidate_id in selected_ids:
        by_id[candidate_id].selection_reason = "nearest_to=" + "|".join(
            sorted(nearest_for[candidate_id])
        )

    remaining: dict[tuple[str, str], list[ContextCandidate]] = defaultdict(list)
    for candidate in candidates:
        if candidate.sample_id not in selected_set:
            remaining[candidate.stratum].append(candidate)
    for members in remaining.values():
        members.sort(
            key=lambda candidate: (
                (
                    candidate.min_ska_distance
                    if candidate.min_ska_distance is not None
                    else float("inf")
                ),
                candidate.sample_id,
            ),
            reverse=True,
        )
    keys = sorted(remaining)
    random.Random(seed).shuffle(keys)
    while len(selected_ids) < max_context:
        progressed = False
        for key in keys:
            if remaining[key]:
                candidate = remaining[key].pop()
                candidate.selection_reason = (
                    f"stratified_background={candidate.country or 'Unknown'}|{candidate.year}"
                )
                selected_ids.append(candidate.sample_id)
                selected_set.add(candidate.sample_id)
                progressed = True
                if len(selected_ids) == max_context:
                    break
        if not progressed:
            break

    selected = [by_id[candidate_id] for candidate_id in selected_ids]
    selected.sort(
        key=lambda candidate: (
            (
                candidate.min_ska_distance
                if candidate.min_ska_distance is not None
                else float("inf")
            ),
            candidate.sample_id,
        )
    )
    covered = {
        focal_id for candidate in selected for focal_id in nearest_for.get(candidate.sample_id, [])
    }
    audit = {
        "screened_candidates": len(candidates),
        "selected_contexts": len(selected),
        "max_context": max_context,
        "nearest_per_focal": nearest_per_focal,
        "focal_samples": len(focal_ids),
        "focal_samples_with_selected_neighbour": len(covered),
        "focal_neighbour_coverage": len(covered) / len(focal_ids) if focal_ids else 0.0,
    }
    return selected, audit


def write_candidate_table(path: Path, candidates: list[ContextCandidate]) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=MANIFEST_FIELDS,
            extrasaction="ignore",
            delimiter="\t",
        )
        writer.writeheader()
        for candidate in candidates:
            row = asdict(candidate)
            for key in ("assembly", "catalogue_path"):
                if row.get(key) and Path(row[key]).is_absolute():
                    row[key] = os.path.relpath(row[key], path.parent)
            writer.writerow(row)
    return path


def write_combined_metadata(
    path: Path, focal: list[Sample], selected: list[ContextCandidate]
) -> Path:
    fields = [
        "sample_id",
        "assembly",
        "collection_date",
        "location",
        "species",
        "lineage",
        "origin",
        "is_reference",
        "patient_id",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for sample in focal:
            writer.writerow(
                {
                    "sample_id": sample.sample_id,
                    "assembly": os.path.relpath(sample.assembly, path.parent),
                    "collection_date": sample.collection_date,
                    "location": sample.location,
                    "species": sample.species,
                    "lineage": sample.lineage,
                    "origin": sample.origin,
                    "is_reference": str(sample.is_reference).lower(),
                    "patient_id": sample.patient_id,
                }
            )
        for candidate in selected:
            writer.writerow(
                {
                    "sample_id": candidate.sample_id,
                    "assembly": os.path.relpath(candidate.assembly, path.parent),
                    "collection_date": candidate.collection_date,
                    "location": candidate.country or "Public_context",
                    "species": candidate.species,
                    "lineage": candidate.lineage,
                    "origin": "context",
                    "is_reference": "false",
                    "patient_id": "",
                }
            )
    return path


def write_audit(path: Path, audit: dict[str, object]) -> Path:
    path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    return path


def prepare_context(
    focal: list[Sample],
    *,
    species: str,
    lineage: str,
    scheme: str,
    st: str,
    output: Path,
    cache_dir: Path,
    metadata_table: Path | None,
    countries: list[str] | None,
    year_from: int | None,
    year_to: int | None,
    host: str | None,
    isolation_source: str | None,
    candidate_pool: int,
    max_context: int,
    nearest_per_focal: int,
    seed: int,
    threads: int,
    source: str,
    dry_run: bool,
    atbfetcher_executable: str = "atbfetcher",
    ska_executable: str = "ska",
    context_source: str = "pathogenwatch",
    catalogue: Path | None = None,
    cglin_export: Path | None = None,
    focal_crosswalk: Path | None = None,
    refresh_catalogue: bool = False,
) -> dict[str, object]:
    """Prepare a frozen, distance-screened public context set for one lineage."""

    if context_source == "pathogenwatch":
        from chronoclade.pathogenwatch_context import prepare_pathogenwatch_context
        if metadata_table is not None or source != "auto":
            raise ContextError("--metadata-table and --source aws/osf require --context-source atb")
        from chronoclade.cglin import CGLINError
        from chronoclade.pathogenwatch import PathogenwatchError
        try:
            return prepare_pathogenwatch_context(
                focal, species=species, lineage=lineage, scheme=scheme, st=st, output=output,
                cache_dir=cache_dir, countries=countries, year_from=year_from, year_to=year_to,
                host=host, isolation_source=isolation_source, candidate_pool=candidate_pool,
                max_context=max_context, nearest_per_focal=nearest_per_focal, seed=seed,
                threads=threads, dry_run=dry_run, ska_executable=ska_executable,
                catalogue=catalogue, cglin_export=cglin_export, focal_crosswalk=focal_crosswalk,
                refresh_catalogue=refresh_catalogue)
        except (CGLINError, PathogenwatchError) as error:
            raise ContextError(str(error)) from None

    if context_source != "atb":
        raise ContextError("--context-source must be pathogenwatch or atb")
    if (catalogue is not None or cglin_export is not None or focal_crosswalk is not None
            or refresh_catalogue):
        raise ContextError("Pathogenwatch catalogue/cgLIN options cannot be used with ATB")
    if not focal:
        raise ContextError("No focal samples were supplied for context preparation")
    output = output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    logs = output / "logs"
    cache_dir = cache_dir.expanduser().resolve()
    resolved_metadata = resolve_context_metadata(metadata_table)
    tool_version = atbfetcher_version(atbfetcher_executable)

    raw_accessions, metadata_candidates = load_same_st_candidates(
        resolved_metadata,
        species=species,
        lineage=lineage,
        scheme=scheme,
        st=st,
    )
    focal_ids = {sample.sample_id for sample in focal}
    raw_accessions = [accession for accession in raw_accessions if accession not in focal_ids]
    metadata_candidates = [
        candidate for candidate in metadata_candidates if candidate.sample_id not in focal_ids
    ]
    (output / "same_st_accessions.txt").write_text(
        "".join(f"{accession}\n" for accession in raw_accessions), encoding="utf-8"
    )

    filtered = filter_candidates(
        metadata_candidates,
        countries=countries,
        year_from=year_from,
        year_to=year_to,
        host=host,
        isolation_source=isolation_source,
    )
    if not filtered:
        raise ContextError("No dated same-ST candidates passed the requested metadata filters")
    pool = stratified_candidate_pool(filtered, limit=candidate_pool, seed=seed)
    write_candidate_table(output / "candidate_pool.tsv", pool)
    accessions_path = write_accessions(output / "candidate_accessions.txt", pool)

    audit: dict[str, object] = {
        "species": species,
        "lineage": lineage,
        "mlst_scheme": scheme,
        "mlst_st": st,
        "context_source": "atb",
        "atbfetcher_version": tool_version,
        "metadata_discovery": "bundled_or_supplied_parquet",
        "context_metadata_snapshot": context_metadata_provenance(resolved_metadata),
        "filters": {
            "countries": countries or [],
            "year_from": year_from,
            "year_to": year_to,
            "host": host,
            "isolation_source": isolation_source,
        },
        "seed": seed,
        "same_st_accessions": len(raw_accessions),
        "dated_hq_metadata_candidates": len(metadata_candidates),
        "metadata_filtered_candidates": len(filtered),
        "candidate_pool_limit": candidate_pool,
        "candidate_pool": len(pool),
        "dry_run": dry_run,
    }
    if dry_run:
        audit["next_step"] = "Rerun without --dry-run to fetch, screen and select assemblies."
        write_audit(output / "context_selection.json", audit)
        return audit

    assemblies = output / "assemblies"
    download_candidates(
        pool,
        accessions_path=accessions_path,
        output=assemblies,
        cache_dir=cache_dir,
        source=source,
        threads=threads,
        executable=atbfetcher_executable,
        log=logs / "atbfetcher_download.log",
    )
    downloaded, missing = attach_downloaded_assemblies(pool, assemblies)
    if not downloaded:
        raise ContextError("atbfetcher completed but no candidate FASTA assemblies were found")

    ska_inputs = write_ska_inputs(output / "ska_inputs.tsv", focal, downloaded)
    distance_path = run_ska_screen(
        inputs=ska_inputs,
        output_dir=output,
        threads=threads,
        executable=ska_executable,
    )
    distance_values = parse_ska_distances(distance_path)
    selected, selection_audit = select_context(
        downloaded,
        [sample.sample_id for sample in focal],
        distance_values,
        max_context=max_context,
        nearest_per_focal=nearest_per_focal,
        seed=seed,
    )
    if not selected:
        raise ContextError("SKA screening did not yield any usable contextual genomes")

    screened = annotate_screening_distances(
        downloaded, [sample.sample_id for sample in focal], distance_values
    )
    write_candidate_table(output / "screened_candidates.tsv", screened)
    manifest = write_candidate_table(output / "context_manifest.tsv", selected)
    combined = write_combined_metadata(output / "combined_metadata.csv", focal, selected)
    audit.update(selection_audit)
    audit.update(
        {
            "downloaded_candidates": len(downloaded),
            "missing_downloads": missing,
            "source": source,
            "outputs": {
                "manifest": str(manifest),
                "combined_metadata": str(combined),
                "ska_distances": str(distance_path),
            },
        }
    )
    write_audit(output / "context_selection.json", audit)
    return audit
