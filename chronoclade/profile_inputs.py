"""Resolve focal genomes and frozen typing before acquiring context assemblies.

Collection membership is a query dataset. Public discovery always omits credentials;
only an explicitly requested collection and analysis/assembly downloads may use them.
"""

from __future__ import annotations

from copy import deepcopy
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import time
from typing import Any, Iterable
from urllib.parse import urlencode, urlsplit
import zipfile

from chronoclade.cglin import annotate_catalogue, download_cglin_export, load_cglin_export
from chronoclade.context import ContextCandidate
from chronoclade.context_refinement import refine_candidate_pool
from chronoclade.metadata import SAFE_IDENTIFIER, Sample
from chronoclade.pathogenwatch import (
    PathogenwatchClient,
    PathogenwatchError,
    content_hash,
    deduplicate_catalogue,
    load_catalogue,
    normalize_country,
    normalize_dates,
)
from chronoclade.pathogenwatch_download import DownloadError, download_assemblies, request_download
from chronoclade.public_typing import annotate_public_typing, load_public_typing
from chronoclade.query_typing import load_query_typing, type_query_assemblies


class ProfileInputError(ValueError):
    """Input identities or frozen profile data cannot be resolved safely."""


class CGMLSTExportIdentityError(ProfileInputError):
    """A response contains a genome outside the exact requested batch."""


def collection_id(value: str) -> str:
    """Accept a short UUID or a Pathogenwatch collection URL with an optional slug."""
    text = value.strip()
    if "://" in text:
        parsed = urlsplit(text)
        if parsed.scheme != "https" or parsed.hostname not in {
            "pathogen.watch",
            "next.pathogen.watch",
        }:
            raise ProfileInputError("Collection URL must use HTTPS on pathogen.watch")
        if parsed.username or parsed.password or parsed.port not in {None, 443}:
            raise ProfileInputError("Collection URL cannot contain credentials or a custom port")
        match = re.fullmatch(r"/collections/([^/]+)/?", parsed.path)
        if not match:
            raise ProfileInputError("Expected a Pathogenwatch /collections/UUID URL")
        text = match[1]
    ident = text.split("-", 1)[0]
    if not re.fullmatch(r"[1-9A-HJ-NP-Za-km-z]{20,22}", ident):
        raise ProfileInputError("Collection identifier must be a short UUID")
    return ident


def _csv(path: Path) -> list[dict]:
    try:
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            if not reader.fieldnames:
                raise ProfileInputError("Metadata requires a CSV header")
            rows = [
                {key: (value or "").strip() for key, value in row.items() if key} for row in reader
            ]
    except OSError as error:
        raise ProfileInputError(f"Cannot read metadata: {path}") from error
    if not rows:
        raise ProfileInputError("Metadata has no samples")
    return rows


def _freeze_file(path: Path, output: Path, label: str) -> dict:
    data = path.read_bytes()
    destination = output / "inputs" / f"{label}{path.suffix}"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(data)
    return {"path": str(destination.resolve()), "sha256": hashlib.sha256(data).hexdigest()}


def _canonical(raw: dict, *, origin: str, species: str = "", base: Path | None = None) -> dict:
    row = deepcopy(raw)
    flat = {
        key.rsplit("/", 1)[-1].casefold(): value
        for key, value in raw.items()
        if key.startswith("metadata/")
    }
    metadata = raw.get("metadata") or {}
    if isinstance(metadata, list):
        metadata = dict(metadata)
    flat.update({str(key).casefold(): value for key, value in metadata.items()})
    source = str(raw.get("source_genome_id") or raw.get("uuid") or "")
    sample = str(raw.get("sample_id") or ("PW_" + source if source else raw.get("name") or ""))
    if not SAFE_IDENTIFIER.fullmatch(sample):
        raise ProfileInputError("Each sample requires a unique safe sample_id or source UUID")
    st = str(raw.get("mlst_st") or raw.get("mlst") or raw.get("typing/MLST/ST") or "")
    lineage = str(raw.get("lineage") or ("ST" + st if st and st != "-" else "Unassigned"))
    country_raw = (
        raw.get("country")
        or flat.get("country")
        or raw.get("location")
        or flat.get("location")
        or ""
    )
    region = raw.get("region") or flat.get("city/region") or flat.get("region") or ""
    date = (
        raw.get("collection_date")
        or flat.get("collection date")
        or flat.get("collection_date")
        or flat.get("date")
        or flat.get("year")
        or ""
    )
    row.update(
        sample_id=sample,
        origin=origin,
        species=str(raw.get("species") or species),
        lineage=lineage,
        mlst_st=st,
        source_genome_id=source,
        numeric_source_id=raw.get("numeric_source_id", raw.get("id")),
        country=normalize_country(country_raw),
        country_raw=raw.get("country_raw", country_raw),
        region=str(region),
        nuts2=str(raw.get("nuts2") or flat.get("nuts2") or ""),
        location=str(raw.get("location") or country_raw or "Unknown"),
        geography_raw=raw.get("geography_raw", {"country": country_raw, "region": region}),
        aliases=list(raw.get("aliases") or []),
    )
    row.update(normalize_dates(raw, {"collection_date": date}, raw))
    # Empty or partial dates remain explicit, never replaced by upload dates.
    if raw.get("date_precision"):
        row.update(
            {
                key: raw[key]
                for key in (
                    "collection_date",
                    "date_precision",
                    "date_start",
                    "date_end",
                    "date_raw",
                )
                if key in raw
            }
        )
    for field in ("sampleAccession", "runAccession", "assemblyAccession", "biosample", "accession"):
        if raw.get(field):
            row["aliases"].append(str(raw[field]))
    for key, value in raw.items():
        if (
            key.startswith("metadata/")
            and any(
                x in key.lower()
                for x in ("run accession", "sample accession", "assembly accession")
            )
            and value
        ):
            row["aliases"].append(str(value))
    row["aliases"] += [sample, source, str(raw.get("name") or "")]
    row["biosample_accessions"] = sorted(
        set(row.get("biosample_accessions") or [])
        | {a for a in row["aliases"] if a.startswith("SAM")}
    )
    row["run_accessions"] = sorted(
        set(row.get("run_accessions") or [])
        | {a for a in row["aliases"] if re.fullmatch(r"[SED]RR[0-9]+", a)}
    )
    row["assembly_accessions"] = sorted(
        set(row.get("assembly_accessions") or [])
        | {a for a in row["aliases"] if a.startswith(("GCA_", "GCF_"))}
    )
    row["aliases"] = sorted(set(value for value in row["aliases"] if value))
    for source_field, collection_field in (
        ("source_length", "assemblyStats/Genome length"),
        ("source_n50", "assemblyStats/N50"),
        ("source_contigs", "assemblyStats/No. contigs"),
    ):
        if collection_field in raw:
            row[source_field] = raw[collection_field]
    assembly = raw.get("assembly")
    if (
        isinstance(assembly, str)
        and assembly
        and (base is not None or Path(assembly).is_absolute())
    ):
        path = Path(assembly).expanduser()
        path = path if path.is_absolute() else base / path
        if not path.is_file():
            raise ProfileInputError(f"Assembly does not exist for {sample}: {path}")
        row["assembly"] = str(path.resolve())
        row["assembly_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    else:
        row.pop("assembly", None)
    return row


def _collection(client: PathogenwatchClient, ident: str, output: Path) -> tuple[list[dict], dict]:
    details = client.request_json(
        "GET", "/api/collections/details", params={"uuid": ident}, authenticated=True
    )
    if not isinstance(details, dict) or details.get("uuid") != ident:
        raise ProfileInputError("Collection details do not match the requested UUID")
    size = details.get("size")
    if type(size) is not int or size < 1:
        raise ProfileInputError("Collection must advertise a positive genome count")
    rows, pages, seen = [], [], set()
    for page_number in range(1, size + 2):
        page = client.request_json(
            "GET",
            "/api/collections/genomes",
            params={"uuid": ident, "page": page_number},
            authenticated=True,
        )
        if (
            not isinstance(page, dict)
            or not isinstance(page.get("genomes"), list)
            or type(page.get("hasMore")) is not bool
        ):
            raise ProfileInputError("Unsupported collection pagination response")
        pages.append(page)
        if page["hasMore"] and not page["genomes"]:
            raise ProfileInputError("Collection pagination stalled")
        for member in page["genomes"]:
            uuid = member.get("uuid")
            if not uuid or uuid in seen:
                raise ProfileInputError("Collection repeats or omits a source UUID")
            seen.add(uuid)
            row = _canonical(member, origin="local", species=details.get("organismName", ""))
            row["organism_id"] = str(details.get("organismId") or "")
            rows.append(row)
        if not page["hasMore"]:
            break
    if len(rows) != size:
        raise ProfileInputError("Collection membership does not reconcile with its advertised size")
    snapshot = {"schema_version": 1, "details": details, "pages": pages}
    snapshot["snapshot_sha256"] = content_hash(snapshot)
    path = output / "collection_snapshot.json"
    path.write_text(json.dumps(snapshot, indent=2) + "\n")
    return rows, {
        "collection_uuid": ident,
        "collection_access": details.get("access", "unknown"),
        "collection_snapshot": str(path.resolve()),
        "collection_snapshot_sha256": snapshot["snapshot_sha256"],
        "downloads": details.get("downloads") or [],
        "member_count": len(rows),
    }


def parse_cgmlst_export(data: bytes, requested: Iterable[str], *, job: str) -> list[dict]:
    """Read the verified long-format allele export; never invent a full locus universe."""
    if data.startswith(b"\x1f\x8b"):
        data = gzip.decompress(data)
    if zipfile.is_zipfile(io.BytesIO(data)):
        archive = zipfile.ZipFile(io.BytesIO(data))
        files = [member for member in archive.infolist() if not member.is_dir()]
        if len(files) != 1 or files[0].file_size > 200_000_000:
            raise ProfileInputError("Expected one bounded cgMLST CSV export")
        data = archive.read(files[0])
    reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig")))
    if not {"Genome ID", "Gene", "Allele ID"}.issubset(reader.fieldnames or []):
        raise ProfileInputError("Unsupported cgMLST export columns")
    expected, profiles, ambiguous = set(requested), {}, {}
    for item in reader:
        ident, locus, allele = item["Genome ID"], item["Gene"], item["Allele ID"].strip()
        if ident not in expected:
            raise CGMLSTExportIdentityError("cgMLST export contains an unrequested genome")
        if not locus:
            raise ProfileInputError("cgMLST export contains an empty locus")
        calls = profiles.setdefault(ident, {})
        missing = allele in {"", "0", "-", "?", "None", "null"}
        if not missing and not (
            allele.isdigit() and int(allele) > 0 or re.fullmatch(r"[a-fA-F0-9]{40}", allele)
        ):
            raise ProfileInputError("cgMLST export has an unsupported allele value")
        allele = "" if missing else str(int(allele)) if allele.isdigit() else allele.lower()
        if locus in calls and calls[locus] != allele:
            ambiguous.setdefault(ident, set()).add(locus)
        calls[locus] = allele
    result = []
    for ident, calls in profiles.items():
        conflicts = ambiguous.get(ident, set())
        result.append(
            {
                "source_genome_id": ident,
                "cgmlst_scheme": "pathogenwatch:" + job,
                "cgmlst_scheme_version": job,
                "cgmlst_loci": sorted(calls),
                "cgmlst_locus_universe_complete": False,
                "cgmlst_profile": {
                    key: value for key, value in calls.items() if value and key not in conflicts
                },
                "cgmlst_novel_alleles": {},
                "cgmlst_ambiguous_loci": sorted(conflicts),
                "cgmlst_source": "pathogenwatch_analysis_export",
                "cgmlst_status": "resolved" if any(calls.values()) else "unassigned",
                "cgmlst_export_sha256": hashlib.sha256(data).hexdigest(),
            }
        )
    return result


def _analysis_exports(
    rows: list[dict],
    downloads: list[dict],
    client: PathogenwatchClient,
    output: Path,
    provenance: dict,
) -> list[dict]:
    if not client.api_key:
        provenance["analysis_export_status"] = "unavailable_no_credentials"
        return rows
    output.mkdir(parents=True, exist_ok=True)
    for advertised in downloads:
        name, job = advertised.get("name", ""), advertised.get("job", "")
        if name not in {"cgmlst", "klebsiella-lincodes"} or not re.fullmatch(
            r"[A-Za-z0-9_-]+", job
        ):
            continue
        eligible = [
            row
            for row in rows
            if str(row.get("numeric_source_id", "")).isdigit()
            and (
                name != "cgmlst"
                or not row.get("cgmlst_profile")
                and not row.get("cgmlst_novel_alleles")
            )
            and (name != "klebsiella-lincodes" or not row.get("cglin_raw"))
        ]
        if not eligible:
            continue
        try:
            if name == "klebsiella-lincodes":
                manifest = download_cglin_export(
                    eligible,
                    output / ("exports_" + job),
                    api_key=client.api_key,
                    job=job,
                    scheme_version=job,
                    base_url=client.base_url,
                )
                rows = annotate_catalogue(rows, manifest["assignments"])
                for row in rows:
                    if row.get("cglin_scheme_version") == job:
                        row["cglin_scope_source"] = "server_advertised_analysis_job"
                        row["cglin_database_version_status"] = "unavailable"
            else:
                assignments = []
                for start in range(0, len(eligible), 100):
                    batch = eligible[start : start + 100]
                    requested_ids = [row["source_genome_id"] for row in batch]
                    for attempt in range(3):
                        data = request_download(
                            "/api/downloads/cgmlst?" + urlencode({"job": job}),
                            api_key=client.api_key,
                            base_url=client.base_url,
                            body={"ids": ",".join(str(row["numeric_source_id"]) for row in batch)},
                        )
                        suffix = f"_retry{attempt}" if attempt else ""
                        destination = output / f"cgmlst_{job}_{start}{suffix}.bin"
                        destination.write_bytes(data)
                        try:
                            parsed = parse_cgmlst_export(data, requested_ids, job=job)
                        except CGMLSTExportIdentityError:
                            # Concurrent upstream jobs can return another selection's
                            # file. Preserve it, reject every call, then retry exactly
                            # this batch; never filter unrelated genomes into a result.
                            provenance.setdefault("analysis_export_identity_retries", []).append(
                                {
                                    "job": job,
                                    "attempt": attempt + 1,
                                    "requested_source_ids": requested_ids,
                                    "rejected_export": str(destination.resolve()),
                                    "rejected_export_sha256": hashlib.sha256(data).hexdigest(),
                                }
                            )
                            if attempt == 2:
                                raise
                            time.sleep(2**attempt)
                            continue
                        assignments.extend(parsed)
                        break
                # No complete universe/call fraction is inferred from the observed export.
                by_id = {row["source_genome_id"]: row for row in assignments}
                rows = [
                    dict(
                        row,
                        **{
                            key: value
                            for key, value in by_id.get(row["source_genome_id"], {}).items()
                            if key.startswith("cgmlst_")
                        },
                    )
                    for row in rows
                ]
            provenance.setdefault("analysis_exports", []).append(
                {
                    "download": name,
                    "job": job,
                    "status": "downloaded",
                    "requested_count": len(eligible),
                    "scope_provenance": "server-advertised analysis job; typing database version and fingerprint unavailable",
                }
            )
        except (PathogenwatchError, DownloadError, ProfileInputError, ValueError) as error:
            provenance.setdefault("analysis_exports", []).append(
                {"download": name, "job": job, "status": "unavailable", "reason": str(error)}
            )
    return rows


def _grouped_analysis_exports(
    rows: list[dict], client: PathogenwatchClient, output: Path, provenance: dict,
    *, download_names: set[str] | None = None,
) -> list[dict]:
    """Discover server-advertised jobs for exact already-resolved public IDs."""
    eligible = [
        row
        for row in rows
        if str(row.get("numeric_source_id", "")).isdigit()
        and (
            (download_names == {"klebsiella-lincodes"} and not row.get("cglin_raw"))
            or (download_names is None and not row.get("cgmlst_profile")
                and not row.get("cgmlst_novel_alleles"))
        )
    ]
    if not eligible:
        return rows
    if not client.api_key:
        provenance["analysis_export_status"] = "unavailable_no_credentials"
        return rows
    output.mkdir(parents=True, exist_ok=True)
    requested = {int(row["numeric_source_id"]): row for row in eligible}
    try:
        groups = client.request_json(
            "POST", "/api/genomes/group", body={"ids": sorted(requested)}, authenticated=True
        )
        if not isinstance(groups, list):
            raise ProfileInputError("Unsupported genome-group export capability response")
        returned: set[int] = set()
        resolved: dict[str, dict] = {}
        for index, group in enumerate(groups):
            ids = group.get("ids")
            if not isinstance(ids, list) or any(
                type(ident) is not int or ident not in requested or ident in returned
                for ident in ids
            ):
                raise ProfileInputError(
                    "Genome-group capabilities contain unrequested/repeated IDs"
                )
            returned.update(ids)
            members = [requested[ident] for ident in ids]
            if any(
                row.get("organism_id") and str(row["organism_id"]) != str(group.get("organismId"))
                for row in members
            ):
                raise ProfileInputError(
                    "Genome-group export organism does not match public identity"
                )
            advertised = group.get("downloads") or []
            if download_names is not None:
                advertised = [item for item in advertised if item.get("name") in download_names]
            if download_names is None and not any(item.get("name") == "cgmlst" for item in advertised):
                provenance.setdefault("profile_export_unavailable", []).append(
                    {
                        "organism_id": group.get("organismId"),
                        "reason": "server_does_not_advertise_cgmlst_job",
                    }
                )
            exported = _analysis_exports(
                members, advertised, client, output / f"organism_{index}", provenance
            )
            resolved.update({row["source_genome_id"]: row for row in exported})
        snapshot = {"groups": groups, "requested_numeric_ids": sorted(requested)}
        snapshot["sha256"] = content_hash(snapshot)
        (output / "export_capabilities.json").write_text(json.dumps(snapshot, indent=2) + "\n")
        if returned != set(requested):
            provenance.setdefault("profile_export_unavailable", []).append(
                {
                    "reason": "source_ids_omitted_by_export_capability_endpoint",
                    "numeric_ids": sorted(set(requested) - returned),
                }
            )
        return [resolved.get(row.get("source_genome_id"), row) for row in rows]
    except (PathogenwatchError, ProfileInputError) as error:
        provenance.setdefault("profile_export_unavailable", []).append({"reason": str(error)})
        return rows


def _accession_input(path: Path) -> tuple[list[str], list[dict]]:
    """Read a one-ID-per-line list or a CSV with an exact accession/source ID."""
    text = path.read_text(encoding="utf-8-sig")
    first = text.splitlines()[0] if text.splitlines() else ""
    overrides: list[dict] = []
    if "," in first or first.strip() in {"accession", "source_genome_id"}:
        overrides = _csv(path)
        identifiers = [
            row.get("accession") or row.get("source_genome_id") or "" for row in overrides
        ]
        if any(not ident for ident in identifiers):
            raise ProfileInputError(
                "Accession CSV requires accession or source_genome_id on every row"
            )
    else:
        identifiers = [
            line.strip()
            for line in text.splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
    if not identifiers or len(identifiers) != len(set(identifiers)):
        raise ProfileInputError("Accession list must contain unique nonempty identifiers")
    return identifiers, overrides


def _join_accession(accession: str, candidates: list[dict]) -> dict:
    matches = [
        row
        for row in candidates
        if accession
        in set(row.get("aliases") or [])
        | set(row.get("run_accessions") or [])
        | set(row.get("assembly_accessions") or [])
        | set(row.get("biosample_accessions") or [])
        | {str(row.get("biosample") or ""), str(row.get("source_genome_id") or "")}
    ]
    if len(matches) != 1:
        raise ProfileInputError(
            f"Accession {accession} has {len(matches)} exact matches; supply an exact source_genome_id or frozen catalogue"
        )
    return deepcopy(matches[0])


def _apply_overrides(rows: list[dict], overrides: list[dict], base: Path) -> list[dict]:
    output = deepcopy(rows)
    for override in overrides:
        identifier = (
            override.get("source_genome_id")
            or override.get("accession")
            or override.get("sample_id")
        )
        matches = [
            row
            for row in output
            if identifier in row.get("aliases", []) or identifier == row["sample_id"]
        ]
        if len(matches) != 1:
            raise ProfileInputError(
                f"Metadata override {identifier} does not match exactly one query"
            )
        target = matches[0]
        conflicts = target.setdefault("metadata_override_conflicts", [])
        for key, value in override.items():
            if value and key not in {"source_genome_id", "accession"}:
                if (
                    target.get(key) is not None
                    and target.get(key) != ""
                    and target.get(key) != value
                ):
                    conflicts.append(
                        {"field": key, "source_value": target.get(key), "manual_value": value}
                    )
                target[key] = value
        if override.get("collection_date"):
            for key in ("date_precision", "date_raw", "date_start", "date_end"):
                target.pop(key, None)
        target["manual_metadata_override"] = True
        target.update(_canonical(target, origin="local", base=base))
    return output


def _sample(row: dict) -> Sample:
    return Sample(
        sample_id=row["sample_id"],
        assembly=Path(row["assembly"]),
        collection_date=row.get("collection_date", ""),
        location=row.get("location") or row.get("country") or "Unknown",
        species=row["species"],
        lineage=row.get("lineage") or "Unassigned",
        origin=row.get("origin", "local"),
        is_reference=str(row.get("is_reference", "")).lower() in {"true", "1", "yes"},
        patient_id=row.get("patient_id", ""),
    )


def materialise_assemblies(
    records: Iterable[dict], *, output: Path, client: PathogenwatchClient | None = None
) -> list[Sample]:
    """Acquire only the records explicitly selected by the caller."""
    rows = deepcopy(list(records))
    missing = [row for row in rows if not row.get("assembly")]
    if missing:
        client = client or PathogenwatchClient()
        if not client.api_key:
            raise ProfileInputError(
                "Assembly acquisition requires a configured Pathogenwatch API key"
            )
        if any(not row.get("source_genome_id") for row in missing):
            raise ProfileInputError("Assembly acquisition requires exact source genome IDs")
        paths, audit = download_assemblies(
            missing,
            output=output / "assemblies",
            cache_dir=output / "assembly_cache",
            api_key=client.api_key,
            base_url=client.base_url,
        )
        (output / "assembly_acquisition.json").write_text(json.dumps(audit, indent=2) + "\n")
        if any(row["source_genome_id"] not in paths for row in missing):
            raise ProfileInputError(
                "Some selected assemblies could not be acquired; see assembly_acquisition.json"
            )
        for row in missing:
            row["assembly"] = str(paths[row["source_genome_id"]])
    return [_sample(row) for row in rows]


def _refined_context_pool(
    rows: list[dict], queries: list[dict], *, limit: int, seed: int
) -> tuple[list[dict], dict]:
    """Refine independent species/ST pools before their bounded allele exports."""
    if not rows or limit == 0:
        return [], {"method": "query-group-refinement", "selected": 0, "lineages": []}
    groups: dict[tuple[str, str], list[dict]] = {}
    for row in rows:
        groups.setdefault((row["species"].casefold(), row["lineage"]), []).append(row)
    quotas = {scope: 0 for scope in groups}
    budget = min(limit, len(rows))
    while budget:
        progressed = False
        for scope in sorted(groups):
            if quotas[scope] < len(groups[scope]):
                quotas[scope] += 1
                budget -= 1
                progressed = True
                if not budget:
                    break
        if not progressed:
            break
    selected, audits = [], []
    for scope in sorted(groups):
        quota = quotas[scope]
        if not quota:
            continue
        members = groups[scope]
        focal = [row for row in queries if (row["species"].casefold(), row["lineage"]) == scope]
        candidates = [
            ContextCandidate(
                sample_id=row["sample_id"],
                species=row["species"],
                lineage=row["lineage"],
                mlst_scheme=row.get("mlst_scheme") or "mlst",
                mlst_st=row.get("mlst_st") or row["lineage"].removeprefix("ST"),
                source_genome_id=row["source_genome_id"],
                collection_date=row.get("collection_date", ""),
                country=row.get("country", "Unknown"),
                date_precision=row.get("date_precision", "missing"),
            )
            for row in members
        ]
        # A one-slot budget must still allow an evidence match to take priority.
        refined, audit = refine_candidate_pool(
            candidates, members, focal, limit=max(2, quota), seed=seed
        )
        refined = refined[:quota]
        chosen = {item.source_genome_id: item for item in refined}
        selected.extend(
            dict(
                row, profile_pool_selection_reason=chosen[row["source_genome_id"]].selection_reason
            )
            for row in members
            if row["source_genome_id"] in chosen
        )
        audit["settings"]["requested_profile_limit"] = quota
        audit["counts"]["selected"] = len(refined)
        audit["counts"]["priority_selected"] = sum(
            item.selection_reason.startswith(("lineage_priority:", "cgmlst_priority:"))
            for item in refined
        )
        audit["counts"]["background_selected"] = len(refined) - audit["counts"]["priority_selected"]
        for item in audit["candidates"].values():
            item["selected"] = item["source_genome_id"] in chosen
        audits.append({"species": scope[0], "lineage": scope[1], "audit": audit})
    return selected, {
        "method": "query-group-refinement-before-profile-export",
        "selected": len(selected),
        "lineages": audits,
    }


def _stable_records(value: Any) -> Any:
    """Ignore acquisition clocks while preserving biological/source dates and hashes."""
    if isinstance(value, dict):
        return {
            key: _stable_records(item)
            for key, item in value.items()
            if key not in {"retrieved_at", "elapsed_seconds", "seconds", "requests", "attempts"}
            and not key.endswith("_retrieved_at")
        }
    if isinstance(value, list):
        return [_stable_records(item) for item in value]
    return value


def _frozen_cglin(rows: list[dict], assignments: list[dict]) -> list[dict]:
    """Apply the explicit export as the sole cgLIN namespace, including missing IDs.

    Live/native typing can still supply allele profiles. Its cluster release must
    not replace the user's frozen cluster assignments or leak into unmatched IDs.
    """
    clean = [
        {key: value for key, value in row.items() if not key.startswith("cglin_") and key != "cgst"}
        for row in rows
    ]
    annotated = annotate_catalogue(clean, assignments)
    for row in annotated:
        row["cglin_scope_source"] = "user_frozen_export"
        row["cglin_database_version_status"] = (
            "unavailable" if row.get("cglin_scheme_version") == "unknown" else "user_declared"
        )
    return annotated


def resolve_profile_inputs(
    metadata: Path | None,
    *,
    collection: str | None = None,
    accessions: Path | None = None,
    species: str | None = None,
    output: Path,
    client: PathogenwatchClient | None = None,
    catalogue: Path | None = None,
    public_typing: Path | None = None,
    query_typing: Path | None = None,
    typing_config: Path | None = None,
    cglin_export: Path | None = None,
    profile_limit: int = 500,
    lin_min_context: int = 20,
    seed: int = 1729,
) -> dict[str, Any]:
    """Freeze input identity, metadata and profile evidence; defer public assemblies.

    JSON contract: ``queries`` and ``context`` are lists of canonical record dicts;
    ``context`` is a bounded profile pool plus explicitly provided comparison rows.
    ``catalogue_rows`` retains the complete deduplicated public metadata separately.
    ``provenance`` contains frozen input paths/SHA256s, collection/catalogue scopes,
    exact identity/deduplication audit and per-dataset profile coverage. Records use
    sample_id, origin (local/context), species, lineage, collection_date with date
    bounds/precision, country/raw geography, region and explicit nuts2. Public
    records carry source_genome_id and numeric_source_id; assembly is set only for
    existing local files. cgmlst_*, cglin_* and hiercc_* fields retain their imported
    scheme, version, database provenance and exact calls. Missing profiles are
    profile_status=missing. No context assembly is downloaded here.

    For the verified Pathogenwatch cgMLST long CSV, cgmlst_loci lists observed
    loci only and cgmlst_locus_universe_complete=False. Neither a full-universe
    completeness fraction nor a typing-database fingerprint is invented. The
    server-advertised job freezes the scheme/version scope. Raw exports are saved.
    """
    if profile_limit < 0:
        raise ProfileInputError("Profile limit must be nonnegative")
    output.mkdir(parents=True, exist_ok=True)
    provenance: dict[str, Any] = {
        "schema_version": 1,
        "inputs": {},
        "profile_limit": profile_limit,
        "seed": seed,
        "public_context_credential_free": True,
    }
    for name, path in [
        ("metadata", metadata),
        ("accessions", accessions),
        ("catalogue", catalogue),
        ("public_typing", public_typing),
        ("query_typing", query_typing),
        ("typing_config", typing_config),
        ("cglin_export", cglin_export),
    ]:
        if path is not None:
            provenance["inputs"][name] = _freeze_file(Path(path), output, name)
    frozen_rows = load_catalogue(catalogue)["rows"] if catalogue else []
    public_assignments = load_public_typing(public_typing) if public_typing else []
    cglin_assignments = load_cglin_export(cglin_export) if cglin_export else []
    if cglin_export:
        provenance["cglin_authority"] = {
            "source": "user_frozen_export",
            "input": provenance["inputs"]["cglin_export"],
            "policy": "Explicit export is authoritative for queries and context; unmatched IDs remain unassigned in this frozen namespace.",
            "scheme_versions": sorted({row["cglin_scheme_version"] for row in cglin_assignments}),
            "database_version_inferred_from_live_jobs": False,
        }
    provided_context: list[dict] = []
    if collection:
        if accessions:
            raise ProfileInputError("Choose a collection or accession list as the query dataset")
        client = client or PathogenwatchClient()
        queries, collection_provenance = _collection(client, collection_id(collection), output)
        provenance.update(collection_provenance)
        if public_assignments:
            queries = annotate_public_typing(queries, public_assignments)
        queries = _analysis_exports(
            queries,
            collection_provenance["downloads"],
            client,
            output / "query_exports",
            provenance,
        )
        if metadata:
            queries = _apply_overrides(queries, _csv(metadata), Path(metadata).parent)
    elif accessions:
        identifiers, accession_metadata = _accession_input(Path(accessions))
        candidates = [
            _canonical(dict(row, sample_id="PW_" + row["source_genome_id"]), origin="local")
            for row in frozen_rows
        ]
        if not candidates:
            client = client or PathogenwatchClient()
            source_ids = [
                ident
                for ident in identifiers
                if re.fullmatch(r"[1-9A-HJ-NP-Za-km-z]{20,22}", ident)
            ]
            for ident in source_ids:
                detail = client.request_json("GET", "/api/genomes/details", params={"id": ident})
                if (
                    not isinstance(detail, dict)
                    or detail.get("uuid") != ident
                    or detail.get("projectAccess", "PUBLIC") != "PUBLIC"
                ):
                    raise ProfileInputError(
                        "Explicit source UUID did not resolve to the exact public genome"
                    )
                canonical = _canonical(detail, origin="local", species=species or "")
                canonical["organism_id"] = str(detail.get("organismId") or "")
                candidates.append(canonical)
            search_ids = [ident for ident in identifiers if ident not in source_ids]
            organisms = client.supported_organisms() if search_ids else []
            organism = next(
                (
                    row
                    for row in organisms
                    if row.get("fullName", "").casefold() == (species or "").casefold()
                ),
                None,
            )
            if search_ids and organism is None:
                raise ProfileInputError("Accession queries require species or a frozen catalogue")
            for accession in search_ids:
                response = client.request_json(
                    "POST",
                    "/api/search/genomes",
                    body={
                        "organismId": str(organism["organismId"]),
                        "searchText": accession,
                        "qc": [True, False],
                    },
                    params={"limit": 100, "sort": "id"},
                )
                count = (
                    response.get("meta", {}).get("count") if isinstance(response, dict) else None
                )
                found = response.get("genomes") if isinstance(response, dict) else None
                if (
                    type(count) is not int
                    or count < 0
                    or not isinstance(found, list)
                    or len(found) != count
                    or count > 100
                ):
                    raise ProfileInputError(
                        "Accession search must be complete and bounded for an exact identity join"
                    )
                for search in response.get("genomes", []):
                    if search.get("projectAccess") != "PUBLIC":
                        raise ProfileInputError("Public accession search returned a private record")
                    detail = client.request_json(
                        "GET", "/api/genomes/details", params={"id": search["uuid"]}
                    )
                    if detail.get("uuid") != search["uuid"]:
                        raise ProfileInputError("Accession detail identity mismatch")
                    row = _canonical({**search, **detail}, origin="local", species=species or "")
                    if row["source_genome_id"] not in {
                        item["source_genome_id"] for item in candidates
                    }:
                        candidates.append(row)
        queries = [_join_accession(accession, candidates) for accession in identifiers]
        if accession_metadata:
            queries = _apply_overrides(queries, accession_metadata, Path(accessions).parent)
        if metadata:
            queries = _apply_overrides(queries, _csv(metadata), Path(metadata).parent)
    elif metadata:
        imported = [
            _canonical(
                row,
                origin="context" if row.get("origin", "").casefold() == "context" else "local",
                species=species or "",
                base=Path(metadata).parent,
            )
            for row in _csv(metadata)
        ]
        queries = [row for row in imported if row["origin"] == "local"]
        provided_context = [
            dict(row, provided_context=True) for row in imported if row["origin"] == "context"
        ]
        provenance["provided_context_count"] = len(provided_context)
    else:
        raise ProfileInputError("Provide assembly metadata, an accession list, or a collection")
    if not queries:
        raise ProfileInputError("At least one focal/local query is required")
    all_inputs = queries + provided_context
    if len({row["sample_id"] for row in all_inputs}) != len(all_inputs):
        raise ProfileInputError("Query sample IDs must be unique")
    if any(not row["species"] for row in all_inputs):
        raise ProfileInputError("Each query needs a species identity")
    if query_typing:
        assignments = {row["sample_id"]: row for row in load_query_typing(query_typing)}
        for row in all_inputs:
            assignment = assignments.get(row["sample_id"])
            if assignment:
                if assignment["species"].casefold() != row["species"].casefold():
                    raise ProfileInputError("Query typing species differs from input metadata")
                if (
                    row.get("assembly_sha256")
                    and row["assembly_sha256"] != assignment["assembly_sha256"]
                ):
                    raise ProfileInputError(
                        "Query typing assembly SHA256 differs from input assembly"
                    )
                row.update(assignment)
    if public_assignments:
        queries = annotate_public_typing(queries, public_assignments)
        provided_context = annotate_public_typing(provided_context, public_assignments)
    if cglin_export:
        queries = _frozen_cglin(queries, cglin_assignments)
        provided_context = _frozen_cglin(provided_context, cglin_assignments)
    queries_with_existing_profiles = {
        row["sample_id"] for row in queries
        if row.get("cgmlst_profile") or row.get("cgmlst_novel_alleles")
    }
    if not collection and any(
        str(row.get("numeric_source_id", "")).isdigit()
        and not row.get("cgmlst_profile")
        and not row.get("cgmlst_novel_alleles")
        for row in queries
    ):
        client = client or PathogenwatchClient()
        queries = _grouped_analysis_exports(queries, client, output / "query_exports", provenance)
    if typing_config:
        need_typing = [
            row
            for row in queries
            if not row.get("cgmlst_profile") and not row.get("cgmlst_novel_alleles")
        ]
        if need_typing:
            samples = materialise_assemblies(
                need_typing, output=output / "query_typing_inputs", client=client
            )
            assignments = type_query_assemblies(
                samples, Path(typing_config), output / "native_query_typing"
            )
            by_sample = {row["sample_id"]: row for row in assignments}
            for row in queries:
                if row["sample_id"] in by_sample:
                    row.update(by_sample[row["sample_id"]])
    if cglin_export:
        queries = _frozen_cglin(queries, cglin_assignments)
    # An available allele profile does not imply the LIN assignment is present.
    # Obtain only missing query LIN codes before defining the contextual pool.
    frozen_typing_only = client is None and catalogue is not None and (
        public_typing is not None or query_typing is not None
    )
    if not cglin_export and not frozen_typing_only and any(
        row["sample_id"] in queries_with_existing_profiles
        and str(row.get("numeric_source_id", "")).isdigit()
        and "klebsiella" in row["species"].casefold()
        and not row.get("cglin_raw")
        for row in queries
    ):
        client = client or PathogenwatchClient()
        if client.api_key:
            queries = _grouped_analysis_exports(
                queries, client, output / "query_lin_exports", provenance,
                download_names={"klebsiella-lincodes"},
            )
    scopes = {(row["species"].casefold(), row["lineage"]) for row in queries}
    context = [
        _canonical(dict(row, sample_id="PW_" + row["source_genome_id"]), origin="context")
        for row in frozen_rows
    ]
    if not catalogue and profile_limit > 0:
        client = client or PathogenwatchClient()
        organisms = client.supported_organisms()
        for index, (species_name, lineage) in enumerate(sorted(scopes)):
            organism = next(
                (row for row in organisms if row.get("fullName", "").casefold() == species_name),
                None,
            )
            st = lineage.removeprefix("ST")
            if organism is None or not st.isdigit():
                provenance.setdefault("context_discovery_skipped", []).append(
                    {
                        "species": species_name,
                        "lineage": lineage,
                        "reason": "species_or_ST_unavailable",
                    }
                )
                continue
            snapshot = client.freeze_catalogue(
                output / f"public_catalogue_{index}.json",
                organism_id=str(organism["organismId"]),
                st=st,
            )
            context.extend(
                _canonical(dict(row, sample_id="PW_" + row["source_genome_id"]), origin="context")
                for row in snapshot["rows"]
            )
            provenance.setdefault("public_catalogues", []).append(snapshot["provenance"])
    catalogue_rows = deepcopy(context)
    if catalogue_rows:
        catalogue_rows, catalogue_audit = deduplicate_catalogue(catalogue_rows)
        provenance["catalogue_deduplication"] = catalogue_audit
    if public_assignments:
        catalogue_rows = annotate_public_typing(catalogue_rows, public_assignments)
    if cglin_export:
        catalogue_rows = _frozen_cglin(catalogue_rows, cglin_assignments)
    # LIN assignment precedes sampling: fetch codes for the whole same-ST pool,
    # without fetching its profiles or assemblies.
    if (profile_limit > 0 and any(row.get("cglin_raw") for row in queries)
            and any(not row.get("cglin_raw") for row in catalogue_rows)):
        client = client or PathogenwatchClient()
        catalogue_rows = _grouped_analysis_exports(
            catalogue_rows, client, output / "catalogue_lin_exports", provenance,
            download_names={"klebsiella-lincodes"},
        )
        if cglin_export:
            catalogue_rows = _frozen_cglin(catalogue_rows, cglin_assignments)
        catalogue_by_id = {row["source_genome_id"]: row for row in catalogue_rows}
        context = [catalogue_by_id.get(row["source_genome_id"], row) for row in context]
    provenance["catalogue_metadata_count"] = len(catalogue_rows)
    excluded_inputs = queries + provided_context
    query_ids = {
        row.get("source_genome_id") for row in excluded_inputs if row.get("source_genome_id")
    }
    query_aliases = set().union(*(set(row.get("aliases") or []) for row in excluded_inputs))
    context = [
        row
        for row in context
        if (row["species"].casefold(), row["lineage"]) in scopes
        and row.get("source_genome_id") not in query_ids
    ]
    if context:
        context, audit = deduplicate_catalogue(context, query_aliases)
        provenance["context_deduplication"] = audit
    if public_assignments:
        context = annotate_public_typing(context, public_assignments)
    if cglin_export:
        context = _frozen_cglin(context, cglin_assignments)
    unbounded_context_count = len(context)
    from chronoclade.adaptive_context import adaptive_cglin_context

    if profile_limit > 0:
        queries, adaptive_context, adaptive_audit = adaptive_cglin_context(
            queries, context, min_context=lin_min_context
        )
    else:
        adaptive_context, adaptive_audit = [], {"datasets": [], "subgroups": []}
    if adaptive_audit["datasets"]:
        adaptive_ids = {r["sample_id"] for r in adaptive_context if r.get("analysis_dataset")}
        fallback_queries = [r for r in queries if not r.get("analysis_dataset")]
        fallback_context = [r for r in context if r["sample_id"] not in adaptive_ids and any(
            (r["species"].casefold(), r["lineage"]) == (q["species"].casefold(), q["lineage"])
            for q in fallback_queries
        )]
        fallback_context, refinement_audit = _refined_context_pool(
            fallback_context, fallback_queries, limit=profile_limit, seed=seed
        )
        context = [r for r in adaptive_context if r["sample_id"] in adaptive_ids] + fallback_context
        provenance["adaptive_context_selection"] = adaptive_audit
    else:
        context, refinement_audit = _refined_context_pool(
            context, queries, limit=profile_limit, seed=seed
        )
    provenance["eligible_public_context_count"] = unbounded_context_count
    provenance["bounded_public_context_count"] = len(context)
    provenance["context_pool_selection"] = refinement_audit
    if not collection and any(
        str(row.get("numeric_source_id", "")).isdigit()
        and not row.get("cgmlst_profile")
        and not row.get("cgmlst_novel_alleles")
        for row in context
    ):
        client = client or PathogenwatchClient()
        context = _grouped_analysis_exports(context, client, output / "context_exports", provenance)
    # A collection's advertised job is usable only for same-species context.
    if collection and context and profile_limit:
        ranked = sorted(
            context,
            key=lambda row: hashlib.sha256(
                f"{seed}:{row['source_genome_id']}".encode()
            ).hexdigest(),
        )
        export_rows = ranked if adaptive_audit["datasets"] else ranked[:profile_limit]
        selected_ids = {row["source_genome_id"] for row in export_rows}
        selected = [
            row
            for row in context
            if row["source_genome_id"] in selected_ids and not row.get("cgmlst_profile")
        ]
        exported = _analysis_exports(
            selected,
            provenance.get("downloads", []),
            client,
            output / "context_exports",
            provenance,
        )
        by_id = {row["source_genome_id"]: row for row in exported}
        context = [by_id.get(row["source_genome_id"], row) for row in context]
    context.extend(dict(row, provided_context=True) for row in provided_context)
    if cglin_export:
        queries = _frozen_cglin(queries, cglin_assignments)
        context = _frozen_cglin(context, cglin_assignments)
        provenance["cglin_authority"]["matched_records"] = {
            label: sum(row.get("cglin_export_record_count", 0) > 0 for row in rows)
            for label, rows in [
                ("queries", queries),
                ("context", context),
                ("catalogue_rows", catalogue_rows),
            ]
        }
    coverage = {}
    for label, rows in [("queries", queries), ("context", context)]:
        for row in rows:
            row["profile_status"] = (
                "available"
                if row.get("cgmlst_profile") or row.get("cgmlst_novel_alleles")
                else "missing"
            )
        coverage[label] = {
            "total": len(rows),
            "profiles_available": sum(row["profile_status"] == "available" for row in rows),
            "profiles_missing": sum(row["profile_status"] == "missing" for row in rows),
        }
    provenance["coverage"] = coverage
    result = {
        "queries": queries,
        "context": context,
        "catalogue_rows": catalogue_rows,
        "provenance": provenance,
    }
    result["provenance"]["resolved_records_sha256"] = content_hash(
        {"queries": queries, "context": context}
    )
    (output / "profile_inputs.json").write_text(json.dumps(result, indent=2) + "\n")
    return result
