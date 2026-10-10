"""Adapt frozen profile input records without retrieval or typing decisions."""

from collections.abc import Mapping, Sequence
import re

from chronoclade.typing_scopes import typing_scope
from chronoclade.metadata_dates import normalize_dates
from chronoclade.sample_labels import sample_labels
from .model import AlleleMatrix, DatasetError, LocusCatalogue, PreparedDataset, SAMPLE_COLUMNS


def samples_from_records(records: Sequence[Mapping]) -> tuple[dict, ...]:
    """Adapt canonical metadata without copying large profiles into the samples table.

    Provider identity resolution and field precedence are deliberately upstream.
    This adapter only selects existing metadata and normalizes date precision.
    Raw date evidence and accession aliases remain available for inspection.
    """
    labels = sample_labels(records)
    result = []
    for record in records:
        ident = record.get("sample_id")
        row = {key: record.get(key) or None for key in SAMPLE_COLUMNS}
        row["label"] = record.get("label") or labels.get(ident, ident)
        origin = record.get("origin")
        row["role"] = record.get("role") or (
            "input" if origin in {"local", "query", "focal", "input"} else origin
        )
        row["assembly_reference"] = record.get("assembly_reference") or record.get("assembly")
        dates = normalize_dates({}, {"collection_date": record.get("collection_date") or ""}, {})
        for key in ("date_start", "date_end", "date_precision"):
            row[key] = record.get(key) or dates[key] or None
        row["date_raw"] = record.get("date_raw", dates["date_raw"])
        for key in (
            "source_genome_id",
            "aliases",
            "run_accessions",
            "biosample_accessions",
            "assembly_accessions",
            "metadata_status",
            "country_raw",
            "geography_raw",
            "accession",
            "run_accession",
            "runAccession",
            "sample_accession",
            "sampleAccession",
            "assembly_accession",
            "assemblyAccession",
            "biosample",
            "host_raw",
            "isolation_source_raw",
            "hospital",
            "patient_id",
            "location",
        ):
            if key in record:
                row[key] = record[key]
        result.append(row)
    return tuple(result)


def _record_lineages(records):
    rows = []
    for record in records:
        for kind, field_name in (("cglin", "cglin_raw"), ("hiercc", "hiercc_codes")):
            evidence = {key: value for key, value in record.items() if key.startswith(kind + "_")}
            if not evidence:
                continue
            scope = typing_scope(record, kind)
            rows.append(
                {
                    "sample_id": record["sample_id"],
                    "kind": kind,
                    "scheme_id": record.get(kind + "_scheme") or None,
                    "scheme_version": scope[1]
                    if scope
                    else record.get(kind + "_scheme_version") or None,
                    "database_version": record.get(kind + "_database_version") or None,
                    "assignment": record.get(field_name) or None,
                    "assignment_method": record.get(kind + "_source") or "frozen_record",
                    "resolution_status": record.get(kind + "_status") or "unverified",
                    "evidence": evidence,
                }
            )
    return tuple(rows)


def _record_crosswalk(records):
    fields = (
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
    )
    return tuple(
        {
            "sample_id": row["sample_id"],
            "resolution_status": "ambiguous"
            if row.get("identity_conflict")
            else "resolved"
            if row.get("identity_resolved")
            else "unresolved",
            **{key: row[key] for key in fields if key in row},
        }
        for row in records
    )


def _record_retrieval(records):
    return tuple(
        {
            "sample_id": row["sample_id"],
            "operation": "frozen_profile",
            "source": row.get("cgmlst_source") or "frozen_record",
            "status": row.get("cgmlst_status") or "unavailable",
            "evidence": {
                key: value
                for key, value in row.items()
                if key.startswith("cgmlst_")
                and key not in {"cgmlst_profile", "cgmlst_novel_alleles", "cgmlst_loci"}
            },
        }
        for row in records
    )


def from_profile_records(
    records: Sequence[Mapping],
    catalogues: Sequence[LocusCatalogue],
    *,
    lineages=None,
    crosswalk=None,
    provenance=None,
    conflicts=(),
    exclusions=(),
    retrieval=None,
    parameters=None,
) -> PreparedDataset:
    """Convert frozen canonical records using explicit, compatible catalogues.

    No typing or provider lookup occurs. Incompatible or unversioned profiles
    raise an error instead of being silently discarded. Missing profiles remain
    as samples; unavailable/failed profile status produces an all-null row.
    Existing complete lineage assignments and accession crosswalks are retained.
    Evidence supplied explicitly overrides automatic frozen-record evidence; no
    provider, identity resolution or source-precedence decision is invented.
    """
    matrices = []
    assigned = set()
    for catalogue in catalogues:
        catalogue.validate()
        members, calls = [], {}
        for row in records:
            if typing_scope(row, "cgmlst") != (
                catalogue.scheme_id,
                catalogue.scheme_version,
            ):
                continue
            database = row.get("cgmlst_database_version")
            if database and catalogue.database_version and database != catalogue.database_version:
                continue
            checksum = row.get("cgmlst_database_sha256") or None
            if checksum is not None:
                if not isinstance(checksum, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", checksum):
                    raise DatasetError("Frozen profile database fingerprint is not a SHA256 digest")
                checksum = checksum.lower()
            # A declared database cannot be silently mixed with a different or
            # unknown frozen fingerprint, even when scheme versions coincide.
            if checksum != catalogue.database_sha256:
                continue
            ident = row["sample_id"]
            if ident in assigned:
                raise DatasetError("A profile record matches multiple supplied catalogues")
            assigned.add(ident)
            members.append(ident)
            regular = row.get("cgmlst_profile") or {}
            novel = row.get("cgmlst_novel_alleles") or {}
            if not isinstance(regular, Mapping) or not isinstance(novel, Mapping):
                raise DatasetError("Frozen profiles and novel alleles must be locus mappings")
            if set(regular) & set(novel):
                raise DatasetError("Called and novel alleles overlap at a locus")
            combined = dict(regular, **novel)
            ambiguous = set(row.get("cgmlst_ambiguous_loci") or ())
            if set(combined) & ambiguous:
                raise DatasetError("An ambiguous locus cannot also have an allele call")
            if row.get("cgmlst_status") in {"conflict", "failed", "unassigned", "unsupported"}:
                combined = {}
            calls[ident] = combined
        if members:
            matrices.append(AlleleMatrix.from_profiles(catalogue, members, calls))
    for row in records:
        if (row.get("cgmlst_profile") or row.get("cgmlst_novel_alleles")) and row[
            "sample_id"
        ] not in assigned:
            raise DatasetError(
                "Frozen profile has no compatible explicitly supplied catalogue or database version/fingerprint"
            )
    samples = samples_from_records(records)
    if provenance is None:
        provenance = tuple(
            {"sample_id": row["sample_id"], "field": key, "source": "frozen_record", "value": value}
            for row in samples
            for key, value in row.items()
            if key not in {"sample_id", "label", "role"} and value is not None
        )
    dataset = PreparedDataset(
        samples=samples,
        profiles=tuple(matrices),
        lineages=_record_lineages(records) if lineages is None else tuple(lineages),
        crosswalk=_record_crosswalk(records) if crosswalk is None else tuple(crosswalk),
        provenance=tuple(provenance),
        conflicts=tuple(conflicts),
        exclusions=tuple(exclusions),
        retrieval=_record_retrieval(records) if retrieval is None else tuple(retrieval),
        parameters=dict(parameters or {}),
    )
    dataset.validate()
    return dataset
