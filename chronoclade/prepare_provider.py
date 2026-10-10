"""Query-only provider integration and complete public profile context retrieval.

Public searches never authenticate. Only exact resolved query membership and
analysis exports use configured credentials. No context assembly is acquired.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

from chronoclade.datasets import DatasetError, LocusCatalogue
from chronoclade.datasets.model import check_json
from chronoclade.pathogenwatch import PathogenwatchClient, normalize_country
from chronoclade.profile_inputs import (
    _accession_input,
    _analysis_exports,
    _apply_overrides,
    _canonical,
    _collection,
    _csv,
    _frozen_cglin,
    _grouped_analysis_exports,
    _join_accession,
    collection_id,
    materialise_assemblies,
)
from chronoclade.typing_scopes import typing_scope


def _save(path, value):
    check_json(value, location=path.name)
    path.write_text(json.dumps(value, indent=2) + "\n")


def _has_profile(row):
    return bool(row.get("cgmlst_profile") or row.get("cgmlst_novel_alleles"))


def _has_lineage(row):
    kind = "cglin" if "klebsiella" in row.get("species", "").casefold() else "hiercc"
    return bool(
        typing_scope(row, kind)
        and row.get("cglin_raw" if kind == "cglin" else "hiercc_codes")
        and row.get(kind + "_status") in {"complete", "resolved", "predicted"}
        and not row.get(kind + "_provisional")
    )


def apply_metadata_overrides(rows, overrides, base):
    """Keep field-specific precedence rather than labelling every field as user supplied."""
    sources = {}
    for row in rows:
        fields = set()
        aliases = set(row.get("aliases") or []) | {row["sample_id"], row.get("source_genome_id")}
        for override in overrides:
            ident = (
                override.get("source_genome_id")
                or override.get("accession")
                or override.get("sample_id")
            )
            if ident in aliases:
                fields.update(
                    k
                    for k, value in override.items()
                    if value and k not in {"source_genome_id", "accession"}
                )
        sources[row.get("source_genome_id") or row["sample_id"]] = fields
    updated = _apply_overrides(rows, overrides, base)
    for row in updated:
        row["prepare_user_metadata_fields"] = sorted(
            sources.get(row.get("source_genome_id") or row["sample_id"], set())
        )
    return updated


def _resolve_accessions(path, species, client):
    identifiers, overrides = _accession_input(path)
    candidates = []
    organisms = None
    for ident in identifiers:
        if re.fullmatch(r"[1-9A-HJ-NP-Za-km-z]{20,22}", ident):
            searches = [{"uuid": ident}]
        else:
            organisms = organisms or client.supported_organisms()
            organism = next(
                (
                    r
                    for r in organisms
                    if r.get("fullName", "").casefold() == (species or "").casefold()
                ),
                None,
            )
            if organism is None:
                raise DatasetError("Accession queries require an exact supported --species")
            response = client.request_json(
                "POST",
                "/api/search/genomes",
                body={
                    "organismId": str(organism["organismId"]),
                    "searchText": ident,
                    "qc": [True, False],
                },
                params={"limit": 100, "sort": "id"},
            )
            count = response.get("meta", {}).get("count")
            searches = response.get("genomes")
            if type(count) is not int or not isinstance(searches, list) or count != len(searches):
                raise DatasetError("Accession search is incomplete; supply an exact source UUID")
        for search in searches:
            detail = client.request_json(
                "GET", "/api/genomes/details", params={"id": search["uuid"]}
            )
            if (
                detail.get("uuid") != search["uuid"]
                or detail.get("projectAccess", "PUBLIC") != "PUBLIC"
            ):
                raise DatasetError("Exact accession detail identity/public-access mismatch")
            row = _canonical({**search, **detail}, origin="local", species=species or "")
            row["organism_id"] = str(detail.get("organismId") or "")
            if row["source_genome_id"] not in {r["source_genome_id"] for r in candidates}:
                candidates.append(row)
    rows = [_join_accession(ident, candidates) for ident in identifiers]
    if overrides:
        rows = apply_metadata_overrides(rows, overrides, path.parent)
    return rows


def enrich_exact_ena(rows, *, fetch=None):
    """Fill missing metadata only from an exact run or BioSample-linked ENA row."""

    def fetch_remote(accession):
        fields = "run_accession,sample_accession,collection_date,country,host,isolation_source"
        url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urlencode(
            {"accession": accession, "result": "read_run", "fields": fields, "format": "json"}
        )
        with urlopen(url, timeout=30) as response:
            return json.load(response)

    fetch = fetch or fetch_remote
    provenance, conflicts, ledger = [], [], []
    for row in rows:
        links = sorted(
            set(row.get("run_accessions") or []) | set(row.get("biosample_accessions") or [])
        )
        if not links:
            ledger.append(
                {
                    "sample_id": row["sample_id"],
                    "operation": "ena_metadata",
                    "status": "unavailable",
                    "reason": "no_exact_linked_accession",
                }
            )
            continue
        found = []
        for accession in links:
            try:
                response = fetch(accession)
                if not isinstance(response, list):
                    raise TypeError("Unsupported ENA response")
                exact = [
                    item
                    for item in response
                    if accession in {item.get("run_accession"), item.get("sample_accession")}
                ]
                if len(exact) != len(response):
                    raise ValueError("ENA response contains an unrelated identity")
                found.extend(exact)
                ledger.append(
                    {
                        "sample_id": row["sample_id"],
                        "operation": "ena_metadata",
                        "accession": accession,
                        "status": "retrieved" if exact else "unavailable",
                        "records": len(exact),
                    }
                )
            except (OSError, ValueError, TypeError):
                ledger.append(
                    {
                        "sample_id": row["sample_id"],
                        "operation": "ena_metadata",
                        "accession": accession,
                        "status": "failed",
                        "reason": "exact_link_metadata_request_failed",
                    }
                )
        for field in ("collection_date", "country", "host", "isolation_source"):
            values = sorted({str(item[field]) for item in found if item.get(field)})
            for value in values:
                selected = len(values) == 1 and row.get(field) in {None, "", "Unknown"}
                provenance.append(
                    {
                        "sample_id": row["sample_id"],
                        "field": field,
                        "source": "ENA_exact_link",
                        "value": value,
                        "selected": selected,
                    }
                )
                if (
                    row.get(field) not in {None, "", "Unknown"}
                    and row[field] != value
                    or len(values) > 1
                ):
                    conflicts.append(
                        {
                            "sample_id": row["sample_id"],
                            "field": field,
                            "source": "ENA_exact_link",
                            "value": value,
                        }
                    )
            if len(values) == 1 and row.get(field) in {None, "", "Unknown"}:
                row[field] = normalize_country(values[0]) if field == "country" else values[0]
                if field == "collection_date":
                    for key in ("date_precision", "date_start", "date_end", "date_raw"):
                        row.pop(key, None)
    return provenance, conflicts, ledger


def catalogues_for_records(rows, explicit=()):
    """Use explicit universes; observed server exports are always incomplete."""
    catalogues = list(explicit)
    scopes = {}
    for row in rows:
        scope = typing_scope(row, "cgmlst")
        if not scope:
            if _has_profile(row):
                raise DatasetError("Provider returned an unscoped called profile")
            continue
        key = (
            *scope,
            row.get("cgmlst_database_version") or None,
            row.get("cgmlst_database_sha256") or None,
        )
        scopes.setdefault(key, []).append(row)
    for (scheme, version, database_version, digest), members in scopes.items():
        if any(
            (c.scheme_id, c.scheme_version, c.database_version, c.database_sha256)
            == (scheme, version, database_version, digest)
            for c in catalogues
        ):
            continue
        declared = [tuple(r["cgmlst_loci"]) for r in members if r.get("cgmlst_loci")]
        loci = sorted(
            set().union(
                *(
                    set(r.get("cgmlst_profile") or {})
                    | set(r.get("cgmlst_novel_alleles") or {})
                    | set(r.get("cgmlst_loci") or [])
                    for r in members
                )
            )
        )
        if loci:
            complete = bool(
                declared
                and all(set(v) == set(loci) for v in declared)
                and all(r.get("cgmlst_locus_universe_complete") is True for r in members)
            )
            catalogues.append(
                LocusCatalogue(
                    scheme,
                    version,
                    tuple(loci),
                    "provider_frozen_evidence",
                    database_version=database_version,
                    complete=complete,
                    database_sha256=digest,
                )
            )
    return tuple(catalogues)


def assign_existing_profiles(rows, config, output):
    """Run the configured assigner on compatible existing profiles, without an assembly."""
    from chronoclade.query_typing import (
        _assignment_result,
        _provenance,
        _run,
        normalise_typing_species,
    )

    organisms = {normalise_typing_species(k): v for k, v in config["organisms"].items()}
    assigned = set()
    for row in rows:
        organism = organisms.get(normalise_typing_species(row["species"]))
        if not organism or not _has_profile(row) or _has_lineage(row):
            continue
        caller, assigner = organism["cgmlst"], organism["assignment"]
        if typing_scope(row, "cgmlst") != (caller["scheme"], caller["scheme_version"]):
            continue
        caller_provenance = _provenance(caller, ("index_dir",))
        digest = caller_provenance["databases"]["index_dir"]["sha256"]
        if row.get("cgmlst_database_sha256") and row["cgmlst_database_sha256"] != digest:
            continue
        loci = assigner["canonical_loci"]
        calls = dict(row.get("cgmlst_profile") or {}, **(row.get("cgmlst_novel_alleles") or {}))
        if set(calls) - set(loci):
            continue
        ambiguous = set(row.get("cgmlst_ambiguous_loci") or [])
        raw = {
            "scheme": caller["scheme"],
            "schemeSize": len(loci),
            "genes": loci,
            "code": "_".join(
                str(calls.get(locus) or "") if locus not in ambiguous else "" for locus in loci
            ),
            "st": "unassigned",
        }
        directory = output / hashlib.sha256(row["sample_id"].encode()).hexdigest()[:24]
        directory.mkdir(parents=True, exist_ok=True)
        if assigner["kind"] == "plincer":
            args = [
                "classify",
                "-",
                "--scheme-file",
                str(assigner["scheme_file"]),
                "--profiles-file",
                str(assigner["profiles_file"]),
                "--alleles-db",
                str(assigner["alleles_db"]),
            ]
            fields = ("scheme_file", "profiles_file", "alleles_db", "loci_file")
        else:
            args = ["assign", "-", "--exact", "--reference-db", str(assigner["reference_db"])]
            fields = ("reference_db", "loci_file")
        assignment = _run(
            assigner["command"] + args, json.dumps(raw).encode(), directory / "assignment.json"
        )
        row.update(_assignment_result(assignment, assigner))
        row["typing_assignment_evidence"] = assignment
        row["typing_assignment_provenance"] = _provenance(assigner, fields)
        kind = "cglin" if assigner["kind"] == "plincer" else "hiercc"
        row[kind + "_database_sha256"] = hashlib.sha256(
            json.dumps(
                {
                    key: info["sha256"]
                    for key, info in row["typing_assignment_provenance"]["databases"].items()
                },
                sort_keys=True,
            ).encode()
        ).hexdigest()
        assigned.add(row["sample_id"])
    return assigned


def resolve_queries(
    input_value,
    work,
    *,
    input_kind,
    species=None,
    metadata=None,
    client=None,
    typing_config=None,
    query_typing=None,
    public_typing=None,
    cglin_export=None,
    enrich_metadata=True,
):
    """Resolve and type inputs, without public context discovery."""
    client = client or PathogenwatchClient()
    work.mkdir(parents=True, exist_ok=True)
    provider = {"input_kind": input_kind, "context_discovery": "not_run"}
    if input_kind == "collection":
        rows, details = _collection(client, collection_id(str(input_value)), work)
        provider.update(details)
        rows = _analysis_exports(
            rows, details["downloads"], client, work / "query_exports", provider
        )
    elif input_kind == "accessions":
        rows = _resolve_accessions(Path(input_value).resolve(), species, client)
    elif input_kind == "assemblies":
        path = Path(input_value).resolve()
        rows = [
            _canonical(r, origin="local", species=species or "", base=path.parent)
            for r in _csv(path)
        ]
        if any(not r.get("assembly") for r in rows):
            raise DatasetError("Assembly metadata requires an existing assembly on every row")
    else:
        raise DatasetError("Unsupported live preparation input kind")
    if metadata:
        rows = apply_metadata_overrides(rows, _csv(Path(metadata)), Path(metadata).resolve().parent)
    if not rows or len({r["sample_id"] for r in rows}) != len(rows):
        raise DatasetError("Inputs require nonempty unique exact sample identities")
    if species:
        for row in rows:
            if row.get("species") and row["species"].casefold() != species.casefold():
                raise DatasetError("Input species differs from --species")
            row["species"] = species
    if any(not r.get("species") for r in rows):
        raise DatasetError("Each input needs --species or a provider species identity")
    rows = _grouped_analysis_exports(rows, client, work / "query_exports_grouped", provider)
    rows = _grouped_analysis_exports(
        rows, client, work / "query_lineages", provider, download_names={"klebsiella-lincodes"}
    )
    if public_typing:
        from chronoclade.public_typing import annotate_public_typing, load_public_typing

        rows = annotate_public_typing(rows, load_public_typing(Path(public_typing)))
    if cglin_export:
        from chronoclade.cglin import load_cglin_export

        rows = _frozen_cglin(rows, load_cglin_export(Path(cglin_export)))
    if query_typing:
        from chronoclade.query_typing import load_query_typing

        assignments = {r["sample_id"]: r for r in load_query_typing(Path(query_typing))}
        if set(assignments) - {r["sample_id"] for r in rows}:
            raise DatasetError("Query typing contains an identity outside the input cohort")
        for row in rows:
            assignment = assignments.get(row["sample_id"])
            if assignment:
                if assignment["species"].casefold() != row["species"].casefold() or (
                    row.get("assembly_sha256")
                    and row["assembly_sha256"] != assignment["assembly_sha256"]
                ):
                    raise DatasetError(
                        "Query typing species/assembly hash differs from the exact input"
                    )
                row.update(assignment)
    operations = {
        "provider_retrieval": "run_exact_queries",
        "context_discovery": "not_run",
        "typing": "not_run",
        "lineage_assignment": "not_run",
        "typing_tool_installation": "not_configured",
        "reference_database": "unavailable_not_configured",
    }
    ledger = []
    if typing_config:
        from chronoclade.query_typing import (
            QueryTypingError,
            load_typing_config,
            type_query_assemblies,
        )

        try:
            config = load_typing_config(Path(typing_config))
        except QueryTypingError:
            operations.update(
                reference_database="unavailable_configured_database_or_tool_not_ready",
                typing_tool_installation="configured_not_ready",
            )
            setup_manifest = Path(typing_config).resolve().parent / "typing_setup.json"
            if setup_manifest.is_file():
                setup = json.loads(setup_manifest.read_text())
                if all(
                    t.get("status") == "installed" for t in setup.get("tools", {}).values()
                ) and setup.get("tools"):
                    operations["typing_tool_installation"] = (
                        "installed_reference_database_unavailable"
                    )
        else:
            operations.update(
                reference_database="ready_validated_configuration",
                typing_tool_installation="configured_available",
            )
            assigned = assign_existing_profiles(rows, config, work / "native_lineage_assignment")
            if assigned:
                operations["lineage_assignment"] = "run_existing_compatible_profiles"
            need = [
                r
                for r in rows
                if r["sample_id"] not in assigned and (not _has_profile(r) or not _has_lineage(r))
            ]
            if need:
                samples = materialise_assemblies(
                    need, output=work / "query_assemblies", client=client
                )
                assignments = type_query_assemblies(
                    samples, Path(typing_config), work / "native_typing"
                )
                by_id = {r["sample_id"]: r for r in assignments}
                for row, sample in zip(need, samples, strict=True):
                    row.update(by_id[row["sample_id"]])
                    row["assembly"] = str(sample.assembly)
                    # Caller supplies its ordered full universe, distinct from an observed export.
                    row["cgmlst_locus_universe_complete"] = True
                operations.update(
                    typing="run_missing_profile_or_lineage_inputs",
                    lineage_assignment="run_missing_profile_or_lineage_inputs",
                )
    for row in rows:
        row.update(role="input", identity_resolved=bool(row.get("source_genome_id")))
        ledger.append(
            {
                "sample_id": row["sample_id"],
                "operation": "prepare_typing",
                "status": "ready" if _has_profile(row) and _has_lineage(row) else "incomplete",
                "reason": None
                if _has_profile(row) and _has_lineage(row)
                else operations["reference_database"],
            }
        )
    provenance = [
        {
            "sample_id": r["sample_id"],
            "field": key,
            "source": "user"
            if key in r.get("prepare_user_metadata_fields", []) or input_kind == "assemblies"
            else "input_provider",
            "value": r[key],
            "selected": True,
        }
        for r in rows
        for key in ("country", "region", "nuts2", "host", "isolation_source", "collection_date")
        if r.get(key)
    ]
    conflicts = [
        {
            "sample_id": r["sample_id"],
            "field": c["field"],
            "source": "provider",
            "value": c["source_value"],
        }
        for r in rows
        for c in r.get("metadata_override_conflicts", [])
    ]
    if enrich_metadata:
        extra, clashes, ena_ledger = enrich_exact_ena(rows)
        provenance.extend(extra)
        conflicts.extend(clashes)
        ledger.extend(ena_ledger)
    _save(
        work / "resolved_queries.json",
        {
            "records": rows,
            "provenance": provenance,
            "conflicts": conflicts,
            "retrieval": ledger,
            "parameters": {"provider_prepare": provider, "prepare_operations": operations},
        },
    )
    return rows, provenance, conflicts, ledger, provider, operations
