import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

import pytest

from chronoclade.pathogenwatch import PathogenwatchClient
from chronoclade.profile_inputs import (
    ProfileInputError,
    collection_id,
    materialise_assemblies,
    parse_cgmlst_export,
    resolve_profile_inputs,
)

UUID = "jX5cwsoyJ1KDquMqssUzAD"


def client_fixture(members, *, api_key="", downloads=None, size=None):
    calls = []

    def transport(method, path, *, params, headers, body):
        calls.append((path, params, headers, body))
        if path == "/api/collections/details":
            return {
                "uuid": UUID,
                "size": len(members) if size is None else size,
                "organismName": "Staphylococcus aureus",
                "organismId": "1280",
                "access": "PUBLIC",
                "downloads": downloads or [],
            }
        if path == "/api/collections/genomes":
            page = params["page"]
            return {"genomes": members[page - 1 : page], "hasMore": page < len(members)}
        if path == "/api/organisms/supported":
            return []
        raise AssertionError(path)

    return PathogenwatchClient(api_key=api_key, transport=transport), calls


@pytest.fixture(autouse=True)
def no_implicit_live_client(monkeypatch):
    """Every network-like interaction in this module must use fixture transports."""
    monkeypatch.setattr(
        "chronoclade.profile_inputs.PathogenwatchClient", lambda: client_fixture([])[0]
    )


def member(ident="x", **extra):
    return {
        "uuid": ident,
        "id": 17,
        "name": "Visible name",
        "species": "Staphylococcus aureus",
        "typing/MLST/ST": "45",
        "metadata/User defined metadata/Collection Date": "2012",
        "metadata/User defined metadata/Country": "GB",
        **extra,
    }


@pytest.mark.parametrize(
    "text",
    [
        UUID,
        f"https://pathogen.watch/collections/{UUID}-readable-slug",
        f"https://next.pathogen.watch/collections/{UUID}",
    ],
)
def test_collection_identifier(text):
    assert collection_id(text) == UUID


@pytest.mark.parametrize(
    "text",
    [
        "https://evil.example/collections/" + UUID,
        "https://pathogen.watch/genomes/" + UUID,
        "https://user:pass@pathogen.watch/collections/" + UUID,
        "bad",
    ],
)
def test_collection_identifier_rejects_untrusted_url(text):
    with pytest.raises(ProfileInputError):
        collection_id(text)


def test_collection_all_members_focal_partial_dates_and_missing_profiles(tmp_path):
    client, calls = client_fixture(
        [member("one"), member("two", **{"metadata/User defined metadata/Collection Date": ""})]
    )
    result = resolve_profile_inputs(None, collection=UUID, output=tmp_path, client=client)
    assert [row["sample_id"] for row in result["queries"]] == ["PW_one", "PW_two"]
    assert all(row["origin"] == "local" for row in result["queries"])
    assert result["queries"][0]["date_precision"] == "year"
    assert result["queries"][1]["date_precision"] == "missing"
    assert result["queries"][0]["country"] == "United Kingdom"
    assert result["queries"][0]["country_raw"] == "GB"
    assert result["provenance"]["coverage"]["queries"]["profiles_missing"] == 2
    assert not any("download" in call[0] for call in calls)
    assert (tmp_path / "collection_snapshot.json").exists()


def test_collection_auth_is_not_used_for_context_discovery(tmp_path):
    client, calls = client_fixture([member()], api_key="fixture-key")
    resolve_profile_inputs(None, collection=UUID, output=tmp_path, client=client)
    assert all(
        call[2].get("X-API-Key") == "fixture-key"
        for call in calls
        if call[0].startswith("/api/collections")
    )
    assert all(
        "X-API-Key" not in call[2] for call in calls if not call[0].startswith("/api/collections")
    )


def test_collection_size_mismatch_is_not_partial_success(tmp_path):
    client, _ = client_fixture([member()], size=2)
    with pytest.raises(ProfileInputError, match="reconcile"):
        resolve_profile_inputs(None, collection=UUID, output=tmp_path, client=client)


def export(rows):
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=["Genome ID", "Gene", "Allele ID"])
    writer.writeheader()
    writer.writerows(rows)
    return gzip.compress(stream.getvalue().encode())


def test_verified_export_shape_exact_ids_and_ambiguous_loci():
    raw = export(
        [
            {"Genome ID": "A", "Gene": "g1", "Allele ID": "2"},
            {"Genome ID": "A", "Gene": "g2", "Allele ID": "3"},
            {"Genome ID": "A", "Gene": "g2", "Allele ID": "4"},
            {"Genome ID": "A", "Gene": "g3", "Allele ID": "0"},
        ]
    )
    row = parse_cgmlst_export(raw, ["A"], job="cgmlst-1280-2")[0]
    assert row["cgmlst_profile"] == {"g1": "2"}
    assert row["cgmlst_ambiguous_loci"] == ["g2"]
    assert row["cgmlst_loci"] == ["g1", "g2", "g3"]
    assert row["cgmlst_locus_universe_complete"] is False
    assert "cgmlst_called_fraction" not in row
    with pytest.raises(ProfileInputError, match="unrequested"):
        parse_cgmlst_export(raw, ["B"], job="cgmlst-1280-2")


def test_collection_exports_are_profiles_not_assemblies(tmp_path, monkeypatch):
    client, calls = client_fixture(
        [member("A")], api_key="fixture", downloads=[{"name": "cgmlst", "job": "cgmlst-1280-2"}]
    )

    def fetch(route, **kwargs):
        assert route == "/api/downloads/cgmlst?job=cgmlst-1280-2"
        assert kwargs["body"] == {"ids": "17"}
        return export([{"Genome ID": "A", "Gene": "g1", "Allele ID": "2"}])

    monkeypatch.setattr("chronoclade.profile_inputs.request_download", fetch)
    result = resolve_profile_inputs(None, collection=UUID, output=tmp_path, client=client)
    assert result["queries"][0]["profile_status"] == "available"
    assert "assembly" not in result["queries"][0]


def test_manual_metadata_conflicts_are_audited(tmp_path):
    path = tmp_path / "metadata.csv"
    path.write_text(
        "source_genome_id,collection_date,country,region,nuts2\nA,2013-02,France,Paris,FR10\n"
    )
    client, _ = client_fixture([member("A")])
    result = resolve_profile_inputs(path, collection=UUID, output=tmp_path / "out", client=client)
    row = result["queries"][0]
    assert row["collection_date"] == "2013-02"
    assert row["country"] == "France"
    assert row["nuts2"] == "FR10"
    assert any(item["field"] == "collection_date" for item in row["metadata_override_conflicts"])


def test_local_assembly_import_hashes_and_never_needs_native_setup_if_typed(tmp_path, monkeypatch):
    assembly = tmp_path / "q.fa"
    assembly.write_text(">q\nATGC\n")
    metadata = tmp_path / "metadata.csv"
    metadata.write_text(
        "sample_id,assembly,species,lineage,collection_date,location\nq,q.fa,Klebsiella pneumoniae,ST147,,GB\n"
    )
    typing = tmp_path / "typing.json"
    typing.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "assignments": [
                    {
                        "sample_id": "q",
                        "species": "Klebsiella pneumoniae",
                        "assembly_sha256": hashlib.sha256(assembly.read_bytes()).hexdigest(),
                        "cgmlst_scheme": "klebsiella_1",
                        "cgmlst_scheme_version": "v1",
                        "cgmlst_loci": ["g1"],
                        "cgmlst_profile": {"g1": "1"},
                    }
                ],
            }
        )
    )
    config = tmp_path / "config.json"
    config.write_text("{}")
    client, _ = client_fixture([])
    monkeypatch.setattr(
        "chronoclade.profile_inputs.type_query_assemblies",
        lambda *a, **k: pytest.fail("already typed"),
    )
    result = resolve_profile_inputs(
        metadata, query_typing=typing, typing_config=config, client=client, output=tmp_path / "out"
    )
    row = result["queries"][0]
    assert row["profile_status"] == "available"
    assert row["assembly"] == str(assembly)
    samples = materialise_assemblies([row], output=tmp_path / "out")
    assert samples[0].sample_id == "q"
    assert (
        result["provenance"]["inputs"]["metadata"]["sha256"]
        == hashlib.sha256(metadata.read_bytes()).hexdigest()
    )


def test_deferred_assembly_acquisition_only_selected_rows(tmp_path, monkeypatch):
    seen = []

    def download(rows, **kwargs):
        seen.extend(row["source_genome_id"] for row in rows)
        path = tmp_path / "a.fa"
        path.write_text(">a\nACGT\n")
        return {"A": path}, [{"status": "downloaded"}]

    monkeypatch.setattr("chronoclade.profile_inputs.download_assemblies", download)
    client, _ = client_fixture([], api_key="fixture")
    records = [
        {
            "sample_id": "PW_A",
            "source_genome_id": "A",
            "origin": "context",
            "species": "E. coli",
            "lineage": "ST1",
        }
    ]
    result = materialise_assemblies(records, output=tmp_path, client=client)
    assert seen == ["A"]
    assert result[0].origin == "context"
    assert "assembly" not in records[0]


def frozen_catalogue(tmp_path):
    def transport(method, path, *, params, headers, body):
        assert "X-API-Key" not in headers
        if path == "/api/organisms/supported":
            return [{"organismId": "573", "fullName": "Klebsiella pneumoniae", "typing": ["MLST"]}]
        if path == "/api/search/genomes":
            return {
                "meta": {"count": 2, "endCursor": "B"},
                "genomes": [
                    {"uuid": "A", "projectAccess": "PUBLIC", "organism": "Klebsiella pneumoniae"},
                    {"uuid": "B", "projectAccess": "PUBLIC", "organism": "Klebsiella pneumoniae"},
                ],
            }
        if path == "/api/genomes/details":
            ident = params["id"]
            return {
                "uuid": ident,
                "id": 1 if ident == "A" else 2,
                "organismId": "573",
                "mlst": "147",
                "sampleAccession": "SAMN123" if ident == "A" else "SAMN456",
                "runAccession": "SRR123" if ident == "A" else "SRR456",
                "qc": True,
                "metadata": {"Country": "GB", "Date": "2012"},
                "name": "Name with spaces",
            }
        raise AssertionError(path)

    client = PathogenwatchClient(api_key="", transport=transport)
    path = tmp_path / "catalogue.json"
    client.freeze_catalogue(path)
    return path


def test_offline_accessions_exact_identity_imported_profiles(tmp_path, monkeypatch):
    catalogue = frozen_catalogue(tmp_path)
    accessions = tmp_path / "accessions.txt"
    accessions.write_text("SRR123\n")
    profiles = tmp_path / "profiles.json"
    profiles.write_text(
        json.dumps(
            [
                {
                    "source_genome_id": ident,
                    "cgmlst_scheme": "klebsiella_1",
                    "cgmlst_scheme_version": "v1",
                    "cgmlst_loci": ["g1", "g2"],
                    "cgmlst_profile": {"g1": "1", "g2": "2"},
                }
                for ident in ["A", "B"]
            ]
        )
    )
    monkeypatch.setattr(
        "chronoclade.profile_inputs.PathogenwatchClient",
        lambda: pytest.fail("offline must not create client"),
    )
    result = resolve_profile_inputs(
        None,
        accessions=accessions,
        catalogue=catalogue,
        public_typing=profiles,
        output=tmp_path / "out",
    )
    assert result["queries"][0]["source_genome_id"] == "A"
    assert result["queries"][0]["origin"] == "local"
    assert result["queries"][0]["profile_status"] == "available"
    assert result["context"][0]["source_genome_id"] == "B"
    assert result["context"][0]["sample_id"] == "PW_B"
    assert "assembly" not in result["context"][0]
    assert {row["source_genome_id"] for row in result["catalogue_rows"]} == {"A", "B"}
    assert all(row["cgmlst_profile"] == {"g1": "1", "g2": "2"} for row in result["catalogue_rows"])


def test_accession_never_substitutes_a_near_match(tmp_path):
    catalogue = frozen_catalogue(tmp_path)
    accessions = tmp_path / "accessions.txt"
    accessions.write_text("SRR12\n")
    with pytest.raises(ProfileInputError, match="0 exact matches"):
        resolve_profile_inputs(
            None,
            accessions=accessions,
            catalogue=catalogue,
            output=tmp_path / "out",
            client=client_fixture([])[0],
        )


def test_failed_export_keeps_metadata_explicitly_missing(tmp_path, monkeypatch):
    from chronoclade.pathogenwatch_download import DownloadError

    client, _ = client_fixture(
        [member("A")], api_key="fixture", downloads=[{"name": "cgmlst", "job": "cgmlst-1280-2"}]
    )

    def fetch(*args, **kwargs):
        raise DownloadError("Pathogenwatch download HTTP 401")

    monkeypatch.setattr("chronoclade.profile_inputs.request_download", fetch)
    result = resolve_profile_inputs(None, collection=UUID, output=tmp_path, client=client)
    assert result["provenance"]["analysis_exports"][0]["status"] == "unavailable"
    assert result["queries"][0]["profile_status"] == "missing"
    assert "assembly" not in result["queries"][0]


def test_explicit_provided_context_stays_context_and_keeps_assembly(tmp_path):
    catalogue = frozen_catalogue(tmp_path)
    for ident in ["q", "ctx"]:
        (tmp_path / f"{ident}.fa").write_text(f">{ident}\nACGT\n")
    metadata = tmp_path / "combined_metadata.csv"
    metadata.write_text(
        "sample_id,assembly,species,lineage,origin,source_genome_id\n"
        "q,q.fa,Klebsiella pneumoniae,ST147,local,A\n"
        "ctx,ctx.fa,Klebsiella pneumoniae,ST999,context,C\n"
    )
    result = resolve_profile_inputs(
        metadata, catalogue=catalogue, output=tmp_path / "out", profile_limit=0
    )
    assert [row["sample_id"] for row in result["queries"]] == ["q"]
    provided = next(row for row in result["context"] if row["sample_id"] == "ctx")
    assert provided["origin"] == "context"
    assert provided["provided_context"] is True
    assert provided["lineage"] == "ST999"
    assert provided["assembly"] == str(tmp_path / "ctx.fa")
    assert result["provenance"]["provided_context_count"] == 1


def test_accession_csv_metadata_is_applied_to_exact_identity(tmp_path):
    catalogue = frozen_catalogue(tmp_path)
    path = tmp_path / "accessions.csv"
    path.write_text("accession,sample_id,collection_date,region\nSRR123,my_query,2015-06,Oxford\n")
    result = resolve_profile_inputs(
        None,
        accessions=path,
        catalogue=catalogue,
        output=tmp_path / "out",
        client=client_fixture([])[0],
    )
    query = result["queries"][0]
    assert query["source_genome_id"] == "A"
    assert query["sample_id"] == "my_query"
    assert query["collection_date"] == "2015-06"
    assert query["date_precision"] == "month"
    assert query["region"] == "Oxford"


def test_imported_collection_profiles_skip_unnecessary_export(tmp_path, monkeypatch):
    client, _ = client_fixture(
        [member("A")], api_key="fixture", downloads=[{"name": "cgmlst", "job": "cgmlst-1280-2"}]
    )
    profiles = tmp_path / "profiles.json"
    profiles.write_text(
        json.dumps(
            [
                {
                    "source_genome_id": "A",
                    "cgmlst_scheme": "s_aureus_core2208",
                    "cgmlst_scheme_version": "v1",
                    "cgmlst_loci": ["g1"],
                    "cgmlst_profile": {"g1": "1"},
                }
            ]
        )
    )
    monkeypatch.setattr(
        "chronoclade.profile_inputs.request_download",
        lambda *a, **k: pytest.fail("already imported"),
    )
    result = resolve_profile_inputs(
        None, collection=UUID, public_typing=profiles, output=tmp_path / "out", client=client
    )
    assert result["queries"][0]["cgmlst_scheme"] == "s_aureus_core2208"
    assert result["queries"][0]["profile_status"] == "available"


def test_zero_profile_limit_skips_public_discovery(tmp_path, monkeypatch):
    client, calls = client_fixture([member("A")])
    monkeypatch.setattr(
        client,
        "supported_organisms",
        lambda: pytest.fail("zero budget must skip catalogue discovery"),
    )
    result = resolve_profile_inputs(
        None, collection=UUID, output=tmp_path, client=client, profile_limit=0
    )
    assert result["context"] == []
    assert result["catalogue_rows"] == []
    assert len(result["queries"]) == 1
    assert all(call[0].startswith("/api/collections/") for call in calls)


def test_context_analysis_pool_is_bounded_full_metadata_is_retained(tmp_path):
    catalogue = frozen_catalogue(tmp_path)
    path = tmp_path / "accessions.txt"
    path.write_text("SRR123\n")
    result = resolve_profile_inputs(
        None, accessions=path, catalogue=catalogue, output=tmp_path / "out", profile_limit=0
    )
    assert result["context"] == []
    assert len(result["catalogue_rows"]) == 2
    assert result["provenance"]["eligible_public_context_count"] == 1
    assert result["provenance"]["bounded_public_context_count"] == 0
    result = resolve_profile_inputs(
        None, accessions=path, catalogue=catalogue, output=tmp_path / "out2", profile_limit=1
    )
    assert len(result["context"]) == 1
    assert len(result["catalogue_rows"]) == 2


def test_group_advertises_jobs_for_accession_query_and_bounded_context(tmp_path, monkeypatch):
    catalogue = frozen_catalogue(tmp_path)
    accessions = tmp_path / "accessions.txt"
    accessions.write_text("SRR123\n")
    config = tmp_path / "not_ready.json"
    config.write_text("{}")
    calls = []

    def transport(method, path, *, body, params, headers):
        assert path == "/api/genomes/group"
        assert headers["X-API-Key"] == "fixture-key"
        calls.append(body["ids"])
        return [
            {
                "organismId": "573",
                "ids": body["ids"],
                "downloads": [{"name": "cgmlst", "job": "cgmlst-573-2"}],
            }
        ]

    client = PathogenwatchClient(api_key="fixture-key", transport=transport)

    def fetch(route, **kwargs):
        assert route == "/api/downloads/cgmlst?job=cgmlst-573-2"
        return export(
            [
                {"Genome ID": "A" if ident == "1" else "B", "Gene": "g1", "Allele ID": "1"}
                for ident in kwargs["body"]["ids"].split(",")
            ]
        )

    monkeypatch.setattr("chronoclade.profile_inputs.request_download", fetch)
    monkeypatch.setattr(
        "chronoclade.profile_inputs.type_query_assemblies",
        lambda *a, **k: pytest.fail("exported query already typed"),
    )
    result = resolve_profile_inputs(
        None,
        accessions=accessions,
        catalogue=catalogue,
        typing_config=config,
        output=tmp_path / "out",
        client=client,
        profile_limit=1,
    )
    assert result["queries"][0]["profile_status"] == "available"
    assert result["context"][0]["profile_status"] == "available"
    assert calls == [[1], [2]]
    assert (tmp_path / "out/query_exports/export_capabilities.json").exists()
    assert (tmp_path / "out/context_exports/export_capabilities.json").exists()


def test_group_unadvertised_job_is_explicit_missing_coverage(tmp_path):
    catalogue = frozen_catalogue(tmp_path)
    accessions = tmp_path / "accessions.txt"
    accessions.write_text("SRR123\n")

    def transport(method, path, *, body, params, headers):
        return [{"organismId": "573", "ids": body["ids"], "downloads": []}]

    client = PathogenwatchClient(api_key="fixture-key", transport=transport)
    result = resolve_profile_inputs(
        None,
        accessions=accessions,
        catalogue=catalogue,
        output=tmp_path / "out",
        client=client,
        profile_limit=0,
    )
    assert result["queries"][0]["profile_status"] == "missing"
    assert (
        result["provenance"]["profile_export_unavailable"][0]["reason"]
        == "server_does_not_advertise_cgmlst_job"
    )


@pytest.mark.parametrize("kind", ["cglin", "hiercc"])
def test_one_slot_context_budget_prioritises_compatible_query_group(kind):
    from chronoclade.profile_inputs import _refined_context_pool

    species = "Klebsiella pneumoniae" if kind == "cglin" else "Escherichia coli"
    query = {
        "sample_id": "q",
        "species": species,
        "lineage": "ST1",
        kind + "_scheme": "scgMLST629_S" if kind == "cglin" else "ecoli_hiercc",
        kind + "_scheme_version": "v1",
    }
    if kind == "cglin":
        query["cglin_raw"] = "1_2_3_4_5_6_7_8_9_10"
    else:
        query["hiercc_codes"] = {"HC1100": "7"}
    near = dict(
        query, sample_id="PW_Z", source_genome_id="Z", country="France", collection_date="2019"
    )
    far = dict(
        query,
        sample_id="PW_A",
        source_genome_id="A",
        country="United Kingdom",
        collection_date="2012",
    )
    if kind == "cglin":
        far["cglin_raw"] = "9_2_3_4_5_6_7_8_9_10"
    else:
        far["hiercc_codes"] = {"HC1100": "99"}
    chosen, audit = _refined_context_pool([far, near], [query], limit=1, seed=1)
    assert [row["source_genome_id"] for row in chosen] == ["Z"]
    assert chosen[0]["profile_pool_selection_reason"] == "lineage_priority:q"
    assert audit["lineages"][0]["audit"]["counts"]["priority_selected"] == 1
    near[kind + "_scheme_version"] = "different-v2"
    _, audit = _refined_context_pool([far, near], [query], limit=1, seed=1)
    assert audit["lineages"][0]["audit"]["counts"]["priority_selected"] == 0
    assert (
        audit["lineages"][0]["audit"]["query_status"]["q"]
        == "no_comparable_evidence_same_st_fallback"
    )


@pytest.mark.parametrize("header,identity", [("accession", "SRR123"), ("source_genome_id", "A")])
def test_single_column_accession_csv_header_is_not_an_identity(tmp_path, header, identity):
    catalogue = frozen_catalogue(tmp_path)
    path = tmp_path / "accessions.csv"
    path.write_text(header + "\n" + identity + "\n")
    result = resolve_profile_inputs(
        None, accessions=path, catalogue=catalogue, output=tmp_path / "out", profile_limit=0
    )
    assert [row["source_genome_id"] for row in result["queries"]] == ["A"]


def test_explicit_source_uuid_uses_exact_details_without_search_or_species_guess(tmp_path):
    ident = "ooy8fFUy3fXqa7WaWX15eS"
    path = tmp_path / "sources.csv"
    path.write_text("source_genome_id\n" + ident + "\n")
    calls = []

    def transport(method, route, *, body, params, headers):
        calls.append(route)
        assert route == "/api/genomes/details"
        assert params == {"id": ident}
        assert "X-API-Key" not in headers
        return {
            "uuid": ident,
            "id": 17,
            "name": "Name with spaces",
            "species": "Staphylococcus aureus",
            "organismId": "1280",
            "mlst": "45",
            "metadata": {"Date": "2012", "Country": "GB"},
        }

    client = PathogenwatchClient(api_key="", transport=transport)
    result = resolve_profile_inputs(
        None, accessions=path, output=tmp_path / "out", client=client, profile_limit=0
    )
    assert calls == ["/api/genomes/details"]
    assert result["queries"][0]["source_genome_id"] == ident
    assert result["queries"][0]["species"] == "Staphylococcus aureus"
    assert result["queries"][0]["lineage"] == "ST45"
    assert result["queries"][0]["collection_date"] == "2012"


def test_query_context_cglin_exports_share_advertised_job_scope_without_database_claim(
    tmp_path, monkeypatch
):
    from chronoclade.cglin import normalise_assignment
    from chronoclade.profile_inputs import _analysis_exports, _refined_context_pool

    job = "lincodes-3390273-2"
    seen = []

    def download(rows, destination, **kwargs):
        seen.append(kwargs["scheme_version"])
        assert kwargs["scheme_version"] == job
        return {
            "assignments": [
                normalise_assignment(
                    {
                        "source_genome_id": row["source_genome_id"],
                        "cglin_raw": "1_2_3_4_5_6_7_8_9_10",
                        "cglin_scheme": "scgMLST629_S",
                        "cglin_scheme_version": kwargs["scheme_version"],
                    }
                )
                for row in rows
            ]
        }

    monkeypatch.setattr("chronoclade.profile_inputs.download_cglin_export", download)
    client = PathogenwatchClient(api_key="fixture")
    common = {"species": "Klebsiella pneumoniae", "lineage": "ST1"}
    query = dict(common, sample_id="q", source_genome_id="A", numeric_source_id=1)
    context = dict(common, sample_id="PW_B", source_genome_id="B", numeric_source_id=2)
    advertised = [{"name": "klebsiella-lincodes", "job": job}]
    provenance = {}
    queries = _analysis_exports([query], advertised, client, tmp_path / "queries", provenance)
    contexts = _analysis_exports([context], advertised, client, tmp_path / "context", provenance)
    chosen, audit = _refined_context_pool(contexts, queries, limit=1, seed=1)
    assert seen == [job, job]
    assert chosen[0]["cglin_scope_source"] == "server_advertised_analysis_job"
    assert chosen[0]["cglin_database_version_status"] == "unavailable"
    assert not chosen[0].get("cglin_database_sha256")
    assert audit["lineages"][0]["audit"]["counts"]["priority_eligible"] == 1


@pytest.mark.parametrize("query_provider", ["export", "native"])
def test_explicit_frozen_cglin_survives_query_typing_and_context_live_exports(
    tmp_path, monkeypatch, query_provider
):
    """Use frozen groups before selection and retain them after allele acquisition."""
    from chronoclade.cglin import normalise_assignment

    assembly = tmp_path / "q.fa"
    assembly.write_text(">q\nACGT\n")
    metadata = tmp_path / "metadata.csv"
    metadata.write_text(
        "sample_id,source_genome_id,numeric_source_id,assembly,species,lineage,collection_date,location\n"
        "q,A,1,q.fa,Klebsiella pneumoniae,ST147,2019,Greece\n"
    )
    frozen = tmp_path / "cglin.csv"
    frozen.write_text(
        "source_genome_id,cglin_raw,cglin_scheme,cglin_scheme_version,cglin_source\n"
        "A,1_2_3_4_5_6_7_8_9_10,scgMLST629_S,unknown,provided_pathogenwatch_export\n"
        "B,9_2_3_4_5_6_7_8_9_10,scgMLST629_S,unknown,provided_pathogenwatch_export\n"
        "C,1_2_3_4_5_6_7_8_9_10,scgMLST629_S,unknown,provided_pathogenwatch_export\n"
    )
    client, _ = client_fixture([], api_key="fixture")
    monkeypatch.setattr(
        client,
        "supported_organisms",
        lambda: [{"fullName": "Klebsiella pneumoniae", "organismId": "573"}],
    )
    monkeypatch.setattr(
        client,
        "freeze_catalogue",
        lambda *a, **k: {
            "rows": [
                {
                    "source_genome_id": ident,
                    "numeric_source_id": number,
                    "species": "Klebsiella pneumoniae",
                    "lineage": "ST147",
                    "country": "Greece",
                    "collection_date": "2019",
                    "aliases": [ident],
                    "qc_pass": True,
                }
                for ident, number in [("B", 2), ("C", 3)]
            ],
            "provenance": {},
        },
    )
    exported = []

    def live_typing(rows, *args, **kwargs):
        result = []
        for row in rows:
            exported.append(row["source_genome_id"])
            assignment = normalise_assignment(
                {
                    "source_genome_id": row["source_genome_id"],
                    "cglin_raw": "8_8_8_8_8_8_8_8_8_8",
                    "cglin_scheme_version": "current-live-release",
                    "cglin_source": "live_export",
                }
            )
            assignment.pop("source_genome_id")
            assignment.update(
                cglin_scope_source="server_advertised_analysis_job",
                cglin_database_sha256="live-database-fingerprint",
            )
            if row["origin"] == "context" or query_provider == "export":
                assignment["cgmlst_profile"] = {"g1": "1"}
            result.append(dict(row, **assignment))
        return result

    monkeypatch.setattr("chronoclade.profile_inputs._grouped_analysis_exports", live_typing)
    typing_config = None
    if query_provider == "native":
        typing_config = tmp_path / "config.json"
        typing_config.write_text("{}")
        monkeypatch.setattr(
            "chronoclade.profile_inputs.type_query_assemblies",
            lambda *a, **k: [
                dict(
                    normalise_assignment(
                        {
                            "source_genome_id": "A",
                            "cglin_raw": "7_7_7_7_7_7_7_7_7_7",
                            "cglin_scheme_version": "native-release",
                        }
                    ),
                    sample_id="q",
                    cgmlst_profile={"g1": "1"},
                )
            ],
        )
    result = resolve_profile_inputs(
        metadata,
        cglin_export=frozen,
        typing_config=typing_config,
        output=tmp_path / "out",
        client=client,
        profile_limit=1,
    )
    assert exported == ["A", "C"]
    assert result["context"][0]["source_genome_id"] == "C"
    assert result["context"][0]["profile_pool_selection_reason"] == "lineage_priority:q"
    for row in result["queries"] + result["context"]:
        assert row["cglin_raw"] == "1_2_3_4_5_6_7_8_9_10"
        assert row["cglin_scheme_version"] == "unknown"
        assert row["cglin_source"] == "provided_pathogenwatch_export"
        assert row["cglin_scope_source"] == "user_frozen_export"
        assert row["cglin_database_version_status"] == "unavailable"
        assert row["cglin_frozen_export_sha256"] == hashlib.sha256(frozen.read_bytes()).hexdigest()
        assert "cglin_database_sha256" not in row
        assert row["cgmlst_profile"] == {"g1": "1"}
        assert row["profile_status"] == "available"
    authority = result["provenance"]["cglin_authority"]
    assert authority["matched_records"] == {"queries": 1, "context": 1, "catalogue_rows": 2}
    assert authority["scheme_versions"] == ["unknown"]
    assert authority["database_version_inferred_from_live_jobs"] is False
    assert {row["source_genome_id"] for row in result["catalogue_rows"]} == {"B", "C"}
    assert all(
        row["cglin_scope_source"] == "user_frozen_export" for row in result["catalogue_rows"]
    )
    assert all(row["cglin_export_record_count"] == 1 for row in result["catalogue_rows"])
    assert all(not row.get("cgmlst_profile") for row in result["catalogue_rows"])


def test_explicit_frozen_cglin_unmatched_ids_cannot_keep_live_namespace():
    from chronoclade.cglin import normalise_assignment
    from chronoclade.profile_inputs import _frozen_cglin

    live = normalise_assignment(
        {
            "source_genome_id": "unmatched",
            "cglin_raw": "1_2_3_4_5_6_7_8_9_10",
            "cglin_scheme_version": "live-version",
            "cgst": "42",
        }
    )
    live.update(cgmlst_profile={"g1": "1"}, cglin_database_sha256="live-fingerprint")
    frozen = normalise_assignment(
        {
            "source_genome_id": "another",
            "cglin_raw": "1_2_3_4_5_6_7_8_9_10",
            "cglin_scheme_version": "unknown",
        }
    )
    row = _frozen_cglin([live], [frozen])[0]
    assert row["cglin_status"] == "missing"
    assert row["cglin_group_5"] == ""
    assert row["cglin_scheme_version"] == "unknown"
    assert row["cglin_scope_source"] == "user_frozen_export"
    assert row["cglin_export_record_count"] == 0
    assert row["cgst"] == ""
    assert "cglin_database_sha256" not in row
    assert row["cgmlst_profile"] == {"g1": "1"}


@pytest.mark.parametrize("mismatch_count", [1, 3])
def test_cgmlst_identity_mismatch_retries_exact_batch_without_accepting_foreign_calls(
    tmp_path, monkeypatch, mismatch_count
):
    from chronoclade.profile_inputs import _analysis_exports

    wrong = export([{"Genome ID": "foreign", "Gene": "g1", "Allele ID": "9"}])
    correct = export([{"Genome ID": "A", "Gene": "g1", "Allele ID": "1"}])
    requests, waits = [], []

    def fetch(route, **kwargs):
        requests.append((route, kwargs["body"]))
        return wrong if len(requests) <= mismatch_count else correct

    monkeypatch.setattr("chronoclade.profile_inputs.request_download", fetch)
    monkeypatch.setattr("chronoclade.profile_inputs.time.sleep", waits.append)
    provenance = {}
    rows = _analysis_exports(
        [{"source_genome_id": "A", "numeric_source_id": 17}],
        [{"name": "cgmlst", "job": "cgmlst-573-2"}],
        PathogenwatchClient(api_key="fixture"),
        tmp_path,
        provenance,
    )
    assert all(
        request == ("/api/downloads/cgmlst?job=cgmlst-573-2", {"ids": "17"}) for request in requests
    )
    rejected = provenance["analysis_export_identity_retries"]
    assert len(rejected) == mismatch_count
    for item in rejected:
        assert item["requested_source_ids"] == ["A"]
        assert item["rejected_export_sha256"] == hashlib.sha256(wrong).hexdigest()
        assert Path(item["rejected_export"]).read_bytes() == wrong
    assert len({item["rejected_export"] for item in rejected}) == mismatch_count
    if mismatch_count == 1:
        assert rows[0]["cgmlst_profile"] == {"g1": "1"}
        assert len(requests) == 2
        assert waits == [1]
        assert provenance["analysis_exports"][-1]["status"] == "downloaded"
    else:
        assert not rows[0].get("cgmlst_profile")
        assert len(requests) == 3
        assert waits == [1, 2]
        assert provenance["analysis_exports"][-1]["status"] == "unavailable"


def test_cgmlst_malformed_locus_is_not_retried_as_transient_identity(tmp_path, monkeypatch):
    from chronoclade.profile_inputs import _analysis_exports

    requests = []

    def fetch(*args, **kwargs):
        requests.append(kwargs["body"])
        return export([{"Genome ID": "A", "Gene": "", "Allele ID": "1"}])

    monkeypatch.setattr("chronoclade.profile_inputs.request_download", fetch)
    provenance = {}
    result = _analysis_exports(
        [{"source_genome_id": "A", "numeric_source_id": 17}],
        [{"name": "cgmlst", "job": "cgmlst-573-2"}],
        PathogenwatchClient(api_key="fixture"),
        tmp_path,
        provenance,
    )
    assert requests == [{"ids": "17"}]
    assert not result[0].get("cgmlst_profile")
    assert "analysis_export_identity_retries" not in provenance
    assert provenance["analysis_exports"][-1]["reason"] == "cgMLST export contains an empty locus"
