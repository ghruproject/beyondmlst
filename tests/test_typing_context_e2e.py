"""Offline E. coli typing-to-context integration through the real CLI/provider."""

import csv
import hashlib
import json

from typer.testing import CliRunner

from chronoclade.cli import app
from chronoclade.metadata import read_metadata
from chronoclade.pathogenwatch import PathogenwatchClient
import chronoclade.pathogenwatch_context as provider
import chronoclade.query_typing as typing

SPECIES = "Escherichia coli"
LOCI = [f"locus{i}" for i in range(10)]


def assignment(ident, cluster="15", allele="1", *, with_profile=True):
    row = dict(
        source_genome_id=ident,
        hiercc_scheme="ecoli_hiercc",
        hiercc_scheme_version="2026-01-01",
        hiercc_codes={"HC1100": cluster, "HC10": "123"},
    )
    if with_profile:
        row.update(
            cgmlst_scheme="ecoli_1",
            cgmlst_scheme_version="2026-01-01",
            cgmlst_loci=LOCI,
            cgmlst_profile={locus: allele for locus in LOCI},
        )
    return row


def fixture(tmp_path, *, lookup=False, with_profile=True):
    records = {
        "self": ("SAMN90001", "UK"),
        "match1": ("SAMN90002", "UK"),
        "match2": ("SAMN90003", "India"),
        "wide1": ("SAMN90004", "France"),
        "wide2": ("SAMN90005", "Nigeria"),
        "wide3": ("SAMN90006", "Japan"),
        "wide4": ("SAMN90007", "USA"),
    }

    def transport(method, path, *, body, params, headers):
        if path.endswith("/supported"):
            return [
                {"organismId": "562", "fullName": SPECIES, "typing": ["MLST"], "other": ["HierCC"]}
            ]
        if path.endswith("/search/genomes"):
            assert body["organismId"] == "562"
            assert body["mlst"] == ["131"]
            return {
                "meta": {"count": len(records), "endCursor": "end"},
                "genomes": [
                    {
                        "uuid": ident,
                        "projectAccess": "PUBLIC",
                        "name": accession,
                        "organism": SPECIES,
                        "location": None,
                    }
                    for ident, (accession, _) in records.items()
                ],
            }
        ident = params["id"]
        accession, country = records[ident]
        return {
            "uuid": ident,
            "id": list(records).index(ident) + 1,
            "organismId": "562",
            "mlst": "131",
            "name": accession,
            "sampleAccession": accession,
            "qc": True,
            "metadata": [["Country", country], ["Date", "2020"]],
        }

    catalogue = tmp_path / "catalogue.json"
    PathogenwatchClient(transport=transport).freeze_catalogue(
        catalogue, organism_id="562", st="131"
    )
    public = tmp_path / "public.json"
    public.write_text(
        json.dumps(
            {
                "assignments": [
                    assignment(
                        ident,
                        "15" if ident.startswith("match") or ident == "self" else "99",
                        "1" if ident.startswith("match") or ident == "self" else "2",
                        with_profile=with_profile,
                    )
                    for ident in records
                ]
            }
        )
    )
    assembly = tmp_path / "query.fasta"
    assembly.write_text(">chromosome\n" + "ACGT" * 100 + "\n")
    sample_id = "SAMN90001" if lookup else "F1"
    focal = tmp_path / "focal.csv"
    with focal.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["sample_id", "assembly", "collection_date", "location", "species", "lineage", "origin"]
        )
        writer.writerow([sample_id, str(assembly), "2025", "UK", SPECIES, "ST131", "local"])
    query = assignment(sample_id, with_profile=with_profile)
    query.pop("source_genome_id")
    query.update(
        sample_id=sample_id,
        species=SPECIES,
        assembly_sha256=hashlib.sha256(assembly.read_bytes()).hexdigest(),
    )
    query_path = tmp_path / "query.json"
    query_path.write_text(json.dumps({"schema_version": 1, "assignments": [query]}))
    return focal, catalogue, public, query_path


def args(paths, output, *, pool=3, maximum=3):
    focal, catalogue, public, _ = paths
    return [
        "prepare-context",
        str(focal),
        "--scheme",
        "ecoli",
        "--catalogue",
        str(catalogue),
        "--public-typing",
        str(public),
        "--output",
        str(output),
        "--candidate-pool",
        str(pool),
        "--max-context",
        str(maximum),
        "--nearest-per-focal",
        "1",
        "--seed",
        "7",
    ]


def tables(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def offline_services(monkeypatch):
    monkeypatch.setattr(provider, "load_api_key", lambda: "offline-test-only")

    def download(rows, *, output, **kwargs):
        output.mkdir(parents=True, exist_ok=True)
        paths, ledger = {}, []
        for row in rows:
            ident = row["source_genome_id"]
            path = output / f"{ident}.fasta"
            path.write_text(">chr\n" + "ACGT" * 100 + "\n")
            paths[ident] = path
            ledger.append({"source_genome_id": ident, "status": "downloaded"})
        return paths, ledger

    def screen(*, inputs, output_dir, **kwargs):
        ids = [line.split("\t")[0] for line in inputs.read_text().splitlines()]
        assert len(set(ids)) == len(ids)
        path = output_dir / "distances.tsv"
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["Sample1", "Sample2", "Distance", "Mismatches (proportion)"])
            for index, ident in enumerate(ids[1:]):
                assert ident != ids[0]
                writer.writerow([ids[0], ident, index + 1, 0.001])
        return path

    monkeypatch.setattr(provider, "download_assemblies", download)
    monkeypatch.setattr(provider, "run_ska_screen", screen)


def test_cli_imported_query_refines_same_st_pool_with_background(tmp_path, monkeypatch):
    paths = fixture(tmp_path)
    offline_services(monkeypatch)
    output = tmp_path / "context"
    result = CliRunner().invoke(app, args(paths, output) + ["--query-typing", str(paths[3])])
    assert result.exit_code == 0, result.output + repr(result.exception)
    audit = json.loads((output / "context_selection.json").read_text())
    assert audit["same_st_accessions"] == 7
    assert audit["query_typing_source"] == "imported"
    refinement = audit["lineage_refinement"]
    assert refinement["counts"]["priority_selected"] == 2
    assert refinement["counts"]["background_selected"] == 1
    assert refinement["settings"]["hiercc_level"] == "HC1100"
    pool = tables(output / "candidate_pool.tsv")
    assert all(row["pool_selection_reason"] for row in pool)
    assert all(row["pool_selection_reason"].startswith("lineage_priority:") for row in pool[:2])
    assert len(tables(output / "context_manifest.tsv")) == 3
    combined = read_metadata(output / "combined_metadata.csv")
    assert combined[0].sample_id == "F1" and len(combined) == 4
    frozen = json.loads((output / "context_catalogue.json").read_text())
    assert frozen["focal_rows"][0]["hiercc_codes"]["HC1100"] == "15"


def test_public_accession_lookup_preserves_query_and_excludes_self(tmp_path, monkeypatch):
    paths = fixture(tmp_path, lookup=True)
    offline_services(monkeypatch)
    output = tmp_path / "lookup"
    result = CliRunner().invoke(app, args(paths, output))
    assert result.exit_code == 0, result.output + repr(result.exception)
    audit = json.loads((output / "context_selection.json").read_text())
    assert audit["query_typing_source"] == "public_identity_lookup"
    assert audit["lineage_refinement"]["counts"]["eligible_same_st"] == 6
    frozen = json.loads((output / "context_catalogue.json").read_text())
    focal = frozen["focal_rows"][0]
    assert focal["sample_id"] == "SAMN90001"
    assert focal["public_typing_matched_source_genome_id"] == "self"
    assert focal["hiercc_codes"]["HC1100"] == "15"
    assert {row["source_genome_id"] for row in tables(output / "candidate_pool.tsv")[:2]} == {
        "match1",
        "match2",
    }
    assert all(row["source_genome_id"] != "self" for row in tables(output / "context_manifest.tsv"))
    assert read_metadata(output / "combined_metadata.csv")[0].sample_id == "SAMN90001"


def test_native_route_types_downloaded_pool_before_bounded_ska_screen(tmp_path, monkeypatch):
    paths = fixture(tmp_path)
    offline_services(monkeypatch)
    calls = []

    def native(samples, config_path, output, **kwargs):
        calls.append([sample.sample_id for sample in samples])
        result = []
        for sample in samples:
            ident = sample.sample_id.removeprefix("PW_")
            match = ident.startswith("match") or ident in {"F1", "self"}
            record = assignment(ident, "15" if match else "99", "1" if match else "2")
            record.pop("source_genome_id")
            record.update(
                sample_id=sample.sample_id,
                species=sample.species,
                assembly_sha256=hashlib.sha256(sample.assembly.read_bytes()).hexdigest(),
            )
            result.append(record)
        return result

    monkeypatch.setattr(typing, "type_query_assemblies", native)
    output = tmp_path / "native"
    config = tmp_path / "config.json"
    config.write_text("{}")
    result = CliRunner().invoke(
        app, args(paths, output, pool=7, maximum=3) + ["--typing-config", str(config)]
    )
    assert result.exit_code == 0, result.output + repr(result.exception)
    assert calls[0] == ["F1"] and len(calls[1]) == 7
    audit = json.loads((output / "context_selection.json").read_text())
    assert audit["query_typing_source"] == "native"
    assert (
        audit["native_public_typing_sha256"]
        == hashlib.sha256((output / "native_public_typing.json").read_bytes()).hexdigest()
    )
    assert audit["native_lineage_refinement"]["counts"]["selected"] == 3
    assert audit["native_lineage_refinement"]["counts"]["priority_selected"] == 2
    assert audit["stage_losses"]["native_lineage_refinement"] == 4
    assert len((output / "ska_inputs.tsv").read_text().splitlines()) == 4
    assert len(tables(output / "context_manifest.tsv")) == 3
    assert all(
        json.loads(row["hiercc_codes"])["HC1100"] for row in tables(output / "context_manifest.tsv")
    )
    native_public = json.loads((output / "native_public_typing.json").read_text())
    assert len(native_public["assignments"]) == 7
    assert {r["source_genome_id"] for r in native_public["assignments"]} == {
        "self",
        "match1",
        "match2",
        "wide1",
        "wide2",
        "wide3",
        "wide4",
    }


def test_imported_typing_assembly_mismatch_fails_before_download(tmp_path, monkeypatch):
    paths = fixture(tmp_path)
    payload = json.loads(paths[3].read_text())
    payload["assignments"][0]["assembly_sha256"] = "0" * 64
    paths[3].write_text(json.dumps(payload))
    monkeypatch.setattr(
        provider,
        "download_assemblies",
        lambda *a, **kw: (_ for _ in ()).throw(AssertionError("download must not run")),
    )
    output = tmp_path / "mismatch"
    result = CliRunner().invoke(app, args(paths, output) + ["--query-typing", str(paths[3])])
    assert result.exit_code != 0
    assert "assembly hash does not match" in result.output
    assert not (output / "context_manifest.tsv").exists()


def test_unavailable_hiercc_level_uses_audited_same_st_fallback(tmp_path):
    paths = fixture(tmp_path, with_profile=False)
    output = tmp_path / "unknown-level"
    result = CliRunner().invoke(
        app,
        args(paths, output)
        + ["--query-typing", str(paths[3]), "--hiercc-level", "HC999999", "--dry-run"],
    )
    assert result.exit_code == 0, result.output + repr(result.exception)
    audit = json.loads((output / "context_selection.json").read_text())
    refinement = audit["lineage_refinement"]
    assert refinement["counts"]["priority_eligible"] == 0
    assert refinement["counts"]["background_selected"] == 3
    assert refinement["query_status"]["F1"] == "no_comparable_evidence_same_st_fallback"
    assert refinement["counts"]["comparison_reasons"]["hiercc_unassigned_or_unversioned"] == 7
    assert all(
        row["pool_selection_reason"] == "balanced_same_st_fallback"
        for row in tables(output / "candidate_pool.tsv")
    )
