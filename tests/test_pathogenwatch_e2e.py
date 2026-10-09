"""Offline provider-to-manifest/report integration; no credentials or native tools."""

import csv
import importlib.util
import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from chronoclade.cli import app
from chronoclade.metadata import read_metadata
from chronoclade.pathogenwatch import PathogenwatchClient
import chronoclade.pathogenwatch_context as provider


def make_fixture(tmp_path):
    records = {
        "a": ("SAMN10001", "India", "2019", True),
        "a_duplicate": ("SAMN10001", "India", "2019", True),
        "b": ("SAMN10002", "United Kingdom", "2020-02", True),
        "c": ("SAMN10003", "", "2021-03-04", True),
        "undated": ("SAMN10004", "France", "", True),
        "qc_failure": ("SAMN10005", "Germany", "2022", False),
        "singleton": ("SAMN10006", "France", "2023", True),
    }

    def transport(method, path, *, body, params, headers):
        assert "X-API-Key" not in headers
        if path.endswith("/supported"):
            return [
                {
                    "organismId": "573",
                    "fullName": "Klebsiella pneumoniae",
                    "typing": ["MLST"],
                    "other": ["LIN Codes"],
                }
            ]
        if path.endswith("/search/genomes"):
            return {
                "meta": {"count": len(records), "endCursor": "end"},
                "genomes": [
                    {
                        "uuid": identifier,
                        "projectAccess": "PUBLIC",
                        "name": sample,
                        "organism": "Klebsiella pneumoniae",
                        "location": None,
                    }
                    for identifier, (sample, *_rest) in records.items()
                ],
            }
        identifier = params["id"]
        sample, country, day, qc = records[identifier]
        return {
            "uuid": identifier,
            "id": list(records).index(identifier) + 1,
            "organismId": "573",
            "mlst": "147",
            "name": sample,
            "qc": qc,
            "sampleAccession": sample,
            "metadata": [["Country", country], ["Date", day]],
        }

    catalogue = tmp_path / "catalogue.json"
    PathogenwatchClient(transport=transport).freeze_catalogue(catalogue)
    export = tmp_path / "cglin.csv"
    with export.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Genome ID", "LIN code", "cgST"])
        for identifier in records:
            code = "0,0,197,0,4,0,93,0,0,0"
            if identifier == "singleton":
                code = "1,1,198,0,4,0,93,0,0,0"
            if identifier == "undated":
                code = "0,0,197,0,4,-,-,-,-,-"
            writer.writerow([identifier, code, "19550"])
    assembly = tmp_path / "focal.fasta"
    assembly.write_text(">chr\n" + "ACGT" * 100 + "\n")
    focal = tmp_path / "focal.csv"
    with focal.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "sample_id",
                "assembly",
                "collection_date",
                "location",
                "species",
                "lineage",
                "origin",
                "is_reference",
            ]
        )
        writer.writerow(
            [
                "F1",
                str(assembly),
                "2018",
                "India",
                "Klebsiella pneumoniae",
                "ST147",
                "local",
                "true",
            ]
        )
    return focal, catalogue, export


def read_table(path, delimiter="\t"):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def mock_download_and_screen(monkeypatch):
    monkeypatch.setattr(provider, "load_api_key", lambda: "offline-fixture-only")

    def download(rows, *, output, **kwargs):
        output.mkdir(parents=True, exist_ok=True)
        paths, ledger = {}, []
        for row in rows:
            identifier = row["source_genome_id"]
            if identifier == "c":
                ledger.append(
                    {
                        "source_genome_id": identifier,
                        "status": "failed",
                        "error": "deliberate fixture partial batch failure",
                    }
                )
                continue
            path = output / f"{identifier}.fasta"
            path.write_text(">chr\n" + "ACGT" * 100 + "\n")
            paths[identifier] = path
            ledger.append(
                {
                    "source_genome_id": identifier,
                    "status": "downloaded",
                    "bytes": path.stat().st_size,
                    "request_count": 1,
                }
            )
        return paths, ledger

    def screen(*, inputs, output_dir, **kwargs):
        rows = [line.split("\t")[0] for line in inputs.read_text().splitlines()]
        path = output_dir / "context_distances.tsv"
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["Sample1", "Sample2", "Distance", "Mismatches (proportion)"])
            for index, identifier in enumerate(rows[1:]):
                writer.writerow([rows[0], identifier, index + 1, ".0001"])
        return path

    monkeypatch.setattr(provider, "download_assemblies", download)
    monkeypatch.setattr(provider, "run_ska_screen", screen)


def prepare_args(focal, catalogue, export, output, *, dry=False):
    command = [
        "prepare-context",
        str(focal),
        "--scheme",
        "klebsiella",
        "--catalogue",
        str(catalogue),
        "--cglin-export",
        str(export),
        "--output",
        str(output),
        "--cache-dir",
        str(output / "cache"),
        "--candidate-pool",
        "8",
        "--max-context",
        "2",
        "--nearest-per-focal",
        "1",
        "--seed",
        "20261009",
    ]
    return command + (["--dry-run"] if dry else [])


def test_cli_offline_discovery_download_selection_manifest_geography_replay(tmp_path, monkeypatch):
    focal, catalogue, export = make_fixture(tmp_path)
    mock_download_and_screen(monkeypatch)
    runner = CliRunner()
    output = tmp_path / "context"
    result = runner.invoke(app, prepare_args(focal, catalogue, export, output))
    assert result.exit_code == 0, result.output + repr(result.exception)
    audit = json.loads((output / "context_selection.json").read_text())
    assert audit["context_source"] == "pathogenwatch"
    assert audit["same_st_accessions"] == 7
    assert audit["candidate_pool"] == 4
    assert audit["selected_contexts"] == 2
    assert audit["stage_losses"]["undated_or_invalid"] == 1
    assert audit["stage_losses"]["download_failure"] == 1
    manifest = read_table(output / "context_manifest.tsv")
    assert len(manifest) == 2
    assert all(
        r["selection_reason"] and r["source_genome_id"] and r["assembly_sha256"] for r in manifest
    )
    assert len(read_metadata(output / "combined_metadata.csv")) == 3
    geography = json.loads(
        (output / "context_geography/country_composition_audit.json").read_text()
    )
    assert geography["sample_units"] == 5
    assert geography["qc_eligible_records"] == 6
    assert any(
        s["unresolved_units"] == 1
        for s in geography["summaries"]
        if s["cohort"] == "public_catalogue" and s["depth"] == 7
    )
    assert "United Kingdom" in (output / "context_geography/index.html").read_text()
    script = Path(__file__).resolve().parents[1] / "scripts/st147_pilot.py"
    spec = importlib.util.spec_from_file_location("st147_pilot", script)
    pilot = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pilot)
    replay = pilot.replay(output, tmp_path / "replayed")
    assert replay["selection_identical"] and replay["country_counts_identical"]
    result2 = runner.invoke(app, prepare_args(focal, catalogue, export, output))
    assert result2.exit_code == 0, result2.output
    assert [(r["source_genome_id"], r["selection_reason"]) for r in manifest] == [
        (r["source_genome_id"], r["selection_reason"])
        for r in read_table(output / "context_manifest.tsv")
    ]
    assert "offline-fixture-only" not in (output / "context_selection.json").read_text()


def test_dry_run_without_download_credentials_and_tampered_snapshot(tmp_path, monkeypatch):
    focal, catalogue, export = make_fixture(tmp_path)
    monkeypatch.setattr(
        provider, "load_api_key", lambda: pytest.fail("dry-run requested credentials")
    )
    output = tmp_path / "dry"
    runner = CliRunner()
    result = runner.invoke(app, prepare_args(focal, catalogue, export, output, dry=True))
    assert result.exit_code == 0, result.output
    assert not (output / "assemblies").exists()
    assert (output / "context_geography/index.html").exists()
    payload = json.loads(catalogue.read_text())
    payload["rows"][0]["country"] = "France"
    catalogue.write_text(json.dumps(payload))
    failed = runner.invoke(app, prepare_args(focal, catalogue, export, output, dry=True))
    assert failed.exit_code != 0 and "hash mismatch" in failed.output


def test_unsupported_species_is_rejected(tmp_path):
    focal, catalogue, export = make_fixture(tmp_path)
    focal.write_text(focal.read_text().replace("Klebsiella pneumoniae", "Escherichia coli"))
    result = CliRunner().invoke(
        app, prepare_args(focal, catalogue, export, tmp_path / "unsupported", dry=True)
    )
    assert result.exit_code != 0
    assert "currently supports Klebsiella pneumoniae" in result.output
    assert "ATB" not in result.output and "--context-source" not in result.output


@pytest.mark.parametrize(
    "option,value",
    [("--source", "aws"), ("--context-source", "atb"), ("--metadata-table", "legacy.parquet")],
)
def test_removed_provider_options_are_rejected(tmp_path, option, value):
    focal, catalogue, export = make_fixture(tmp_path)
    output = tmp_path / "removed-option"
    result = CliRunner().invoke(
        app, prepare_args(focal, catalogue, export, output, dry=True) + [option, value]
    )
    assert result.exit_code != 0
    assert "No such option" in result.output
    assert not output.exists()


def test_context_help_exposes_only_pathogenwatch_route():
    result = CliRunner().invoke(app, ["prepare-context", "--help"])
    assert result.exit_code == 0, result.output
    assert "--catalogue" in result.output
    for removed in ("--source", "--context-source", "--metadata-table", "atbfetcher"):
        assert removed not in result.output


def test_interval_date_candidate_retains_bounds_and_filters_explicitly(tmp_path):
    from chronoclade.context import filter_candidates, stratified_candidate_pool

    normalized = {
        "source_genome_id": "interval",
        "country": "India",
        "biosample": "SAMN1",
        "dated_cohort_eligible": True,
        "date_precision": "interval",
        "date_start": "2018-01-01",
        "date_end": "2019-12-31",
        "collection_date": "2018/2019",
    }
    candidate = provider._candidate(
        normalized,
        species="Klebsiella pneumoniae",
        lineage="ST147",
        scheme="klebsiella",
        st="147",
        snapshot="fixture",
        catalogue_path=tmp_path / "frozen.json",
    )
    assert candidate.collection_date.startswith("[2018.")
    assert candidate.date_start == "2018-01-01" and candidate.date_end == "2019-12-31"
    assert candidate.year == "2018"
    assert filter_candidates([candidate], year_from=2018, year_to=2019) == [candidate]
    assert filter_candidates([candidate], year_to=2017) == []
    assert stratified_candidate_pool([candidate], limit=8, seed=7) == [candidate]
