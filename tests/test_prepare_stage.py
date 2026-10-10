"""Independent prepare preserves portable evidence without running analysis/providers."""

import json
import shutil
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

import pytest
import typer
from typer.testing import CliRunner

from chronoclade.datasets import DatasetError, LocusCatalogue, load_dataset
from chronoclade.prepare_cli import prepare_command
from chronoclade.prepare_stage import run_prepare


@pytest.fixture
def frozen(tmp_path):
    records = [
        {
            "sample_id": "query",
            "origin": "local",
            "species": "Klebsiella pneumoniae",
            "collection_date": "2024-02",
            "country": "Greece",
            "cgmlst_scheme": "example:cgmlst",
            "cgmlst_scheme_version": "v1",
            "cgmlst_status": "resolved",
            "cgmlst_profile": {"a": "001"},
            "cgmlst_novel_alleles": {"b": "A" * 40},
            "cglin_scheme": "example:cglin",
            "cglin_scheme_version": "lin-v1",
            "cglin_status": "resolved",
            "cglin_raw": "1.2.3.4.5.6.7",
            "identity_conflict": True,
            "identity_candidates": ["SAMN1", "SAMN2"],
        },
        {
            "sample_id": "missing",
            "role": "input",
            "species": "Klebsiella pneumoniae",
            "cgmlst_scheme": "example:cgmlst",
            "cgmlst_scheme_version": "v1",
            "cgmlst_status": "unassigned",
            "collection_date": "invalid-date",
        },
    ]
    catalogue = LocusCatalogue(
        "example:cgmlst", "v1", ("a", "b", "never_called"), "frozen-schema.json"
    )
    input_path, catalogues = tmp_path / "profiles.json", tmp_path / "catalogues.json"
    input_path.write_text(json.dumps(records))
    catalogues.write_text(json.dumps([asdict(catalogue)]))
    return input_path, catalogues, records


def test_profiles_publish_full_catalogue_missing_samples_and_lineage(frozen, tmp_path):
    source, catalogues, _ = frozen
    result = run_prepare(source, tmp_path / "prepared", catalogues=catalogues)
    dataset = load_dataset(result.dataset_manifest)
    assert dataset.sample_ids == ("query", "missing")
    assert dataset.profiles[0].catalogue.loci == ("a", "b", "never_called")
    assert dataset.profiles[0].profile("query") == {
        "a": "001",
        "b": "A" * 40,
        "never_called": None,
    }
    assert dataset.profiles[0].profile("missing") == {
        "a": None,
        "b": None,
        "never_called": None,
    }
    assert dataset.samples[0]["date_end"] == "2024-02-29"
    assert dataset.samples[1]["date_precision"] == "invalid"
    assert dataset.lineages[0]["assignment"] == "1.2.3.4.5.6.7"
    assert dataset.crosswalk[0]["resolution_status"] == "ambiguous"
    audit = json.loads(result.audit_path.read_text())
    assert audit["readiness"] == "incomplete"
    assert audit["samples"][0]["readiness"] == "ready"
    assert audit["samples"][1]["profile_ready"] is False
    assert audit["operations"]["reference_database"] == "not_checked_not_required_for_import"
    assert "Archivo" in result.report_path.read_text()
    assert "No profile" not in result.report_path.read_text()
    assert "no profile calling was performed" in result.report_path.read_text()
    stage = json.loads(result.stage_manifest.read_text())
    assert stage["status"] == "complete"
    assert stage["sample_ids"] == ["query", "missing"]
    assert "sources/input.json" in {item["path"] for item in stage["artifacts"]}


def test_explicit_metadata_provenance_conflicts_and_retrieval_preserved(frozen, tmp_path):
    source, catalogues, records = frozen
    evidence = {
        "records": records,
        "provenance": [
            {
                "sample_id": "query",
                "field": "country",
                "source": "user",
                "value": "Greece",
                "selected": True,
            }
        ],
        "conflicts": [
            {"sample_id": "query", "field": "country", "source": "provider", "value": "UK"}
        ],
        "retrieval": [
            {
                "sample_id": "missing",
                "operation": "lineage_assignment",
                "status": "unavailable",
                "reason": "reference_database_unavailable",
            }
        ],
        "parameters": {"frozen_export": "experiment-v1"},
    }
    source.write_text(json.dumps(evidence))
    result = run_prepare(source, tmp_path / "prepared", catalogues=catalogues)
    dataset = load_dataset(result.dataset_manifest)
    for name in ("provenance", "conflicts", "retrieval"):
        assert list(getattr(dataset, name)) == evidence[name]
    assert dataset.parameters["frozen_export"] == "experiment-v1"


def test_local_assembly_is_copied_and_remains_portable(frozen, tmp_path):
    source, catalogues, records = frozen
    assembly = tmp_path / "query.fasta"
    assembly.write_text(">query\nACGT\n")
    records[0]["assembly"] = "query.fasta"
    source.write_text(json.dumps(records))
    result = run_prepare(source, tmp_path / "prepared", catalogues=catalogues)
    dataset = load_dataset(result.dataset_manifest)
    reference = dataset.samples[0]["assembly_reference"]
    assert reference.startswith("assemblies/")
    assert not Path(reference).is_absolute()
    assert (result.dataset_manifest.parent / reference).read_bytes() == assembly.read_bytes()
    assert dataset.parameters["prepare"]["assemblies"][0]["path"] == reference
    moved = tmp_path / "moved"
    shutil.move(result.dataset_manifest.parent, moved)
    assembly.unlink()
    relocated = load_dataset(moved / "dataset.json")
    assert (moved / relocated.samples[0]["assembly_reference"]).read_text() == ">query\nACGT\n"
    copied = run_prepare(moved / "dataset.json", tmp_path / "imported")
    assert (copied.dataset_manifest.parent / reference).read_text() == ">query\nACGT\n"
    imported = load_dataset(copied.dataset_manifest)
    assert imported.parameters["prepare"]["source_dataset_id"] == relocated.dataset_id


def test_existing_bundle_import_validates_and_preserves_evidence(frozen, tmp_path):
    source, catalogues, _ = frozen
    original = run_prepare(source, tmp_path / "first", catalogues=catalogues)
    copied = run_prepare(original.dataset_manifest.parent, tmp_path / "second")
    first, second = load_dataset(original.dataset_manifest), load_dataset(copied.dataset_manifest)
    assert first.samples == second.samples
    assert first.lineages == second.lineages
    assert first.crosswalk == second.crosswalk
    (original.dataset_manifest.parent / "samples.csv").write_text("tampered")
    with pytest.raises(DatasetError, match="checksum/size"):
        run_prepare(original.dataset_manifest, tmp_path / "corrupt-import")
    assert not (tmp_path / "corrupt-import").exists()


def test_complete_lineage_is_a_separate_check(frozen, tmp_path):
    source, catalogues, records = frozen
    records[0]["cglin_status"] = "unassigned"
    records[0]["cglin_raw"] = None
    source.write_text(json.dumps(records))
    result = run_prepare(source, tmp_path / "prepared", catalogues=catalogues)
    audit = json.loads(result.audit_path.read_text())
    assert audit["samples"][0]["profile_ready"] is True
    assert audit["samples"][0]["lineage_ready"] is False
    assert audit["samples"][0]["readiness"] == "incomplete"


@pytest.mark.parametrize(
    "status,provisional,ready",
    [
        ("complete", False, True),
        ("partial", False, False),
        ("provisional", True, False),
        ("complete", True, False),
    ],
)
def test_authentic_cglin_status_readiness(frozen, tmp_path, status, provisional, ready):
    source, catalogues, records = frozen
    records[0].update(
        cglin_scheme="scgMLST629_S",
        cglin_raw="0,0,107,0,0,0,0,0,5,0",
        cglin_status=status,
        cglin_provisional=provisional,
    )
    source.write_text(json.dumps(records))
    result = run_prepare(source, tmp_path / "prepared", catalogues=catalogues)
    audit = json.loads(result.audit_path.read_text())
    assert audit["samples"][0]["lineage_ready"] is ready
    assert load_dataset(result.dataset_manifest).lineages[0]["resolution_status"] == status


def test_malformed_lineage_evidence_fails_with_field_specific_error(frozen, tmp_path):
    source, catalogues, records = frozen
    source.write_text(
        json.dumps(
            {
                "records": records,
                "lineages": [
                    {
                        "sample_id": "query",
                        "kind": "cglin",
                        "scheme_id": "test",
                        "scheme_version": "v1",
                        "database_version": None,
                        "assignment": "0,0,107,0,0,0,0,0,5,0",
                        "assignment_method": "frozen",
                        "resolution_status": "complete",
                        "evidence": ["invalid-shape"],
                    }
                ],
            }
        )
    )
    with pytest.raises(DatasetError, match="lineage evidence must be an object"):
        run_prepare(source, tmp_path / "prepared", catalogues=catalogues)
    assert not (tmp_path / "prepared").exists()


def test_unknown_lineage_version_is_not_ready(frozen, tmp_path):
    source, catalogues, records = frozen
    records[0]["cglin_scheme_version"] = "unknown"
    source.write_text(json.dumps(records))
    result = run_prepare(source, tmp_path / "prepared", catalogues=catalogues)
    audit = json.loads(result.audit_path.read_text())
    assert audit["samples"][0]["lineage_ready"] is False


def test_hiercc_assignment_and_supplied_context_are_preserved(frozen, tmp_path):
    source, catalogues, records = frozen
    records[0] = {key: value for key, value in records[0].items() if not key.startswith("cglin_")}
    records[0].update(
        {
            "species": "Escherichia coli",
            "origin": "context",
            "hiercc_scheme": "example:hiercc",
            "hiercc_scheme_version": "v1",
            "hiercc_codes": {"HC0": "1", "HC10": "23", "HC1100": "55"},
            "hiercc_status": "resolved",
        }
    )
    source.write_text(json.dumps(records))
    result = run_prepare(source, tmp_path / "prepared", catalogues=catalogues)
    dataset = load_dataset(result.dataset_manifest)
    assert dataset.samples[0]["role"] == "context"
    assert dataset.lineages[0]["assignment"] == records[0]["hiercc_codes"]
    audit = json.loads(result.audit_path.read_text())
    assert audit["context_ids"] == ["query"]
    assert audit["samples"][0]["required_lineage_kind"] == "hiercc"
    assert audit["samples"][0]["lineage_ready"] is True


def test_bundled_assembly_tampering_fails_import(frozen, tmp_path):
    source, catalogues, records = frozen
    (tmp_path / "query.fasta").write_text(">query\nACGT\n")
    records[0]["assembly"] = "query.fasta"
    source.write_text(json.dumps(records))
    result = run_prepare(source, tmp_path / "prepared", catalogues=catalogues)
    reference = load_dataset(result.dataset_manifest).samples[0]["assembly_reference"]
    (result.dataset_manifest.parent / reference).write_text("tampered")
    with pytest.raises(DatasetError, match="checksum/size"):
        run_prepare(result.dataset_manifest, tmp_path / "tampered-import")
    assert not (tmp_path / "tampered-import").exists()


def test_unknown_species_does_not_claim_lineage_ready(frozen, tmp_path):
    source, catalogues, records = frozen
    records[0]["species"] = "Example species"
    source.write_text(json.dumps(records))
    result = run_prepare(source, tmp_path / "prepared", catalogues=catalogues)
    audit = json.loads(result.audit_path.read_text())
    assert audit["samples"][0]["required_lineage_kind"] is None
    assert audit["samples"][0]["lineage_ready"] is False
    assert "not assessed" in audit["samples"][0]["reasons"][0]


def test_incomplete_catalogue_does_not_claim_ready(frozen, tmp_path):
    source, catalogues, _ = frozen
    values = json.loads(catalogues.read_text())
    values[0]["complete"] = False
    catalogues.write_text(json.dumps({"catalogues": values}))
    result = run_prepare(source, tmp_path / "prepared", catalogues=catalogues)
    audit = json.loads(result.audit_path.read_text())
    assert audit["samples"][0]["readiness"] == "incomplete"
    assert "explicitly incomplete" in audit["samples"][0]["reasons"][0]


def test_conflicting_known_database_fingerprint_fails(frozen, tmp_path):
    source, catalogues, records = frozen
    records[0]["cgmlst_database_sha256"] = "a" * 64
    source.write_text(json.dumps(records))
    with pytest.raises(DatasetError, match="fingerprint"):
        run_prepare(source, tmp_path / "prepared", catalogues=catalogues)
    assert not (tmp_path / "prepared").exists()


def test_catalogue_completeness_cannot_be_inferred(frozen, tmp_path):
    source, catalogues, _ = frozen
    values = json.loads(catalogues.read_text())
    del values[0]["complete"]
    catalogues.write_text(json.dumps(values))
    with pytest.raises(DatasetError, match="explicit complete"):
        run_prepare(source, tmp_path / "prepared", catalogues=catalogues)


def test_credentials_fail_before_snapshot_publication(frozen, tmp_path):
    source, catalogues, records = frozen
    records[0]["api_key"] = "not-a-real-key"
    source.write_text(json.dumps(records))
    with pytest.raises(DatasetError, match="credentials"):
        run_prepare(source, tmp_path / "prepared", catalogues=catalogues)
    assert not (tmp_path / "prepared").exists()


def test_output_is_immutable(frozen, tmp_path):
    source, catalogues, _ = frozen
    result = run_prepare(source, tmp_path / "prepared", catalogues=catalogues)
    original = result.dataset_manifest.read_bytes()
    with pytest.raises(DatasetError, match="already exists"):
        run_prepare(source, tmp_path / "prepared", catalogues=catalogues)
    assert result.dataset_manifest.read_bytes() == original


def test_failure_during_report_leaves_no_partial_output(frozen, tmp_path, monkeypatch):
    from chronoclade import prepare_stage

    source, catalogues, _ = frozen

    def fail(audit):
        raise OSError("report publication interrupted")

    monkeypatch.setattr(prepare_stage, "_report", fail)
    with pytest.raises(OSError, match="interrupted"):
        run_prepare(source, tmp_path / "prepared", catalogues=catalogues)
    assert not (tmp_path / "prepared").exists()
    assert not list(tmp_path.glob(".prepared.prepare-*"))


def test_species_fill_and_conflicts(frozen, tmp_path):
    source, catalogues, records = frozen
    records[0].pop("species")
    source.write_text(json.dumps(records))
    result = run_prepare(
        source, tmp_path / "filled", catalogues=catalogues, species="Klebsiella pneumoniae"
    )
    assert load_dataset(result.dataset_manifest).samples[0]["species"] == "Klebsiella pneumoniae"
    with pytest.raises(DatasetError, match="differs"):
        run_prepare(source, tmp_path / "conflict", catalogues=catalogues, species="E. coli")


def test_report_escapes_input_labels(frozen, tmp_path):
    source, catalogues, records = frozen
    records[0]["label"] = "<script>bad()</script>"
    source.write_text(json.dumps(records))
    result = run_prepare(source, tmp_path / "prepared", catalogues=catalogues)
    assert "&lt;script&gt;bad()&lt;/script&gt;" in result.report_path.read_text()
    assert "<script>bad()" not in result.report_path.read_text()


def test_cli_success_and_actionable_missing_catalogue(frozen, tmp_path):
    source, catalogues, _ = frozen
    app = typer.Typer()
    app.command("prepare")(prepare_command)
    runner = CliRunner()
    success = runner.invoke(
        app, [str(source), "--catalogues", str(catalogues), "--out", str(tmp_path / "prepared")]
    )
    assert success.exit_code == 0, success.output
    assert "Dataset manifest:" in success.output
    failure = runner.invoke(app, [str(source), "--out", str(tmp_path / "missing")])
    assert failure.exit_code == 2
    assert "explicit locus catalogue" in failure.output


def test_prepare_import_has_no_scientific_stage_or_provider_dependency():
    check = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import sys; import chronoclade.prepare_stage; "
                "forbidden = {'torch', 'chronoclade.profile_inputs', 'chronoclade.workflow', "
                "'chronoclade.profile_analysis', 'chronoclade.pathogenwatch'}; "
                "assert not forbidden.intersection(sys.modules), forbidden.intersection(sys.modules)"
            ),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert check.returncode == 0, check.stderr


@pytest.mark.parametrize("input_kind", ["collection", "accessions", "assemblies"])
def test_live_modes_reject_wrong_input_contract(frozen, tmp_path, input_kind):
    source, _, _ = frozen
    with pytest.raises(DatasetError, match="short UUID|Metadata has no samples|Assembly metadata"):
        run_prepare(source, tmp_path / "prepared", input_kind=input_kind)
