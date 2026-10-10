"""Frozen independent profile stages retain partitions, full trees and exact selections."""

import json

import pytest
from typer.testing import CliRunner

from chronoclade.cli import app
from chronoclade.cgmlst.partitions import partition_records, profile_records
from chronoclade.cgmlst.workflow import run_cgmlst
from chronoclade.datasets import LocusCatalogue, from_profile_records, write_dataset
from chronoclade.selection_manifest import load_selection_manifest


def fixture_dataset(tmp_path):
    records = []
    for ident, role, code, allele in [
        ("q", "local", "0,0,1,0,1,0,0,0,0,0", "1"),
        ("c1", "context", "0,0,1,0,1,0,0,0,0,0", "2"),
        ("c2", "context", "0,0,1,0,1,1,0,0,0,0", "3"),
        ("outside", "context", "0,0,1,0,2,0,0,0,0,0", "4"),
        ("unresolved", "local", None, "5"),
    ]:
        records.append(
            {
                "sample_id": ident,
                "origin": role,
                "species": "Klebsiella pneumoniae",
                "mlst_st": "39",
                "collection_date": "2020",
                "country": "Greece",
                "cgmlst_scheme": "test",
                "cgmlst_scheme_version": "1",
                "cgmlst_profile": {"a": allele, "b": "1"},
                "cgmlst_status": "resolved",
                "cglin_scheme": "scgMLST629_S",
                "cglin_scheme_version": "1",
                "cglin_raw": code,
                "cglin_status": "complete" if code else "missing",
            }
        )
    dataset = from_profile_records(
        records, [LocusCatalogue("test", "1", ("a", "b"), "frozen fixture")]
    )
    return dataset, write_dataset(dataset, tmp_path / "prepared")


def test_partitions_explicit_depth_full_membership_and_unresolved_audit(tmp_path):
    dataset, _ = fixture_dataset(tmp_path)
    records = profile_records(dataset)
    blocks, audit = partition_records(records)
    assert [row["sample_id"] for row in blocks[0]["records"]] == ["c1", "c2", "q"]
    assert audit["unmatched_context_ids"] == ["outside"]
    assert audit["unresolved_inputs"][0]["sample_id"] == "unresolved"
    assert {row["level"] for row in audit["available_level_counts"]} == {5, 6, 7}
    narrowed, _ = partition_records(records, lin_level=6)
    assert [row["sample_id"] for row in narrowed[0]["records"]] == ["c1", "q"]
    records[1]["cglin_database_sha256"] = "a" * 64
    scoped, _ = partition_records(records)
    assert "c1" not in {row["sample_id"] for row in scoped[0]["records"]}


def test_real_stage_has_full_tree_and_validated_selections_without_assembly_retrieval(
    tmp_path, monkeypatch
):
    _, manifest = fixture_dataset(tmp_path)

    def forbidden(*args, **kwargs):
        raise AssertionError("Independent frozen profile analysis must not acquire assemblies")

    monkeypatch.setattr("chronoclade.profile_inputs.materialise_assemblies", forbidden)
    result = run_cgmlst(
        manifest,
        tmp_path / "cgmlst",
        replicates=2,
        context_size=1,
        nearest_per_query=1,
        bootstrap_replicates=0,
    )
    record = json.loads(result.manifest_path.read_text())
    assert record["status"] == "partial"
    block = record["blocks"][0]
    assert block["input_count"] == 1 and block["context_count"] == 2
    directory = result.manifest_path.parent
    ensemble = json.loads((directory / block["ensemble"]).read_text())
    # Inspect the emitted exact saved selection via the public loader.
    for path in (directory / block["ensemble"]).parent.glob("selections/selection-*.json"):
        selection = load_selection_manifest(path)
        assert "q" in selection["selected_sample_ids"]
    analysis = json.loads((directory / block["block_id"] / "profile_analysis.json").read_text())
    assert analysis["records_count"] == 3
    assert ensemble
    newick_files = list((directory / block["block_id"]).glob("*.nwk"))
    if not newick_files:
        newick_files = list((directory / block["block_id"]).glob("*.newick"))
    assert newick_files
    tree = newick_files[0].read_text()
    assert all(ident in tree for ident in ("q", "c1", "c2"))
    assert result.report_path.is_file()
    index_html = result.report_path.read_text()
    assert "Available LIN depths" in index_html and "ST39" in index_html
    block_html = (directory / block["report"]).read_text()
    assert "cgMLST profile report" in block_html
    with pytest.raises(ValueError, match="empty"):
        run_cgmlst(manifest, tmp_path / "cgmlst")


def test_cgmlst_cli_errors_are_clear_and_do_not_publish_completion(tmp_path):
    result = CliRunner().invoke(
        app, ["cgmlst", str(tmp_path / "missing.json"), "--out", str(tmp_path / "out")]
    )
    assert result.exit_code == 2
    assert not (tmp_path / "out/cgmlst.json").exists()


def test_lineage_evidence_cannot_change_validated_identity_or_paths(tmp_path):
    from dataclasses import replace

    dataset, _ = fixture_dataset(tmp_path)
    lineage = dict(dataset.lineages[0])
    lineage["evidence"] = dict(
        lineage["evidence"], sample_id="escape", role="context", cgmlst_profile={"a": "999"}
    )
    changed = replace(dataset, lineages=(lineage, *dataset.lineages[1:]))
    rows = profile_records(changed)
    assert rows[0]["sample_id"] == "q" and rows[0]["role"] == "input"
    assert rows[0]["cgmlst_profile"]["a"] == "1"
    rows[0]["mlst_st"] = "39/../../escape"
    blocks, _ = partition_records(rows)
    assert all("/" not in block["block_id"] and ".." not in block["block_id"] for block in blocks)
    duplicated = replace(dataset, lineages=(*dataset.lineages, dataset.lineages[0]))
    with pytest.raises(ValueError, match="Repeated lineage"):
        profile_records(duplicated)


def test_cli_stage_help_is_lazy():
    import subprocess
    import sys

    script = (
        "import sys; from chronoclade.cli import app; from typer.testing import CliRunner; "
        "assert CliRunner().invoke(app,['--help']).exit_code==0; "
        "assert 'chronoclade.workflow' not in sys.modules; "
        "assert 'chronoclade.lineage' not in sys.modules; "
        "assert 'torch' not in sys.modules; assert 'esm' not in sys.modules"
    )
    result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize(
    "update,message", [({"kind": None}, "kind"), ({"evidence": ["bad"]}, "mapping")]
)
def test_malformed_generic_lineage_evidence_has_clear_error(tmp_path, update, message):
    from dataclasses import replace

    dataset, _ = fixture_dataset(tmp_path)
    changed = dict(dataset.lineages[0], **update)
    with pytest.raises(ValueError, match=message):
        profile_records(replace(dataset, lineages=(changed, *dataset.lineages[1:])))


def test_unknown_lineage_namespace_is_not_a_shared_partition(tmp_path):
    dataset, _ = fixture_dataset(tmp_path)
    rows = profile_records(dataset)
    for row in rows:
        row["cglin_scheme_version"] = "unknown"
    blocks, audit = partition_records(rows)
    assert blocks == []
    assert len(audit["unresolved_inputs"]) == 2


def test_distinct_mlst_namespaces_do_not_share_a_st_block(tmp_path):
    dataset, _ = fixture_dataset(tmp_path)
    rows = profile_records(dataset)
    rows[0]["mlst_scheme"] = "primary-scheme"
    rows[1]["mlst_scheme"] = "another-scheme"
    blocks, audit = partition_records(rows)
    assert "c1" in audit["unmatched_context_ids"]
    assert all(row["sample_id"] != "c1" for row in blocks[0]["records"])


def test_missing_mlst_metadata_can_compare_with_one_explicit_namespace(tmp_path):
    dataset, _ = fixture_dataset(tmp_path)
    rows = profile_records(dataset)
    rows[1]["mlst_scheme"] = "primary-scheme"
    blocks, audit = partition_records(rows)
    assert {r["sample_id"] for r in blocks[0]["records"]} == {"q", "c1", "c2"}
    assert next(r for r in blocks[0]["records"] if r["sample_id"] == "q")["mlst_scheme"] is None
    assert audit["ambiguous_mlst_namespace_ids"] == []
