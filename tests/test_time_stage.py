"""Dating reads saved evidence, preserves products and reports unavailable runs."""

import json

import pytest

from chronoclade.errors import WorkflowError
from chronoclade.time_stage.workflow import run_time
from chronoclade.time_stage.comparison import compare_time_runs
from chronoclade.tree_stage.workflow import run_tree
from test_tree_stage import selection_fixture, native_stubs as _tree_native_stubs


@pytest.fixture
def native_stubs(monkeypatch):
    return _tree_native_stubs.__wrapped__(monkeypatch)


def test_unsupported_time_needs_no_native_tools_or_tree_rebuild(
    tmp_path, native_stubs, monkeypatch
):
    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    tree = run_tree(selection, output=tmp_path / "tree")
    monkeypatch.setattr(
        "chronoclade.time_stage.workflow.require_tools",
        lambda names: (_ for _ in ()).throw(AssertionError("no tools for unsupported evidence")),
    )
    result = run_time(tmp_path / "tree/tree.json", output=tmp_path / "time")
    assert result["dating_status"] == "unsupported" and result["dated_tree"] is None
    assert run_time(tmp_path / "tree/tree.json", output=tmp_path / "time") == result
    assert native_stubs == ["corrected"]
    assert json.loads((tmp_path / "tree/tree.json").read_text()) == tree
    with pytest.raises(WorkflowError, match="three saved eligible"):
        run_time(tmp_path / "tree/tree.json", output=tmp_path / "override", allow_unsupported=True)
    assert not (tmp_path / "override/time.json").exists()


def test_ensemble_comparison_keeps_unsupported_and_absent_runs(tmp_path, native_stubs):
    selection, ensemble = selection_fixture(tmp_path, dates=("", "", ""))
    run_tree(selection, output=tmp_path / "tree")
    run_time(tmp_path / "tree/tree.json", output=tmp_path / "time")
    comparison = compare_time_runs(
        ensemble, [tmp_path / "time/time.json"], output=tmp_path / "comparison"
    )
    assert comparison["requested_runs"] == 2 and comparison["dated_runs"] == 0
    assert [run["status"] for run in comparison["runs"]] == ["unsupported", "unavailable"]
    assert comparison["sensitivity"]["rate"]["assessable_runs"] == 0
    assert "not confidence intervals" in comparison["interpretation"]
    assert comparison["target_pair_agreement"] == []


def test_time_rejects_fabricated_gate_in_rehashed_manifest(tmp_path, native_stubs):
    from chronoclade.stage_artifacts import load_result

    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    run_tree(selection, output=tmp_path / "tree")
    run_time(tmp_path / "tree/tree.json", output=tmp_path / "time")
    manifest = tmp_path / "time/time.json"
    edited = json.loads(manifest.read_text())
    edited["temporal_assessment"]["supported"] = True
    manifest.write_text(json.dumps(edited))
    with pytest.raises(WorkflowError, match="differs from the saved tree"):
        load_result(manifest, "time")


@pytest.mark.integration
def test_native_treetime_dating_saved_tree_without_phylogeny_rebuild(
    tmp_path, native_stubs, monkeypatch
):
    import shutil
    from chronoclade.stage_artifacts import load_result

    if not shutil.which("treetime"):
        pytest.skip("Native TreeTime unavailable")
    import chronoclade.lineage as lineage

    original = lineage._run_core_phylogeny

    def synthetic_tree(files, members, *args):
        result = original(files, members, *args)
        newick = (
            "("
            + ",".join(
                f"{sample.sample_id}:{(index + 1) * 0.001}" for index, sample in enumerate(members)
            )
            + ");"
        )
        files.tree.write_text(newick)
        files.starting_tree.write_text(newick)
        return result

    monkeypatch.setattr(lineage, "_run_core_phylogeny", synthetic_tree)
    selection, _ = selection_fixture(tmp_path)
    run_tree(selection, output=tmp_path / "tree", assess_temporal=False)
    result = run_time(tmp_path / "tree/tree.json", output=tmp_path / "time", allow_unsupported=True)
    assert result["dating_status"] == "dated" and result["dating_override"]
    assert result["clock_confidence"]["rate"] > 0
    assert load_result(tmp_path / "time/time.json", "time") == result
    assert native_stubs == ["corrected"]
    report = tmp_path / "time" / result["report"]["path"]
    assert "data-cc-network-library" in report.read_text()


def test_comparison_retains_failed_job_reason(tmp_path, native_stubs):
    selection, ensemble = selection_fixture(tmp_path, dates=("", "", ""))
    run_tree(selection, output=tmp_path / "tree")
    with pytest.raises(WorkflowError):
        run_time(tmp_path / "tree/tree.json", output=tmp_path / "time", allow_unsupported=True)
    failed = next((tmp_path / "time/jobs").glob("*/job.json"))
    result = compare_time_runs(ensemble, [failed], output=tmp_path / "compare")
    assert result["runs"][0]["status"] == "failed"
    assert "three saved eligible" in result["runs"][0]["reason"]


def test_new_tree_invalidates_time_but_preserves_previous_time_evidence(tmp_path, native_stubs):
    from chronoclade.stage_artifacts import load_result

    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    run_tree(selection, output=tmp_path / "tree")
    first = run_time(tmp_path / "tree/tree.json", output=tmp_path / "time")
    # The source snapshot remains valid while the root tree.json pointer advances.
    run_tree(selection, output=tmp_path / "tree", force=True)
    assert load_result(tmp_path / "time/time.json", "time") == first
    second = run_time(tmp_path / "tree/tree.json", output=tmp_path / "time")
    assert second["fingerprint"] != first["fingerprint"]
    assert second["source"]["tree"]["sha256"] != first["source"]["tree"]["sha256"]
