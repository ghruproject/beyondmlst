"""Stage sequencing and selected-assembly boundaries, without native executables."""

import csv
import json
import zipfile
from pathlib import Path

import pytest

from chronoclade.errors import WorkflowError
from chronoclade.staged_workflow import run_staged_workflow, select_assembly_context


def row(identifier, **kwargs):
    return {
        "sample_id": identifier,
        "species": "Klebsiella pneumoniae",
        "lineage": "ST147",
        "origin": "context",
        "collection_date": "2020",
        "country": "UK",
        **kwargs,
    }


def pair_table(path, triples):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id_1", "sample_id_2", "distance"])
        writer.writerows(triples)
    return {"paths": {"pairwise_distances": str(path)}, "genetic_groups": []}


def test_selection_fair_neighbours_then_metadata_and_diversity(tmp_path):
    queries = [row("q1", origin="local"), row("q2", origin="local")]
    contexts = [
        row(f"c{i}", collection_date=str(2000 + i), country=f"Country{i}") for i in range(1, 9)
    ]
    triples = [("q1", f"c{i}", i / 100) for i in range(1, 9)]
    triples += [("q2", f"c{i}", abs(4 - i) / 100) for i in range(1, 9)]
    triples += [(f"c{i}", f"c{j}", abs(i - j) / 100) for i in range(1, 9) for j in range(i + 1, 9)]
    analysis = pair_table(tmp_path / "pairs.csv", triples)
    selected, audit = select_assembly_context(queries, contexts, analysis, size=6, seed=11)
    assert len(selected) == 6
    assert [r["sample_id"] for r in selected[:2]] == ["c1", "c4"]
    assert any(
        r["selection_reason"] == "genetic_group_time_region_representative" for r in selected
    )
    assert any(r["selection_reason"] == "cgmlst_diversity_representative" for r in selected)
    assert audit["queries_without_selected_comparable_neighbour"] == []
    reversed_selection, reversed_audit = select_assembly_context(
        list(reversed(queries)), list(reversed(contexts)), analysis, size=6, seed=11
    )
    assert selected == reversed_selection
    assert audit == reversed_audit
    assert all("selection_reason" not in r for r in contexts)


def test_pins_aliases_budget_and_background_without_profiles(tmp_path):
    contexts = [row("a", accession="ACC1"), row("b", source_genome_id="UUID2"), row("c")]
    selected, audit = select_assembly_context(
        [row("q")], contexts, {}, size=2, include=["ACC1", "a", "UUID2"]
    )
    assert [r["sample_id"] for r in selected] == ["a", "b"]
    assert all(r["selection_reason"] == "user_requested" for r in selected)
    assert audit["queries_without_selected_comparable_neighbour"] == ["q"]
    with pytest.raises(WorkflowError, match="exceed"):
        select_assembly_context([], contexts, {}, size=1, include=["a", "b"])
    with pytest.raises(WorkflowError, match="uniquely"):
        select_assembly_context([], contexts, {}, include=["missing"])
    assert select_assembly_context([], contexts, {}, size=0)[0] == []
    selected, _ = select_assembly_context([], contexts, {}, size=50)
    assert len(selected) == 3


@pytest.fixture
def stage_mocks(monkeypatch):
    calls = {"resolve": [], "analysis": [], "materialise": [], "workflow": []}
    inputs = {
        "queries": [row("q", origin="local")],
        "context": [row("c", accession="ERR123")],
        "provenance": {"source": "frozen fixture"},
    }

    def resolve(*args, **kwargs):
        calls["resolve"].append(kwargs)
        return inputs

    def analyse(records, *, output, **kwargs):
        calls["analysis"].append(records)
        return pair_table(output / "pairs.csv", [("q", "c", 0.01)])

    def report(analysis, *, directory, **kwargs):
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / "report.html"
        path.write_text("fast profiles only")
        return path

    def index(output, stages):
        (output / "index.html").write_text(json.dumps(stages))

    def materialise(records, **kwargs):
        calls["materialise"].append(records)
        return records

    def workflow(samples, *, mode, output, **kwargs):
        calls["workflow"].append(mode)
        directory = output / "Klebsiella_pneumoniae__ST147"
        directory.mkdir(parents=True, exist_ok=True)
        for name in ("report.html", "report.json", "supporting_results.zip"):
            (directory / name).write_text(mode)
        (directory / "genetic_tree.svg").write_text(mode)
        (directory / "nearest_neighbours.csv").write_text(mode)
        return {
            "lineages": [
                {
                    "species": "Klebsiella pneumoniae",
                    "lineage": "ST147",
                    "status": "completed",
                    "analysis_mode": mode,
                    "outputs": {"html_report": str(directory / "report.html")},
                }
            ]
        }

    monkeypatch.setattr("chronoclade.profile_inputs.resolve_profile_inputs", resolve)
    monkeypatch.setattr("chronoclade.profile_inputs.materialise_assemblies", materialise)
    monkeypatch.setattr("chronoclade.profile_analysis.analyse_profiles", analyse)
    monkeypatch.setattr("chronoclade.profile_report.write_profile_report", report)
    monkeypatch.setattr("chronoclade.profile_report.write_stage_index", index)
    monkeypatch.setattr("chronoclade.workflow.run_workflow", workflow)
    return calls, inputs


def test_fast_never_materialises_or_invokes_native_workflow(tmp_path, stage_mocks):
    calls, _ = stage_mocks
    result = run_staged_workflow(None, collection="collection-id", output=tmp_path, mode="fast")
    assert calls["materialise"] == calls["workflow"] == []
    assert len(calls["analysis"]) == 1
    assert set(result["stages"]) == {"fast"}
    assert Path(result["stages"]["fast"]["report"]).read_text() == "fast profiles only"


def test_full_only_corrected_and_downloads_selected_records(tmp_path, stage_mocks):
    calls, _ = stage_mocks
    result = run_staged_workflow(
        None, collection="id", output=tmp_path, mode="full", context_size=1
    )
    assert calls["workflow"] == ["corrected"]
    assert [r["sample_id"] for r in calls["materialise"][0]] == ["q", "c"]
    assert json.loads(
        (tmp_path / "assembly" / "Klebsiella_pneumoniae__ST147" / "sample_labels.json").read_text()
    ) == {"q": "q", "c": "ERR123"}
    assert set(result["stages"]) == {"fast", "full"}
    directory = tmp_path / "assembly" / "Klebsiella_pneumoniae__ST147"
    assert (directory / "report.full.html").read_text() == "corrected"
    assert not (directory / "report.finish.html").exists()


def test_finish_sequences_and_preserves_independent_reports(tmp_path, stage_mocks):
    calls, _ = stage_mocks
    result = run_staged_workflow(
        None, collection="id", output=tmp_path, mode="finish", context_size=1
    )
    assert calls["workflow"] == ["corrected", "full"]
    assert len(calls["materialise"]) == 1
    assert set(result["stages"]) == {"fast", "full", "finish"}
    directory = tmp_path / "assembly" / "Klebsiella_pneumoniae__ST147"
    assert (directory / "report.full.html").read_text() == "corrected"
    assert (directory / "report.finish.html").read_text() == "full"
    with zipfile.ZipFile(directory / "supporting_results.full.zip") as archive:
        assert json.loads(archive.read("report.json"))["analysis_mode"] == "corrected"
    with zipfile.ZipFile(directory / "supporting_results.finish.zip") as archive:
        assert json.loads(archive.read("report.json"))["analysis_mode"] == "full"
    assert json.loads((tmp_path / "full-summary.json").read_text())["lineages"][0]["outputs"][
        "html_report"
    ].endswith("stages/full/report.html")
    assert json.loads((tmp_path / "finish-summary.json").read_text())["lineages"][0]["outputs"][
        "html_report"
    ].endswith("stages/finish/report.html")


def test_dry_run_has_no_fetch_typing_analysis_or_materialisation(tmp_path, stage_mocks):
    calls, _ = stage_mocks
    result = run_staged_workflow(
        None, collection="id", output=tmp_path, mode="finish", dry_run=True
    )
    assert result["stages"] == ["fast", "full", "finish"]
    assert all(not values for values in calls.values())


def test_changed_selection_invalidates_stale_dated_outputs(tmp_path, stage_mocks):
    calls, _ = stage_mocks
    run_staged_workflow(None, collection="id", output=tmp_path, mode="finish", context_size=1)
    directory = tmp_path / "assembly" / "Klebsiella_pneumoniae__ST147"
    dated = directory / "timetree"
    dated.mkdir()
    (dated / "stale.nwk").write_text("stale")
    result = run_staged_workflow(
        None, collection="id", output=tmp_path, mode="fast", context_size=0
    )
    assert set(result["stages"]) == {"fast"}
    assert not (directory / "report.finish.html").exists()
    assert not (tmp_path / "finish.html").exists()
    assert not dated.exists()


def test_invalid_pin_rejected_before_any_assembly_download(tmp_path, stage_mocks):
    calls, _ = stage_mocks
    with pytest.raises(WorkflowError, match="absent|pool|Requested"):
        run_staged_workflow(
            None, collection="id", output=tmp_path, mode="full", include_genomes=["missing"]
        )
    assert calls["materialise"] == []


def test_pin_outside_query_lineage_is_not_silently_ignored(tmp_path, stage_mocks):
    calls, inputs = stage_mocks
    inputs["context"].append(row("other", lineage="ST1"))
    with pytest.raises(WorkflowError, match="lineage|select|Requested"):
        run_staged_workflow(
            None, collection="id", output=tmp_path, mode="full", include_genomes=["other"]
        )
    assert calls["materialise"] == []


def test_stage_snapshots_keep_figures_and_tables_immutable(tmp_path, stage_mocks):
    calls, _ = stage_mocks
    result = run_staged_workflow(
        None, collection="id", output=tmp_path, mode="finish", context_size=1
    )
    full_directory = Path(
        json.loads((tmp_path / "full-summary.json").read_text())["lineages"][0]["outputs"][
            "html_report"
        ]
    ).parent
    finished_directory = Path(result["lineages"][0]["outputs"]["html_report"]).parent
    assert full_directory != finished_directory
    assert json.loads((full_directory / "sample_labels.json").read_text()) == {
        "q": "q",
        "c": "ERR123",
    }
    assert json.loads((finished_directory / "sample_labels.json").read_text()) == {
        "q": "q",
        "c": "ERR123",
    }
    assert (full_directory / "genetic_tree.svg").read_text() == "corrected"
    assert (full_directory / "nearest_neighbours.csv").read_text() == "corrected"
    assert (finished_directory / "genetic_tree.svg").read_text() == "full"
    assert (finished_directory / "nearest_neighbours.csv").read_text() == "full"
    assert not (full_directory / "stages").exists()


def test_corrected_snapshot_bundle_excludes_stale_dated_evidence(tmp_path, stage_mocks):
    from chronoclade.staged_workflow import _snapshot_stage

    directory = tmp_path / "lineage"
    directory.mkdir()
    (directory / "report.html").write_text('<img src="genetic_tree.svg">')
    (directory / "report.json").write_text("{}")
    (directory / "genetic_tree.svg").write_text("corrected")
    (directory / "date_randomisation.csv").write_text("stale dated result")
    (directory / "node_dates.csv").write_text("stale inferred dates")
    (directory / "timetree").mkdir()
    (directory / "timetree" / "timetree.nexus").write_text("dated")
    (directory / "cohort_nj.nwk").write_text("(a,b);")
    from chronoclade.report import write_supporting_bundle

    write_supporting_bundle(directory)
    record = {
        "outputs": {
            "html_report": str(directory / "report.html"),
            "tree_figure": str(directory / "genetic_tree.svg"),
        }
    }
    _snapshot_stage(record, "full")
    snapshot = Path(record["outputs"]["html_report"]).parent
    assert record["outputs"]["html_report"].endswith("stages/full/report.html")
    assert json.loads((snapshot / "report.json").read_text())["outputs"]["tree_figure"] == str(
        snapshot / "genetic_tree.svg"
    )
    with zipfile.ZipFile(snapshot / "supporting_results.zip") as archive:
        assert "cohort_nj.nwk" in archive.namelist()
        assert "genetic_tree.svg" in archive.namelist()
        assert not any(
            "date_randomisation" in name or "node_dates" in name or "timetree" in name
            for name in archive.namelist()
        )


def test_unsupported_finish_snapshot_excludes_stale_dated_products(tmp_path, stage_mocks):
    from chronoclade.staged_workflow import _snapshot_stage
    from chronoclade.report import write_supporting_bundle

    directory = tmp_path / "lineage"
    directory.mkdir()
    (directory / "report.html").write_text("Date test did not support dating")
    (directory / "report.json").write_text("{}")
    (directory / "temporal_signal.json").write_text('{"supported":false}')
    (directory / "date_randomisation.csv").write_text("real permutation evidence")
    (directory / "node_dates.csv").write_text("stale inferred dates")
    (directory / "timetree_with_confidence.svg").write_text("stale dated figure")
    (directory / "timetree").mkdir()
    (directory / "timetree" / "timetree.nexus").write_text("stale dated tree")
    write_supporting_bundle(directory)
    record = {
        "temporal_signal_supported": False,
        "outputs": {"html_report": str(directory / "report.html")},
    }
    _snapshot_stage(record, "finish")
    snapshot = Path(record["outputs"]["html_report"]).parent
    with zipfile.ZipFile(snapshot / "supporting_results.zip") as archive:
        assert "temporal_signal.json" in archive.namelist()
        assert "date_randomisation.csv" in archive.namelist()
        assert not any("node_dates" in name or "timetree" in name for name in archive.namelist())


@pytest.mark.parametrize("stale_default", [False, True])
def test_snapshot_uses_advertised_corrected_report_and_normalises_entry(tmp_path, stale_default):
    from chronoclade.staged_workflow import _snapshot_stage

    directory = tmp_path / "lineage"
    directory.mkdir()
    report = directory / "corrected_report.html"
    report.write_text('<img src="genetic_tree.svg"><a href="../../fast/report.html">Fast</a>')
    (directory / "genetic_tree.svg").write_text("corrected relationships")
    if stale_default:
        (directory / "report.html").write_text("stale dated report")
    record = {
        "analysis_mode": "corrected",
        "outputs": {"html_report": str(report)},
    }
    _snapshot_stage(record, "full")
    entry = directory / "stages" / "full" / "report.html"
    assert record["outputs"]["html_report"] == str(entry)
    assert 'src="genetic_tree.svg"' in entry.read_text()
    assert 'href="../../../../fast/report.html"' in entry.read_text()
    assert "stale dated report" not in entry.read_text()
    saved = json.loads((entry.parent / "report.json").read_text())
    assert saved["outputs"]["html_report"] == str(entry)
    assert 'src="stages/full/genetic_tree.svg"' in (directory / "report.full.html").read_text()
    with zipfile.ZipFile(entry.parent / "supporting_results.zip") as archive:
        assert archive.read("genetic_tree.svg").decode() == "corrected relationships"
        assert json.loads(archive.read("report.json"))["outputs"]["html_report"] == str(entry)


def test_cg_datasets_remain_separate_through_assembly_selection(tmp_path, stage_mocks, monkeypatch):
    calls, inputs = stage_mocks
    inputs["queries"][0]["analysis_dataset"] = "ST147_CG147"
    inputs["context"][0]["analysis_dataset"] = "ST147_CG147"
    inputs["queries"].append(row("q2", origin="local", analysis_dataset="ST147_CG2"))
    inputs["context"].append(row("c2", analysis_dataset="ST147_CG2"))

    def analyse(records, *, output, **kwargs):
        queries = [r for r in records if r["origin"] == "local"]
        public = [r for r in records if r["origin"] == "context"]
        return pair_table(output / "pairs.csv", [(q["sample_id"], c["sample_id"], 0.01)
                                                for q in queries for c in public])

    monkeypatch.setattr("chronoclade.profile_analysis.analyse_profiles", analyse)
    result = run_staged_workflow(None, collection="id", output=tmp_path, mode="full", context_size=1)
    assert len(result["fast_datasets"]) == 2
    assert {frozenset(r["sample_id"] for r in group) for group in calls["materialise"]} == {
        frozenset({"q", "c"}), frozenset({"q2", "c2"})}
    assert {r["lineage"] for group in calls["materialise"] for r in group} == {
        "ST147_CG147", "ST147_CG2"}
