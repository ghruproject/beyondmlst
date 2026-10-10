"""Exact selection, immutable source evidence and independent tree jobs."""

import json
from io import StringIO
from pathlib import Path

import pytest
from Bio import Phylo

from chronoclade.datasets import PreparedDataset, samples_from_records, write_dataset
from chronoclade.errors import WorkflowError
from chronoclade.selection_manifest import write_selection_ensemble
from chronoclade.stage_artifacts import load_tree_result
from chronoclade.tree_stage.workflow import run_tree


def selection_fixture(tmp_path, *, dates=("2018", "2019", "2020")):
    rows = []
    for index, ident in enumerate(("q", "c", "d")):
        assembly = tmp_path / (ident + ".fasta")
        assembly.write_text(">contig\n" + "A" * 2000 + "\n")
        rows.append(
            {
                "sample_id": ident,
                "origin": "local" if ident == "q" else "context",
                "species": "Klebsiella pneumoniae",
                "collection_date": dates[index],
                "country": "Greece" if ident != "d" else "UK",
                "assembly_reference": str(assembly),
            }
        )
    dataset = PreparedDataset(samples=samples_from_records(rows))
    manifest = write_dataset(dataset, tmp_path / "prepared")
    distance = tmp_path / "distances.csv"
    distance.write_text(
        "sample_id_1,sample_id_2,distance,allele_differences,shared_called_loci\nq,c,0.1,1,10\nq,d,0.2,2,10\nc,d,0.1,1,10\n"
    )
    ensemble = write_selection_ensemble(
        rows[:1],
        rows[1:],
        {"paths": {"pairwise_distances": str(distance)}},
        dataset_manifest=manifest,
        output=tmp_path / "selections",
        size=2,
        replicates=2,
    )
    return ensemble.parent / "selections/selection-001.json", ensemble


@pytest.fixture
def native_stubs(monkeypatch):
    import chronoclade.lineage as lineage

    monkeypatch.setattr("chronoclade.tree_stage.workflow.require_tools", lambda names: None)
    monkeypatch.setattr(
        "chronoclade.tree_stage.workflow.tools",
        lambda names: {name: "frozen-test-tool" for name in names},
    )
    calls = []

    def core(files, members, reference, threads, force, mode):
        calls.append(mode)
        ids = [s.sample_id for s in members]
        newick = "(" + ",".join(f"{ident}:0.001" for ident in ids) + ");\n"
        files.starting_tree.write_text(newick)
        files.tree.write_text(newick)
        sequences = "".join(f">{ident}\n" + "A" * 2000 + "\n" for ident in ids)
        files.alignment.write_text(sequences)
        files.filtered_alignment.write_text(sequences)
        files.importations.write_text("Node\tBeg\tEnd\n")
        return files.tree, files.filtered_alignment, None

    monkeypatch.setattr(lineage, "_run_core_phylogeny", core)
    return calls


def test_tree_exact_ids_hashes_and_independent_resume(tmp_path, native_stubs):
    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    result = run_tree(selection, output=tmp_path / "tree")
    assert result["selected_sample_ids"] == ["q", "c", "d"]
    assert result["undated_sample_ids"] == ["q", "c", "d"]
    assert not result["temporal_assessment"]["supported"]
    assert load_tree_result(tmp_path / "tree/tree.json") == result
    assert run_tree(selection, output=tmp_path / "tree") == result
    assert native_stubs == ["corrected"]
    (tmp_path / "q.fasta").write_text(">contig\n" + "C" * 2000 + "\n")
    changed = run_tree(selection, output=tmp_path / "tree")
    assert changed["fingerprint"] != result["fingerprint"]
    assert len(native_stubs) == 2
    original_job = tmp_path / "tree" / result["corrected_tree"]["path"]
    assert original_job.is_file()
    report = tmp_path / "tree" / changed["report"]["path"]
    assert "data-cc-network-library" in report.read_text()
    assert "Date randomisations" not in report.read_text()


def test_failed_replacement_preserves_last_manifest(tmp_path, native_stubs, monkeypatch):
    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    saved = run_tree(selection, output=tmp_path / "tree")

    def broken(*args, **kwargs):
        raise WorkflowError("native failure")

    monkeypatch.setattr("chronoclade.lineage._run_core_phylogeny", broken)
    with pytest.raises(WorkflowError, match="native failure"):
        run_tree(selection, output=tmp_path / "tree", force=True)
    assert load_tree_result(tmp_path / "tree/tree.json") == saved
    assert any(
        json.loads(path.read_text())["status"] == "failed"
        for path in (tmp_path / "tree/jobs").glob("*/job.json")
    )


def test_acquisition_failure_has_audit_without_replacement(tmp_path):
    from chronoclade.tree_stage.workflow import acquire_assemblies

    row = {"sample_id": "missing", "source_genome_id": "exact-id"}
    with pytest.raises(WorkflowError, match="no replacements"):
        acquire_assemblies(
            [row],
            tmp_path / "dataset.json",
            tmp_path,
            api_key="",
            base_url="https://pathogen.watch",
            threads=1,
        )
    assert json.loads((tmp_path / "assembly_audit.json").read_text())[0]["sample_id"] == "missing"


def test_tree_rejects_changed_artifact_and_source(tmp_path, native_stubs):
    selection, _ = selection_fixture(tmp_path, dates=("", "", ""))
    result = run_tree(selection, output=tmp_path / "tree")
    artifact = tmp_path / "tree" / result["corrected_tree"]["path"]
    artifact.write_text("(q:1,c:1,outsider:1);")
    with pytest.raises(WorkflowError, match="sha256"):
        load_tree_result(tmp_path / "tree/tree.json")


def test_clustered_dates_preserve_exact_partial_intervals():
    import random
    from chronoclade.metadata import Sample
    from chronoclade.tree_stage.temporal import date_clusters, permute_cluster_dates

    tree = Phylo.read(StringIO("((a:1,b:1):1,c:1,d:1);"), "newick")
    samples = [
        Sample(i, Path("."), date, "", "sp", "lin", "local")
        for i, date in zip("abcd", ("2018", "2018", "2019-02", "2020-03-04"))
    ]
    groups = date_clusters(tree, samples)
    assert [group["sample_ids"] for group in groups] == [["a", "b"], ["c"], ["d"]]
    shuffled = permute_cluster_dates(groups, random.Random(3))
    assert shuffled["a"] == shuffled["b"]
    assert set(shuffled.values()) == {"2018", "2019-02", "2020-03-04"}


def test_clustered_assessment_requires_complete_null_refits(tmp_path, monkeypatch):
    from chronoclade.metadata import Sample
    from chronoclade.tree_stage.temporal import run_clustered_assessment

    tree = tmp_path / "tree.nwk"
    tree.write_text("(a:0.001,b:0.002,c:0.003);")
    clock = tmp_path / "clock.txt"
    clock.write_text("--rate: 0.001\n--r^2: 0.95\n")
    samples = [
        Sample(i, Path("."), value, "", "species", "lineage", "local")
        for i, value in zip("abc", ("2018", "2019-02", "2020-03-04"))
    ]
    monkeypatch.setattr(
        "chronoclade.tree_stage.temporal._run_randomised_clock",
        lambda **kwargs: {"rate": 0.0001, "r_squared": 0.1},
    )
    args = dict(
        tree=tree,
        sequence_length=2000,
        samples=samples,
        observed_clock=clock,
        randomisations=20,
        randomisation_jobs=2,
        seed=1,
        output=tmp_path / "assessment.json",
        p_value_threshold=0.05,
    )
    result = run_clustered_assessment(**args)
    assert result["supported"] and result["temporal"]["cluster_count"] == 3
    assert result["temporal"]["p_value_r_squared"] == 1 / 21
    monkeypatch.setattr(
        "chronoclade.tree_stage.temporal._run_randomised_clock", lambda **kwargs: None
    )
    incomplete = run_clustered_assessment(**args)
    assert not incomplete["supported"] and incomplete["code"] == "not_assessed"
