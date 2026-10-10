"""Genome exploration preserves locus identity, fixed denominators and selections."""
from dataclasses import replace
import json
import shutil

import numpy as np
import pytest

from chronoclade.datasets import write_dataset
from chronoclade.esm2.genome_analysis import analyse_genomes, nearest_ties
from chronoclade.esm2.proteins import EmbeddingError
from chronoclade.esm2.workflow import run_esm2
from chronoclade.selection_manifest import load_selection_ensemble
from test_esm2_sample_distances import dataset, artifacts


def cohort():
    prepared = dataset()
    samples = tuple(dict(row, role="input" if row["sample_id"] == "a" else "context",
                         species="Klebsiella pneumoniae", mlst_st="147", mlst_scheme="test:mlst")
                    for row in prepared.samples)
    lineages = tuple({"sample_id": row["sample_id"], "kind": "cglin", "scheme_id": "test:lin",
                      "scheme_version": "v1", "database_version": "db1", "assignment_method": "frozen-fixture", "assignment": "1,2,3,4,5,6,7", "resolution_status": "resolved"}
                     for row in samples)
    return replace(prepared, samples=samples, lineages=lineages)


def test_locus_identity_and_chord_geometry(tmp_path):
    paths = artifacts(tmp_path, [("a", "x", "px"), ("a", "y", "py"),
                                 ("b", "x", "py"), ("b", "y", "px"),
                                 ("c", "x", "px"), ("c", "y", "py")])
    summary, arrays, pairs, _ = analyse_genomes(dataset(), *paths)
    assert arrays["distances"][0, 1] == pytest.approx(1)
    assert arrays["distances"][0, 2] == 0
    chord_squared = np.sum((arrays["coordinates"][0] - arrays["coordinates"][1]) ** 2)
    assert chord_squared == pytest.approx(2)
    assert arrays["mapped_mask"].all()
    assert summary["groups"]["by_sample"]["a"] == summary["groups"]["by_sample"]["c"]
    assert summary["biological_comparison"]["synonymous_protein_identity"]
    assert pairs[1]["protein_identical_changed_alleles"] == 2


def test_fixed_panel_retains_missing_mask_and_excluded_context(tmp_path):
    paths = artifacts(tmp_path, [("a", "x", "px"), ("a", "y", "px"), ("b", "x", "py"),
                                 ("c", "x", "nx"), ("c", "y", "nx")])
    summary, arrays, _, _ = analyse_genomes(dataset(), *paths, panel_loci=["x", "y"])
    assert summary["sample_ids"] == ["a", "c"]
    assert arrays["mapped_mask"].tolist() == [[True, True], [True, False], [True, True]]
    assert summary["exclusions"][0]["missing_panel_loci"] == ["y"]
    assert summary["panel"]["count"] == 2


def test_all_nearest_boundary_ties_retained():
    matrix = np.array([[0, 1, 1], [1, 0, 0], [1, 0, 0]], dtype=float)
    assert len(nearest_ties(matrix, ["a", "b", "c"], 1)[0]["neighbours"]) == 2


def test_full_saved_vector_stage_selection_and_portability(tmp_path):
    prepared = cohort()
    dataset_path = write_dataset(prepared, tmp_path / "prepared")
    embeddings, mapping = artifacts(tmp_path)
    result = run_esm2(dataset_path, tmp_path / "esm2", embeddings_manifest=embeddings,
                      mapping_csv=mapping, selection_runs=2, selection_size=1,
                      nearest_per_input=1, baseline=True)
    manifest = json.loads(result.manifest_path.read_text())
    assert manifest["status"] == "complete"
    block = manifest["blocks"][0]
    ensemble = load_selection_ensemble(result.manifest_path.parent / block["ensemble"])
    assert ensemble["source"]["method"] == "esm2"
    assert ensemble["target_ids"] == ["a"]
    assert (tmp_path / "esm2" / block["report"]).is_file()
    assert list((tmp_path / "esm2").rglob("protein_ordination.html"))
    assert list((tmp_path / "esm2").rglob("temporal_diagnostics.json"))
    baseline_path = next((tmp_path / "esm2").rglob("profile_analysis.json"))
    baseline = json.loads(baseline_path.read_text())
    assert ".pending-" not in baseline_path.read_text()
    assert (baseline_path.parent / baseline["cohorts"][0]["tree_path"]).is_file()
    moved = tmp_path / "moved"
    moved.mkdir()
    shutil.move(tmp_path / "prepared", moved / "prepared")
    shutil.move(tmp_path / "esm2", moved / "esm2")
    assert load_selection_ensemble(moved / "esm2" / block["ensemble"])["source"]["method"] == "esm2"


def test_missing_input_blocks_selection_instead_of_dropping_input(tmp_path):
    dataset_path = write_dataset(cohort(), tmp_path / "prepared")
    paths = artifacts(tmp_path, [("a", "x", "px"), ("b", "x", "py"), ("b", "y", "py"),
                                 ("c", "x", "nx"), ("c", "y", "nx")])
    result = run_esm2(dataset_path, tmp_path / "esm2", embeddings_manifest=paths[0],
                      mapping_csv=paths[1], panel_loci=["x", "y"], baseline=False)
    manifest = json.loads(result.manifest_path.read_text())
    assert manifest["status"] == "partial"
    assert manifest["blocks"][0]["ensemble"] is None
    assert manifest["blocks"][0]["missing_mandatory_ids"] == ["a"]


def test_no_fake_allele_id_sequence_path(tmp_path):
    dataset_path = write_dataset(cohort(), tmp_path / "prepared")
    with pytest.raises(EmbeddingError, match="allele catalogue"):
        run_esm2(dataset_path, tmp_path / "esm2")
    assert not (tmp_path / "esm2").exists()
