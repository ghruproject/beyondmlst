"""Scale correctness checks on categorical counts, storage, ordination and NJ."""

import copy
import json
from pathlib import Path
import shutil
import numpy as np
import pytest
from Bio import Phylo

from chronoclade.cgmlst.scale_distance import prepare_records, write_distances, categorical_blocks
from chronoclade.cgmlst.scale_tree import IndexedTree, leading_pcoa, rapidnj_tree
from chronoclade.cgmlst.scale_setup import DEFAULT_RAPIDNJ
from chronoclade.cgmlst.scale import analyse_profiles_scale
from chronoclade.matrix_distances import BinaryDistanceEvidence
from chronoclade.profile_analysis import _pair


def record(name, calls, **kwargs):
    loci = [f"l{i}" for i in range(len(calls))]
    return dict(
        sample_id=name,
        role="context",
        origin="context",
        species="Klebsiella pneumoniae",
        cgmlst_scheme="scheme",
        cgmlst_scheme_version="1",
        cgmlst_loci=loci,
        cgmlst_profile=dict(zip(loci, calls)),
        **kwargs,
    )


def test_exact_categorical_and_missing_evidence_matches_legacy(tmp_path):
    rows = [
        record("a", [1, 2, 3, 4]),
        record("b", [900000, 2, 0, 4]),
        record("c", [1, 0, 3, 9]),
        record("d", [0, 0, 3, 4]),
        record("e", [1, 2, 3, 4], cgmlst_database_sha256="a" * 64),
        record("f", [1, 2, 3, 4], cgmlst_database_sha256="b" * 64),
    ]
    rows = prepare_records(rows)
    evidence, cohorts, exclusions = write_distances(rows, tmp_path, min_overlap=0.6, batch_size=2)
    restored = BinaryDistanceEvidence(evidence.manifest_path)
    assert restored.sample_ids == tuple(row["sample_id"] for row in rows)
    for i, left in enumerate(rows):
        for right in rows[i + 1 :]:
            legacy, reason = _pair(left, right, 0.6)
            observed = restored.get(left["sample_id"], right["sample_id"])
            if reason:
                assert observed is None
                assert (left["sample_id"], right["sample_id"]) not in restored.value_map()
            else:
                for key in legacy:
                    assert observed[key] == legacy[key]
                assert (
                    restored.value_map()[left["sample_id"], right["sample_id"]]
                    == legacy["distance"]
                )
                assert (
                    restored.value_map(raw=True)[left["sample_id"], right["sample_id"]]
                    == legacy["allele_differences"]
                )
    assert exclusions == [{"sample_id": "d", "reason": "insufficient_called_overlap"}]
    for cohort in cohorts:
        assert np.isfinite(evidence.arrays["distance"][np.ix_(cohort, cohort)]).all()
    assert sum(map(len, cohorts)) == 5


def test_storage_rejects_tampering_and_imputed_missing_counts(tmp_path):
    evidence, _, _ = write_distances(
        [record("a", [1, 1]), record("b", [1, 2])], tmp_path, min_overlap=1
    )
    path = tmp_path / "distance.npy"
    changed = np.load(path, mmap_mode="r+")
    changed[0, 1] = 0
    changed.flush()
    with pytest.raises(ValueError, match="checksum"):
        BinaryDistanceEvidence(evidence.manifest_path)
    from chronoclade.artifacts import file_sha256

    manifest = json.loads(evidence.manifest_path.read_text())
    manifest["paths"]["distance"]["sha256"] = file_sha256(path)
    evidence.manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="disagree|Asymmetric"):
        BinaryDistanceEvidence(evidence.manifest_path)


def test_pcoa_leading_coordinates_match_dense_geometry(tmp_path):
    rng = np.random.default_rng(17)
    calls = rng.integers(1, 7, (35, 30)).astype(float)
    matrix = np.zeros((35, 35))
    for a, b, raw, shared in categorical_blocks(calls, batch_size=8):
        matrix[a, b] = raw / shared
        matrix[b, a] = (raw / shared).T
    coords, audit = leading_pcoa(matrix, seed=10, batch_size=7, tolerance=1e-10)
    centered = np.eye(35) - np.ones((35, 35)) / 35
    gram = -0.5 * centered @ matrix**2 @ centered
    values, vectors = np.linalg.eigh(gram)
    expected = vectors[:, -2:] @ np.diag(np.maximum(values[-2:], 0)) @ vectors[:, -2:].T
    assert np.allclose(coords @ coords.T, expected, atol=1e-9)
    assert audit["pcoa_full_spectrum_computed"] is False
    assert audit["pcoa_positive_axis_fraction"] is None


def test_indexed_exact_tree_distances_survive_reroot_and_deepcopy():
    from io import StringIO

    legacy = Phylo.read(StringIO("((a:0.1,b:0.2):0.3,(c:0.4,d:0.5):0.6);"), "newick")
    indexed = IndexedTree.from_tree(copy.deepcopy(legacy))
    for a in ("a", "b", "c", "d"):
        for b in ("a", "b", "c", "d"):
            assert indexed.distance(a, b) == pytest.approx(legacy.distance(a, b))
    indexed.root_with_outgroup("c")
    indexed = copy.deepcopy(indexed)
    legacy.root_with_outgroup("c")
    for a in ("a", "b", "c", "d"):
        assert indexed.distance(a) == pytest.approx(legacy.distance(a))


def test_scale_rejects_requested_bootstrap_before_output(tmp_path):
    with pytest.raises(ValueError, match="explicit bootstrap_replicates=0"):
        analyse_profiles_scale(
            [record("a", [1, 2])], output=tmp_path / "out", bootstrap_replicates=30
        )
    assert not (tmp_path / "out").exists()


def test_rapidnj_missing_never_falls_back(tmp_path):
    with pytest.raises(RuntimeError, match="unavailable"):
        rapidnj_tree(
            np.array([[0, 0.5], [0.5, 0]]),
            ["a", "b"],
            tmp_path / "nj.nwk",
            executable="no-such-rapidnj",
        )
    assert not (tmp_path / "nj.nwk").exists()


@pytest.mark.integration
def test_real_rapidnj_all_tips_and_additive_distances(tmp_path):
    binary = shutil.which("rapidnj") or (str(DEFAULT_RAPIDNJ) if DEFAULT_RAPIDNJ.exists() else None)
    if not binary:
        pytest.skip("RapidNJ source backend is not installed")
    names = ["weird name:1", "b", "c", "d"]
    matrix = np.array(
        [[0, 0.3, 0.9, 1.0], [0.3, 0, 1.0, 1.1], [0.9, 1.0, 0, 0.9], [1.0, 1.1, 0.9, 0]]
    )
    tree, audit = rapidnj_tree(matrix, names, tmp_path / "nj.nwk", executable=binary)
    assert {tip.name for tip in tree.get_terminals()} == set(names)
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            assert tree.distance(a, b) == pytest.approx(matrix[i, j], abs=2e-6)
    assert audit["backend"] == "RapidNJ"
    assert not (tmp_path / "nj.phy").exists()


@pytest.mark.integration
def test_scale_end_to_end_ids_ties_network_and_tree_view(tmp_path):
    binary = shutil.which("rapidnj") or (str(DEFAULT_RAPIDNJ) if DEFAULT_RAPIDNJ.exists() else None)
    if not binary:
        pytest.skip("RapidNJ source backend is not installed")
    rows = [
        record("q", [1, 1, 1, 1], country="UK", collection_date="2020"),
        record("a", [2, 1, 1, 1], country="India", collection_date="2021"),
        record("b", [999, 1, 1, 1], country="India", collection_date="2022"),
    ]
    rows[0].update(role="input", origin="local")
    result = analyse_profiles_scale(rows, output=tmp_path, rapidnj_executable=binary, min_overlap=1)
    assert len(result["nearest_neighbours"]) == 2
    assert all(item["tied_neighbours"] == 2 for item in result["nearest_neighbours"])
    assert result["coverage"]["available_profiles"] == 3
    assert result["cohorts"][0]["tree_display_sample_ids"] == ["a", "b", "q"]
    network = result["location_network"][0]
    assert network["sample_count"] == 3
    assert "all optimal histories" in network["interpretation"]
    view = Path(result["cohorts"][0]["tree_figure"]).read_text()
    assert "Collapse all" in view and "nodes=new Map" in view
    assert result["root_to_tip"][0]["root"] == "q"
    assert result["root_to_tip"][0]["root_sensitivity"]
    assert not (tmp_path / "pairwise_distances.csv").exists()
    assert len(Path(result["cohorts"][0]["pcoa_csv"]).read_text().splitlines()) == 4


def test_nearest_uses_raw_mismatch_ties_with_variable_denominators(tmp_path):
    from chronoclade.cgmlst.scale import _nearest

    rows = [record("q", [1, 1, 1, 1]), record("a", [2, 0, 0, 1]), record("b", [2, 2, 1, 1])]
    rows[0].update(role="input", origin="local")
    rows = prepare_records(rows)
    evidence, cohorts, _ = write_distances(rows, tmp_path, min_overlap=0.5)
    nearest = _nearest(rows, evidence, cohorts, set())
    assert [item["context_id"] for item in nearest] == ["a"]
    assert nearest[0]["allele_differences"] == 1
    assert nearest[0]["shared_called_loci"] == 2
    assert nearest[0]["distance"] == 0.5


def test_requested_csv_export_keeps_unavailable_comparisons(tmp_path):
    import csv

    rows = [record("a", [1, 1, 1]), record("b", [2, 1, 1]), record("c", [1, 0, 0])]
    evidence, _, _ = write_distances(rows, tmp_path / "binary", min_overlap=0.9)
    path = evidence.export_pairs(tmp_path / "requested.csv", pairs=[("a", "b"), ("a", "c")])
    with path.open() as handle:
        exported = list(csv.DictReader(handle))
    assert len(exported) == 2
    assert exported[0]["status"] == "comparable"
    assert exported[1]["status"] == "unavailable"
    assert exported[1]["distance"] == ""
    assert exported[1]["shared_called_loci"] == "1"
    with pytest.raises(ValueError, match="explicit"):
        evidence.export_pairs(tmp_path / "all.csv")


@pytest.mark.parametrize("key,value", [("sample_ids", None), ("scope", []), ("paths", None)])
def test_binary_malformed_types_fail_as_validation_errors(tmp_path, key, value):
    evidence, _, _ = write_distances([record("a", [1, 1]), record("b", [1, 2])], tmp_path)
    data = json.loads(evidence.manifest_path.read_text())
    data[key] = value
    evidence.manifest_path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        BinaryDistanceEvidence(evidence.manifest_path)
