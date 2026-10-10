"""Fixed-panel ESM2 diagnostics must preserve identity, missingness and date uncertainty."""

import csv
from dataclasses import replace
import hashlib
import json

import numpy as np
import pytest

from chronoclade.datasets import AlleleMatrix, LocusCatalogue, PreparedDataset, samples_from_records
from chronoclade.esm2.proteins import EmbeddingError, file_sha256
from chronoclade.esm2.runtime import POOLING
from chronoclade.esm2.sample_distances import analyze_embedding_dates
from chronoclade.esm2.storage import write_npz


def dataset(dates=None, profiles=None, order=None):
    dates = dates or {"a": "2023-01-01", "b": "2022-01-01", "c": "2021-01-01"}
    profiles = profiles or {
        "a": {"x": "1", "y": "1"},
        "b": {"x": "2", "y": "2"},
        "c": {"x": "3", "y": "3"},
    }
    order = order or list(dates)
    samples = samples_from_records(
        [
            {"sample_id": ident, "collection_date": dates[ident], "origin": "local"}
            for ident in order
        ]
    )
    catalogue = LocusCatalogue("test:cgmlst", "v1", ("x", "y"), "frozen-test", "db1")
    matrix = AlleleMatrix.from_profiles(catalogue, order, profiles)
    return PreparedDataset(samples, (matrix,))


def artifacts(tmp_path, mappings=None, directions=None):
    """Small vectors with a real ESM2 manifest contract, without loading weights."""
    directions = directions or {"px": (1, 0), "py": (0, 1), "nx": (-1, 0)}
    mappings = mappings or [
        ("a", "x", "px"),
        ("a", "y", "px"),
        ("b", "x", "py"),
        ("b", "y", "py"),
        ("c", "x", "nx"),
        ("c", "y", "nx"),
    ]
    ids = {name: hashlib.sha256(name.encode()).hexdigest() for name in directions}
    vectors = np.zeros((len(ids), 320), dtype=np.float32)
    for index, direction in enumerate(directions.values()):
        vectors[index, : len(direction)] = direction
    vectors_path = tmp_path / "vectors.npz"
    write_npz(vectors_path, protein_ids=np.asarray(list(ids.values())), embeddings=vectors)
    records = [
        {"record_id": f"{sample}.{locus}", "protein_id": ids[protein]}
        for sample, locus, protein in mappings
    ]
    manifest = {
        "schema": "chronoclade.esm2.protein-embeddings",
        "schema_version": 1,
        "status": "complete",
        "parameters": {"requested_model": "8M"},
        "actual_device": "cpu",
        "protein_ids": list(ids.values()),
        "mapping": records,
        "provenance": {
            "model": "esm2_t6_8M_UR50D",
            "model_sha256": "a" * 64,
            "dimension": 320,
            "representation_layer": 6,
            "actual_device": "cpu",
            "precision": "float32",
            "pooling": POOLING,
        },
        "artifacts": {
            "embeddings": {
                "path": vectors_path.name,
                "shape": list(vectors.shape),
                "dtype": "float32",
                "sha256": file_sha256(vectors_path),
            }
        },
    }
    path = tmp_path / "embeddings.json"
    path.write_text(json.dumps(manifest))
    mapping_csv = tmp_path / "mapping.csv"
    with mapping_csv.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["sample_id", "locus", "record_id"])
        writer.writerows((sample, locus, f"{sample}.{locus}") for sample, locus, _ in mappings)
    return path, mapping_csv


def mutate_manifest(path, update):
    manifest = json.loads(path.read_text())
    update(manifest)
    path.write_text(json.dumps(manifest))


def test_real_sample_fixed_locus_distances_negative_slope_and_residuals(tmp_path):
    result = analyze_embedding_dates(dataset(), *artifacts(tmp_path))
    assert result["reference"]["sample_id"] == "a"
    assert result["reference"]["date_independent"]
    assert result["panel"]["loci"] == ["x", "y"]
    assert [row["distance"] for row in result["samples"]] == pytest.approx([0, 1, 2])
    assert result["regression"]["slope"] == pytest.approx(-1)
    assert result["regression"]["r_squared"] == pytest.approx(1)
    assert all(abs(row["residual"]) < 1e-8 for row in result["points"])
    assert result["units"]["slope"] == "mean locus cosine distance/year"
    json.dumps(result, allow_nan=False)


def test_reference_and_distances_are_invariant_to_metadata_dates_and_input_order(tmp_path):
    paths = artifacts(tmp_path)
    first = analyze_embedding_dates(dataset(), *paths)
    second = analyze_embedding_dates(
        dataset({"a": "", "b": "1999", "c": "2003-02"}, order=["c", "b", "a"]), *paths
    )
    assert first["reference"] == second["reference"]
    assert [(r["sample_id"], r["distance"]) for r in first["samples"]] == [
        (r["sample_id"], r["distance"]) for r in second["samples"]
    ]
    assert len(second["points"]) == 2
    assert second["samples"][0]["distance"] == 0
    assert second["samples"][0]["date_status"] == "missing"
    assert second["samples"][1]["date_interval"]["precision"] == "year"
    assert second["samples"][2]["date_interval"]["end"] == "2003-02-28"
    assert second["regression"]["status"] == "insufficient_data"


def test_canonical_intervals_override_raw_collection_date(tmp_path):
    prepared = dataset()
    samples = list(prepared.samples)
    samples[1] = {
        **samples[1],
        "collection_date": "1900",
        "date_start": "2020-02-01",
        "date_end": "2022-03-01",
        "date_precision": "interval",
    }
    result = analyze_embedding_dates(
        replace(prepared, samples=tuple(samples)), *artifacts(tmp_path)
    )
    point = next(row for row in result["points"] if row["sample_id"] == "b")
    assert point["date_start"] == "2020-02-01"
    assert point["date_end"] == "2022-03-01"
    assert point["date_lower"] > 2020
    assert point["date_upper"] > 2022


def test_synonymous_alleles_can_share_protein_and_constant_distance_is_preserved(tmp_path):
    mappings = [(sample, locus, "px") for sample in ("a", "b", "c") for locus in ("x", "y")]
    result = analyze_embedding_dates(dataset(), *artifacts(tmp_path, mappings))
    assert [row["distance"] for row in result["samples"]] == [0, 0, 0]
    assert result["regression"]["slope"] == 0
    assert result["regression"]["pearson_r"] is None


def test_locus_distances_are_averaged_instead_of_cosine_of_genome_average(tmp_path):
    mappings = [
        ("a", "x", "px"),
        ("a", "y", "py"),
        ("b", "x", "py"),
        ("b", "y", "px"),
        ("c", "x", "px"),
        ("c", "y", "py"),
    ]
    result = analyze_embedding_dates(dataset(), *artifacts(tmp_path, mappings))
    assert result["samples"][1]["distance"] == pytest.approx(1)
    assert result["samples"][2]["distance"] == 0


def test_sparse_mapping_has_one_fixed_intersection_and_audits_wholly_unmapped(tmp_path):
    prepared = dataset(
        {"a": "2020", "b": "2021", "c": "", "d": "2022"},
        profiles={"a": {"x": "1", "y": "1"}, "b": {"x": "2"}, "c": {"x": "3"}, "d": {}},
    )
    paths = artifacts(tmp_path, [("a", "x", "px"), ("a", "y", "px"), ("b", "x", "py")])
    result = analyze_embedding_dates(prepared, *paths)
    assert result["panel"]["loci"] == ["x"]
    assert result["counts"]["total_samples"] == 4
    assert result["counts"]["eligible_samples"] == 2
    assert result["samples"][2]["reason"] == "no_mapped_proteins"
    assert result["samples"][3]["reason"] == "missing_profile"
    explicit = analyze_embedding_dates(prepared, *paths, panel_loci=["x", "y"])
    assert explicit["samples"][1]["distance"] is None
    assert explicit["samples"][1]["reason"] == "incomplete_fixed_panel"
    assert explicit["samples"][1]["missing_panel_loci"] == ["y"]
    with pytest.raises(EmbeddingError, match="Explicit reference"):
        analyze_embedding_dates(prepared, *paths, panel_loci=["x", "y"], reference_sample_id="b")


def test_no_shared_loci_fails_instead_of_pairwise_missingness(tmp_path):
    paths = artifacts(tmp_path, [("a", "x", "px"), ("b", "y", "py")])
    with pytest.raises(EmbeddingError, match="at least one"):
        analyze_embedding_dates(dataset(), *paths)


@pytest.mark.parametrize(
    "change,match",
    [
        (("unknown", "x", "a.x"), "unknown sample"),
        (("a", "other-scheme", "a.x"), "outside the selected scheme"),
        (("a", "x", "not-a-record"), "unknown retained record"),
        (("a", "x", "a.x"), "Duplicate mapping cell"),
    ],
)
def test_unknown_identity_wrong_scheme_and_duplicate_mapping_fail(tmp_path, change, match):
    paths = artifacts(tmp_path)
    with paths[1].open("a", newline="") as stream:
        csv.writer(stream).writerow(change)
    with pytest.raises(EmbeddingError, match=match):
        analyze_embedding_dates(dataset(), *paths)


def test_mapped_missing_allele_and_inconsistent_same_allele_protein_fail(tmp_path):
    paths = artifacts(tmp_path)
    prepared = dataset(profiles={"a": {"x": "1", "y": "1"}, "b": {}, "c": {}})
    with pytest.raises(EmbeddingError, match="has no allele call"):
        analyze_embedding_dates(prepared, *paths)
    prepared = dataset(profiles={sample: {"x": "1", "y": "1"} for sample in ("a", "b", "c")})
    with pytest.raises(EmbeddingError, match="inconsistent protein IDs"):
        analyze_embedding_dates(prepared, *paths)


def test_exactly_one_scheme_required(tmp_path):
    paths = artifacts(tmp_path)
    prepared = dataset()
    with pytest.raises(EmbeddingError, match="exactly one"):
        analyze_embedding_dates(replace(prepared, profiles=()), *paths)
    other = replace(
        prepared.profiles[0], catalogue=replace(prepared.profiles[0].catalogue, scheme_id="other")
    )
    with pytest.raises(EmbeddingError, match="exactly one"):
        analyze_embedding_dates(replace(prepared, profiles=(*prepared.profiles, other)), *paths)


@pytest.mark.parametrize(
    "kind",
    [
        "checksum",
        "shape",
        "non-finite",
        "zero",
        "wrong-model",
        "model-hash",
        "pooling",
        "device",
        "id-order",
        "dtype",
    ],
)
def test_saved_artifact_validation_fails_closed(tmp_path, kind):
    paths = artifacts(tmp_path)
    manifest = json.loads(paths[0].read_text())
    vectors_path = tmp_path / manifest["artifacts"]["embeddings"]["path"]
    if kind in {"shape", "non-finite", "zero", "id-order", "dtype"}:
        with np.load(vectors_path, allow_pickle=False) as arrays:
            vectors, ids = arrays["embeddings"].copy(), arrays["protein_ids"].copy()
        if kind == "shape":
            vectors = vectors[:, :2]
        elif kind == "non-finite":
            vectors[0, 0] = np.nan
        elif kind == "zero":
            vectors[0] = 0
        elif kind == "id-order":
            ids = ids[::-1]
        else:
            vectors = vectors.astype(np.float64)
        write_npz(vectors_path, embeddings=vectors, protein_ids=ids)
        manifest["artifacts"]["embeddings"]["sha256"] = file_sha256(vectors_path)
    elif kind == "checksum":
        vectors_path.write_bytes(vectors_path.read_bytes() + b"altered")
    elif kind == "wrong-model":
        manifest["provenance"]["model"] = "esm2_t12_35M_UR50D"
    elif kind == "model-hash":
        manifest["provenance"]["model_sha256"] = "unknown"
    elif kind == "pooling":
        manifest["provenance"]["pooling"] = "mean-with-padding"
    else:
        manifest["provenance"]["actual_device"] = "mps"
    paths[0].write_text(json.dumps(manifest))
    with pytest.raises(EmbeddingError):
        analyze_embedding_dates(dataset(), *paths)


def test_explicit_reference_supported_without_ancestor_inference(tmp_path):
    result = analyze_embedding_dates(dataset(), *artifacts(tmp_path), reference_sample_id="c")
    assert result["reference"]["selection"] == "explicit"
    assert [row["distance"] for row in result["samples"]] == [2, 1, 0]


@pytest.mark.parametrize("panel", [[], ["outside"], ["x", "x"], "x"])
def test_invalid_fixed_panels_fail(tmp_path, panel):
    with pytest.raises(EmbeddingError):
        analyze_embedding_dates(dataset(), *artifacts(tmp_path), panel_loci=panel)


def test_invalid_and_future_dates_retain_descriptive_sample_distances(tmp_path):
    prepared = dataset({"a": "invalid", "b": "9999-01-01", "c": ""})
    result = analyze_embedding_dates(prepared, *artifacts(tmp_path))
    assert [row["date_status"] for row in result["samples"]] == ["invalid", "future", "missing"]
    assert [row["distance"] for row in result["samples"]] == [0, 1, 2]
    assert result["points"] == []
    assert len(result["exclusions"]) == 3


def test_missing_matrix_row_is_retained_and_never_imputed(tmp_path):
    prepared = dataset()
    matrix = AlleleMatrix.from_profiles(
        prepared.profiles[0].catalogue, ["a", "b"], {"a": {"x": "1"}, "b": {"x": "2"}}
    )
    paths = artifacts(tmp_path, [("a", "x", "px"), ("b", "x", "py")])
    result = analyze_embedding_dates(replace(prepared, profiles=(matrix,)), *paths)
    assert result["samples"][2]["distance"] is None
    assert result["samples"][2]["reason"] == "missing_profile"
    assert result["counts"]["total_samples"] == 3


def test_manifest_record_aliases_must_remain_unique(tmp_path):
    paths = artifacts(tmp_path)
    mutate_manifest(paths[0], lambda manifest: manifest["mapping"].append(manifest["mapping"][0]))
    with pytest.raises(EmbeddingError, match="repeated record"):
        analyze_embedding_dates(dataset(), *paths)


def test_mapping_column_contract_and_blank_cells_fail(tmp_path):
    paths = artifacts(tmp_path)
    paths[1].write_text("sample,locus,record_id\na,x,a.x\n")
    with pytest.raises(EmbeddingError, match="requires exactly"):
        analyze_embedding_dates(dataset(), *paths)
    paths[1].write_text("sample_id,locus,record_id\na,x,\n")
    with pytest.raises(EmbeddingError, match="missing or malformed"):
        analyze_embedding_dates(dataset(), *paths)


def test_mixed_species_in_one_scheme_require_separate_eligible_cohorts(tmp_path):
    prepared = dataset()
    prepared = replace(
        prepared,
        samples=tuple(
            {
                **row,
                "species": "Klebsiella pneumoniae"
                if row["sample_id"] != "c"
                else "Escherichia coli",
            }
            for row in prepared.samples
        ),
    )
    with pytest.raises(EmbeddingError, match="one species group"):
        analyze_embedding_dates(prepared, *artifacts(tmp_path))


def test_species_normalization_matches_baseline_and_unmapped_species_is_audited(tmp_path):
    prepared = dataset()
    prepared = replace(
        prepared,
        samples=tuple(
            {
                **row,
                "species": {
                    "a": "Klebsiella pneumoniae",
                    "b": "klebsiella_pneumoniae",
                    "c": "Escherichia coli",
                }[row["sample_id"]],
            }
            for row in prepared.samples
        ),
    )
    paths = artifacts(tmp_path, [("a", "x", "px"), ("b", "x", "py")])
    result = analyze_embedding_dates(prepared, *paths)
    assert result["counts"]["eligible_samples"] == 2
    assert result["samples"][2]["reason"] == "no_mapped_proteins"
