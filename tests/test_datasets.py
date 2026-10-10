"""Prepared bundles must preserve evidence and fail closed on incomplete/corrupt data."""

import csv
from dataclasses import replace
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from chronoclade.datasets import (
    AlleleMatrix,
    DatasetError,
    LocusCatalogue,
    PreparedDataset,
    from_profile_records,
    load_dataset,
    samples_from_records,
    write_dataset,
)
from chronoclade.datasets import storage


@pytest.fixture
def catalogue():
    return LocusCatalogue(
        "example:cgmlst",
        "v1",
        ("a", "b", "never_called"),
        "official-schema.tsv",
        database_version="db1",
    )


@pytest.fixture
def records():
    return [
        {
            "sample_id": "query",
            "origin": "local",
            "species": "Klebsiella pneumoniae",
            "mlst_st": "39",
            "collection_date": "2024-02",
            "country": "Greece",
            "cgmlst_scheme": "example:cgmlst",
            "cgmlst_scheme_version": "v1",
            "cgmlst_profile": {"a": "001"},
            "cgmlst_novel_alleles": {"b": "A" * 40},
            "cgmlst_status": "resolved",
            "cglin_scheme": "example:cglin",
            "cglin_scheme_version": "lin-v1",
            "cglin_raw": "1.2.3.4.5.6.7",
            "cglin_status": "resolved",
            "cglin_source": "frozen_export",
        },
        {
            "sample_id": "PW_source",
            "source_genome_id": "source",
            "origin": "context",
            "species": "Klebsiella pneumoniae",
            "runAccession": "ERR123",
            "biosample": "SAMN123",
            "collection_date": "invalid-date",
            "cgmlst_scheme": "example:cgmlst",
            "cgmlst_scheme_version": "v1",
            "cgmlst_profile": {"a": "0"},
            "cgmlst_novel_alleles": {},
            "hiercc_scheme": "example:hiercc",
            "hiercc_scheme_version": "hier-v1",
            "hiercc_codes": {"HC0": "1", "HC10": "20"},
            "hiercc_status": "resolved",
            "identity_conflict": True,
            "identity_candidates": ["SAMN123", "SAMN456"],
        },
    ]


def _seal_manifest(path, manifest):
    identity = {
        key: value for key, value in manifest.items() if key not in {"dataset_id", "created_at"}
    }
    text = json.dumps(
        identity, ensure_ascii=False, sort_keys=True, allow_nan=False, separators=(",", ":")
    )
    manifest["dataset_id"] = "sha256:" + hashlib.sha256(text.encode()).hexdigest()
    path.write_text(json.dumps(manifest), encoding="utf-8")


def _reseal_artifact(root, descriptor):
    path = root / descriptor["path"]
    descriptor["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    descriptor["bytes"] = path.stat().st_size


def test_round_trip_preserves_categories_dates_and_evidence(tmp_path, catalogue, records):
    provenance = (
        {
            "sample_id": "query",
            "field": "country",
            "source": "user",
            "value": "Greece",
            "selected": True,
        },
    )
    conflicts = (
        {
            "sample_id": "query",
            "field": "country",
            "source": "provider",
            "value": "France",
            "selected": False,
        },
    )
    dataset = from_profile_records(
        records,
        [catalogue],
        provenance=provenance,
        conflicts=conflicts,
        parameters={"typing_config_ref": "local-config", "minimum_coverage": 0.9},
    )
    dataset = replace(dataset, exclusions=({"sample_id": "PW_source", "reason": "missing_date"},))
    manifest_path = write_dataset(dataset, tmp_path / "prepared")
    restored = load_dataset(manifest_path)
    assert restored.samples == dataset.samples
    assert restored.lineages == dataset.lineages
    assert restored.crosswalk == dataset.crosswalk
    assert restored.provenance == provenance
    assert restored.conflicts == conflicts
    assert restored.exclusions == dataset.exclusions
    assert restored.retrieval == dataset.retrieval
    assert restored.parameters == dataset.parameters
    assert restored.sample_ids == ("query", "PW_source")
    assert restored.samples[1]["label"] == "ERR123"
    assert restored.samples[0]["date_end"] == "2024-02-29"
    assert restored.samples[1]["date_precision"] == "invalid"
    assert restored.samples[1]["date_raw"] == "invalid-date"
    assert restored.crosswalk[1]["resolution_status"] == "ambiguous"
    assert restored.crosswalk[1]["identity_candidates"] == ["SAMN123", "SAMN456"]
    assert restored.lineages[0]["assignment"] == "1.2.3.4.5.6.7"
    assert restored.lineages[1]["assignment"] == {"HC0": "1", "HC10": "20"}
    matrix = restored.profiles[0]
    assert isinstance(matrix.codes, np.memmap)
    assert not matrix.codes.flags.writeable
    assert matrix.codes.dtype == np.uint32
    assert matrix.profile("query") == {"a": "001", "b": "A" * 40, "never_called": None}
    assert matrix.profile("PW_source") == dict.fromkeys(catalogue.loci)
    with (manifest_path.parent / "profiles/0000/profiles.csv").open(newline="") as stream:
        exported = list(csv.DictReader(stream))
    assert exported[0]["never_called"] == "0"
    assert exported[0]["b"] == "A" * 40
    manifest = json.loads(manifest_path.read_text())
    assert manifest["status"] == "complete"
    assert manifest["schema_version"] == 1
    assert manifest["software"]["version"]
    assert "samples" in manifest["tables"]
    assert "profiles" in manifest
    assert not Path(manifest["profiles"][0]["codes"]["path"]).is_absolute()
    assert (
        len(manifest_path.read_bytes())
        < (manifest_path.parent / "profiles/0000/alleles.npy").stat().st_size + 10_000
    )


def test_dataset_identity_is_independent_of_destination(tmp_path, catalogue, records):
    dataset = from_profile_records(records, [catalogue])
    first = write_dataset(dataset, tmp_path / "one")
    second = write_dataset(dataset, tmp_path / "two")
    assert load_dataset(first).dataset_id == load_dataset(second).dataset_id


def test_observed_catalogue_remains_explicitly_incomplete(tmp_path, catalogue, records):
    incomplete = replace(catalogue, complete=False, source="observed-provider-export.csv")
    dataset = from_profile_records(records, [incomplete])
    manifest = write_dataset(dataset, tmp_path / "prepared")
    assert load_dataset(manifest).profiles[0].catalogue.complete is False
    assert load_dataset(manifest).profiles[0].catalogue.loci == catalogue.loci


def test_nulls_are_not_allele_matches(catalogue):
    matrix = AlleleMatrix.from_profiles(catalogue, ["a", "b"], {"a": {"a": 0}, "b": {"a": None}})
    shared = [
        (x, y)
        for x, y in zip(matrix.profile("a").values(), matrix.profile("b").values())
        if x is not None and y is not None
    ]
    assert shared == []
    assert matrix.categories == ()


@pytest.mark.parametrize("missing", [None, 0, "0", "", "-", "?", "null", "None", "unknown"])
def test_export_missing_markers_are_internal_nulls(catalogue, missing):
    matrix = AlleleMatrix.from_profiles(catalogue, ["sample"], {"sample": {"a": missing}})
    assert matrix.profile("sample")["a"] is None
    assert matrix.categories == ()


def test_caller_evidence_can_explicitly_replace_automatic_tables(catalogue, records):
    dataset = from_profile_records(
        records, [catalogue], lineages=(), crosswalk=(), provenance=(), retrieval=()
    )
    assert dataset.lineages == dataset.crosswalk == dataset.provenance == dataset.retrieval == ()


@pytest.mark.parametrize(
    "change,match",
    [
        ({"cgmlst_scheme_version": "other"}, "compatible"),
        ({"cgmlst_scheme_version": None}, "compatible"),
        ({"cgmlst_database_version": "other"}, "database version"),
        ({"cgmlst_profile": {"outside_catalogue": "1"}}, "unsupported loci"),
        ({"cgmlst_profile": {"b": "1"}}, "overlap"),
        ({"cgmlst_ambiguous_loci": ["a"]}, "ambiguous locus"),
    ],
)
def test_adapter_rejects_incompatible_or_conflicting_calls(catalogue, records, change, match):
    records[0].update(change)
    with pytest.raises(DatasetError, match=match):
        from_profile_records(records, [catalogue])


def test_failed_profile_is_null_with_retained_status(catalogue, records):
    records[0]["cgmlst_status"] = "failed"
    dataset = from_profile_records(records, [catalogue])
    assert dataset.profiles[0].profile("query") == dict.fromkeys(catalogue.loci)
    assert dataset.retrieval[0]["status"] == "failed"


def test_different_database_hashes_require_distinct_explicit_catalogues(
    catalogue, records, tmp_path
):
    records[0]["cgmlst_database_sha256"] = "a" * 64
    records[1]["cgmlst_database_sha256"] = "b" * 64
    first = replace(catalogue, database_sha256="a" * 64)
    second = replace(catalogue, database_sha256="b" * 64)
    with pytest.raises(DatasetError, match="fingerprint"):
        from_profile_records(records, [catalogue])
    with pytest.raises(DatasetError, match="fingerprint"):
        from_profile_records(records, [first])
    dataset = from_profile_records(records, [first, second])
    assert len(dataset.profiles) == 2
    assert dataset.profiles[0].sample_ids == ("query",)
    assert dataset.profiles[1].sample_ids == ("PW_source",)
    restored = load_dataset(write_dataset(dataset, tmp_path / "prepared"))
    assert restored.profiles[0].catalogue.database_sha256 == "a" * 64
    assert restored.profiles[1].catalogue.database_sha256 == "b" * 64
    assert restored.retrieval[1]["evidence"]["cgmlst_database_sha256"] == "b" * 64


def test_known_and_unknown_database_hashes_cannot_be_merged(catalogue, records):
    records[0]["cgmlst_database_sha256"] = "a" * 64
    known = replace(catalogue, database_sha256="a" * 64)
    with pytest.raises(DatasetError, match="fingerprint"):
        from_profile_records(records, [catalogue])
    dataset = from_profile_records(records, [known, catalogue])
    assert len(dataset.profiles) == 2


def test_database_fingerprint_can_define_unversioned_profile_scope(catalogue, records):
    records[0]["cgmlst_scheme_version"] = None
    records[0]["cgmlst_database_sha256"] = "a" * 64
    frozen = replace(catalogue, scheme_version="sha256:" + "a" * 64, database_sha256="a" * 64)
    dataset = from_profile_records(records, [frozen, catalogue])
    assert dataset.profiles[0].catalogue.scheme_version == "sha256:" + "a" * 64


@pytest.mark.parametrize(
    "change,match",
    [
        ({"sample_id": "unsafe/id"}, "safe"),
        (
            {"date_start": "2024-02-29", "date_end": "2024-02-01", "date_precision": "interval"},
            "inconsistent",
        ),
        (
            {"date_start": "2024-02-01", "date_end": "2024-02-15", "date_precision": "month"},
            "full month",
        ),
        ({"date_precision": "guessed"}, "precision"),
        ({"role": "other"}, "role"),
    ],
)
def test_sample_validation(catalogue, records, change, match):
    samples = list(samples_from_records(records))
    samples[0].update(change)
    with pytest.raises(DatasetError, match=match):
        PreparedDataset(tuple(samples)).validate()


def test_missing_dates_and_locations_are_not_invented():
    records = [{"sample_id": "local", "origin": "local"}]
    dataset = from_profile_records(records, [])
    sample = dataset.samples[0]
    assert sample["date_precision"] == "missing"
    assert sample["date_start"] is sample["date_end"] is sample["country"] is sample["host"] is None


def test_unversioned_lineage_is_retained_without_inventing_database_version(catalogue, records):
    del records[0]["cglin_scheme_version"]
    lineage = from_profile_records(records, [catalogue]).lineages[0]
    assert lineage["scheme_version"] is None
    assert lineage["database_version"] is None
    assert lineage["assignment"] == records[0]["cglin_raw"]


def test_export_fingerprint_is_not_labelled_as_database_version(catalogue, records):
    del records[0]["cglin_scheme_version"]
    records[0]["cglin_frozen_export_sha256"] = "b" * 64
    lineage = from_profile_records(records, [catalogue]).lineages[0]
    assert lineage["scheme_version"] == "exportsha256:" + "b" * 64
    assert lineage["database_version"] is None


@pytest.mark.parametrize("key", ["api_key", "authorization", "access_token", "password"])
def test_credentials_are_not_serialised(tmp_path, catalogue, records, key):
    dataset = from_profile_records(records, [catalogue])
    unsafe = replace(dataset, parameters={"provider": {key: "do-not-write"}})
    with pytest.raises(DatasetError, match="credentials"):
        write_dataset(unsafe, tmp_path / "prepared")
    assert not (tmp_path / "prepared").exists()


def test_existing_bundle_is_not_replaced(tmp_path, catalogue, records):
    dataset = from_profile_records(records, [catalogue])
    manifest = write_dataset(dataset, tmp_path / "prepared")
    original = manifest.read_bytes()
    with pytest.raises(DatasetError, match="already exists"):
        write_dataset(dataset, manifest.parent)
    assert manifest.read_bytes() == original


def test_failed_write_never_publishes_completed_manifest(tmp_path, catalogue, records, monkeypatch):
    dataset = from_profile_records(records, [catalogue])

    def fail(path, value):
        if path.name == "dataset.json":
            raise OSError("simulated interruption before manifest")
        path.write_text(json.dumps(value))

    monkeypatch.setattr(storage, "_write_json", fail)
    with pytest.raises(OSError, match="interruption"):
        write_dataset(dataset, tmp_path / "prepared")
    assert not (tmp_path / "prepared").exists()
    assert list(tmp_path.glob(".prepared.pending-*")) == []


@pytest.mark.parametrize(
    "artifact", ["samples.npy", "profiles/0000/alleles.npy", "lineages.jsonl", "samples.csv"]
)
def test_corrupt_artifacts_are_rejected(tmp_path, catalogue, records, artifact):
    manifest = write_dataset(from_profile_records(records, [catalogue]), tmp_path / "prepared")
    with (manifest.parent / artifact).open("ab") as stream:
        stream.write(b"corrupted")
    with pytest.raises(DatasetError, match="checksum/size"):
        load_dataset(manifest)


@pytest.mark.parametrize(
    "change,match",
    [
        ({"status": "running"}, "not complete"),
        ({"schema_version": 999}, "schema version"),
        ({"sample_ids": ["wrong"]}, "sample IDs"),
    ],
)
def test_invalid_manifest_is_rejected(tmp_path, catalogue, records, change, match):
    path = write_dataset(from_profile_records(records, [catalogue]), tmp_path / "prepared")
    manifest = json.loads(path.read_text())
    manifest.update(change)
    _seal_manifest(path, manifest)
    with pytest.raises(DatasetError, match=match):
        load_dataset(path)


@pytest.mark.parametrize(
    "unsafe", ["../outside.npy", "/tmp/outside.npy", "profiles/../samples.npy", "a\\b.npy"]
)
def test_artifact_path_traversal_is_rejected(tmp_path, catalogue, records, unsafe):
    path = write_dataset(from_profile_records(records, [catalogue]), tmp_path / "prepared")
    manifest = json.loads(path.read_text())
    manifest["tables"]["samples"]["data"]["path"] = unsafe
    _seal_manifest(path, manifest)
    with pytest.raises(DatasetError, match="path|relative"):
        load_dataset(path)


def test_symlink_outside_bundle_is_rejected(tmp_path, catalogue, records):
    path = write_dataset(from_profile_records(records, [catalogue]), tmp_path / "prepared")
    samples = path.parent / "samples.npy"
    outside = tmp_path / "outside.npy"
    samples.rename(outside)
    samples.symlink_to(outside)
    with pytest.raises(DatasetError, match="outside"):
        load_dataset(path)


def test_resealed_invalid_code_array_still_fails_validation(tmp_path, catalogue, records):
    path = write_dataset(from_profile_records(records, [catalogue]), tmp_path / "prepared")
    manifest = json.loads(path.read_text())
    descriptor = manifest["profiles"][0]["codes"]
    codes = np.load(path.parent / descriptor["path"])
    codes[0, 0] = 999
    np.save(path.parent / descriptor["path"], codes, allow_pickle=False)
    _reseal_artifact(path.parent, descriptor)
    _seal_manifest(path, manifest)
    with pytest.raises(DatasetError, match="unknown category"):
        load_dataset(path)


def test_manifest_tampering_without_resealing_is_rejected(tmp_path, catalogue, records):
    path = write_dataset(from_profile_records(records, [catalogue]), tmp_path / "prepared")
    manifest = json.loads(path.read_text())
    manifest["parameters"]["changed"] = True
    path.write_text(json.dumps(manifest))
    with pytest.raises(DatasetError, match="identity"):
        load_dataset(path)


def test_reader_cannot_open_writable_memory_maps(tmp_path, catalogue, records):
    path = write_dataset(from_profile_records(records, [catalogue]), tmp_path / "prepared")
    with pytest.raises(DatasetError, match="read-only"):
        load_dataset(path, mmap_mode="r+")
    assert not isinstance(load_dataset(path, mmap_mode=None).profiles[0].codes, np.memmap)
