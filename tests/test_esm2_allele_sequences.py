"""Real DNA is required; missing cells and synonymous identities stay explicit."""

import csv
from dataclasses import replace
import hashlib
import json
import subprocess
import sys

import pytest

from chronoclade.datasets import AlleleMatrix, LocusCatalogue, PreparedDataset, samples_from_records
from chronoclade.esm2.allele_sequences import (
    CATALOGUE_SCHEMA,
    SequenceMappingError,
    prepare_allele_sequences,
)
from chronoclade.esm2.proteins import file_sha256, read_proteins


@pytest.fixture
def dataset():
    catalogue = LocusCatalogue(
        "test:cgmlst",
        "v1",
        ("a", "b"),
        "frozen-schema",
        database_version="db1",
        database_sha256="a" * 64,
    )
    samples = samples_from_records(
        [
            {"sample_id": ident, "species": "Klebsiella pneumoniae", "origin": "local"}
            for ident in ("first", "synonym", "missing", "without_profile")
        ]
    )
    matrix = AlleleMatrix.from_profiles(
        catalogue,
        ("first", "synonym", "missing"),
        {
            "first": {"a": "001", "b": "1"},
            "synonym": {"a": "002", "b": "1"},
            "missing": {"a": "novel-sha1-identity"},
        },
    )
    return PreparedDataset(samples, (matrix,))


def frozen(tmp_path, dataset, rows=None, *, scope_change=None):
    csv_path = tmp_path / "alleles.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("locus", "allele", "dna"))
        writer.writerows(
            rows
            if rows is not None
            else [
                ("a", "001", "ATGAAATAA"),
                ("a", "002", "ATGAAGTAA"),  # Explicit synthetic synonymous codon control.
                ("b", "1", "GTGGCTTAG"),  # Table-11 alternative initiator -> methionine.
            ]
        )
    cat = dataset.profiles[0].catalogue
    scope = {
        "species": "Klebsiella pneumoniae",
        "scheme_id": cat.scheme_id,
        "scheme_version": cat.scheme_version,
        "database_version": cat.database_version,
        "database_sha256": cat.database_sha256,
    }
    scope.update(scope_change or {})
    manifest = {
        "schema": CATALOGUE_SCHEMA,
        "schema_version": 1,
        "status": "complete",
        "scope": scope,
        "artifact": {"path": csv_path.name, "sha256": file_sha256(csv_path)},
    }
    path = tmp_path / "catalogue.json"
    path.write_text(json.dumps(manifest))
    return path


def test_synonymous_alleles_deduplicate_proteins_without_losing_identity(tmp_path, dataset):
    result = prepare_allele_sequences(dataset, frozen(tmp_path, dataset), tmp_path / "out")
    inputs = read_proteins(result.fasta_path)
    assert len(inputs.proteins) == 2
    assert inputs.records == 3
    alleles = {(row["locus"], row["allele"]): row for row in result.manifest["alleles"]}
    first, synonym = alleles["a", "001"], alleles["a", "002"]
    assert first["record_id"] != synonym["record_id"]
    assert first["dna_sha256"] != synonym["dna_sha256"]
    assert first["protein_id"] == synonym["protein_id"] == hashlib.sha256(b"MK").hexdigest()
    assert alleles["b", "1"]["initiation_to_methionine"] is True
    assert alleles["b", "1"]["terminal_stop_removed"] is True
    assert {p.sequence for p in inputs.proteins} == {"MK", "MA"}
    with result.mapping_path.open(newline="") as stream:
        mapping = list(csv.DictReader(stream))
    assert set(mapping[0]) == {"sample_id", "locus", "record_id"}
    assert len(mapping) == 4
    record_proteins = {r["record_id"]: r["protein_id"] for r in inputs.mapping}
    assert all(row["record_id"] in record_proteins for row in mapping)
    assert result.manifest["status"] == "complete"
    assert result.manifest["coverage"] == "partial"
    for key, path in (("proteins", result.fasta_path), ("sample_loci", result.mapping_path)):
        assert result.manifest["artifacts"][key]["sha256"] == file_sha256(path)
    assert json.loads(result.manifest_path.read_text()) == result.manifest


def test_all_samples_and_missing_masks_are_audited(tmp_path, dataset):
    result = prepare_allele_sequences(dataset, frozen(tmp_path, dataset), tmp_path / "out")
    assert [row["sample_id"] for row in result.manifest["samples"]] == list(dataset.sample_ids)
    cells = {(r["sample_id"], r["locus"]): r for r in result.manifest["cells"]}
    assert len(cells) == 8
    assert cells["missing", "a"]["allele"] == "novel-sha1-identity"
    assert cells["missing", "a"]["status"] == "missing_sequence"
    assert cells["missing", "b"]["status"] == "missing_call"
    assert cells["without_profile", "a"]["status"] == "missing_profile"
    assert all(
        "protein_id" not in cells[key]
        for key in (("missing", "a"), ("missing", "b"), ("without_profile", "a"))
    )


@pytest.mark.parametrize(
    "dna,reason",
    [
        ("ATGAAATA", "out_of_frame"),
        ("ATGNNNTAA", "ambiguous_or_unsupported_dna"),
        ("ATG---TAA", "ambiguous_or_unsupported_dna"),
        ("ATGAA?TAA", "ambiguous_or_unsupported_dna"),
        ("ATGTGATAA", "internal_stop"),
        ("GCTAAATAA", "invalid_start_codon"),
        ("ATGAAA", "missing_terminal_stop"),
        ("ATGTAA", "protein_exceeds_max_length"),
    ],
)
def test_invalid_cds_is_excluded_with_precise_reason(tmp_path, dataset, dna, reason):
    rows = [("a", "001", dna), ("a", "002", "ATGAAGTAA")]
    kwargs = {"max_length": 1} if reason == "protein_exceeds_max_length" else {}
    # A two-residue sequence tests overlength without changing frame or stops.
    if reason == "protein_exceeds_max_length":
        rows[0] = ("a", "001", "ATGAAATAA")
    result = prepare_allele_sequences(
        dataset, frozen(tmp_path, dataset, rows), tmp_path / "out", **kwargs
    )
    first = next(r for r in result.manifest["alleles"] if r["allele"] == "001")
    assert first["status"] == "invalid_cds"
    assert first["reason"] == reason
    assert "protein_id" not in first
    assert first["record_id"] not in result.fasta_path.read_text()


def test_stopless_input_requires_explicit_policy(tmp_path, dataset):
    manifest = frozen(tmp_path, dataset, [("a", "001", "gtgaaa")])
    result = prepare_allele_sequences(
        dataset, manifest, tmp_path / "out", terminal_stop="allow_absent"
    )
    record = result.manifest["alleles"][0]
    assert record["status"] == "mapped"
    assert record["terminal_stop_removed"] is False
    assert record["initiation_to_methionine"] is True
    assert read_proteins(result.fasta_path).proteins[0].sequence == "MK"


@pytest.mark.parametrize("options", [{"require_complete": True}, {"invalid_policy": "error"}])
def test_failed_strict_run_saves_audit(tmp_path, dataset, options):
    manifest = frozen(tmp_path, dataset, [("a", "001", "ATGNNNTAA")])
    with pytest.raises(SequenceMappingError, match="saved audit") as failure:
        prepare_allele_sequences(dataset, manifest, tmp_path / "out", **options)
    audit = json.loads(failure.value.manifest_path.read_text())
    assert audit["status"] == "failed"
    assert len(audit["samples"]) == 4
    assert audit["counts"]["invalid_alleles"] == 1


@pytest.mark.parametrize(
    "change,match",
    [
        ({"species": "Escherichia coli"}, "species"),
        ({"scheme_id": "other"}, "scheme_id"),
        ({"scheme_version": "v2"}, "scheme_version"),
        ({"database_version": "db2"}, "database_version"),
        ({"database_sha256": "b" * 64}, "database_sha256"),
        ({"database_sha256": None}, "known database SHA256"),
    ],
)
def test_incompatible_catalogue_fails_before_output(tmp_path, dataset, change, match):
    manifest = frozen(tmp_path, dataset, scope_change=change)
    with pytest.raises(SequenceMappingError, match=match):
        prepare_allele_sequences(dataset, manifest, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_unknown_profile_fingerprint_is_not_accepted(tmp_path, dataset):
    matrix = replace(
        dataset.profiles[0], catalogue=replace(dataset.profiles[0].catalogue, database_sha256=None)
    )
    unknown = replace(dataset, profiles=(matrix,))
    manifest = frozen(tmp_path, dataset)
    with pytest.raises(SequenceMappingError, match="database_sha256"):
        prepare_allele_sequences(unknown, manifest, tmp_path / "out")


@pytest.mark.parametrize("species", [None, "Escherichia coli"])
def test_missing_or_mixed_sample_species_fails(tmp_path, dataset, species):
    manifest = frozen(tmp_path, dataset)
    samples = tuple(
        {**row, "species": species} if i == 3 else row for i, row in enumerate(dataset.samples)
    )
    with pytest.raises(SequenceMappingError, match="species"):
        prepare_allele_sequences(replace(dataset, samples=samples), manifest, tmp_path / "out")


@pytest.mark.parametrize(
    "rows,match",
    [
        ([("a", "001", "ATGAAATAA"), ("a", "001", "ATGAAGTAA")], "repeats"),
        ([("foreign", "1", "ATGAAATAA")], "outside"),
        ([("a", "0", "ATGAAATAA")], "missing marker"),
    ],
)
def test_catalogue_identity_conflicts_are_errors(tmp_path, dataset, rows, match):
    with pytest.raises(SequenceMappingError, match=match):
        prepare_allele_sequences(dataset, frozen(tmp_path, dataset, rows), tmp_path / "out")


def test_frozen_checksum_is_verified(tmp_path, dataset):
    manifest = frozen(tmp_path, dataset)
    (tmp_path / "alleles.csv").write_text("locus,allele,dna\na,001,ATGAAATAA\n")
    with pytest.raises(SequenceMappingError, match="SHA256 mismatch"):
        prepare_allele_sequences(dataset, manifest, tmp_path / "out")


def test_output_directory_cannot_overwrite_an_existing_bundle(tmp_path, dataset):
    manifest = frozen(tmp_path, dataset)
    prepare_allele_sequences(dataset, manifest, tmp_path / "out")
    old = (tmp_path / "out/sequence_mapping.json").read_bytes()
    with pytest.raises(SequenceMappingError, match="new or empty"):
        prepare_allele_sequences(dataset, manifest, tmp_path / "out")
    assert (tmp_path / "out/sequence_mapping.json").read_bytes() == old


def test_module_import_does_not_import_torch():
    process = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; import chronoclade.esm2.allele_sequences; "
            "assert 'torch' not in sys.modules; assert 'esm' not in sys.modules",
        ],
        capture_output=True,
        text=True,
    )
    assert process.returncode == 0, process.stderr


def test_complete_mapping_is_consumable_by_existing_sample_distance_reader(tmp_path, dataset):
    from chronoclade.esm2.sample_distances import _read_mapping

    complete = replace(dataset, samples=dataset.samples[:2])
    matrix = AlleleMatrix.from_profiles(
        dataset.profiles[0].catalogue,
        complete.sample_ids,
        {ident: dataset.profiles[0].profile(ident) for ident in complete.sample_ids},
    )
    complete = replace(complete, profiles=(matrix,))
    result = prepare_allele_sequences(
        complete,
        frozen(tmp_path, complete),
        tmp_path / "out",
        require_complete=True,
        invalid_policy="error",
    )
    inputs = read_proteins(result.fasta_path)
    record_proteins = {row["record_id"]: row["protein_id"] for row in inputs.mapping}
    cells, profiles, count = _read_mapping(result.mapping_path, complete, matrix, record_proteins)
    assert count == 4
    assert profiles["first"]["a"] == "001"
    assert cells["first", "a"] == cells["synonym", "a"]
    assert result.manifest["coverage"] == "complete"
    assert all(sample["complete"] for sample in result.manifest["samples"])


def test_leading_zeroes_are_not_coerced_to_ordinary_allele_numbers(tmp_path, dataset):
    manifest = frozen(tmp_path, dataset, [("a", "1", "ATGAAATAA")])
    result = prepare_allele_sequences(dataset, manifest, tmp_path / "out")
    allele = next(row for row in result.manifest["alleles"] if row["allele"] == "001")
    assert allele["status"] == "missing_sequence"
    assert result.manifest["counts"]["mapped_cells"] == 0
    assert result.fasta_path.read_bytes() == b""
