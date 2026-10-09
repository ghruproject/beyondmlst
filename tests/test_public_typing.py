import hashlib
import json
from copy import deepcopy

import pytest

from chronoclade.public_typing import (
    PublicTypingError,
    annotate_public_typing,
    load_public_typing,
    resolve_focal_typing,
)


def row(ident="exact-PW-id", **extras):
    return dict(
        source_genome_id=ident,
        cgmlst_scheme="ecoli_1",
        cgmlst_scheme_version="2026-01-01",
        cgmlst_loci=["a", "b"],
        cgmlst_profile={"a": "1", "b": "2"},
        **extras,
    )


def load(tmp_path, rows, envelope=False):
    path = tmp_path / "typing.json"
    path.write_text(json.dumps({"assignments": rows} if envelope else rows))
    return load_public_typing(path), path


def test_list_and_envelope_preserve_database_provenance(tmp_path):
    raw = row(cgmlst_database_sha256="a" * 64)
    assignments, path = load(tmp_path, [raw], True)
    result = assignments[0]
    assert result["public_typing_export_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert result["cgmlst_database_sha256"] == "a" * 64
    assert result["cgmlst_status"] == "resolved"
    assignments, _ = load(tmp_path, [row()])
    assert "cgmlst_database_sha256" not in assignments[0]


@pytest.mark.parametrize(
    "changes",
    [
        {"source_genome_id": ""},
        {"source_genome_id": " name "},
        {"cgmlst_scheme_version": ""},
        {"cgmlst_loci": ["a", "a"]},
        {"cgmlst_profile": {"c": "1"}},
        {"cgmlst_profile": {"a": "junk"}},
        {"cgmlst_database_sha256": "bad"},
    ],
)
def test_malformed_records_fail(tmp_path, changes):
    raw = row()
    raw.update(changes)
    with pytest.raises(PublicTypingError):
        load(tmp_path, [raw])


def test_missing_and_novel_calls_are_separate(tmp_path):
    raw = row()
    raw["cgmlst_profile"] = {"a": "0", "b": "A" * 40}
    assignments, _ = load(tmp_path, [raw])
    assert assignments[0]["cgmlst_profile"] == {"b": "a" * 40}
    assert assignments[0]["cgmlst_called_fraction"] == 0.5


def test_hiercc_normalises_levels_and_zero_cluster(tmp_path):
    raw = dict(
        source_genome_id="x",
        hiercc_scheme="ecoli_hiercc",
        hiercc_scheme_version="2026",
        hiercc_codes={"d1100": 15, "HC10": 0},
    )
    assignments, _ = load(tmp_path, [raw])
    assert assignments[0]["hiercc_codes"] == {"HC1100": "15", "HC10": "0"}
    raw["hiercc_codes"]["d10"] = 0
    with pytest.raises(PublicTypingError, match="Duplicate"):
        load(tmp_path, [raw])
    raw["hiercc_codes"] = {"HC10": "-"}
    with pytest.raises(PublicTypingError):
        load(tmp_path, [raw])


def test_exact_id_join_never_uses_name_accession_or_order(tmp_path):
    assignments, _ = load(tmp_path, [row("A")])
    catalogue = [
        {"source_genome_id": "B", "name": "A", "biosample": "A"},
        {"source_genome_id": "A"},
    ]
    annotated = annotate_public_typing(catalogue, assignments)
    assert "cgmlst_profile" not in annotated[0]
    assert annotated[1]["cgmlst_profile"] == {"a": "1", "b": "2"}
    assert catalogue[1] == {"source_genome_id": "A"}


def test_deduplicated_alias_conflict_stays_unassigned(tmp_path):
    a, b = row("A"), row("B")
    b["cgmlst_profile"] = {"a": "2", "b": "2"}
    assignments, _ = load(tmp_path, [a, b])
    result = annotate_public_typing(
        [{"source_genome_id": "A", "source_genome_ids": ["A", "B"]}], assignments
    )[0]
    assert result["cgmlst_status"] == "conflict"
    assert result["cgmlst_profile"] == {}
    assert result["public_typing_record_count"] == 2
    assert len(result["cgmlst_conflicting_assignments"]) == 2


def test_repeated_identical_records_do_not_conflict(tmp_path):
    assignments, _ = load(tmp_path, [row("A"), row("A")])
    result = annotate_public_typing([{"source_genome_id": "A"}], assignments)[0]
    assert result["cgmlst_status"] == "resolved"
    assert result["public_typing_record_count"] == 2


def test_existing_different_typing_is_a_conflict(tmp_path):
    assignments, _ = load(tmp_path, [row("A")])
    old = row("A")
    old["cgmlst_profile"] = {"a": "9", "b": "2"}
    result = annotate_public_typing([old], assignments)[0]
    assert result["cgmlst_status"] == "conflict"


def test_focal_join_requires_previously_verified_identity(tmp_path):
    assignments, _ = load(tmp_path, [row("A")])
    catalogue = annotate_public_typing([{"source_genome_id": "A"}], assignments)
    focal = [
        {"sample_id": "A"},
        {"sample_id": "query", "cglin_matched_source_genome_id": "A", "country": "UK"},
    ]
    original = deepcopy(focal)
    result = resolve_focal_typing(focal, catalogue)
    assert "cgmlst_profile" not in result[0]
    assert result[1]["sample_id"] == "query"
    assert result[1]["country"] == "UK"
    assert result[1]["public_typing_matched_source_genome_id"] == "A"
    assert focal == original


def test_optional_cglin_uses_existing_normaliser(tmp_path):
    raw = dict(
        source_genome_id="A",
        cglin_scheme="scgMLST629_S",
        cglin_scheme_version="v1",
        cglin_raw="1_2_3_4_5_6_7_8_9_10",
    )
    assignments, _ = load(tmp_path, [raw])
    result = annotate_public_typing([{"source_genome_id": "A"}], assignments)[0]
    assert result["cglin_status"] == "complete"
    assert result["cglin_group_7"]


@pytest.mark.parametrize("payload", [{}, [], ["bad"], {"rows": []}])
def test_invalid_envelope(tmp_path, payload):
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(payload))
    with pytest.raises(PublicTypingError):
        load_public_typing(path)


def test_native_novel_calls_and_bad_loci_are_validated(tmp_path):
    raw = row(cgmlst_novel_alleles={"a": "A" * 40})
    raw["cgmlst_profile"] = {"b": "1"}
    assignments, _ = load(tmp_path, [raw])
    assert assignments[0]["cgmlst_called_fraction"] == 1
    assert assignments[0]["cgmlst_novel_alleles"] == {"a": "a" * 40}
    raw["cgmlst_novel_alleles"] = {"bad": "a" * 40}
    with pytest.raises(PublicTypingError):
        load(tmp_path, [raw])


def test_existing_focal_typing_conflict_is_not_overwritten(tmp_path):
    assignments, _ = load(tmp_path, [row("A")])
    focal = row("query", cglin_matched_source_genome_id="A", sample_id="Q")
    focal["cgmlst_profile"] = {"a": "9", "b": "2"}
    result = resolve_focal_typing([focal], assignments)[0]
    assert result["sample_id"] == "Q"
    assert result["cgmlst_status"] == "conflict"
