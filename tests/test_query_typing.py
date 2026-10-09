"""Native Pathogenwatch adapters preserve unresolved alleles and lineage levels."""

import json

import pytest

from chronoclade.query_typing import (
    QueryTypingError,
    _assignment_result,
    _database_hash,
    load_query_typing,
    parse_cgmlst_result,
)


def profile(**updates):
    return {
        "scheme": "ecoli_1",
        "genes": ["a", "b", "c", "d"],
        "code": "1__" + "a" * 40 + "_3",
        "schemeSize": 4,
        "alleles": {"d": [{"id": 3}, {"id": 4}]},
        **updates,
    }


def test_known_novel_missing_ambiguous_are_separate():
    row = parse_cgmlst_result(profile(), "ecoli_1", "2026-01-01")
    assert row["cgmlst_profile"] == {"a": "1"}
    assert row["cgmlst_novel_alleles"] == {"c": "a" * 40}
    assert row["cgmlst_missing_loci"] == ["b"]
    assert row["cgmlst_ambiguous_loci"] == ["d"]
    assert row["cgmlst_loci"] == ["a", "b", "c", "d"]
    assert row["cgmlst_called_fraction"] == 0.5
    assert row["cgmlst_known_fraction"] == 0.25


def test_modern_hit_object_is_not_mistaken_for_duplicate():
    row = parse_cgmlst_result(
        profile(alleles={"a": {"id": 1, "contig": "c", "start": 1}}), "ecoli_1", "2026-01-01"
    )
    assert row["cgmlst_profile"] == {"a": "1", "d": "3"}


@pytest.mark.parametrize(
    "updates",
    [
        {"scheme": "other"},
        {"genes": ["a", "a", "c", "d"]},
        {"code": "1_2"},
        {"code": "1__nonsense_3"},
        {"schemeSize": 5},
    ],
)
def test_invalid_profile_rejected(updates):
    with pytest.raises(QueryTypingError):
        parse_cgmlst_result(profile(**updates), "ecoli_1", "2026-01-01")


def test_plincer_disagreeing_best_matches_keep_only_shared_prefix():
    config = {"kind": "plincer", "scheme": "scgMLST629_S", "scheme_version": "2026-01-01"}
    row = _assignment_result(
        {
            "LINcode": [1, 2, 3, 4, 5, "*", "*", "*", "*", "*"],
            "cgST": "*new",
            "matches": [{"LINcode": [1, 2, 3, 9, 5, 6, 7, 8, 9, 10]}],
        },
        config,
    )
    assert row["cglin_raw"] == "1_2_3_*_*_*_*_*_*_*"
    assert row["cglin_resolved_depth"] == 3
    assert row["cglin_provisional"]
    assert row["cglin_group_5"] == ""


def test_hclink_only_returns_native_supported_levels():
    row = _assignment_result(
        {"hierCC": [["d200", "42"], ["d1100", "9"]], "hierccDistance": 122},
        {"kind": "hclink", "scheme": "Enterobase_ecoli", "scheme_version": "2026-01-01"},
    )
    assert row["hiercc_codes"] == {"HC200": "42", "HC1100": "9"}
    assert "HC0" not in row["hiercc_codes"]


def test_no_hiercc_match_stays_unassigned():
    row = _assignment_result(
        {"hierCC": []},
        {"kind": "hclink", "scheme": "Enterobase_ecoli", "scheme_version": "2026-01-01"},
    )
    assert row["hiercc_status"] == "unassigned"


def test_database_hash_depends_on_file_names_and_bytes(tmp_path):
    (tmp_path / "one").write_text("one")
    original = _database_hash(tmp_path)
    (tmp_path / "one").rename(tmp_path / "two")
    assert _database_hash(tmp_path) != original


@pytest.mark.parametrize(
    "assignments",
    [
        [{"sample_id": "../unsafe", "assembly_sha256": "a" * 64}],
        [{"sample_id": "a", "assembly_sha256": "unknown"}],
        [{"sample_id": "a", "assembly_sha256": "a" * 64}] * 2,
    ],
)
def test_import_requires_safe_unique_ids_and_assembly_identity(tmp_path, assignments):
    path = tmp_path / "typing.json"
    path.write_text(json.dumps({"schema_version": 1, "assignments": assignments}))
    with pytest.raises(QueryTypingError):
        load_query_typing(path)


def test_hclink_blank_levels_are_not_copied_from_neighbour():
    row = _assignment_result(
        {"hierCC": [["HC0", ""], ["HC200", "42"]]},
        {"kind": "hclink", "scheme": "Enterobase_ecoli", "scheme_version": "2026-01-01"},
    )
    assert row["hiercc_codes"] == {"HC200": "42"}


@pytest.mark.parametrize(
    "update",
    [
        {"cgmlst_profile": {"a": "novel"}},
        {"cgmlst_loci": ["a", "a"]},
        {"cgmlst_profile": {"other": "1"}},
        {"hiercc_codes": {"HC200": ""}},
        {"cglin_raw": "1_?_3"},
    ],
)
def test_import_rejects_unsafe_assignments(tmp_path, update):
    path = tmp_path / "typing.json"
    row = {
        "sample_id": "query",
        "species": "Escherichia coli",
        "assembly_sha256": "a" * 64,
        "cgmlst_scheme": "ecoli_1",
        "cgmlst_scheme_version": "2026-01-01",
        "cgmlst_profile": {"a": "1"},
        "cgmlst_loci": ["a"],
        **update,
    }
    path.write_text(json.dumps({"schema_version": 1, "assignments": [row]}))
    with pytest.raises(QueryTypingError):
        load_query_typing(path)


def test_hclink_database_date_is_not_allowed_to_precede_caller(tmp_path):
    import sys
    from chronoclade.query_typing import load_typing_config

    (tmp_path / "loci.json").write_text(json.dumps({"genes": ["a", "b"]}))
    index = tmp_path / "index"
    index.mkdir()
    db = tmp_path / "db"
    db.mkdir()
    (db / "metadata.json").write_text(json.dumps({"datestamp": "2026-01-01 12:00:00"}))
    config = {
        "organisms": {
            "Escherichia coli": {
                "cgmlst": {
                    "command": [sys.executable],
                    "scheme": "ecoli_1",
                    "scheme_version": "2026-01-02",
                    "tool_version": "8.0.0",
                    "index_dir": "index",
                },
                "assignment": {
                    "command": [sys.executable],
                    "kind": "hclink",
                    "scheme": "Enterobase_ecoli",
                    "scheme_version": "2026-01-01",
                    "tool_version": "4.0.1",
                    "reference_db": "db",
                    "loci_file": "loci.json",
                },
            }
        }
    }
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    with pytest.raises(QueryTypingError, match="at least as recent"):
        load_typing_config(path)
    config["organisms"]["Escherichia coli"]["cgmlst"]["scheme_version"] = "2026-01-01"
    path.write_text(json.dumps(config))
    assert (
        load_typing_config(path)["organisms"]["Escherichia coli"]["assignment"]["reference_db"]
        == db
    )


def test_not_ready_config_reports_setup_instructions(tmp_path):
    from chronoclade.query_typing import load_typing_config

    path = tmp_path / "config.json"
    path.write_text(json.dumps({"ready": False, "organisms": {}}))
    with pytest.raises(QueryTypingError, match="setup-typing"):
        load_typing_config(path)


def test_locus_order_must_match_assignment_database():
    from chronoclade.query_typing import _check_locus_order

    _check_locus_order(["a_S", "b_S"], ["a_S", "b_S"])
    with pytest.raises(QueryTypingError, match="locus order"):
        _check_locus_order(["b_S", "a_S"], ["a_S", "b_S"])


def test_species_case_and_underscore_aliases_use_canonical_config_and_cache(tmp_path, monkeypatch):
    import hashlib
    from chronoclade import query_typing
    from chronoclade.metadata import Sample

    assembly = tmp_path / "query.fa"
    assembly.write_text(">query\nACGT\n")
    provenance = {"cgmlst": {}, "assignment": {}}
    monkeypatch.setattr(
        query_typing,
        "load_typing_config",
        lambda path: {
            "organisms": {
                "Klebsiella pneumoniae": {"cgmlst": {}, "assignment": {"kind": "plincer"}}
            }
        },
    )
    monkeypatch.setattr(query_typing, "_provenance", lambda tool, fields: {})
    reused = {
        "sample_id": "query",
        "species": "KLEBSIELLA PNEUMONIAE",
        "assembly_sha256": hashlib.sha256(assembly.read_bytes()).hexdigest(),
        "typing_provenance": provenance,
    }
    sample = Sample("query", assembly, "2026", "UK", "Klebsiella_pneumoniae", "ST1", "query")
    rows = query_typing.type_query_assemblies(
        [sample],
        tmp_path / "config.json",
        tmp_path / "output",
        existing_assignments={"query": reused},
    )
    assert rows[0]["typing_reused"] is True
    assert rows[0]["species"] == "Klebsiella_pneumoniae"
    assert query_typing.load_query_typing(tmp_path / "output/query_typing.json") == rows


def test_import_rejects_independent_unsupported_taxon(tmp_path):
    from chronoclade.query_typing import normalise_typing_species

    with pytest.raises(QueryTypingError, match="Unsupported"):
        normalise_typing_species("Klebsiella oxytoca")
    path = tmp_path / "typing.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "assignments": [
                    {
                        "sample_id": "query",
                        "species": "Salmonella_enterica",
                        "assembly_sha256": "a" * 64,
                    }
                ],
            }
        )
    )
    with pytest.raises(QueryTypingError, match="Unsupported"):
        load_query_typing(path)
