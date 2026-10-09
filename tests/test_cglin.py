import gzip
import hashlib
import json

import pytest

from chronoclade.cglin import (
    CGLINError,
    DEFAULT_SCHEME,
    annotate_catalogue,
    download_cglin_export,
    group_key,
    load_cglin_export,
    normalise_assignment,
    parse_code,
    resolve_focal_assignments,
)


def assignment(ident, code="0,0,197,0,4,0,93,0,0,0", **kwargs):
    return normalise_assignment({"Genome ID": ident, "LIN code": code, **kwargs})


@pytest.mark.parametrize(
    "raw,depth,status",
    [
        ("0,0,197,0,4,0,93,0,0,0", 10, "complete"),
        ("0_0_197_0_4_0_93", 7, "partial"),
        ("0,0,197,0,4,-,-,-,-,-", 5, "partial"),
        ("0,0,197,0,4,*f26e", 5, "partial"),
        ("", 0, "missing"),
        (None, 0, "missing"),
        ("0,0,bad,0,4", 0, "malformed"),
        ("0,0,-,0,4", 0, "malformed"),
        ("0,0,0,0,0,0,0,0,0,0,0", 0, "malformed"),
        ("0,-1,0,0,0", 0, "malformed"),
    ],
)
def test_parse_actual_prefix_only(raw, depth, status):
    components, observed = parse_code(raw)
    assert len(components) == depth
    assert observed == status


def test_full_prefix_version_scoping_and_provisional_status():
    first = assignment("a")
    second = assignment("b", "1,1,198,0,4,0,93,0,0,0")
    newer = assignment("c", version="v2")
    assert first["cglin_group_7"] != second["cglin_group_7"]
    assert first["cglin_group_7"] != newer["cglin_group_7"]
    assert json.loads(first["cglin_group_5"])[2] == [0, 0, 197, 0, 4]
    provisional = assignment("d", "0,0,197,0,4,-,-,-,-,-", cgST="*f26e")
    assert provisional["cglin_status"] == "provisional"
    assert provisional["cglin_code_status"] == "partial"
    assert provisional["cglin_status_5"] == "provisional"
    assert provisional["cglin_status_6"] == "partial"
    assert provisional["cglin_group_6"] == ""
    assert first["cgst"] == ""
    with pytest.raises(CGLINError):
        group_key(DEFAULT_SCHEME, "v1", [0], 0)


def test_actual_export_import_hash_and_duplicate_conflict(tmp_path):
    data = (
        '\ufeff"Genome ID","Genome Name","cgST","Closest cgST","LIN code"\n'
        '"a","name","19550","19550","0,0,197,0,4,0,93,0,0,0"\n'
        '"a","name","19550","19550","0,0,197,0,4,0,93,0,0,0"\n'
    ).encode()
    path = tmp_path / "export.csv"
    path.write_bytes(data)
    exported = load_cglin_export(path, retrieved_at="2026-10-09T12:00:00Z")
    assert exported[0]["cglin_export_sha256"] == hashlib.sha256(data).hexdigest()
    assert exported[0]["cglin_raw"] == "0,0,197,0,4,0,93,0,0,0"
    annotated = annotate_catalogue(
        [{"source_genome_id": "a"}, {"source_genome_id": "missing"}], exported
    )
    assert annotated[0]["cglin_export_record_count"] == 2
    assert annotated[0]["cglin_status_7"] == "resolved"
    assert annotated[1]["cglin_status_5"] == "missing"
    assert annotated[1]["cglin_group_5"] == ""
    conflict = annotate_catalogue(
        [{"source_genome_id": "a"}], exported + [assignment("a", "1,1,197,0,4")]
    )[0]
    assert conflict["cglin_status"] == "conflict"
    assert not conflict["cglin_group_5"]
    assert len(conflict["cglin_conflicting_assignments"]) == 3


def test_deduplicated_sample_checks_assignments_of_all_alias_ids():
    row = {"source_genome_id": "a", "source_genome_ids": ["a", "b"]}
    joined = annotate_catalogue([row], [assignment("b")])[0]
    assert joined["cglin_group_7"]
    conflicted = annotate_catalogue([row], [assignment("a"), assignment("b", "1,0,197,0,4")])[0]
    assert conflicted["cglin_status"] == "conflict"


def test_focal_join_strong_accessions_ambiguity_conflict_no_st_inference():
    public = annotate_catalogue(
        [
            {"source_genome_id": "a", "biosample": "SAMN100", "run_accessions": ["SRR100"]},
            {"source_genome_id": "b", "biosample": "SAMN100", "run_accessions": ["SRR200"]},
            {"source_genome_id": "c", "biosample": "SAMN300"},
        ],
        [assignment("a"), assignment("b"), assignment("c")],
    )
    focal = [
        {"sample_id": "amb", "biosample": "SAMN100"},
        {"sample_id": "match", "run_accessions": ["SRR200"]},
        {"sample_id": "ST147", "lineage": "ST147"},
        {"sample_id": "conf", "biosample": "SAMN300"},
    ]
    rows, audit = resolve_focal_assignments(
        focal, public, [{"sample_id": "conf", "source_genome_id": "a"}]
    )
    assert [row["cglin_join_status"] for row in rows] == [
        "ambiguous",
        "matched",
        "missing",
        "conflict",
    ]
    assert rows[0]["cglin_group_5"] == ""
    assert rows[1]["cglin_matched_source_genome_id"] == "b"
    assert rows[2]["cglin_group_7"] == ""
    assert audit[0]["candidate_source_genome_ids"] == ["a", "b"]


def test_focal_preserves_valid_partial_code_and_reports_comparison_conflict():
    public = annotate_catalogue(
        [{"source_genome_id": "a", "biosample": "SAMN100"}], [assignment("a")]
    )
    focal = {
        "sample_id": "f",
        "biosample": "SAMN100",
        "cglin_raw": "0,0,197,0,4",
        "cglin_scheme": DEFAULT_SCHEME,
        "cglin_scheme_version": "unknown",
        "cglin_export_sha256": "frozen",
    }
    rows, _ = resolve_focal_assignments([focal], public)
    assert rows[0]["cglin_join_status"] == "matched"
    assert rows[0]["cglin_raw"] == focal["cglin_raw"]
    assert rows[0]["cglin_status_6"] == "partial"
    focal["cglin_raw"] = "1,0,197,0,4"
    rows, _ = resolve_focal_assignments([focal], public)
    assert rows[0]["cglin_join_status"] == "conflict"
    assert rows[0]["cglin_raw"] == "1,0,197,0,4"
    rows, _ = resolve_focal_assignments([{"sample_id": "unproven", "cglin_raw": "0,0,197,0,4"}], [])
    assert rows[0]["cglin_status"] == "missing_provenance"
    assert not rows[0]["cglin_group_5"]


def test_export_replay_and_all_missing_statuses(tmp_path):
    path = tmp_path / "rows.json"
    path.write_text(json.dumps([{"source_genome_id": "a", "cglin_raw": "0,0,197,0,4"}]))
    assert load_cglin_export(path) == load_cglin_export(path)
    rows = annotate_catalogue(
        [{"source_genome_id": "bad"}, {"source_genome_id": "none"}], [assignment("bad", "bad")]
    )
    for depth in (5, 6, 7):
        assert [r[f"cglin_status_{depth}"] for r in rows] == ["malformed", "missing"]
        assert all(not r[f"cglin_group_{depth}"] for r in rows)


@pytest.mark.parametrize("content", ["", "name,code\nfoo,0_0\n", '{"wrong":[]}'])
def test_export_requires_source_id_contract(tmp_path, content):
    path = tmp_path / "bad.csv"
    path.write_text(content)
    with pytest.raises(CGLINError):
        load_cglin_export(path)


def test_live_contract_numeric_requests_uuid_join_gzip_partial_batch(tmp_path):
    def transport(numeric_ids, job):
        assert numeric_ids == ["174724", "174727"]
        assert job == "lincodes-3390273-2"
        return gzip.compress(
            b'"Genome ID","LIN code","cgST"\n"uuid-a","0,0,197,0,4,0,93,0,0,0","19550"\n'
        )

    manifest = download_cglin_export(
        [
            {"source_genome_id": "uuid-a", "numeric_source_id": 174724},
            {"source_genome_id": "uuid-b", "numeric_source_id": 174727},
        ],
        tmp_path,
        api_key="secret-fixture",
        transport=transport,
    )
    assert not manifest["complete"]
    assert manifest["missing_source_genome_ids"] == ["uuid-b"]
    assert "secret-fixture" not in (tmp_path / "cglin_export.json").read_text()
    assert (tmp_path / "batch_0000.csv").read_bytes().startswith(b'"Genome ID"')


@pytest.mark.parametrize(
    "data",
    [
        b'"Genome ID","LIN code"\n"other","0,0"\n',
        b'"Genome ID","LIN code"\n"a","0,0"\n"a","0,0"\n',
        b"<html>failure</html>",
    ],
)
def test_download_rejects_unrequested_duplicate_and_corrupt_exports(tmp_path, data):
    with pytest.raises(CGLINError):
        download_cglin_export(
            [{"source_genome_id": "a", "numeric_source_id": 1}],
            tmp_path,
            api_key="secret",
            transport=lambda *args: data,
        )
    assert not (tmp_path / "batch_0000.csv").exists()
    assert not list(tmp_path.glob("*.tmp"))


def test_retry_exhaustion_and_recovery(tmp_path):
    from urllib.error import HTTPError

    calls = []

    def flaky(*args):
        calls.append(args)
        if len(calls) < 2:
            raise HTTPError("https://pathogen.watch", 503, "temporarily unavailable", {}, None)
        return b'"Genome ID","LIN code"\n"a","0,0"\n'

    result = download_cglin_export(
        [{"source_genome_id": "a", "numeric_source_id": 1}],
        tmp_path,
        api_key="secret",
        transport=flaky,
        sleep=lambda _: None,
    )
    assert result["batches"][0]["attempts"] == 2

    def unavailable(*args):
        raise HTTPError("https://pathogen.watch", 401, "secret", {}, None)

    with pytest.raises(CGLINError, match="HTTP 401"):
        download_cglin_export(
            [{"source_genome_id": "a", "numeric_source_id": 1}],
            tmp_path,
            api_key="secret",
            transport=unavailable,
        )
