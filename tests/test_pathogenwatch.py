import json
from urllib.error import HTTPError, URLError

import pytest

from chronoclade.pathogenwatch import (
    PathogenwatchClient,
    PathogenwatchError,
    catalogue_summary,
    deduplicate_catalogue,
    load_api_key,
    load_catalogue,
    normalize_country,
    normalize_record,
)


def search_row(uuid="a", **kwargs):
    return dict(
        uuid=uuid,
        projectAccess="PUBLIC",
        name="name",
        organism="Klebsiella pneumoniae",
        location=None,
        **kwargs,
    )


def detail(uuid="a", **kwargs):
    return dict(
        uuid=uuid,
        id=1,
        organismId="573",
        mlst="147",
        name="SAMN12345",
        qc=True,
        metadata=[["Country", "IN"], ["Date", "2019"]],
        sampleAccession="SAMN12345",
        runAccession="SRR12345",
        studyAccession="PRJNA12345",
        **kwargs,
    )


def normalized(uuid="a", **kwargs):
    return normalize_record(search_row(uuid), detail(uuid, **kwargs), organism_id="573", st="147")


def test_normalization_nested_country_qc_accessions():
    row = normalized()
    assert row["country"] == "India"
    assert row["country_provenance"] == "metadata.Country"
    assert row["biosample"] == "SAMN12345"
    assert row["run_accessions"] == ["SRR12345"]
    assert row["study_accessions"] == ["PRJNA12345"]
    assert row["date_precision"] == "year"
    assert row["collection_date"] == "2019"
    assert row["date_start"] == "2019-01-01"
    assert row["date_end"] == "2019-12-31"
    assert row["source_qc"] is True
    assert "completeness" not in row


@pytest.mark.parametrize(
    "raw,start,end,precision",
    [
        ("2018", "2018-01-01", "2018-12-31", "year"),
        ("2018-02", "2018-02-01", "2018-02-28", "month"),
        ("2020-02-29", "2020-02-29", "2020-02-29", "day"),
        ("2018/2019", "2018-01-01", "2019-12-31", "interval"),
        ("2019-02-31", "", "", "invalid"),
        ("2019/2018", "", "", "invalid"),
        ("", "", "", "missing"),
    ],
)
def test_date_precision(raw, start, end, precision):
    record = detail()
    record["metadata"] = {"Date": raw}
    row = normalize_record(search_row(), record, organism_id="573", st="147")
    assert (row["date_start"], row["date_end"], row["date_precision"]) == (start, end, precision)


def test_country_conflict_unknown_and_qc_failure():
    record = detail()
    record["qc"] = False
    row = normalize_record(dict(search_row(), location="GB"), record, organism_id="573", st="147")
    assert row["country"] == "India" and row["country_conflict"]
    assert row["qc_pass"] is False
    assert normalize_country("United States of America") == normalize_country("US")
    assert normalize_country("United Kingdom:Oxford") == normalize_country("UK")
    record["metadata"] = {}
    row = normalize_record(search_row(), record, organism_id="573", st="147")
    assert row["country"] == "Unknown" and row["date_precision"] == "missing"


def test_private_and_mismatched_details_rejected():
    with pytest.raises(PathogenwatchError, match="non-public"):
        normalize_record(
            dict(search_row(), projectAccess="PRIVATE"), detail(), organism_id="573", st="147"
        )
    with pytest.raises(PathogenwatchError, match="UUID"):
        normalize_record(search_row(), detail("b"), organism_id="573", st="147")
    with pytest.raises(PathogenwatchError, match="organism/ST"):
        normalize_record(search_row(), detail(), organism_id="573", st="1")


def test_deterministic_dedup_and_focal_aliases():
    a, b = normalized("a"), normalized("b")
    a["qc_pass"] = False
    b["country"] = "France"
    c = normalized("c")
    c.update(biosample="", biosample_accessions=[], run_accessions=[], assembly_accessions=[])
    rows, audit = deduplicate_catalogue([c, b, a])
    assert [r["source_genome_id"] for r in rows] == ["b", "c"]
    assert rows[0]["source_genome_ids"] == ["a", "b"] and rows[0]["raw_genome_count"] == 2
    assert rows[1]["identity_resolved"] is False
    assert audit["sample_units"] == 2 and audit["raw_records"] == 3
    assert audit["decisions"][0]["conflicts"]["country"] == ["France", "India"]
    assert deduplicate_catalogue([a, c, b]) == (rows, audit)
    kept, excluded = deduplicate_catalogue([a, b, c], ["SRR12345"])
    assert [r["source_genome_id"] for r in kept] == ["c"]
    assert excluded["focal_excluded_units"] == 1
    # A study is shared across many samples; never deduplicate/exclude via study.
    assert len(deduplicate_catalogue([a, b, c], ["PRJNA12345"])[0]) == 2


def fixture_transport(method, path, *, body, params, headers):
    assert "X-API-Key" not in headers
    if path.endswith("supported"):
        return [
            {
                "organismId": "573",
                "fullName": "Klebsiella pneumoniae",
                "typing": ["MLST"],
                "other": ["LIN Codes"],
            }
        ]
    if path.endswith("genomes"):
        assert body == {"organismId": "573", "mlst": ["147"], "qc": [True, False]}
        uuid = "b" if params.get("after") else "a"
        return {"meta": {"count": 2, "endCursor": uuid}, "genomes": [search_row(uuid)]}
    return detail(params["id"])


def test_freeze_multiple_pages_offline_hash_replay(tmp_path):
    client = PathogenwatchClient(api_key="secret-value", transport=fixture_transport)
    path = tmp_path / "catalogue.json"
    frozen = client.freeze_catalogue(path, page_size=1)
    assert frozen["provenance"]["raw_record_count"] == 2
    assert frozen["provenance"]["sample_unit_count"] == 1
    assert frozen["provenance"]["complete"] is True
    assert load_catalogue(path) == frozen
    assert "secret-value" not in path.read_text()
    assert catalogue_summary(frozen["rows"])["known_country"] == 2
    changed = json.loads(path.read_text())
    changed["rows"][0]["country"] = "France"
    path.write_text(json.dumps(changed))
    with pytest.raises(PathogenwatchError, match="hash mismatch"):
        load_catalogue(path)


@pytest.mark.parametrize(
    "mode", ["private", "count_change", "repeat", "empty", "no_cursor", "no_count"]
)
def test_incomplete_search_never_complete(mode):
    calls = 0

    def transport(method, path, *, body, params, headers):
        nonlocal calls
        calls += 1
        row = search_row("a" if calls == 1 or mode == "repeat" else "b")
        if mode == "private":
            row["projectAccess"] = "PRIVATE"
        meta = {"count": 2, "endCursor": "a"}
        if mode == "count_change" and calls == 2:
            meta["count"] = 3
        if mode == "no_count":
            meta.pop("count")
        if mode == "no_cursor":
            meta.pop("endCursor")
        return {"meta": meta, "genomes": [] if mode == "empty" else [row]}

    with pytest.raises(PathogenwatchError):
        PathogenwatchClient(transport=transport).search("573", "147")


def test_retry_transient_failure_access_denial_and_redaction():
    calls = []

    def transport(method, path, **kwargs):
        calls.append(1)
        if len(calls) == 1:
            raise HTTPError("url", 503, "temporary", {}, None)
        return {"ok": True}

    client = PathogenwatchClient(transport=transport, sleep=lambda _: None)
    assert client.request_json("GET", "/api/example") == {"ok": True} and len(calls) == 2

    def denied(*args, **kwargs):
        raise HTTPError("secret-url", 403, "secret-token", {}, None)

    with pytest.raises(PathogenwatchError, match="HTTP 403") as error:
        PathogenwatchClient(transport=denied).request_json("GET", "/api/example")
    assert "secret" not in str(error.value)

    def failure(*args, **kwargs):
        raise URLError("secret")

    with pytest.raises(PathogenwatchError, match="after retries"):
        PathogenwatchClient(transport=failure, sleep=lambda _: None).request_json(
            "GET", "/api/example"
        )


def test_credentials_permissions_and_environment(tmp_path, monkeypatch):
    monkeypatch.delenv("PATHOGENWATCH_API_KEY", raising=False)
    path = tmp_path / "key.json"
    path.write_text('{"api_key":"secret"}')
    path.chmod(0o644)
    with pytest.raises(PathogenwatchError, match="0600"):
        load_api_key(path)
    path.chmod(0o600)
    assert load_api_key(path) == "secret"
    monkeypatch.setenv("PATHOGENWATCH_API_KEY", "env-secret")
    assert load_api_key(path) == "env-secret"


def test_source_bounds_preserve_year_and_month_precision():
    for start, end, precision in [
        ("2019-01-01", "2019-12-31", "year"),
        ("2020-02-01", "2020-02-29", "month"),
    ]:
        record = detail()
        record.update(metadata={}, startDate=start + "T00:00:00Z", endDate=end + "T00:00:00Z")
        row = normalize_record(search_row(), record, organism_id="573", st="147")
        assert row["date_precision"] == precision


def test_real_search_decorated_country_labels_are_aliases():
    for label, code in [
        ("DE - Germany", "DE"),
        ("TR - Turkey", "TR"),
        ("US - United States of America", "US"),
    ]:
        assert normalize_country(label) == normalize_country(code)
    assert normalize_country("Turkey") == normalize_country("TR")
    assert normalize_country("india") == normalize_country("IN")
    assert "PRJNA12345" not in normalized()["aliases"]


def test_malformed_metadata_pairs_fail_explicitly():
    record = detail()
    record["metadata"] = [["Country", "GB", "extra"]]
    with pytest.raises(PathogenwatchError, match="metadata pairs"):
        normalize_record(search_row(), record, organism_id="573", st="147")


def test_real_country_aliases_do_not_split_groups():
    assert normalize_country("Myanmar") == normalize_country("MM")
    assert normalize_country("Myanmar (Burma)") == normalize_country("MM")
    assert normalize_country("UAE") == normalize_country("AE")


def test_ena_biosample_accessions_are_recognised_and_deduplicated():
    records = []
    for ident in ("a", "b"):
        raw = detail(ident)
        raw.update(name="SAMEA8263140", sampleAccession="SAMEA8263140", runAccession=None)
        row = normalize_record(search_row(ident), raw, organism_id="573", st="147")
        assert row["biosample"] == "SAMEA8263140"
        records.append(row)
    units, audit = deduplicate_catalogue(records, ["SAMEA8263140"])
    assert units == [] and audit["sample_units"] == 1 and audit["focal_excluded_units"] == 1
