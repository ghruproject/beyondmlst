import json
from collections import defaultdict
from pathlib import Path

import pytest

from chronoclade.context_geography import (
    CAUTION,
    country_colour,
    generate_context_geography,
    geography_tables,
)


def row(identifier, country, prefix=(1, 2, 3, 4, 5, 6, 7), biosample=None, **kwargs):
    result = {
        "source": "pathogenwatch",
        "source_genome_id": identifier,
        "country": country,
        "qc_pass": True,
        "biosample": biosample,
        "cglin_scheme": "KpSC",
        "cglin_scheme_version": "v1",
        "cglin_raw": ".".join(map(str, prefix)),
        "cglin_status": "complete" if len(prefix) == 7 else "partial",
    }
    for depth in (5, 6, 7):
        result[f"cglin_group_{depth}"] = (
            json.dumps(["KpSC", "v1", list(prefix[:depth])]) if len(prefix) >= depth else ""
        )
        result[f"cglin_status_{depth}"] = "resolved" if len(prefix) >= depth else "partial"
    result.update(kwargs)
    return result


@pytest.fixture
def catalogue():
    return [
        row("a", "India", biosample="SAM1"),
        row("a2", "India", biosample="SAM1"),
        row("b", "United Kingdom", biosample="SAM2"),
        row("c", "Unknown", biosample="SAM3"),
        row("d", "France", (8, 2, 3, 4, 5, 6, 7), biosample="SAM4"),
        row("partial", "India", (1, 2, 3, 4, 5), biosample="SAM5"),
        row("qc_fail", "Germany", qc_pass=False),
    ]


def test_hand_counts_and_denominators(catalogue):
    result = geography_tables(
        catalogue,
        selected_source_ids=["a2", "d"],
        focal_rows=[row("focal", "India", biosample="SAM1")],
    )
    assert result["raw_catalogue_records"] == 7
    assert result["qc_eligible_records"] == 6
    assert result["sample_units"] == 5
    assert len(result["duplicate_audit"]) == 6
    assert len(result["focal_overlap"]) == 1
    groups = defaultdict(list)
    for item in result["rows"]:
        groups[
            (item["cohort"], item["prefix_depth"], item["prefix_key"], item["group_label"])
        ].append(item)
    for members in groups.values():
        assert sum(r["count"] for r in members) == members[0]["denominator"]
        assert sum(r["percentage"] for r in members) == pytest.approx(100)
        assert members[0]["known_country_n"] + members[0]["unknown_n"] == members[0]["denominator"]
    public = [
        r for r in result["rows"] if r["cohort"] == "public_catalogue" and r["prefix_depth"] == 7
    ]
    group = [r for r in public if r["prefix_key"] == catalogue[0]["cglin_group_7"]]
    assert {r["country"]: r["count"] for r in group} == {
        "India": 1,
        "United Kingdom": 1,
        "Unknown": 1,
    }
    assert sum(r["raw_record_count"] for r in public) == 6
    assert len({r["prefix_key"] for r in public if r["assignment_category"] == "lineage"}) == 2
    assert [r for r in public if r["assignment_category"] == "assignment_coverage"][0][
        "country"
    ] == "India"
    summaries = {(s["cohort"], s["depth"]): s for s in result["summaries"]}
    assert summaries["public_catalogue", 5]["assigned_units"] == 5
    assert summaries["public_catalogue", 7]["assigned_units"] == 4
    assert summaries["selected_context", 7]["total_units"] == 2
    assert summaries["focal_survey", 7]["total_units"] == 1
    assert summaries["public_catalogue", 7]["singleton_groups"] == 1


def test_deduplicated_provider_rows():
    frozen = row(
        "a",
        "India",
        sample_unit_id="biosample:SAM1",
        raw_genome_count=3,
        source_genome_ids=["a", "b", "c"],
    )
    result = geography_tables([frozen], selected_source_ids=["c"])
    assert result["raw_catalogue_records"] == 3
    assert result["qc_eligible_records"] == 3
    assert result["sample_units"] == 1
    assert all(r["raw_record_count"] == 3 for r in result["rows"])


def test_assembly_chain_and_conflict():
    records = [
        row("a", "India", biosample="SAM1", assembly_accessions=["GCA_1"]),
        row("b", "France", assembly_accessions=["GCA_1"]),
    ]
    result = geography_tables(records)
    assert result["sample_units"] == 1
    assert all(r["country"] == "Unknown" for r in result["rows"])
    assert all(r["country_conflict"] for r in result["duplicate_audit"])


def test_exports_offline_and_stable(tmp_path, catalogue):
    result = generate_context_geography(
        catalogue,
        tmp_path,
        depths=(7,),
        selected_source_ids=["a", "d"],
        scope={"snapshot": "2026-10-09", "description": "KpSC ST147 public"},
    )
    assert "United Kingdom" in (tmp_path / "country_composition.csv").read_text()
    assert "Unresolved" in result["report_html"]
    assert CAUTION.replace(">", "&gt;") in result["report_html"]
    assert "<svg" in result["report_html"]
    assert "selected_context" in result["report_html"]
    assert "N=3; known=2; Unknown=1" in result["report_html"]
    from PIL import Image

    for path in map(Path, result["outputs"]):
        assert path.exists()
        if path.suffix == ".png":
            with Image.open(path) as image:
                assert image.width >= 1500 and image.height >= 600
    assert country_colour("India") == country_colour("India")
    assert country_colour("Unknown") != country_colour("Other")


def test_empty_missing_unsupported(tmp_path):
    empty = generate_context_geography([], tmp_path / "empty")
    assert "No eligible" in empty["report_html"]
    assert not list((tmp_path / "empty").glob("*.png"))
    missing = geography_tables([{"source_genome_id": "x", "country": "Unknown", "qc_pass": True}])
    assert all(r["assignment_category"] == "assignment_coverage" for r in missing["rows"])
    assert all(s["assigned_units"] == 0 for s in missing["summaries"])
    with pytest.raises(ValueError, match="5, 6 and 7"):
        geography_tables([], depths=(8,))


def test_duplicate_assignment_conflict_never_inherits_representative_lineage():
    result = geography_tables([row('a', 'India', biosample='SAM1'),
        row('b', 'India', biosample='SAM1', cgst='different')])
    assert result['sample_units'] == 1
    assert all(r['assignment_category'] == 'assignment_coverage' for r in result['rows'])
    assert all(r['assignment_status'] == 'conflict' for r in result['rows'])
    assert all(r['cglin_conflict'] for r in result['duplicate_audit'])
