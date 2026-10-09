import csv
import json
from pathlib import Path

import pytest

from chronoclade.country_network import build_country_network, country_network_report_html


def inputs(tmp_path: Path, *, root=(0.95, 0.05), inner=(0.95, 0.05), locations=None):
    location = tmp_path / "location"
    location.mkdir()
    (location / "GTR.txt").write_text(
        "Character to attribute mapping:\n  A: Italy\n  B: India\n  C: ?\n\n"
    )
    (location / "annotated_tree.nexus").write_text(
        '#NEXUS\nBegin Trees;\n Tree tree1=(a:1[&location="India"],'
        '(b:1[&location="India"],c:1[&location="Italy"])inner:1[&location="Italy"])'
        'root:0[&location="Italy"];\nEnd;\n'
    )
    (location / "confidence.csv").write_text(
        "#name, A, B\nroot, "
        + ",".join(map(str, root))
        + "\ninner, "
        + ",".join(map(str, inner))
        + "\na,0,1\nb,0,1\nc,1,0\n"
    )
    with (tmp_path / "metadata.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample_id", "location"])
        writer.writerows(locations or [("a", "India"), ("b", "India"), ("c", "Italy")])
    return tmp_path


def test_branch_aggregation_and_separate_confidences(tmp_path):
    result = build_country_network(inputs(tmp_path))
    assert result["status"] == "available"
    assert result["nodes"] == [
        {"country": "India", "observed_tips": 2},
        {"country": "Italy", "observed_tips": 1},
    ]
    assert result["edges"] == [
        {
            "source": "Italy",
            "target": "India",
            "branches": 2,
            "high_endpoint_branches": 2,
            "uncertain_branches": 0,
        }
    ]
    assert len(result["branches"]) == 4
    edge = result["branches"][0]
    assert edge["parent_marginal_confidence"] == 0.95
    assert edge["child_marginal_confidence"] == 1
    assert "support" not in edge
    for extension in ("json", "svg", "png"):
        assert (tmp_path / f"country_network.{extension}").is_file()
    assert "not joint support" in country_network_report_html(tmp_path)


def test_ties_and_weak_ancestors(tmp_path):
    result = build_country_network(inputs(tmp_path, root=(0.5, 0.5), inner=(0.6, 0.4)))
    assert result["edges"][0]["uncertain_branches"] == 2
    assert result["edges"][0]["high_endpoint_branches"] == 0
    assert result["branches"][0]["parent_tied"]
    assert result["branches"][0]["parent_candidates"] == "India | Italy"


def test_tied_maximum_even_with_low_threshold(tmp_path):
    result = build_country_network(inputs(tmp_path, root=(0.5, 0.5)), confidence_threshold=0.4)
    assert not result["branches"][0]["both_endpoints_high_confidence"]


def test_missing_reconstruction_writes_placeholder(tmp_path):
    result = build_country_network(tmp_path)
    assert result["status"] == "unavailable"
    assert "Missing input" in result["reason"]
    assert "unavailable" in country_network_report_html(tmp_path)
    assert (tmp_path / "country_network.svg").is_file()
    assert list(csv.DictReader((tmp_path / "country_network_edges.csv").open())) == []


@pytest.mark.parametrize("value", ["nan", "-.1", "1.1", ".2"])
def test_invalid_probabilities_rejected(tmp_path, value):
    inputs(tmp_path)
    path = tmp_path / "location/confidence.csv"
    path.write_text(path.read_text().replace("root, 0.95,0.05", f"root, {value},0.05"))
    result = build_country_network(tmp_path)
    assert result["status"] == "unavailable"
    assert result["edges"] == []


def test_metadata_annotation_disagreement(tmp_path):
    result = build_country_network(
        inputs(tmp_path, locations=[("a", "Italy"), ("b", "India"), ("c", "Italy")])
    )
    assert result["status"] == "unavailable"
    assert "disagrees with metadata" in result["reason"]


def test_exact_tip_join_skips_unmatched(tmp_path):
    result = build_country_network(
        inputs(tmp_path, locations=[("a", "India"), ("c", "Italy"), ("extra", "Italy")])
    )
    assert result["status"] == "partial"
    assert result["edges"][0]["branches"] == 1
    assert len(result["excluded_tips"]) == 2
    assert next(x for x in result["branches"] if x["child"] == "b")["included"] is False


def test_unknown_metadata_equivalent_to_missing_annotation(tmp_path):
    inputs(tmp_path, locations=[("a", "Unknown country"), ("b", "India"), ("c", "Italy")])
    tree = tmp_path / "location/annotated_tree.nexus"
    tree.write_text(tree.read_text().replace('a:1[&location="India"]', 'a:1[&location="?"]'))
    result = build_country_network(tmp_path)
    assert result["status"] == "partial"
    assert result["unknown_tips"] == 1
    assert result["edges"][0]["branches"] == 1
    assert not result["branches"][0]["included"]


def test_no_dating_dependency_and_reproducible_counts(tmp_path):
    inputs(tmp_path)
    first = build_country_network(tmp_path, temporal_supported=False)
    second = build_country_network(tmp_path, temporal_supported=True)
    assert first["edges"] == second["edges"]
    assert first["branches"] == second["branches"]
    assert json.loads((tmp_path / "country_network.json").read_text())["temporal_supported"] is True


def test_duplicate_metadata_is_unavailable(tmp_path):
    result = build_country_network(inputs(tmp_path, locations=[("a", "India"), ("a", "India")]))
    assert result["status"] == "unavailable"
    assert "duplicate metadata" in result["reason"]


def test_unmapped_confidence_is_unavailable(tmp_path):
    inputs(tmp_path)
    confidence = tmp_path / "location/confidence.csv"
    confidence.write_text(confidence.read_text().replace("#name, A, B", "#name, A, Z"))
    assert "Unmapped" in build_country_network(tmp_path)["reason"]


def test_single_location_has_explicit_reason(tmp_path):
    (tmp_path / "metadata.csv").write_text("sample_id,location\na,Italy\nb,Italy\n")
    result = build_country_network(tmp_path)
    assert result["status"] == "unavailable"
    assert (
        result["unavailable_reason"]
        == "Fewer than two known locations; no between-location reconstruction was performed"
    )


def test_malformed_tree_is_unavailable(tmp_path):
    inputs(tmp_path)
    (tmp_path / "location/annotated_tree.nexus").write_text("Tree tree1=((a,b)root;")
    result = build_country_network(tmp_path)
    assert result["status"] == "unavailable"


def test_unknown_observed_tip_with_imputed_known_state(tmp_path):
    inputs(tmp_path, locations=[("a", "Unknown"), ("b", "India"), ("c", "Italy")])
    result = build_country_network(tmp_path)
    assert result["status"] == "partial"
    assert result["unknown_tips"] == 1
    assert result["nodes"] == [
        {"country": "India", "observed_tips": 1},
        {"country": "Italy", "observed_tips": 1},
    ]
    assert result["edges"][0]["branches"] == 1
    assert result["branches"][0]["reason"] == "Missing observed location"
    assert not result["branches"][0]["included"]


@pytest.mark.parametrize("length", ["-1", "1e999"])
def test_invalid_branch_lengths_are_unavailable(tmp_path, length):
    inputs(tmp_path)
    tree = tmp_path / "location/annotated_tree.nexus"
    tree.write_text(tree.read_text().replace("a:1[", f"a:{length}["))
    result = build_country_network(tmp_path)
    assert result["status"] == "unavailable"
    assert "finite and nonnegative" in result["reason"]


def test_missing_node_id_is_unavailable(tmp_path):
    inputs(tmp_path)
    tree = tmp_path / "location/annotated_tree.nexus"
    tree.write_text(tree.read_text().replace(")inner:1[", "):1["))
    result = build_country_network(tmp_path)
    assert result["status"] == "unavailable"
    assert "Missing tree node ID" in result["reason"]


def test_short_metadata_row_is_explicitly_unavailable(tmp_path):
    inputs(tmp_path)
    (tmp_path / "metadata.csv").write_text("sample_id,location\na,India\nb\nc,Italy\n")
    result = build_country_network(tmp_path)
    assert result["status"] == "unavailable"
    assert "Malformed metadata row" in result["reason"]


def test_empty_confidence_table_is_explicitly_unavailable(tmp_path):
    inputs(tmp_path)
    (tmp_path / "location/confidence.csv").write_text("")
    result = build_country_network(tmp_path)
    assert result["status"] == "unavailable"
    assert result["reason"] == "Location confidence table is empty"
