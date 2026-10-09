import csv
import json

import pytest

from chronoclade.metadata import Sample
from chronoclade.neighbourhood import build_neighbourhood_evidence, neighbourhood_report_html


def inputs(tmp_path, context_count=3):
    assembly = tmp_path / "assembly.fa"
    assembly.write_text(">x\nACGT\n")
    samples = [Sample("F", assembly, "2020", "Ward A", "Klebsiella pneumoniae", "ST147", "local")]
    samples += [
        Sample(f"C{i}", assembly, "2019", "India", "Klebsiella pneumoniae", "ST147", "context")
        for i in range(1, context_count + 1)
    ]
    tree = tmp_path / "tree.newick"
    tree.write_text(
        "((F:0.001,C1:0.009):0.001,(C2:0.001,C3:0.001):0.001);"
        if context_count == 3
        else "(" + ",".join(f"{s.sample_id}:0.001" for s in samples) + ");"
    )
    return samples, tree


def pairs(path, values):
    with (path / "clonal_pairwise_distances.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_1", "sample_2", "clonal_snps", "callable_sites"])
        for ident, snps, sites in values:
            writer.writerow(["F", ident, snps, sites])


def test_exact_snp_ties_preserved_and_tree_ranking_independent(tmp_path):
    samples, tree = inputs(tmp_path)
    pairs(tmp_path, [("C1", 2, 100), ("C2", 2, 80), ("C3", 9, 100)])
    result = build_neighbourhood_evidence(tree=tree, samples=samples, output=tmp_path, top_n=1)
    summary = result["focal_summaries"][0]
    assert summary["nearest_by_clonal_snps"] == ["C1", "C2"]
    assert summary["nearest_by_tree_distance"] == ["C2", "C3"]
    assert summary["rankings_agree"] is False
    assert [r["context_sample"] for r in result["rows"] if r["metric"] == "clonal_snps"] == [
        "C1",
        "C2",
    ]
    assert all(r["tied_at_distance"] == 2 for r in result["rows"])
    assert (
        next(r for r in result["rows"] if r["context_sample"] == "C2")["snp_proportion"] == 2 / 80
    )
    assert (tmp_path / "genetic_tree.svg").exists() and (tmp_path / "genetic_tree.png").exists()
    assert result["source_tree_sha256"] and result["source_snp_table_sha256"]
    text = neighbourhood_report_html(tmp_path)
    assert "closest in all public data" in text and "Ward A" not in text
    assert "Recorded location" in text
    assert all(r["focal_country"] == "" and r["focal_location"] == "Ward A" for r in result["rows"])
    json.dumps(result, allow_nan=False)


def test_zero_or_missing_callable_not_ranked_as_zero_snps(tmp_path):
    samples, tree = inputs(tmp_path)
    pairs(tmp_path, [("C1", 0, 0), ("C2", 0, ""), ("C3", 10, 100)])
    result = build_neighbourhood_evidence(tree=tree, samples=samples, output=tmp_path)
    summary = result["focal_summaries"][0]
    assert summary["nearest_by_clonal_snps"] == ["C3"]
    assert len(result["invalid_snp_pairs"]) == 2
    assert summary["missing_snp_context_ids"] == ["C1", "C2"]
    assert summary["tree_comparisons"] == 3


def test_no_context_still_has_genetic_tree_and_explicit_status(tmp_path):
    samples, tree = inputs(tmp_path, context_count=0)
    result = build_neighbourhood_evidence(tree=tree, samples=samples, output=tmp_path)
    assert result["rows"] == [] and result["focal_summaries"][0]["status"] == "no_context"
    assert "No context genomes were analysed" in neighbourhood_report_html(tmp_path)
    assert result["tree_pages"][0]["tip_ids"] == ["F"]


@pytest.mark.parametrize(
    "newick",
    [
        "(F:-0.1,C1:0.1,C2:0.2,C3:0.3);",
        "(F,C1:0.1,C2:0.2,C3:0.3);",
        "(F:NaN,C1:0.1,C2:0.2,C3:0.3);",
    ],
)
def test_invalid_tree_keeps_snp_evidence_without_fabricating_branches(tmp_path, newick):
    samples, tree = inputs(tmp_path)
    tree.write_text(newick)
    pairs(tmp_path, [("C1", 1, 100)])
    (tmp_path / "genetic_tree.svg").write_text("stale")
    result = build_neighbourhood_evidence(tree=tree, samples=samples, output=tmp_path)
    assert result["tree_status"] == "unavailable" and result["tree_pages"] == []
    assert result["focal_summaries"][0]["nearest_by_clonal_snps"] == ["C1"]
    assert not (tmp_path / "genetic_tree.svg").exists()


def test_missing_tree_tips_remain_explicit_and_do_not_remove_snp_neighbour(tmp_path):
    samples, tree = inputs(tmp_path)
    tree.write_text("(F:0.1,C2:0.1);")
    pairs(tmp_path, [("C1", 1, 100), ("C2", 10, 100)])
    result = build_neighbourhood_evidence(tree=tree, samples=samples, output=tmp_path)
    assert result["missing_tree_tip_ids"] == ["C1", "C3"]
    assert result["focal_summaries"][0]["nearest_by_clonal_snps"] == ["C1"]
    assert result["focal_summaries"][0]["nearest_by_tree_distance"] == ["C2"]
    assert "missing from the tree" in neighbourhood_report_html(tmp_path)


def test_pagination_accounts_for_every_tip_and_labels_are_escaped(tmp_path):
    samples, tree = inputs(tmp_path)
    samples[1] = Sample(
        "C1",
        samples[1].assembly,
        "2019",
        "<script>alert(1)</script>",
        "Klebsiella pneumoniae",
        "ST147",
        "context",
    )
    pairs(tmp_path, [("C1", 1, 100), ("C2", 2, 100), ("C3", 3, 100)])
    result = build_neighbourhood_evidence(
        tree=tree, samples=samples, output=tmp_path, max_tips_per_page=2
    )
    assert len(result["tree_pages"]) == 2
    shown = [i for page in result["tree_pages"] for i in page["tip_ids"]]
    assert sorted(shown) == sorted(s.sample_id for s in samples)
    assert "Page 1/2" in (tmp_path / "genetic_tree.svg").read_text()
    html = neighbourhood_report_html(tmp_path)
    assert "<script>alert" not in html and "&lt;script&gt;" in html
    assert "&lt;script&gt;" in (tmp_path / "genetic_tree.svg").read_text()


def test_report_shows_tree_and_concise_summary_before_collapsed_rankings(tmp_path):
    samples, tree = inputs(tmp_path)
    pairs(tmp_path, [("C1", 1, 100), ("C2", 2, 100), ("C3", 3, 100)])
    build_neighbourhood_evidence(tree=tree, samples=samples, output=tmp_path)
    text = neighbourhood_report_html(tmp_path)
    assert (
        text.index("<figure>")
        < text.index("Closest relatives of")
        < text.index("All ranked relatives")
    )
    assert "<h2>" not in text
    assert "<details><summary>All ranked relatives and comparison measures" in text
    assert "<details open" not in text


def test_native_renderers_use_accessions_without_changing_audit_ids(tmp_path):
    samples, tree = inputs(tmp_path)
    pairs(tmp_path, [("C1", 1, 100), ("C2", 2, 100), ("C3", 3, 100)])
    labels = {"F": "F", "C1": "ERR123", "C2": "SAMN456", "C3": "C3"}
    (tmp_path / "sample_labels.json").write_text(json.dumps(labels))

    result = build_neighbourhood_evidence(tree=tree, samples=samples, output=tmp_path)
    text = neighbourhood_report_html(tmp_path)

    assert "Closest relatives of F" in text
    assert "ERR123" in text and ">C1<" not in text
    assert "ERR123" in (tmp_path / "genetic_tree.svg").read_text()
    assert result["rows"][0]["context_sample"] == "C1"
    assert sorted(result["tree_pages"][0]["tip_ids"]) == ["C1", "C2", "C3", "F"]
