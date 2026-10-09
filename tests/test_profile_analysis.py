from pathlib import Path

import pytest

from chronoclade.profile_analysis import analyse_profiles, _date_interval


def record(name, calls=None, **kwargs):
    return {
        "sample_id": name,
        "origin": "context",
        "species": "Klebsiella pneumoniae",
        "cgmlst_scheme": "test",
        "cgmlst_scheme_version": "2026-01-01",
        "cgmlst_loci": ["a", "b", "c", "d"],
        "cgmlst_profile": dict(zip(["a", "b", "c", "d"], calls or [1, 1, 1, 1])),
        **kwargs,
    }


def test_accessions_label_tree_but_preserve_distance_and_tree_identity(tmp_path):
    rows = [record("q", origin="local", country="Greece", collection_date="2019-05"),
            record("PW_public", run_accessions=["SRR32641190"],
                   biosample_accessions=["SAMN46159676"])]
    result = analyse_profiles(rows, output=tmp_path, bootstrap_replicates=0)
    assert result["sample_labels"]["PW_public"] == "SRR32641190"
    assert result["nearest_neighbours"][0]["context_id"] == "PW_public"
    tree = result["cohorts"][0]
    assert "PW_public" in Path(tree["tree_path"]).read_text()
    assert "SRR32641190" in Path(tree["tree_figure"]).read_text()
    svg = Path(tree["tree_figure"]).read_text()
    assert "q | Greece | 2019-05" in svg
    assert "SRR32641190 | Country unknown | Date unknown" in svg
    assert "Input genomes" in svg and "Public comparisons" in svg
    assert "#2166ac" in svg and "#666666" in svg
    assert "Nearest public relatives (ties included)" in svg
    assert "#e66101" in svg
    assert "Inner" not in svg
    assert "SAMN46159676" in Path(result["paths"]["sample_labels"]).read_text()


def test_categorical_distance_ties_and_separate_temporal_outputs(tmp_path):
    rows = [
        record("q", origin="local", collection_date="2020", country="UK"),
        record("a", [99999, 1, 1, 1], collection_date="2021-02", country="UK"),
        record("b", [2, 1, 1, 1], collection_date="2021-02-02", country="France"),
    ]
    result = analyse_profiles(
        rows, output=tmp_path, distance_threshold=0.25, bootstrap_replicates=5
    )
    assert len(result["nearest_neighbours"]) == 2
    assert all(row["allele_differences"] == 1 for row in result["nearest_neighbours"])
    assert all(row["tied_neighbours"] == 2 for row in result["nearest_neighbours"])
    svg = Path(result["cohorts"][0]["tree_figure"]).read_text()
    assert "Nearest public relatives (ties included)" in svg
    # The two nearest public tips and the legend use the orange star outline.
    assert svg.count("stroke: #e66101") >= 3
    assert result["temporal_persistence"][0]["observed_years"] == [2020, 2021]
    assert result["time_place_concentration"][0]["largest_cell_fraction_of_annotated"] == 1 / 3
    assert result["root_to_tip"][0]["points"][2]["precision"] == "year"
    assert result["cohorts"][0]["tree_units"].startswith("fraction")
    assert all(edge["uncertain"] for edge in result["location_network"][0]["edges"])
    for path in result["paths"].values():
        assert Path(path).exists()


def test_missing_loci_not_matches_and_scheme_isolation(tmp_path):
    rows = [
        record("q", origin="local"),
        record("missing", [1, 1, 0, 0]),
        record("version", cgmlst_scheme_version="other"),
        record("untyped", cgmlst_profile={}),
    ]
    result = analyse_profiles(rows, output=tmp_path, bootstrap_replicates=0)
    assert result["nearest_neighbours"][0]["status"] == "no_comparable_context"
    assert len(result["cohorts"]) == 2
    assert {r["sample_id"] for r in result["exclusions"]} == {"missing", "untyped"}
    assert result["coverage"]["available_profiles"] == 2


def test_novel_hashes_are_categorical_and_database_conflicts(tmp_path):
    q = record(
        "q",
        origin="local",
        cgmlst_profile={"a": 1, "b": 1, "c": 1},
        cgmlst_novel_alleles={"d": "a" * 40},
        cgmlst_database_sha256="a" * 64,
    )
    c = record(
        "c",
        cgmlst_profile={"a": 1, "b": 1, "c": 1},
        cgmlst_novel_alleles={"d": "b" * 40},
        cgmlst_database_sha256="a" * 64,
    )
    bad = record("bad", cgmlst_database_sha256="b" * 64)
    result = analyse_profiles([q, c, bad], output=tmp_path, bootstrap_replicates=0)
    assert result["nearest_neighbours"][0]["context_id"] == "c"
    assert result["nearest_neighbours"][0]["allele_differences"] == 1
    assert len(result["cohorts"]) == 2


def test_complete_linkage_does_not_chain(tmp_path):
    result = analyse_profiles(
        [record("a"), record("b", [2, 1, 1, 1]), record("c", [2, 2, 1, 1])],
        output=tmp_path,
        distance_threshold=0.25,
        bootstrap_replicates=8,
    )
    assert sorted(len(g["sample_ids"]) for g in result["genetic_groups"]) == [1, 2]
    assert all(g["valid_bootstrap_replicates"] == 8 for g in result["genetic_groups"])


def test_reproducible_bootstrap(tmp_path):
    rows = [record("a"), record("b", [2, 1, 1, 1]), record("c", [2, 2, 1, 1])]
    a = analyse_profiles(rows, output=tmp_path / "a", seed=7, bootstrap_replicates=10)
    b = analyse_profiles(
        list(reversed(rows)), output=tmp_path / "b", seed=7, bootstrap_replicates=10
    )
    assert a["genetic_groups"] == b["genetic_groups"]


def test_date_precision_and_invalid_dates():
    assert _date_interval({"collection_date": "2020-02"})["end"] == "2020-02-29"
    assert _date_interval({"collection_date": "2020-02-29"})["precision"] == "day"
    assert _date_interval({"collection_date": "2021-02-29"}) is None


def test_empty_and_input_guard(tmp_path):
    result = analyse_profiles([], output=tmp_path)
    assert result["cohorts"] == []
    with pytest.raises(ValueError, match="unique"):
        analyse_profiles([record("a"), record("a")], output=tmp_path)
    with pytest.raises(ValueError, match="limited"):
        analyse_profiles([{}] * 1501, output=tmp_path)


def test_observed_locus_universe_warning_and_missing_location(tmp_path):
    rows = [
        record("a", country="UK", cgmlst_locus_universe_complete=False),
        record("b", country="UK", cgmlst_locus_universe_complete=False),
        record("c", country="", cgmlst_locus_universe_complete=False),
    ]
    result = analyse_profiles(rows, output=tmp_path, bootstrap_replicates=0)
    assert "canonical completeness unknown" in result["cohorts"][0]["coverage_denominator"]
    assert result["location_network"][0]["edges"] == []


def test_species_not_mixed_into_same_tree(tmp_path):
    result = analyse_profiles(
        [record("a"), record("b", species="Escherichia coli")],
        output=tmp_path,
        bootstrap_replicates=0,
    )
    assert len(result["cohorts"]) == 2
    assert result["coverage"]["unavailable_pairs_by_reason"] == {"different_species": 1}


def test_sparse_export_observed_union_is_explicit_and_input_unmodified(tmp_path):
    a = record("a", cgmlst_locus_universe_complete=False)
    b = record(
        "b",
        cgmlst_loci=["a", "b", "c"],
        cgmlst_profile={"a": 1, "b": 1, "c": 1},
        cgmlst_locus_universe_complete=False,
    )
    result = analyse_profiles([a, b], output=tmp_path, min_overlap=0.75, bootstrap_replicates=0)
    assert len(result["cohorts"]) == 1
    assert result["coverage"]["comparable_pairs"] == 1
    assert b["cgmlst_loci"] == ["a", "b", "c"]
    assert "observed export locus union" in result["cohorts"][0]["coverage_denominator"]


def test_one_year_does_not_claim_persistence(tmp_path):
    result = analyse_profiles(
        [record("a", collection_date="2020"), record("b", collection_date="2020-02")],
        output=tmp_path,
        bootstrap_replicates=3,
    )
    assert result["bootstrap_replicates"] == 3
    assert all(
        "cannot be assessed" in group["interpretation"] for group in result["temporal_persistence"]
    )
    assert Path(result["cohorts"][0]["tree_figure"]).exists()


def test_impossible_future_dates_retained_but_excluded_from_time_analysis(tmp_path):
    rows = [
        record("q", origin="local", collection_date="2020", country="UK"),
        record("future", collection_date="2386", country="France"),
    ]
    result = analyse_profiles(rows, output=tmp_path, bootstrap_replicates=0)
    assert rows[1]["collection_date"] == "2386"
    assert result["date_exclusions"] == [
        {"sample_id": "future", "collection_date": "2386", "reason": "future_collection_date"}
    ]
    assert result["temporal_persistence"][0]["observed_years"] == [2020]
    assert result["time_place_concentration"][0]["records_with_time_and_place"] == 1
    assert [point["sample_id"] for point in result["root_to_tip"][0]["points"]] == ["q"]
    neighbour = result["nearest_neighbours"][0]
    assert neighbour["country"] == "France"
    assert neighbour["collection_date"] == "2386"
    assert neighbour["year"] is None
    assert neighbour["date_status"] == "future_collection_date"
    assert neighbour["cohort_id"] == "cohort_1"
    assert any("future collection dates" in warning for warning in result["cohorts"][0]["warnings"])


def test_tree_limit_retains_nearest_ties_without_limiting_distance_search(tmp_path):
    rows = [record("q", origin="local")]
    rows += [record(f"near{i}", [2, 1, 1, 1]) for i in range(3)]
    rows += [record(f"far{i}", [2, 2, 2, 2]) for i in range(4)]
    result = analyse_profiles(rows, output=tmp_path, tree_limit=3, bootstrap_replicates=0)
    assert result["coverage"]["comparable_pairs"] == 28
    assert {r["context_id"] for r in result["nearest_neighbours"]} == {"near0", "near1", "near2"}
    cohort = result["cohorts"][0]
    assert set(cohort["tree_display_sample_ids"]) == {"q", "near0", "near1", "near2"}
    assert len(cohort["sample_ids"]) == 8
    assert "far3" in Path(cohort["tree_path"]).read_text()
    assert "far3" not in Path(cohort["tree_figure"]).read_text()
