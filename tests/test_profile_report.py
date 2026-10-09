import numpy as np
from pathlib import Path

from chronoclade.profile_report import (
    write_corrected_report,
    write_fast_group_index,
    write_profile_report,
    write_stage_index,
)


def test_report_uses_accession_without_changing_analysis_or_style(tmp_path):
    analysis = {"nearest_neighbours": [{"query_id": "q", "context_id": "PW_public",
                "status": "matched", "distance": 0}],
                "sample_labels": {"q": "q", "PW_public": "SRR32641190"}}
    html = write_profile_report(analysis, directory=tmp_path).read_text()
    assert "SRR32641190" in html
    assert "PW_public" not in html
    assert analysis["nearest_neighbours"][0]["context_id"] == "PW_public"
    assert "Archivo" in html


def test_report_orders_summary_figures_and_groups_before_neighbours(tmp_path: Path) -> None:
    (tmp_path / "figures").mkdir()
    (tmp_path / "figures" / "pcoa.svg").write_text(
        '<svg xmlns="http://www.w3.org/2000/svg"></svg>', encoding="utf-8"
    )
    (tmp_path / "figures" / "cohort.nwk").write_text("(Q1:0.1,P1:0.1);\n", encoding="utf-8")
    (tmp_path / "figures" / "cohort_nj.svg").write_text(
        '<svg xmlns="http://www.w3.org/2000/svg"></svg>', encoding="utf-8"
    )
    analysis = {
        "schema_version": 1,
        "records_count": 3,
        "coverage": {"available_profiles": 2, "comparable_pairs": 1},
        "exclusions": [{"sample_id": "Q3", "reason": "no_profile"}],
        "geography": [{"country": "Kenya", "count": 2, "denominator": 3}],
        "paths": {"pcoa": "figures/pcoa.svg", "nj_tree": "figures/missing.svg"},
        "genetic_groups": [{"group_id": "G1", "sample_count": 2}],
        "temporal_persistence": [{"group_id": "G1", "years": [2020, 2022]}],
        "time_place_concentration": [{"group_id": "G2", "time_window": "2023", "location": "Kenya"}],
        "cohorts": [{"cohort_id": "cohort_1", "sample_ids": ["Q1", "P1"],
                     "tree_units": "allele differences", "tree_path": "figures/cohort.nwk",
                     "tree_figure": "figures/cohort_nj.svg",
                     "pcoa_figure": "figures/pcoa.svg", "country_figure": "figures/pcoa.svg",
                     "warnings": []}],
        "root_to_tip": [{"cohort_id": "cohort_1", "root": "Q1", "midpoint_date_slope": None,
                         "midpoint_date_pearson_r": None,
                         "interpretation": "Exploratory arbitrary-root diagnostic; no gate."}],
        "location_network": [{"cohort_id": "cohort_1", "method": "alternate roots", "tested_roots": 2,
                              "ambiguous_edge_reconstructions": 1,
                              "interpretation": "Root-dependent possibilities, not transmission.",
                              "edges": [{"source": "Kenya", "target": "Uganda",
                                         "roots_with_possible_change": 1, "roots_tested": 2,
                                         "root_fraction": 0.5, "uncertain": True}]}],
        "nearest_neighbours": [{"query_id": "Q1", "context_id": "P1", "status": "matched", "distance": 2}],
    }

    output = write_profile_report(
        analysis,
        directory=tmp_path,
        stage="fast",
        provenance={"coverage": {"queries": {"total": 1}, "context": {"total": 8, "profiles_available": 8}}},
    )
    html = output.read_text(encoding="utf-8")

    assert output.name == "profile_report.html"
    assert html.index("id=\"summary\"") < html.index("id=\"geography\"")
    assert html.index("id=\"geography\"") < html.index("id=\"nearest\"")
    assert html.index("id=\"nearest\"") < html.index("id=\"figures\"")
    assert html.index("id=\"nearest\"") < html.index("id=\"groups\"")
    assert 'id="persistence"' not in html
    assert "concentration in time and place" in html
    assert 'src="figures/pcoa.svg"' in html
    assert 'src="figures/cohort_nj.svg"' in html
    assert 'href="figures/cohort.nwk"' in html
    assert "alternate roots" in html
    assert "8 profiles available / 8 records" in html
    assert "missing.svg" not in html
    assert "profile_report" not in html or 'href="profile_report' not in html


def test_zero_profiles_keeps_coverage_and_says_findings_are_unavailable(tmp_path: Path) -> None:
    output = write_profile_report(
        {"schema_version": 1, "records_count": 4, "coverage": {"available_profiles": 0},
         "exclusions": [{"sample_id": f"Q{i}", "reason": "no_profile"} for i in range(4)]},
        directory=tmp_path,
    )
    html = output.read_text(encoding="utf-8")

    assert "No usable profiles were available" in html
    assert "available profiles" in html
    assert "Excluded sample IDs" in html
    assert "No groups were reported" in html
    assert "Nearest-relative results were not generated" in html


def test_report_escapes_user_text_and_neighbour_output_links_only_existing_files(tmp_path: Path) -> None:
    neighbour = tmp_path / "nearest.tsv"
    neighbour.write_text("query\tnearest\n", encoding="utf-8")
    output = write_profile_report(
        {"records_count": 1, "title": "<script>alert(1)</script>",
         "nearest_neighbours": "nearest.tsv",
         "paths": {"location_network": "absent.svg"}},
        directory=tmp_path,
        stage="full",
    )
    html = output.read_text(encoding="utf-8")

    assert "<script>alert(1)</script>" not in html
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in html
    assert 'href="nearest.tsv"' in html
    assert "absent.svg" not in html
    assert "rooted tree" in html


def test_stage_index_links_existing_reports_only(tmp_path: Path) -> None:
    (tmp_path / "fast").mkdir()
    (tmp_path / "fast" / "profile_report.html").write_text("fast", encoding="utf-8")
    index = write_stage_index(
        tmp_path,
        {"fast": {"report": "fast/profile_report.html", "status": "complete"},
         "finish": "finish/profile_report.html"},
    )

    html = index.read_text(encoding="utf-8")
    assert 'href="fast/profile_report.html"' in html
    assert "finish/profile_report.html" not in html


def test_corrected_report_marks_dating_not_assessed_and_keeps_guardrails(tmp_path: Path) -> None:
    output = write_corrected_report(
        {
            "species": "Klebsiella pneumoniae",
            "lineage": "ST23",
            "sample_count": 12,
            "context": {"local_samples": 5, "context_samples": 7},
            "temporal_status": "not_assessed",
            "temporal_signal_reason": "Fewer than three usable distinct collection dates; temporal signal was not assessed.",
        },
        directory=tmp_path,
    )
    html = output.read_text(encoding="utf-8")

    assert "Recombination-adjusted analysis" in html
    assert "Fewer than three usable distinct collection dates" in html
    assert "No temporal test or calendar tree is implied" in html
    assert "not proof of transmission" in html
    assert "country network reconstruction is not available" in html


def test_report_separates_query_context_coverage_and_explains_unknown_support(tmp_path: Path) -> None:
    output = write_profile_report(
        {
            "schema_version": 1,
            "records_count": 4,
            "coverage": {"available_profiles": 2, "comparable_pairs": 1},
            "exclusions": [
                {"sample_id": "Q2", "reason": "profile_unavailable"},
                {"sample_id": "C2", "reason": "profile_unavailable"},
            ],
            "bootstrap_replicates": 10,
            "cohorts": [{
                "cohort_id": "cohort_1", "sample_ids": ["Q1", "C1"],
                "tree_units": "cgMLST mismatches", "coverage_denominator":
                "observed export locus union; canonical completeness unknown",
                "warnings": ["Canonical completeness unknown."],
            }],
            "genetic_groups": [{
                "group_id": "G1", "sample_ids": ["Q1", "C1"],
                "method": "complete linkage", "valid_bootstrap_replicates": 0,
                "minimum_pair_bootstrap_coassignment": None, "singleton": False,
            }],
        },
        directory=tmp_path,
        provenance={"analysis_export_status": "unavailable_no_credentials", "coverage": {
            "queries": {"total": 2, "profiles_available": 1, "profiles_missing": 1},
            "context": {"total": 2, "profiles_available": 1, "profiles_missing": 1},
        }},
    )
    html = output.read_text(encoding="utf-8")

    assert "Query dataset" in html and "Public context" in html
    assert "no API credentials were configured" in html
    assert "Records without profiles" in html
    assert "combine query and public context records" in html
    assert "observed export locus union; canonical completeness unknown" in html
    assert "support is unknown because no valid replicates were available" in html
    assert "not evidence of instability" in html


def test_single_profile_cohort_does_not_display_a_pcoa_as_an_ordination(tmp_path: Path) -> None:
    (tmp_path / "single.svg").write_text(
        '<svg xmlns="http://www.w3.org/2000/svg"></svg>', encoding="utf-8"
    )
    output = write_profile_report(
        {"records_count": 1, "coverage": {"available_profiles": 1},
         "cohorts": [{"cohort_id": "cohort_1", "sample_ids": ["Q1"],
                      "pcoa_figure": "single.svg"}]},
        directory=tmp_path,
    )
    html = output.read_text(encoding="utf-8")

    assert "Profile PCoA is not available" in html
    assert "single.svg" not in html


def test_geography_distinguishes_input_origin_and_full_catalogue(tmp_path: Path) -> None:
    output = write_profile_report(
        {
            "records_count": 3,
            "coverage": {"available_profiles": 2},
            "metadata_geography": [
                {"origin": "local", "country": "Kenya", "region": "Nairobi",
                 "nuts2": "KE01", "count": 1},
                {"origin": "context", "country": "Kenya", "region": "Nairobi",
                 "nuts2": "KE01", "count": 1},
                {"origin": "local", "country": "Unknown", "region": "Unknown",
                 "nuts2": "Unknown", "count": 1},
            ],
            "public_catalogue_geography": [
                {"origin": "context", "country": "Kenya", "region": "Nairobi",
                 "nuts2": "KE01", "count": 27},
            ],
            "genetic_groups": [{"group_id": "G1", "sample_ids": ["Q1", "C1"],
                                "method": "complete linkage", "valid_bootstrap_replicates": 10,
                                "minimum_pair_bootstrap_coassignment": 0.8}],
        },
        directory=tmp_path,
        provenance={
            "coverage": {
                "queries": {"total": 2, "profiles_available": 1, "profiles_missing": 1},
                "context": {"total": 1, "profiles_available": 1, "profiles_missing": 0},
            },
            "catalogue_deduplication": {"raw_records": 40, "retained_units": 30},
            "catalogue_metadata_count": 30,
            "eligible_public_context_count": 10,
            "bounded_public_context_count": 1,
            "profile_limit": 500,
        },
    )
    html = output.read_text(encoding="utf-8")

    assert "Your input genomes" in html and "Public comparison genomes" in html
    assert "2 input genomes · 1 public comparisons" in html
    assert "your 2 input genomes and compares them with 1 additional public genomes" in html
    assert html.index("<h3>Your input genomes") < html.index("<h3>Public comparison genomes")
    assert "Frozen public catalogue" in html
    assert "Public catalogue records before deduplication" in html
    assert "Public context records in profile analysis" in html
    assert "27" in html
    groups_section = html.split('id="groups"', 1)[1].split('id="concentration"', 1)[0]
    assert "descriptive grouping" in groups_section
    assert "describes recurrence across years" not in groups_section


def test_report_renders_numpy_scalar_metrics(tmp_path: Path) -> None:
    output = write_profile_report(
        {"genetic_groups": [{"group_id": "G1", "score": np.float64(0.75)}]},
        directory=tmp_path,
    )
    html = output.read_text(encoding="utf-8")

    assert "0.75" in html


def test_query_only_neighbours_are_summarised_and_expandable(tmp_path: Path) -> None:
    output = write_profile_report(
        {
            "records_count": 23,
            "coverage": {"available_profiles": 23},
            "nearest_neighbours": [
                {"query_id": f"Q{i}", "status": "no_comparable_context"} for i in range(23)
            ],
        },
        directory=tmp_path,
        provenance={"coverage": {"context": {"total": 0}}},
    )
    html = output.read_text(encoding="utf-8")

    assert "No public context profiles were included in this analysis." in html
    assert "Query-level status (23 queries)" in html
    assert html.count('id="Q0"') == 0


def test_nearest_relatives_display_engine_comparison_group_with_legacy_fallback(tmp_path: Path) -> None:
    output = write_profile_report(
        {
            "nearest_neighbours": [
                {"query_id": "Q1", "context_id": "C1", "cohort_id": "cohort_1", "cohort": "obsolete"},
                {"query_id": "Q2", "context_id": "C2", "cohort": "legacy_cohort"},
            ],
        },
        directory=tmp_path,
    )
    html = output.read_text(encoding="utf-8")
    nearest = html.split('id="nearest"', 1)[1].split('id="figures"', 1)[0]

    assert "Comparison group" in nearest
    assert "cohort_1" in nearest
    assert "legacy_cohort" in nearest
    assert "obsolete" not in nearest


def test_rep_section_removed_even_with_historical_persistence_results(tmp_path: Path) -> None:
    analysis = {
        "records_count": 2,
        "coverage": {"available_profiles": 2},
        "temporal_persistence": [{"group_id": "G1", "observed_years": [2020, 2024],
                                  "interpretation": "Persisting strain established."}],
        "time_place_concentration": [{"group_id": "G1", "cells": [
            {"location": "Greece", "year": 2020, "count": 2}]}],
    }
    html = write_profile_report(analysis, directory=tmp_path).read_text()
    assert 'id="persistence"' not in html
    assert "REP-inspired" not in html
    assert "rep-strains" not in html
    assert "Persisting strain established" not in html
    assert 'id="concentration"' in html
    assert 'id="root-to-tip"' in html
    assert analysis["temporal_persistence"][0]["observed_years"] == [2020, 2024]


def test_large_geography_and_network_tables_are_collapsed_and_network_figure_is_first(
    tmp_path: Path,
) -> None:
    (tmp_path / "network.svg").write_text(
        '<svg xmlns="http://www.w3.org/2000/svg"></svg>', encoding="utf-8"
    )
    output = write_profile_report(
        {
            "records_count": 1,
            "coverage": {"available_profiles": 1},
            "public_catalogue_geography": [
                {"country": "Country", "region": f"R{i}", "count": 1} for i in range(87)
            ],
            "cohorts": [{"cohort_id": "C1", "network_figure": "network.svg"}],
            "location_network": [
                {
                    "cohort_id": "C1",
                    "method": "alternate roots",
                    "tested_roots": 9,
                    "edges": [
                        {
                            "source": f"R{i}",
                            "target": f"R{i + 1}",
                            "roots_with_possible_change": 3,
                            "roots_tested": 9,
                            "root_fraction": 1 / 3,
                        }
                        for i in range(72)
                    ],
                }
            ],
        },
        directory=tmp_path,
    )
    html = output.read_text(encoding="utf-8")
    catalogue = html.split('id="geography"', 1)[1].split('id="nearest"', 1)[0]
    network = html.split('id="network"', 1)[1].split('id="cohorts"', 1)[0]

    assert "Frozen public catalogue: 87 records" in catalogue
    assert 'class="coverage-details"><summary>Frozen public catalogue' in catalogue
    assert "Possible location changes (72 edges)" in network
    assert network.index('src="network.svg"') < network.index("Possible location changes")
    assert "Fraction of tested roots with possible change" in network
    assert "not probabilities, confidence scores, or evidence of transmission" in network


def test_neighbour_table_prioritises_formatted_normalized_distance(tmp_path: Path) -> None:
    output = write_profile_report(
        {
            "records_count": 1,
            "coverage": {"available_profiles": 1},
            "nearest_neighbours": [
                {
                    "query_id": "Q1",
                    "context_id": "C1",
                    "status": "matched",
                    "distance": 7 / 628,
                    "allele_differences": 7,
                    "shared_called_loci": 628,
                    "call_overlap": 0.99,
                    "country": "Kenya",
                    "year": 2024,
                }
            ],
        },
        directory=tmp_path,
    )
    html = output.read_text(encoding="utf-8")
    nearest = html.split('id="nearest"', 1)[1].split('id="figures"', 1)[0]

    assert "Normalized distance" in nearest
    assert "0.0111" in nearest
    assert "7 / 628" in nearest
    assert "Kenya" in nearest and "2024" in nearest
    assert nearest.index("Normalized distance") < nearest.index("Allele differences / shared loci")


def test_adaptive_context_report_separates_dataset_pool_from_tree_and_keeps_one_figure_set(tmp_path):
    provenance = {
        "coverage": {"queries": {"total": 9}, "context": {"total": 252}},
        "adaptive_context_selection": {
            "method": "adaptive-cglin-context", "min_context": 20,
            "datasets": [{
                "analysis_dataset": "ST39_CG39", "mlst_st": "39", "clonal_group": "CG39",
                "input_count": 9, "public_cg_pool_count": 787, "selected_context_count": 252,
                "subgroups": [{
                    "input_prefix": [0, 0, 39, 39, 0], "input_count": 7,
                    "context_counts": {"5": 162, "6": 155, "7": 146},
                    "selected_level": 7, "selected_context_count": 146, "limited_context": False,
                }, {
                    "input_prefix": [0, 0, 39, 39, 9], "input_count": 1,
                    "context_counts": {"5": 4, "6": 1, "7": 1},
                    "selected_level": 5, "selected_context_count": 4, "limited_context": True,
                }],
                "selected_context_geography": [
                    {"country": "Unknown", "region": "Unknown", "nuts2": "Unknown", "count": 2, "denominator": 252}
                ],
                "selected_context_years": [{"year": "Unknown", "count": 3, "denominator": 252}],
            }],
        },
    }
    html = write_profile_report({"records_count": 261}, directory=tmp_path, provenance=provenance).read_text()
    assert html.index("Input datasets and their public context") < html.index('id="geography"')
    assert "CG39" in html and "787" in html and "252" in html
    assert "Level 5 context" in html and "Level 7 context" in html
    assert "Limited context" in html and "Minimum met" in html
    assert "actual cgMLST allele differences" in html
    assert "jointly compared loci and distance ties" in html
    assert "before tree subsampling or profile exclusions" in html
    assert "Full public clonal-group pool" in html
    assert "separate figures are not generated for every LIN level" in html
    assert html.count('id="figures"') == 1
    assert "REP-inspired" not in html
    assert "Archivo" in html


def test_adaptive_context_report_escapes_prefix_and_dataset_text(tmp_path):
    html = write_profile_report({}, directory=tmp_path, provenance={
        "adaptive_context_selection": {"min_context": 20, "datasets": [{
            "analysis_dataset": "<script>", "clonal_group": "<script>",
            "subgroups": [{"input_prefix": ["<script>"], "context_counts": {}}],
            "selected_context_years": [{"year": "<script>", "count": 1, "denominator": 1}],
        }]}
    }).read_text()
    assert "<script>" not in html
    assert "&lt;script&gt;" in html


def test_fast_group_index_links_only_existing_contained_reports_and_escapes(tmp_path):
    group = tmp_path / 'ST39_CG39'
    group.mkdir()
    (group / 'profile_report.html').write_text('report')
    outside = tmp_path.parent / 'outside_report.html'
    outside.write_text('outside')
    html = write_fast_group_index(tmp_path, [
        {"species": "<script>", "lineage": "ST39_CG39", "input_count": 9,
         "context_count": 252, "report": "ST39_CG39/profile_report.html"},
        {"lineage": "unsafe", "report": outside},
        {"lineage": "absent", "report": "missing/profile_report.html"},
    ], provenance={"adaptive_context_selection": {"datasets": [{
        "analysis_dataset": "ST39_CG39", "mlst_st": "39", "clonal_group": "CG39",
        "public_cg_pool_count": 787,
    }]}}).read_text()
    assert 'href="ST39_CG39/profile_report.html"' in html
    assert 'outside_report.html' not in html
    assert 'missing/profile_report.html' not in html
    assert '<script>' not in html and '&lt;script&gt;' in html
    assert '787' in html and '252' in html and 'CG39' in html
    assert 'Archivo' in html
    assert '<img' not in html
