import numpy as np
from pathlib import Path

from chronoclade.profile_report import (
    write_corrected_report,
    write_profile_report,
    write_stage_index,
)


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
    assert html.index("id=\"persistence\"") < html.index("id=\"concentration\"")
    assert "observation across collection years" in html
    assert "concentration in a particular time and place" in html
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

    assert "Analysed input metadata" in html
    assert "Query set" in html and "Public context" in html
    assert "Frozen public catalogue" in html
    assert "Public catalogue records before deduplication" in html
    assert "Public context records in profile analysis" in html
    assert "27" in html
    groups_section = html.split('id="groups"', 1)[1].split('id="persistence"', 1)[0]
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


def test_persistence_summary_does_not_claim_recurrence_for_one_year(tmp_path: Path) -> None:
    output = write_profile_report(
        {
            "records_count": 2,
            "coverage": {"available_profiles": 2},
            "temporal_persistence": [
                {
                    "group_id": "G1",
                    "observed_years": [2024],
                    "observed_year_count": 1,
                    "interpretation": "Observed across collection years; gaps do not establish continuous persistence.",
                },
                {"group_id": "G2", "observed_years": [], "observed_year_count": 0},
            ],
            "time_place_concentration": [
                {"group_id": "C1", "cells": [{"location": "Kenya", "year": 2024, "count": 2}]}
            ],
        },
        directory=tmp_path,
    )
    html = output.read_text(encoding="utf-8")
    persistence = html.split('id="persistence"', 1)[1].split('id="concentration"', 1)[0]

    assert "Observed in one collection year (2024); persistence cannot be assessed." in persistence
    assert "year-to-year persistence cannot be assessed" in persistence
    assert "Observed across collection years; gaps do not establish continuous persistence." not in persistence


def test_rep_framework_does_not_classify_sampled_genomes_as_cdc_rep_strains(tmp_path: Path) -> None:
    output = write_profile_report(
        {
            "species": "Klebsiella pneumoniae",
            "records_count": 3,
            "coverage": {"available_profiles": 3},
            "temporal_persistence": [{
                "group_id": "G1", "observed_years": [2020, 2022, 2024],
                "interpretation": "Persisting strain established by three sampled years.",
            }],
        },
        directory=tmp_path,
    )
    html = output.read_text(encoding="utf-8")
    persistence = html.split('id="persistence"', 1)[1].split('id="concentration"', 1)[0]

    assert "Recurrence and persistence in sampled genomes (REP-inspired)" in persistence
    assert "https://www.cdc.gov/foodborne-outbreaks/php/rep-strains/index.html" in persistence
    assert "Reoccurring" in persistence and "Emerging" in persistence and "Persisting" in persistence
    assert persistence.count("Not assessed") == 3
    assert "surveillance denominators" in persistence
    assert "outbreaks or quiet periods" in persistence
    assert "consistent illness or uninterrupted circulation" in persistence
    assert "does not assign official CDC REP designations, including to Klebsiella" in persistence
    assert "Observed across 3 collection years" in persistence
    assert "Persisting strain established" not in persistence
    assert "Stable-group results" not in html
    assert "Genetic group stability refers to locus-bootstrap co-assignment" in html


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
