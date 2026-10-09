import zipfile
from pathlib import Path

from chronoclade.report import (
    assess_temporal_signal,
    write_fast_lineage_report,
    write_lineage_report,
    write_supporting_bundle,
    write_summary_report,
)


def temporal_result(*, rate: float = 1e-6, p_value: float = 0.01) -> dict[str, object]:
    return {
        "observed": {"rate": rate, "r_squared": 0.65},
        "requested_randomisations": 4,
        "successful_randomisations": 4,
        "p_value_r_squared": p_value,
        "randomised": [
            {"rate": 2e-7, "r_squared": 0.10},
            {"rate": 4e-7, "r_squared": 0.15},
            {"rate": 6e-7, "r_squared": 0.20},
            {"rate": 8e-7, "r_squared": 0.25},
        ],
    }


def full_tree_temporal_result() -> dict[str, object]:
    return {
        "method": "full_tree",
        "observed": {
            "rate": 1e-6,
            "r_squared": 0.5,
            "rate_std": 1e-7,
            "rate_lower_95": 8.04e-7,
            "rate_upper_95": 1.196e-6,
        },
        "requested_randomisations": 2,
        "successful_randomisations": 2,
        "criterion": "cr2",
        "cr1_passed": True,
        "cr2_passed": True,
        "cr1_overlapping_randomisations": 0,
        "cr2_overlapping_randomisations": 0,
        "randomised": [
            {
                "rate": 2e-7,
                "r_squared": 0.1,
                "rate_std": 2e-8,
                "rate_lower_95": 1.608e-7,
                "rate_upper_95": 2.392e-7,
            },
            {
                "rate": 4e-7,
                "r_squared": 0.2,
                "rate_std": 2e-8,
                "rate_lower_95": 3.608e-7,
                "rate_upper_95": 4.392e-7,
            },
        ],
    }


def test_temporal_assessment_reports_supported_and_unsupported_results() -> None:
    supported = assess_temporal_signal(temporal_result(), p_value_threshold=0.05)
    unsupported = assess_temporal_signal(temporal_result(p_value=0.25), p_value_threshold=0.05)

    assert supported["code"] == "supported"
    assert supported["supported"] is True
    assert unsupported["code"] == "not_supported"
    assert unsupported["supported"] is False


def test_lineage_report_contains_visuals_verdict_and_guardrail(tmp_path: Path) -> None:
    (tmp_path / "clock").mkdir()
    (tmp_path / "timetree").mkdir()
    (tmp_path / "clock" / "root_to_tip_regression.svg").write_text(
        "<svg xmlns='http://www.w3.org/2000/svg'/>", encoding="utf-8"
    )
    (tmp_path / "timetree" / "timetree.svg").write_text(
        "<svg xmlns='http://www.w3.org/2000/svg'/>", encoding="utf-8"
    )
    (tmp_path / "timetree" / "molecular_clock.txt").write_text(
        "Root-Tip-Regression:\n --rate:\t1.2e-06 +/- 8e-08 (one std-dev)\n"
        " --chi^2:\t12.4\n --r^2:\t0.71\n",
        encoding="utf-8",
    )
    (tmp_path / "timetree" / "dates.tsv").write_text(
        "#node\tdate\tnumeric date\tlower bound\tupper bound\n"
        "NODE_0000000\t1999-01-01\t1999.0\t1997.5\t2000.5\n"
        "S1\t2020-01-01\t2020.0\t2020.0\t2020.0\n",
        encoding="utf-8",
    )
    (tmp_path / "clock" / "rtt.csv").write_text(
        "#Dates inferred from root-to-tip regression.\n"
        "name, date, root-to-tip distance, clock-deviation\n"
        "S1, 2020, 0.001, 0.0\nS2, 2021, 0.002, 0.0\n",
        encoding="utf-8",
    )
    (tmp_path / "clonal_snp_heatmap.svg").write_text(
        "<svg xmlns='http://www.w3.org/2000/svg'/>", encoding="utf-8"
    )
    (tmp_path / "recombination_map.svg").write_text(
        "<svg xmlns='http://www.w3.org/2000/svg'/>", encoding="utf-8"
    )
    report = {
        "species": "E_coli",
        "lineage": "ST131",
        "sample_count": 42,
        "distinct_dates": 12,
        "first_collection_date": "2020-01",
        "last_collection_date": "2024-06",
        "complete_alignment_sites": 12345,
        "recombination_masking": {
            "inferred_importation_intervals": 8,
            "branches_with_inferred_importation": 3,
            "recombination_removed_sites": 2400,
            "recombination_removed_percent": 1.2,
            "incomplete_only_removed_sites": 400,
            "incomplete_only_removed_percent": 0.2,
            "retained_clonal_sites": 197200,
            "retained_clonal_percent": 98.6,
            "boundary_crossing_intervals": 0,
            "longest_inferred_intervals": [
                {
                    "node": "NODE_12",
                    "alignment_start": 1200,
                    "alignment_end": 3599,
                    "length": 2400,
                    "reference_segments": ["chromosome:1200-3599"],
                    "crosses_reference_record_boundary": False,
                }
            ],
        },
        "temporal_signal": temporal_result(),
        "context": {
            "local_samples": 30,
            "context_samples": 12,
            "context_locations": ["Philippines", "United Kingdom"],
            "manifest_available": True,
            "interpretation": "Contextual genomes are available.",
            "nearest_screening_contexts": [
                {
                    "sample_id": "SAMN1",
                    "country": "Philippines",
                    "collection_date": "2022",
                    "nearest_focal": "F1",
                    "min_ska_distance": 4.0,
                }
            ],
        },
        "public_health": {
            "scenario": {
                "code": "mixed",
                "label": "Consistent with local persistence plus additional introductions",
                "confidence": "moderate",
                "decision_rule": "At least two separated focal-only groups are present.",
                "reasons": ["Two candidate local groups were observed."],
                "recommended_follow_up": ["Review the groups separately."],
                "guardrail": "This does not establish direct transmission.",
                "evidence": [
                    {
                        "id": "clonal_distances",
                        "status": "supported",
                        "finding": "Corrected distances were available.",
                    }
                ],
            },
            "distance_summary": {
                "categories": {
                    "focal_focal": {
                        "comparisons": 3,
                        "minimum": 1,
                        "median": 4,
                        "maximum": 20,
                    }
                }
            },
            "topology": {
                "interpretation": "Candidate groups are review aids.",
                "groups": [
                    {
                        "group_id": "LG1",
                        "sample_count": 2,
                        "first_collection_date": "2022",
                        "last_collection_date": "2024",
                        "within_group_clonal_snps": {
                            "comparisons": 1,
                            "minimum": 1,
                            "median": 1,
                            "maximum": 1,
                        },
                        "nearest_context": {"sample_id": "SAMN1", "clonal_snps": 4},
                    }
                ],
            },
            "patient_sensitivity": {
                "patient_metadata_complete": True,
                "excluded_repeated_patient_samples": ["F2"],
                "interpretation": "The earliest isolate per patient was retained.",
            },
        },
    }

    output = write_lineage_report(report, directory=tmp_path, p_value_threshold=0.05)
    text = output.read_text(encoding="utf-8")

    assert "Root-to-tip permutation screen passed" in text
    assert ".decision.proceed { color:var(--green); }" in text
    assert 'class="decision proceed">READY TO TEST' in text
    assert 'class="decision proceed">CONTINUE TO TEST' in text
    assert 'class="decision proceed">PROCEED' in text
    assert 'class="overall supported"' in text
    assert ".overall.supported { color:var(--green); }" in text
    assert "Consistent with local persistence plus additional introductions" in text
    assert "What does the genomic evidence support?" in text
    assert "What should happen next?" in text
    assert "Does divergence increase with sampling time?" in text
    assert "Is this a coherent, dated lineage dataset?" in text
    assert text.index("Is this a coherent, dated lineage dataset?") < text.index(
        "Does divergence increase with sampling time?"
    )
    assert "workflow's configured screening rule for time scaling" in text
    assert "without preserving genetic or outbreak clusters" in text
    assert text.index("Does divergence increase with sampling time?") < text.index(
        "Is the observed fit stronger than shuffled dates?"
    )
    assert text.index("Is the observed fit stronger than shuffled dates?") < text.index(
        "Estimate the dated phylogeny"
    )
    assert text.index("Estimate the dated phylogeny") < text.index(
        "What pattern is consistent with these genomes?"
    )
    assert "Recombination-filtered genomic distances" in text
    assert "Candidate focal groups" in text
    assert "Patient-level sensitivity" in text
    assert "Which parts of the alignment were inferred as recombinant?" in text
    assert "does not mean that every sampled genome acquired that segment" in text
    assert "Recombination removed" in text
    assert "Incomplete-data removed" in text
    assert "recombination_map.svg" in text
    assert "Longest inferred intervals" in text
    assert "1,200–3,599" in text
    assert "chromosome:1200-3599" in text
    assert "clock/root_to_tip_regression.svg" in text
    assert "date_randomisation.svg" in text
    assert "Time-scaled phylogeny" in text
    assert "Public contextual genomes" in text
    assert "SAMN1" in text
    assert "SKA distances" in text
    assert "does not prove direct transmission" in text
    assert (tmp_path / "date_randomisation.svg").is_file()
    assert (tmp_path / "date_randomisation.png").is_file()
    assert (tmp_path / "date_randomisation.csv").is_file()
    assert (tmp_path / "root_to_tip.png").is_file()
    assert (tmp_path / "timetree_confidence.csv").is_file()
    assert (tmp_path / "node_dates.csv").is_file()
    assert (tmp_path / "supporting_results.zip").is_file()
    assert "Median interval width" in text
    assert "Widest interval" in text
    assert "<i>E coli</i><span>ST131</span>" in text
    assert "E_coli</i>" not in text
    assert "Final dated-tree fit" in text
    assert "Open the full-resolution dated phylogeny" in text
    assert 'role="region" aria-label="Evidence ledger"' in text
    assert "prefers-reduced-motion:reduce" in text
    assert "IntersectionObserver" not in text
    assert 'class="stage active"' not in text
    assert ".stage.active" not in text
    assert text.index(
        "Root-to-tip permutation screen passed", text.index("date_randomisation.svg")
    ) < text.index("Evidence and downloads", text.index("date_randomisation.svg"))


def test_lineage_report_warns_when_root_interval_is_disproportionately_wide(
    tmp_path: Path,
) -> None:
    (tmp_path / "clock").mkdir()
    (tmp_path / "clock" / "root_to_tip_regression.svg").write_text(
        "<svg xmlns='http://www.w3.org/2000/svg'/>", encoding="utf-8"
    )
    (tmp_path / "timetree").mkdir()
    (tmp_path / "timetree" / "timetree.svg").write_text("<svg/>", encoding="utf-8")
    (tmp_path / "timetree" / "molecular_clock.txt").write_text(
        " --rate:\t4.4e-07 +/- 7e-08 (one std-dev)\n --r^2:\t0.11\n",
        encoding="utf-8",
    )
    (tmp_path / "timetree" / "dates.tsv").write_text(
        "#node\tdate\tnumeric date\tlower bound\tupper bound\n"
        "NODE_0000000\t1980-01-01\t1980\t1775\t2002\n"
        "NODE_1\t1995-01-01\t1995\t1985\t2000\n"
        "NODE_2\t1996-01-01\t1996\t1988\t2001\n"
        "NODE_3\t1997-01-01\t1997\t1990\t2002\n",
        encoding="utf-8",
    )
    report = {
        "species": "Escherichia_coli",
        "lineage": "ST131",
        "sample_count": 96,
        "distinct_dates": 15,
        "temporal_signal": temporal_result(),
        "context": {},
        "public_health": {},
    }

    output = write_lineage_report(report, directory=tmp_path, p_value_threshold=0.05)
    text = output.read_text(encoding="utf-8")

    assert "REVIEW ROOT DATE" in text
    assert "Root date is poorly constrained" in text
    assert "spans 227.0 years" in text
    assert "Not applicable" not in text


def test_unsupported_report_omits_dated_tree_visual(tmp_path: Path) -> None:
    report = {
        "species": "Klebsiella_pneumoniae",
        "lineage": "ST15",
        "sample_count": 20,
        "distinct_dates": 5,
        "temporal_signal": temporal_result(p_value=0.25),
    }

    output = write_lineage_report(report, directory=tmp_path, p_value_threshold=0.05)
    text = output.read_text(encoding="utf-8")

    assert "Root-to-tip permutation screen failed" in text
    assert 'class="decision stop">DO NOT TIME-SCALE' in text
    assert 'class="overall not_supported"' in text
    assert "Time-scaled phylogeny" not in text
    assert "cannot distinguish" in text


def test_supporting_bundle_is_reproducible(tmp_path: Path) -> None:
    (tmp_path / "result.csv").write_text("sample,value\nS1,1\n", encoding="utf-8")
    (tmp_path / "phipack_recombination_regions.tsv").write_text("obsolete\n", encoding="utf-8")

    first = write_supporting_bundle(tmp_path).read_bytes()
    (tmp_path / "result.csv").touch()
    second = write_supporting_bundle(tmp_path).read_bytes()

    assert first == second
    with zipfile.ZipFile(tmp_path / "supporting_results.zip") as archive:
        assert "phipack_recombination_regions.tsv" not in archive.namelist()


def test_full_tree_report_explains_cr2_without_root_to_tip_p_value(tmp_path: Path) -> None:
    report = {
        "species": "E_coli",
        "lineage": "ST131",
        "sample_count": 20,
        "distinct_dates": 5,
        "temporal_signal": full_tree_temporal_result(),
        "context": {},
        "public_health": {},
    }

    text = write_lineage_report(report, directory=tmp_path, p_value_threshold=0.05).read_text(
        encoding="utf-8"
    )

    assert "Full TreeTime randomisation screen passed" in text
    assert "strict CR2 rule" in text
    assert "CR2 overlaps" in text
    assert "Empirical p" not in text


def test_fast_report_separates_recombination_root_to_tip_and_permutation_stages(
    tmp_path: Path,
) -> None:
    (tmp_path / "clock").mkdir()
    (tmp_path / "clock" / "root_to_tip_regression.svg").write_text(
        "<svg xmlns='http://www.w3.org/2000/svg'/>", encoding="utf-8"
    )
    report = {
        "species": "E_coli",
        "lineage": "ST131",
        "sample_count": 20,
        "distinct_dates": 5,
        "temporal_signal": temporal_result(),
        "recombination": {
            "tested_core_blocks": 3,
            "significant_blocks": 1,
            "recombination_detected": True,
            "minimum_p_value": 0.001,
            "boundary_rule": "No block crossed a reference-contig join.",
            "localisation_limit": (
                "A significant block does not localise a tract and no alignment sites were masked."
            ),
        },
    }

    text = write_fast_lineage_report(report, directory=tmp_path, p_value_threshold=0.05).read_text(
        encoding="utf-8"
    )

    assert "FAST SCREEN" in text
    assert "Is there a fast PHI signal that warrants recombination correction?" in text
    assert "Does genetic divergence increase with sampling time?" in text
    assert "Is the observed fit stronger than shuffled dates?" in text
    assert "no alignment sites were masked" in text
    assert "no biological meaning" in text
    assert "RUN THE FULL ANALYSIS" in text
    assert 'class="decision proceed">PHI SCREEN POSITIVE' in text
    assert 'class="status proceed">RUN THE FULL ANALYSIS' in text
    assert "clock/root_to_tip_regression.svg" in text
    assert "dated phylogeny" in text
    assert "Time-scaled phylogeny" not in text
    assert "Interpret" not in text
    assert "Not applicable" not in text
    assert "phipack.filtered.fasta" not in text
    assert "Masked regions" not in text


def test_summary_report_links_lineage_reports(tmp_path: Path) -> None:
    metrics = temporal_result()
    summary = {
        "lineages": [
            {
                "species": "E_coli",
                "lineage": "ST131",
                "slug": "E_coli__ST131",
                "sample_count": 42,
                "temporal_status": "supported",
                "temporal_signal": metrics,
            }
        ]
    }

    output = write_summary_report(summary, output=tmp_path)
    text = output.read_text(encoding="utf-8")

    assert "E_coli__ST131/report.html" in text
    assert "Supported" in text


def test_embedded_country_svg_keeps_local_glyph_references(tmp_path):
    from chronoclade.report import _context_geography_visual

    directory = tmp_path / "context_geography"
    directory.mkdir()
    (directory / "fragment.html").write_text(
        '<a href="country_composition.csv">Data</a>'
        '<svg><defs><path id="glyph" /></defs><use xlink:href="#glyph" /></svg>'
    )
    embedded = _context_geography_visual(tmp_path)
    assert 'href="context_geography/country_composition.csv"' in embedded
    assert 'xlink:href="#glyph"' in embedded
    assert "context_geography/#glyph" not in embedded
