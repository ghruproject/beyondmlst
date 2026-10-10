from __future__ import annotations

import csv
from html.parser import HTMLParser
import json
from pathlib import Path
import subprocess
import sys

import pytest

from chronoclade.esm2.temporal_report import write_temporal_report
from chronoclade.location_network.colours import country_palette
from chronoclade.report_components.styles import report_styles
from chronoclade.temporal_diagnostics import date_regression, sample_date_interval


class ReportLinks(HTMLParser):
    def __init__(self):
        super().__init__()
        self.files = []
        self.anchors = []
        self.ids = []

    def handle_starttag(self, tag, attrs):
        values = dict(attrs)
        if "id" in values:
            self.ids.append(values["id"])
        for key in ("href", "src"):
            value = values.get(key, "")
            if value.startswith("#"):
                self.anchors.append(value[1:])
            elif value:
                self.files.append(value)


def diagnostics(*, dates=None, constant=False):
    dates = dates if dates is not None else ["2005", "2008-03", "2011-06-02", None]
    samples = [
        {
            "sample_id": f"stable:{index}",
            "label": f"Isolate {index}",
            "role": "input" if index % 2 else "context",
            "country": "Greece" if index < 3 else "India",
            "collection_date": date,
            "distance": 0.1 if constant else index * 0.025,
        }
        for index, date in enumerate(dates, start=1)
    ]
    regression = date_regression(
        samples,
        distance_field="distance",
        distance_units="mean locus cosine distance",
        label="protein",
    )
    embedding = {
        "regression": regression,
        "samples": [
            {**sample, "date_interval": sample_date_interval(sample)} for sample in samples
        ],
        "reference": {"sample_id": "stable:1", "selection": "explicit", "date_independent": True},
        "panel": {"loci": ["locus_A", "locus_B"], "count": 2, "selection": "explicit-fixed"},
        "counts": {
            "total_samples": 5,
            "mapped_samples": 4,
            "eligible_samples": 4,
            "dated_points": len(regression["points"]),
            "excluded_samples": 1,
        },
        "provenance": {
            "model": {
                "model": "esm2_t6_8M_UR50D",
                "representation_layer": 6,
                "runtime_versions": {"torch": "2.6.0", "fair-esm": "2.0.0"},
            }
        },
        "exclusions": regression["exclusions"],
    }
    cg_samples = [{**sample, "root_to_tip": sample["distance"] * 2} for sample in samples]
    cg_regression = date_regression(
        cg_samples,
        distance_field="root_to_tip",
        distance_units="cgMLST mismatch fraction",
        label="cgmlst",
    )
    cgmlst = {
        "cohorts": [
            {
                "cohort_id": "cohort_1",
                "scheme_id": "test-scheme",
                "scheme_version": "1",
                "root": "stable:1",
                "root_selection": "explicit_reference",
                "n_profiles": 4,
                "n_callable_loci": 2,
                "scheme_loci": 3,
                "negative_branch_count": 1,
                "distance_units": "cgMLST mismatch fraction",
                "diagnostic": cg_regression,
                "samples": [
                    {**sample, "date_interval": sample_date_interval(sample)}
                    for sample in cg_samples
                ],
            }
        ],
        "exclusions": [],
    }
    return embedding, cgmlst


def test_temporal_report_uses_existing_style_and_portable_evidence(tmp_path: Path):
    embedding, cgmlst = diagnostics()
    output = write_temporal_report(embedding, cgmlst, tmp_path)
    html = output.read_text()
    assert "<style>" + report_styles() + "</style>" in html
    assert "cgMLST mismatch fraction/year" in html
    assert "mean locus cosine distance/year" in html
    assert "esm2_t6_8M_UR50D" in html and "2.6.0" in html
    assert "Protein embeddings cannot distinguish synonymous DNA changes" in html
    assert "2008-03-01 to 2008-03-31" in html
    assert "missing_collection_date" in html
    parser = ReportLinks()
    parser.feed(html)
    assert set(parser.anchors) <= set(parser.ids)
    assert all((tmp_path / path).is_file() for path in parser.files)
    assert all((tmp_path / path).stat().st_size for path in parser.files)
    saved = json.loads((tmp_path / "temporal_diagnostics.json").read_text())
    assert saved["location_colours"] == country_palette(["Greece", "India"])
    assert saved["embedding"] == embedding
    with (tmp_path / "temporal_points.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 8  # Undated samples retained for both measurements.
    assert len({row["panel"] for row in rows}) == 2
    assert sum(bool(row.get("residual")) for row in rows) == 6
    protein_svg = (tmp_path / "protein_date_distance.svg").read_text()
    assert "Collection date (decimal year; interval midpoint)" in protein_svg
    assert "Samples without an eligible date (no date assigned)" in protein_svg
    for colour in saved["location_colours"].values():
        assert colour in protein_svg
        assert colour in (tmp_path / "cgmlst_date_distance_1.svg").read_text()


@pytest.mark.parametrize("dates", [[None, None, None], ["2005", "2005", "2008"]])
def test_temporal_report_retains_points_when_fit_is_unavailable(tmp_path: Path, dates):
    embedding, cgmlst = diagnostics(dates=dates)
    output = write_temporal_report(embedding, cgmlst, tmp_path)
    html = output.read_text()
    assert "Fit unavailable" in html
    assert "fewer_than_three_distinct_dated_points" in html
    with (tmp_path / "protein_date_distance_points.csv").open(newline="") as handle:
        assert len(list(csv.DictReader(handle))) == len(dates)
    assert (tmp_path / "protein_date_distance.png").is_file()


def test_temporal_report_marks_constant_distances_and_unavailable_cohorts(tmp_path: Path):
    embedding, cgmlst = diagnostics(constant=True)
    cgmlst["cohorts"].append(
        {
            "cohort_id": "no-common-loci",
            "diagnostic": None,
            "status": "unavailable",
            "reason": "no_common_callable_loci",
            "samples": [{"sample_id": "missing:1", "root_to_tip": None}],
        }
    )
    html = write_temporal_report(embedding, cgmlst, tmp_path).read_text()
    assert "Constant distance" in html
    assert "Distances are constant" in html
    assert "Undefined" in html
    assert "no_common_callable_loci" in html
    assert "no-common-loci" in html
    assert (tmp_path / "cgmlst_date_distance_2.svg").is_file()


def test_temporal_report_escapes_labels_and_can_report_no_baseline(tmp_path: Path):
    embedding, _ = diagnostics()
    embedding["samples"][0]["label"] = '<script>alert("x")</script>'
    embedding["regression"]["points"][0]["label"] = '<script>alert("x")</script>'
    html = write_temporal_report(embedding, {"status": "unavailable"}, tmp_path).read_text()
    assert '<script>alert("x")</script>' not in html
    assert "&lt;script&gt;" in html
    assert "No comparable cgMLST cohort was available." in html


def test_temporal_report_import_does_not_load_matplotlib_or_model_runtime():
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; import chronoclade.esm2.temporal_report; "
            "assert 'matplotlib' not in sys.modules; assert 'torch' not in sys.modules; "
            "assert 'esm' not in sys.modules",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
