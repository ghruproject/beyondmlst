"""Shared diagnostics retain date uncertainty and categorical distance meaning."""

from dataclasses import replace
import json

import numpy as np
import pytest

from chronoclade.datasets import AlleleMatrix, LocusCatalogue, PreparedDataset, samples_from_records
from chronoclade.temporal_diagnostics import (
    cgmlst_root_to_tip,
    date_regression,
    fixed_callable_distances,
    sample_date_interval,
)


def regress(points):
    return date_regression(
        points,
        distance_field="distance",
        distance_units="example units",
        label="Example diagnostic",
    )


def points(distances, dates=("2000-01-01", "2001-01-01", "2002-01-01")):
    return [
        {
            "sample_id": f"sample-{index}",
            "collection_date": day,
            "distance": distance,
            "role": "input",
            "country": "UK",
        }
        for index, (day, distance) in enumerate(zip(dates, distances, strict=True))
    ]


def dataset(calls, dates=None, *, complete=True):
    ids = tuple(calls)
    dates = dates or [f"200{index}" for index in range(len(ids))]
    records = [
        {
            "sample_id": ident,
            "label": "Readable " + ident,
            "origin": "local" if index == 0 else "context",
            "species": "E. coli",
            "collection_date": dates[index],
            "country": "UK",
        }
        for index, ident in enumerate(ids)
    ]
    catalogue = LocusCatalogue(
        "example:cgmlst", "v1", ("a", "b", "c"), "frozen.tsv", complete=complete
    )
    matrix = AlleleMatrix.from_profiles(catalogue, ids, calls)
    return PreparedDataset(samples_from_records(records), profiles=(matrix,))


def test_negative_slope_is_retained_and_predictions_are_centred():
    result = regress(points([3, 2, 1]))
    assert result["status"] == "fitted"
    assert result["slope"] == pytest.approx(-1)
    assert result["reference_year"] == 2001
    assert result["centered_intercept"] == 2
    assert result["pearson_r"] == pytest.approx(-1)
    assert result["r_squared"] == pytest.approx(1)
    assert [point["residual"] for point in result["points"]] == pytest.approx([0, 0, 0])
    assert result["points"][0]["role"] == "input"
    assert result["points"][0]["country"] == "UK"
    assert result["slope_units"] == "example units per year"
    assert "not assumed independent" in result["interpretation"]
    assert "p_value" not in result and "confidence_interval" not in result
    json.dumps(result, allow_nan=False)


def test_constant_distances_have_zero_slope_and_null_undefined_statistics():
    result = regress(points([0.25, 0.25, 0.25]))
    assert result["status"] == "fitted"
    assert result["slope"] == 0
    assert result["pearson_r"] is None
    assert result["r_squared"] is None
    assert result["intercept"] == 0.25
    assert [point["predicted"] for point in result["points"]] == [0.25] * 3
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize(
    "dates",
    [
        ("2000", "2000", "2000"),
        ("2000", "2000", "2001"),
    ],
)
def test_three_distinct_dates_are_required(dates):
    result = regress(points([1, 2, 3], dates))
    assert result["status"] == "insufficient_data"
    assert result["slope"] is None
    assert len(result["points"]) == 3


def test_canonical_intervals_take_precedence_and_raw_dates_remain_available():
    point = {
        "sample_id": "interval",
        "collection_date": "unusable raw date",
        "date_start": "2000-02-01",
        "date_end": "2002-12-31",
        "date_precision": "interval",
        "distance": 1,
    }
    interval = sample_date_interval(point)
    assert interval["status"] == "valid"
    assert interval["precision"] == "interval"
    assert interval["start"] == "2000-02-01"
    assert interval["end"] == "2002-12-31"
    assert interval["midpoint_year"] == (interval["year_min"] + interval["year_max"]) / 2
    result = regress([point])
    assert result["points"][0]["collection_date"] == "unusable raw date"
    assert result["points"][0]["date_precision"] == "interval"
    assert sample_date_interval({"collection_date": "2000-02"})["end"] == "2000-02-29"
    assert sample_date_interval({"collection_date": "2000/2001"})["precision"] == "interval"


def test_date_and_distance_exclusions_are_audited_and_finite_json():
    rows = points([1, 2, 3]) + [
        {"sample_id": "missing", "distance": 1},
        {"sample_id": "invalid", "collection_date": "2001-02-29", "distance": 1},
        {"sample_id": "future", "collection_date": "9999", "distance": 1},
        {"sample_id": "not-finite", "collection_date": "2000", "distance": float("nan")},
        {"sample_id": "no-distance", "collection_date": "2000"},
    ]
    result = regress(rows)
    assert result["status"] == "fitted"
    assert {row["sample_id"]: row["reason"] for row in result["exclusions"]} == {
        "missing": "missing_collection_date",
        "invalid": "invalid_collection_date",
        "future": "future_collection_date",
        "not-finite": "nonfinite_or_missing_distance",
        "no-distance": "nonfinite_or_missing_distance",
    }
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize(
    "point",
    [
        {"date_start": "2001-01-01", "date_end": "2000-12-31", "date_precision": "interval"},
        {"date_start": "2000-01-01", "date_end": "2000-01-02", "date_precision": "day"},
        {"date_start": "2000-02-02", "date_end": "2000-02-29", "date_precision": "month"},
        {"date_start": "2000-01-01", "date_end": "2000-12-30", "date_precision": "year"},
        {"collection_date": "2000-01-01garbage"},
        {"collection_date": "2001/2000"},
    ],
)
def test_inconsistent_or_malformed_intervals_are_invalid(point):
    assert sample_date_interval(point)["status"] == "invalid"


def test_large_offsets_still_fit_a_centred_line():
    result = regress(points([10**12, 10**12 + 1, 10**12 + 2]))
    assert result["slope"] == pytest.approx(1)
    assert [point["residual"] for point in result["points"]] == pytest.approx([0, 0, 0])
    json.dumps(result, allow_nan=False)


def test_fixed_loci_use_categories_and_a_shared_denominator():
    data = dataset(
        {
            "PW_full_source_z": {"a": "99999", "b": "same", "c": "same"},
            "PW_full_source_a": {"a": "2", "b": "same", "c": None},
            "PW_full_source_c": {"a": "hash", "b": "different", "c": "same"},
        }
    )
    ids, loci, distances = fixed_callable_distances(data.profiles[0])
    assert ids == ("PW_full_source_a", "PW_full_source_c", "PW_full_source_z")
    assert loci == ("a", "b")
    np.testing.assert_array_equal(distances, [[0, 1, 0.5], [1, 0, 1], [0.5, 1, 0]])
    result = cgmlst_root_to_tip(data)
    cohort = result["cohorts"][0]
    assert cohort["root"] == "PW_full_source_a"
    assert cohort["n_callable_loci"] == 2
    assert cohort["scheme_loci"] == 3
    assert set(cohort["sample_ids"]) == set(data.sample_ids)
    assert all(ident in cohort["tree_newick"] for ident in data.sample_ids)
    assert cohort["diagnostic"]["slope_units"] == "allele-distance units on 2 fixed callable loci per year"
    assert (
        next(row for row in cohort["samples"] if row["sample_id"] == cohort["root"])["root_to_tip"]
        == 0
    )
    json.dumps(result, allow_nan=False)


def test_guide_tree_and_root_do_not_depend_on_dates_or_row_order():
    calls = {"z": {"a": "1", "b": "1"}, "a": {"a": "2", "b": "1"}, "c": {"a": "2", "b": "3"}}
    first = cgmlst_root_to_tip(dataset(calls))["cohorts"][0]
    reversed_calls = dict(reversed(list(calls.items())))
    second = cgmlst_root_to_tip(dataset(reversed_calls, ["2019", "2017", "2018"]))["cohorts"][0]
    assert first["tree_newick"] == second["tree_newick"]
    assert first["root"] == second["root"] == "a"
    assert {row["sample_id"]: row["root_to_tip"] for row in first["samples"]} == {
        row["sample_id"]: row["root_to_tip"] for row in second["samples"]
    }


def test_explicit_root_retains_full_id_and_unknown_reference_fails():
    data = dataset({"a": {"a": "1"}, "z": {"a": "2"}})
    cohort = cgmlst_root_to_tip(data, reference_sample_id="z")["cohorts"][0]
    assert cohort["root"] == "z" and cohort["root_selection"] == "explicit_reference"
    with pytest.raises(ValueError, match="Unknown reference"):
        cgmlst_root_to_tip(data, reference_sample_id="missing")


def test_undated_profiles_keep_distances_and_separate_date_exclusions():
    data = dataset({"a": {"a": "1"}, "b": {"a": "2"}, "c": {"a": "3"}}, ["2000", "", "invalid"])
    cohort = cgmlst_root_to_tip(data)["cohorts"][0]
    assert len(cohort["samples"]) == 3
    assert all(row["root_to_tip"] is not None for row in cohort["samples"])
    assert len(cohort["diagnostic"]["points"]) == 1
    assert len(cohort["diagnostic"]["exclusions"]) == 2
    assert cohort["samples"][1]["date_interval"]["status"] == "missing"


def test_no_common_loci_is_audited_without_imputed_distances():
    data = dataset({"a": {"a": "1"}, "b": {"b": "2"}, "null": {}}, complete=False)
    result = cgmlst_root_to_tip(data)
    assert result["status"] == "unavailable"
    assert result["cohorts"][0]["reason"] == "no_common_callable_loci"
    assert result["cohorts"][0]["catalogue_complete"] is False
    assert result["cohorts"][0]["tree_newick"] is None
    assert result["exclusions"] == [
        {
            "sample_id": "null",
            "scheme_id": "example:cgmlst",
            "reason": "profile_has_no_callable_loci",
        }
    ]


def test_scheme_and_species_scopes_remain_separate():
    data = dataset({"a": {"a": "1"}, "b": {"a": "2"}, "c": {"a": "3"}})
    other = replace(
        data.profiles[0], catalogue=replace(data.profiles[0].catalogue, scheme_version="v2")
    )
    data = replace(
        data,
        profiles=(data.profiles[0], other),
        samples=tuple(
            {**sample, "species": "Other"} if sample["sample_id"] == "c" else sample
            for sample in data.samples
        ),
    )
    result = cgmlst_root_to_tip(data)
    assert len(result["cohorts"]) == 4
    assert all(
        len({sample["species"] for sample in cohort["samples"]}) == 1
        for cohort in result["cohorts"]
    )
    explicit = cgmlst_root_to_tip(data, reference_sample_id="a")
    assert sum(cohort["reason"] == "reference_not_in_cohort" for cohort in explicit["cohorts"]) == 2


def test_dataset_without_profiles_is_honestly_unavailable():
    data = replace(dataset({"a": {"a": "1"}}), profiles=())
    result = cgmlst_root_to_tip(data)
    assert result["status"] == "unavailable" and result["cohorts"] == []
    assert result["exclusions"] == [{"sample_id": "a", "reason": "profile_unavailable"}]


def test_stable_sample_id_cannot_be_confused_with_generated_internal_node():
    data = dataset({"Inner1": {"a": "1"}, "b": {"a": "2"}, "c": {"a": "3"}})
    cohort = cgmlst_root_to_tip(data)["cohorts"][0]
    assert cohort["root"] == "Inner1"
    by_id = {row["sample_id"]: row["root_to_tip"] for row in cohort["samples"]}
    assert by_id == {"Inner1": 0, "b": 1, "c": 1}
    assert cohort["tree_newick"].count("Inner1") == 1
