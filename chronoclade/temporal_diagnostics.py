"""Descriptive date diagnostics shared by embedding and categorical-profile routes.

Date interval midpoints are an explicit approximation. Neither a regression nor
an arbitrary, date-independent NJ root establishes temporal signal or a clock.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import io
import math

import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from chronoclade.datasets.model import AlleleMatrix, PreparedDataset
from chronoclade.metadata_dates import sample_date_interval


_NOTE = (
    "Exploratory regression using collection-date interval midpoints; not a temporal-signal "
    "gate, molecular clock, substitution rate, or significance test. Shared ancestry and "
    "sampling structure can produce association; observations are not assumed independent."
)


def date_regression(
    points: Sequence[Mapping],
    *,
    distance_field: str,
    distance_units: str,
    label: str,
) -> dict:
    """Fit a centred ordinary least-squares line with an explicit date audit.

    Three distinct interval midpoints are required. Undefined correlation/R²
    for constant distances are null; a fitted constant slope is zero. There are
    no inferential intervals, p-values, independent-tip claims, or pass/fail gate.
    """
    accepted, exclusions = [], []
    for point in points:
        interval = sample_date_interval(point)
        reason = interval["reason"]
        if reason is None:
            try:
                value = float(point[distance_field])
                if not math.isfinite(value):
                    raise ValueError("nonfinite distance")
            except (KeyError, ValueError, TypeError, OverflowError):
                reason = "nonfinite_or_missing_distance"
        if reason:
            exclusions.append(
                {
                    "sample_id": point.get("sample_id"),
                    "collection_date": point.get("collection_date"),
                    "date_start": point.get("date_start"),
                    "date_end": point.get("date_end"),
                    "date_precision": point.get("date_precision") or interval["precision"],
                    "reason": reason,
                }
            )
        else:
            accepted.append({**point, **interval, distance_field: value})
    years = np.asarray([point["midpoint_year"] for point in accepted], dtype=float)
    values = np.asarray([point[distance_field] for point in accepted], dtype=float)
    distinct = len(np.unique(years))
    result = {
        "label": label,
        "distance_field": distance_field,
        "distance_units": distance_units,
        "slope_units": f"{distance_units} per year",
        "status": "insufficient_data",
        "reason": "fewer_than_three_distinct_dated_points",
        "n_points": len(accepted),
        "n_distinct_dates": distinct,
        "slope": None,
        "intercept": None,
        "centered_intercept": None,
        "reference_year": None,
        "pearson_r": None,
        "r_squared": None,
        "points": accepted,
        "exclusions": exclusions,
        "interpretation": _NOTE,
    }
    if len(accepted) < 3 or distinct < 3:
        return result
    reference = float(np.mean(years))
    # Taking the mean after scaling also handles finite values near float limits.
    value_scale = float(np.max(np.abs(values)))
    center = float(np.mean(values / value_scale) * value_scale) if value_scale else 0.0
    x, y = years - reference, values - center
    # Scale both centred axes to avoid overflowing sums of squares.
    x_scale, y_scale = float(np.max(np.abs(x))), float(np.max(np.abs(y)))
    scaled_x = x / x_scale
    if y_scale == 0:
        slope, correlation, r_squared = 0.0, None, None
    else:
        scaled_y = y / y_scale
        xx = float(scaled_x @ scaled_x)
        yy = float(scaled_y @ scaled_y)
        xy = float(scaled_x @ scaled_y)
        slope = (xy / xx) * (y_scale / x_scale)
        correlation = float(np.clip(xy / math.sqrt(xx * yy), -1, 1))
        r_squared = correlation**2
    predicted = center + slope * x
    intercept = center - slope * reference
    if (
        not math.isfinite(slope)
        or not math.isfinite(intercept)
        or not np.all(np.isfinite(predicted))
    ):
        result["reason"] = "numerically_unrepresentable_regression"
        return result
    result.update(
        status="fitted",
        reason=None,
        slope=slope,
        intercept=intercept,
        centered_intercept=center,
        reference_year=reference,
        pearson_r=correlation,
        r_squared=r_squared,
    )
    for point, prediction in zip(accepted, predicted, strict=True):
        point["predicted"] = float(prediction)
        residual = point[distance_field] - point["predicted"]
        point["residual"] = residual if math.isfinite(residual) else None
    return result


def fixed_callable_distances(matrix: AlleleMatrix, *, sample_ids=None):
    """Return sorted IDs, a fixed common callable locus set, and mismatch fractions.

    Categorical vocabulary codes are compared only for equality. Missing calls
    never count as matches and the denominator never changes between pairs.
    """
    ids = tuple(sorted(matrix.sample_ids if sample_ids is None else sample_ids))
    positions = {ident: index for index, ident in enumerate(matrix.sample_ids)}
    codes = matrix.codes[[positions[ident] for ident in ids]]
    common = np.all(codes != 0, axis=0) if ids else np.zeros(len(matrix.catalogue.loci), bool)
    loci = tuple(
        locus for locus, included in zip(matrix.catalogue.loci, common, strict=True) if included
    )
    distances = np.zeros((len(ids), len(ids)), dtype=float)
    if loci:
        called = codes[:, common]
        for index in range(len(ids)):
            distances[index, :index] = np.mean(called[:index] != called[index], axis=1)
            distances[:index, index] = distances[index, :index]
    return ids, loci, distances


def cgmlst_root_to_tip(dataset: PreparedDataset, *, reference_sample_id=None) -> dict:
    """Build independent guide NJ cohorts and descriptive root-to-tip regressions.

    Profiles are separated by scheme/database and species. The root is an
    explicit reference or lexicographically first full stable ID, never selected
    to strengthen date association. Negative NJ branches are retained and counted.
    """
    dataset.validate()
    samples = {sample["sample_id"]: sample for sample in dataset.samples}
    if reference_sample_id is not None and reference_sample_id not in samples:
        raise ValueError(f"Unknown reference sample ID: {reference_sample_id}")
    cohorts, exclusions, represented = [], [], set()
    for matrix in sorted(
        dataset.profiles,
        key=lambda item: (
            item.catalogue.scheme_id,
            item.catalogue.scheme_version,
            item.catalogue.database_version or "",
            item.catalogue.database_sha256 or "",
        ),
    ):
        species_groups = {}
        for index, ident in enumerate(matrix.sample_ids):
            if not np.any(matrix.codes[index]):
                exclusions.append(
                    {
                        "sample_id": ident,
                        "scheme_id": matrix.catalogue.scheme_id,
                        "reason": "profile_has_no_callable_loci",
                    }
                )
                continue
            species = str(samples[ident].get("species") or "").casefold().replace("_", " ")
            species_groups.setdefault(species, []).append(ident)
            represented.add(ident)
        for species, group in sorted(species_groups.items()):
            ids, loci, distances = fixed_callable_distances(matrix, sample_ids=group)
            cohort = {
                "cohort_id": f"cohort_{len(cohorts) + 1}",
                "scheme_id": matrix.catalogue.scheme_id,
                "scheme_version": matrix.catalogue.scheme_version,
                "database_version": matrix.catalogue.database_version,
                "database_sha256": matrix.catalogue.database_sha256,
                "species": species,
                "sample_ids": list(ids),
                "n_profiles": len(ids),
                "n_callable_loci": len(loci),
                "callable_loci": list(loci),
                "scheme_loci": len(matrix.catalogue.loci),
                "catalogue_complete": matrix.catalogue.complete,
                "distance_units": "cgMLST mismatch fraction",
                "status": "unavailable",
                "reason": None,
                "root": reference_sample_id or ids[0],
                "root_selection": "explicit_reference"
                if reference_sample_id
                else "lexicographic_sample_id",
                "tree_newick": None,
                "negative_branch_count": 0,
                "diagnostic": None,
                "samples": [
                    {
                        **samples[ident],
                        "root_to_tip": None,
                        "date_interval": sample_date_interval(samples[ident]),
                        "n_callable_loci": len(loci),
                    }
                    for ident in ids
                ],
            }
            cohorts.append(cohort)
            if reference_sample_id is not None and reference_sample_id not in ids:
                cohort["reason"] = "reference_not_in_cohort"
                continue
            if len(ids) < 2:
                cohort["reason"] = "fewer_than_two_profiles"
                continue
            if not loci:
                cohort["reason"] = "no_common_callable_loci"
                continue
            tree = DistanceTreeConstructor().nj(
                DistanceMatrix(
                    list(ids), [distances[index, : index + 1].tolist() for index in range(len(ids))]
                )
            )
            tips = {tip.name: tip for tip in tree.get_terminals()}
            # Stable IDs may themselves be named "Inner1"; target actual tips,
            # never the generated internal-node name matched first by Bio.Phylo.
            tree.root_with_outgroup(tips[cohort["root"]])
            for node in tree.get_nonterminals():
                node.name = None
            cohort["negative_branch_count"] = sum(
                (node.branch_length or 0) < 0 for node in tree.find_clades()
            )
            buffer = io.StringIO()
            Phylo.write(tree, buffer, "newick", format_branch_length="%.17g")
            cohort["tree_newick"] = buffer.getvalue().strip()
            points = [
                {
                    **samples[ident],
                    "root_to_tip": float(tree.distance(tips[ident])),
                    "n_callable_loci": len(loci),
                }
                for ident in ids
            ]
            cohort["samples"] = [
                {**point, "date_interval": sample_date_interval(point)} for point in points
            ]
            diagnostic = date_regression(
                points,
                distance_field="root_to_tip",
                distance_units=cohort["distance_units"],
                label="cgMLST guide NJ root-to-tip",
            )
            cohort.update(status="available", diagnostic=diagnostic)
    for ident in sorted(set(samples) - represented):
        if not any(item["sample_id"] == ident for item in exclusions):
            exclusions.append({"sample_id": ident, "reason": "profile_unavailable"})
    return {
        "method": "cgmlst_nj_root_to_tip",
        "status": "available"
        if any(cohort["status"] == "available" for cohort in cohorts)
        else "unavailable",
        "cohorts": cohorts,
        "exclusions": exclusions,
        "interpretation": (
            "Date-independent arbitrary tip-root NJ guide tree on a fixed set of loci callable "
            "in every included profile. Distances are categorical mismatch fractions, not SNP "
            "distances or substitutions per site. Negative NJ branches are retained. " + _NOTE
        ),
    }
