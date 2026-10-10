"""Clustered date permutation on a frozen corrected biological tree.

Murray et al. 2016 (doi:10.1111/2041-210X.12466): maximal monophyletic
clusters sharing a sampling date are exchangeability units. The phylogeny is
inferred without dates. Raw date strings (including partial-date precision) are
permuted between clusters; no midpoint or invented date is introduced.
"""

from concurrent.futures import ThreadPoolExecutor
import math
import random
import tempfile
from pathlib import Path

import numpy as np
from Bio import Phylo

from chronoclade.artifacts import write_json
from chronoclade.temporal import _write_dates, _run_randomised_clock, parse_clock
from chronoclade.temporal_report import assess_temporal_signal


def date_clusters(tree, samples):
    """Find disjoint maximal date-homogeneous clades; undated tips are omitted."""
    by_id = {sample.sample_id: sample for sample in samples}
    clusters = []

    def visit(node):
        members = [tip.name for tip in node.get_terminals() if tip.name in by_id]
        dates = {by_id[ident].collection_date for ident in members}
        if members and len(dates) == 1:
            clusters.append({"sample_ids": sorted(members), "collection_date": dates.pop()})
        else:
            for child in node.clades:
                visit(child)

    visit(tree.root)
    return clusters


def permute_cluster_dates(clusters, rng):
    dates = [cluster["collection_date"] for cluster in clusters]
    rng.shuffle(dates)
    return {
        ident: value
        for cluster, value in zip(clusters, dates, strict=True)
        for ident in cluster["sample_ids"]
    }


def _distance_date_association(tree, samples):
    """Descriptive pairwise distance/date correlation, explicitly not a p-value."""
    from chronoclade.metadata_dates import date_interval

    dated, values = [], {}
    for sample in samples:
        interval = date_interval({"collection_date": sample.collection_date})
        if interval:
            values[sample.sample_id] = (interval["year_min"] + interval["year_max"]) / 2
            dated.append(sample.sample_id)
    genetic, temporal = [], []
    for index, left in enumerate(dated):
        for right in dated[index + 1 :]:
            genetic.append(tree.distance(left, right))
            temporal.append(abs(values[left] - values[right]))
    r = None
    if len(genetic) > 1 and np.std(genetic) and np.std(temporal):
        number = float(np.corrcoef(genetic, temporal)[0, 1])
        r = number if math.isfinite(number) else None
    return {
        "method": "descriptive Pearson correlation of pairwise patristic distances and absolute date-midpoint gaps",
        "pairs": len(genetic),
        "correlation": r,
        "inference": "Pairs share tips and ancestry; this descriptive value has no independent-pair p-value or gate threshold.",
    }


def run_clustered_assessment(
    *,
    tree,
    sequence_length,
    samples,
    observed_clock,
    randomisations,
    randomisation_jobs,
    seed,
    output,
    p_value_threshold,
):
    biological_tree = Phylo.read(tree, "newick")
    clusters = date_clusters(biological_tree, samples)
    observed = parse_clock(observed_clock)
    rng = random.Random(seed)
    assignments = [permute_cluster_dates(clusters, rng) for _ in range(randomisations)]
    with tempfile.TemporaryDirectory(
        prefix="clustered-permutations-", dir=Path(output).parent
    ) as tmp:
        temporary = Path(tmp)
        jobs = []
        for index, assigned in enumerate(assignments):
            dates_path = temporary / f"dates-{index}.csv"
            _write_dates(dates_path, samples, [assigned[s.sample_id] for s in samples])
            jobs.append((dates_path, temporary / f"run-{index}"))

        def run(job):
            return _run_randomised_clock(
                tree=tree, sequence_length=sequence_length, date_path=job[0], run_dir=job[1]
            )

        with ThreadPoolExecutor(
            max_workers=max(1, min(randomisation_jobs, randomisations))
        ) as pool:
            null = [value for value in pool.map(run, jobs) if value is not None]
    temporal = {
        "method": "root_to_tip",
        "permutation_policy": "maximal_monophyletic_same_date_clusters",
        "observed": observed,
        "requested_randomisations": randomisations,
        "successful_randomisations": len(null),
        "seed": seed,
        "randomised": null,
        "p_value_r_squared": (1 + sum(r["r_squared"] >= observed["r_squared"] for r in null))
        / (len(null) + 1),
        "clusters": clusters,
        "cluster_count": len(clusters),
        "eligible_sample_ids": [s.sample_id for s in samples],
        "root_policy": "least-squares root refitted for observed and every permutation",
        "partial_date_policy": "original year/month/day strings retained; TreeTime interprets partial-date intervals",
        "cluster_weighting": "One exchangeable date per cluster; assigned date expands to all cluster tips. Tip weights in each clock fit are unchanged; tip-level date multiplicities need not be preserved.",
        "confounding": _distance_date_association(biological_tree, samples),
        "citation": "https://doi.org/10.1111/2041-210X.12466",
    }
    assessment = assess_temporal_signal(temporal, p_value_threshold=p_value_threshold)
    if len(clusters) < 3 or len({c["collection_date"] for c in clusters}) < 3:
        assessment = {
            "code": "not_assessed",
            "supported": False,
            "reason": "Fewer than three independent date clusters or distinct cluster dates.",
        }
    assessment.update(
        temporal=temporal,
        p_value_threshold=p_value_threshold,
        limitation="Clustered root-to-tip permutations are a screen on a fixed topology; passing does not establish model adequacy or remove all sampling bias.",
    )
    write_json(output, assessment)
    return assessment
