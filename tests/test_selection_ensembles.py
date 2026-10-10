"""Mandatory identity, interchangeable representatives and portable tree inputs."""

import csv
import json
import shutil

import pytest

from chronoclade.datasets import PreparedDataset, samples_from_records, write_dataset
from chronoclade.errors import WorkflowError
from chronoclade.selection_manifest import (
    load_selection_ensemble,
    load_selection_manifest,
    write_selection_ensemble,
)
from chronoclade.selections import (
    DistanceEvidence,
    select_assembly_context,
    select_context_ensemble,
    selected_context_records,
)


def row(identifier, **kwargs):
    return {
        "sample_id": identifier,
        "origin": "context",
        "country": "UK",
        "collection_date": "2020",
        **kwargs,
    }


def analysis(tmp_path, records):
    path = tmp_path / "distances.csv"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "sample_id_1",
                "sample_id_2",
                "distance",
                "allele_differences",
                "shared_called_loci",
            ],
        )
        writer.writeheader()
        writer.writerows(records)
    return {"paths": {"pairwise_distances": str(path)}}


def pair(a, b, raw, normalized):
    return {
        "sample_id_1": a,
        "sample_id_2": b,
        "allele_differences": raw,
        "distance": normalized,
        "shared_called_loci": 100,
    }


def test_raw_allele_boundary_ties_union_with_pins_overrun(tmp_path):
    queries = [row("q1"), row("q2")]
    context = [row(i) for i in ("a", "b", "c", "d", "pin", "far")]
    distances = analysis(
        tmp_path,
        [
            pair("q1", "a", 2, 0.02),
            pair("q1", "b", 2, 0.20),
            pair("q1", "c", 2, 0.01),
            pair("q1", "far", 3, 0.005),
            pair("q2", "d", 1, 0.01),
        ],
    )
    selections, ensemble = select_context_ensemble(
        queries,
        context,
        distances,
        size=1,
        nearest_per_query=1,
        include=["pin"],
        replicates=3,
    )
    first = selections[0]
    assert first["selected_context_ids"] == ["pin", "a", "d", "b", "c"]
    assert first["budget_overrun"] == 4
    assert first["mandatory_sample_ids"] == first["selected_sample_ids"]
    assert first["nearest_neighbours"][0]["boundary_tie_ids"] == ["a", "b", "c"]
    assert first["nearest_neighbours"][0]["tie_expanded"]
    assert ensemble["distinct_selections"] == 1
    assert selections[1]["duplicate_of"] == "selection-001"
    assert selections[1]["alternative_status"] == "duplicate"
    assert selections[1]["variation_status"] == "no_nonmandatory_choices"
    assert "far" not in first["selected_context_ids"]


def test_zero_budget_keeps_inputs_nearest_and_pins(tmp_path):
    distances = analysis(tmp_path, [pair("q", "a", 1, 0.01)])
    contexts, audit = select_assembly_context(
        [row("q")],
        [row("a"), row("pin")],
        distances,
        size=0,
        include=["pin"],
    )
    assert {r["sample_id"] for r in contexts} == {"a", "pin"}
    assert audit["budget_overrun"] == 2


def test_pinned_nearest_preserves_both_reasons(tmp_path):
    distances = analysis(tmp_path, [pair("q", "a", 1, 0.01)])
    selections, _ = select_context_ensemble([row("q")], [row("a")], distances, include=["a"])
    assert selections[0]["decisions"][0]["reasons"] == [
        "user_requested",
        "cgmlst_neighbour:q:rank=1",
    ]


def test_alternatives_change_only_optional_ids_and_freeze_quotas():
    context = [row(f"a{i}") for i in range(12)]
    distances = {"distance_evidence": DistanceEvidence((pair("q", "a0", 0, 0),))}
    selections, ensemble = select_context_ensemble(
        [row("q")],
        context,
        distances,
        size=4,
        replicates=5,
        seed=73,
    )
    repeated = select_context_ensemble(
        [row("q")], list(reversed(context)), distances, size=4, replicates=5, seed=73
    )
    assert (selections, ensemble) == repeated
    assert ensemble["distinct_selections"] == 5
    for selection in selections:
        assert selection["mandatory_sample_ids"] == ["q", "a0"]
        assert selection["quotas"] == selections[0]["quotas"]
        assert selection["strata_coverage"] == [
            {"stratum": ["unassigned", "2020", "UK", "Unknown"], "available": 12, "selected": 4}
        ]
    overlap = ensemble["pairwise_overlap"][0]
    left, right = (set(s["selected_sample_ids"]) for s in selections[:2])
    assert overlap["shared_ids"] == sorted(left & right)
    assert overlap["jaccard"] == len(left & right) / len(left | right)


def test_missing_dates_are_retained_and_reported():
    selections, _ = select_context_ensemble([], [row("x", collection_date=None)], {}, size=1)
    assert selections[0]["selected_context_ids"] == ["x"]
    assert selections[0]["missing_metadata"]["collection_date"] == ["x"]
    assert selections[0]["strata_coverage"][0]["stratum"][1] == "Unknown"


def test_distance_adapter_labels_esm2_without_cgmlst_policy():
    evidence = DistanceEvidence(
        ({"sample_id_1": "q", "sample_id_2": "x", "distance": 0.5},),
        method="esm2",
        distance_definition="cosine distance",
    )
    selections, _ = select_context_ensemble([row("q")], [row("x")], {"distance_evidence": evidence})
    assert selections[0]["source_method"] == "esm2"
    assert selections[0]["nearest_ranking"] == "supplied_distance"
    assert selections[0]["decisions"][0]["reason"] == "esm2_neighbour:q:rank=1"


def test_incomplete_raw_evidence_never_mixes_raw_and_normalized_ranks(tmp_path):
    distances = analysis(tmp_path, [pair("q", "a", 2, 0.2), pair("q", "b", None, 0.001)])
    selections, _ = select_context_ensemble(
        [row("q")], [row("a"), row("b")], distances, size=0, nearest_per_query=1
    )
    assert selections[0]["required_nearest_context_ids"] == ["a"]


def test_selection_pool_and_evidence_fail_closed():
    with pytest.raises(WorkflowError, match="disjoint"):
        select_context_ensemble([row("a")], [row("a")], {})
    with pytest.raises(WorkflowError, match="uniquely"):
        select_context_ensemble(
            [], [row("a", accession="alias"), row("b", accession="alias")], {}, include=["alias"]
        )
    with pytest.raises(WorkflowError, match="Conflicting"):
        select_context_ensemble(
            [row("q")],
            [row("a")],
            {
                "distance_evidence": DistanceEvidence(
                    (pair("q", "a", 1, 0.01), pair("a", "q", 2, 0.02))
                )
            },
        )


@pytest.fixture
def bundle(tmp_path):
    queries, context = [row("q", origin="local")], [row("a"), row("b")]
    dataset = PreparedDataset(samples=samples_from_records(queries + context))
    dataset_manifest = write_dataset(dataset, tmp_path / "prepared")
    distances = analysis(tmp_path, [pair("q", "a", 1, 0.01), pair("q", "b", 2, 0.02)])
    manifest = write_selection_ensemble(
        queries,
        context,
        distances,
        dataset_manifest=dataset_manifest,
        output=tmp_path / "runs",
        replicates=2,
        size=1,
        nearest_per_query=1,
    )
    return manifest, queries, context, distances


def test_manifest_roundtrip_relative_portability_and_exact_tree_ids(bundle, tmp_path):
    manifest, queries, context, _ = bundle
    ensemble = load_selection_ensemble(manifest)
    selection_path = manifest.parent / ensemble["selections"][0]["path"]
    selection = load_selection_manifest(selection_path)
    assert selection["selected_sample_ids"] == ["q", "a"]
    assert selected_context_records(selection, queries, context)[0]["sample_id"] == "a"
    moved = tmp_path / "moved"
    moved.mkdir()
    shutil.copytree(tmp_path / "prepared", moved / "prepared")
    shutil.copytree(tmp_path / "runs", moved / "runs")
    shutil.rmtree(tmp_path / "prepared")
    shutil.rmtree(tmp_path / "runs")
    assert load_selection_ensemble(moved / "runs" / "ensemble.json") == ensemble
    assert selection["source"]["distance_evidence"]["sha256"]
    assert selection["source"]["dataset"]["sha256"]


def test_manifest_tampered_distances_and_dataset_rejected(bundle):
    manifest, _, _, _ = bundle
    evidence = manifest.parent / "distance_evidence.json"
    evidence.write_text("{}")
    with pytest.raises(WorkflowError, match="source.distance_evidence.sha256"):
        load_selection_ensemble(manifest)


def test_manifest_mandatory_and_order_validation(bundle):
    manifest, _, _, _ = bundle
    ensemble = load_selection_ensemble(manifest)
    selection_path = manifest.parent / ensemble["selections"][0]["path"]
    selection = load_selection_manifest(selection_path)
    selection["selected_sample_ids"] = ["a", "q"]
    selection_path.write_text(json.dumps(selection))
    with pytest.raises(WorkflowError, match="selected_sample_ids"):
        load_selection_manifest(selection_path)
    with pytest.raises(WorkflowError, match="selections\\[0\\].sha256"):
        load_selection_ensemble(manifest)


def test_selection_writer_does_not_overwrite_or_publish_unknown_pool(bundle):
    manifest, queries, context, distances = bundle
    dataset_manifest = manifest.parent.parent / "prepared" / "dataset.json"
    with pytest.raises(WorkflowError, match="already exists"):
        write_selection_ensemble(
            queries, context, distances, dataset_manifest=dataset_manifest, output=manifest.parent
        )
    output = manifest.parent.parent / "invalid"
    with pytest.raises(WorkflowError, match="outside source.dataset"):
        write_selection_ensemble(
            queries,
            context + [row("unknown")],
            distances,
            dataset_manifest=dataset_manifest,
            output=output,
        )
    assert not output.exists()


def test_partition_exact_block_identity_enforced(bundle):
    from chronoclade.artifacts import file_sha256, write_json

    manifest, queries, context, distances = bundle
    dataset_manifest = manifest.parent.parent / "prepared" / "dataset.json"
    partition_path = manifest.parent.parent / "partition.json"
    partition = {
        "schema": "chronoclade.cgmlst.partitions",
        "schema_version": 1,
        "dataset_id": json.loads(dataset_manifest.read_text())["dataset_id"],
        "dataset_sha256": file_sha256(dataset_manifest),
        "blocks": [{"block_id": "block-1", "sample_ids": ["q", "a", "b"]}],
    }
    write_json(partition_path, partition)
    output = manifest.parent.parent / "partitioned_selection"
    ensemble_path = write_selection_ensemble(
        queries,
        context,
        dict(distances, partition={"block_id": "block-1"}),
        dataset_manifest=dataset_manifest,
        partition_manifest=partition_path,
        output=output,
    )
    ensemble = load_selection_ensemble(ensemble_path)
    assert ensemble["source"]["partition_block_id"] == "block-1"
    subset_output = manifest.parent.parent / "wrong_subset"
    with pytest.raises(WorkflowError, match="exact source.partition"):
        write_selection_ensemble(
            queries,
            context[:1],
            dict(distances, partition={"block_id": "block-1"}),
            dataset_manifest=dataset_manifest,
            partition_manifest=partition_path,
            output=subset_output,
        )
    assert not subset_output.exists()


def test_distance_ids_validated_before_any_publication(bundle):
    manifest, queries, context, _ = bundle
    evidence = DistanceEvidence((pair("q", "outside", 1, 0.01),))
    output = manifest.parent.parent / "unknown_evidence"
    with pytest.raises(WorkflowError, match="Distance evidence.*outside"):
        write_selection_ensemble(
            queries,
            context,
            {"distance_evidence": evidence},
            dataset_manifest=manifest.parent.parent / "prepared" / "dataset.json",
            output=output,
        )
    assert not output.exists()


def test_overlap_evidence_and_dataset_hash_are_validated(bundle):
    manifest, _, _, _ = bundle
    ensemble = load_selection_ensemble(manifest)
    ensemble["pairwise_overlap"][0]["intersection"] = 0
    manifest.write_text(json.dumps(ensemble))
    with pytest.raises(WorkflowError, match="pairwise_overlap"):
        load_selection_ensemble(manifest)
    dataset_manifest = manifest.parent.parent / "prepared" / "dataset.json"
    dataset_manifest.write_text(dataset_manifest.read_text() + "\n")
    with pytest.raises(WorkflowError, match="source.dataset.sha256"):
        load_selection_ensemble(manifest)


def test_selection_date_strata_use_validated_intervals_not_raw_prefixes():
    context = [row('canonical',collection_date=None,date_start='2018-01-01',date_end='2020-12-31',date_precision='interval'),
               row('invalid',collection_date='09/01/2019'), row('future',collection_date='2099')]
    selected,_ = select_context_ensemble([],context,{},size=3)
    audit=selected[0]
    assert {tuple(row['stratum']) for row in audit['strata_coverage']} == {
        ('unassigned','2018..2020','UK','Unknown'), ('unassigned','Unknown','UK','Unknown')}
    assert audit['missing_metadata']['collection_date'] == ['future','invalid']
