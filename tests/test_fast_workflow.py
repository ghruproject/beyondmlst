import csv
import json
from pathlib import Path

from chronoclade.adaptive_context import adaptive_cglin_context
from chronoclade.fast_workflow import route_datasets, run_fast_datasets


def record(name, cg, origin="context", calls=None):
    return dict(sample_id=name, source_genome_id=name, species="Klebsiella pneumoniae",
                lineage="ST39", mlst_st="39", origin=origin, country="Greece",
                collection_date="2019", cgmlst_scheme="test", cgmlst_scheme_version="1",
                cgmlst_loci=["a", "b", "c", "d"],
                cgmlst_profile=dict(zip(["a", "b", "c", "d"], calls or [1, 1, 1, 1])),
                cglin_scheme="scgMLST629_S", cglin_scheme_version="1",
                cglin_raw=f"0,0,107,{cg},0,0,0,0,0,0", cglin_status="complete",
                cglin_export_row={"Clonal Group": "39" if cg == 0 else "10192"})


def test_fast_partitions_cg_before_neighbour_search_and_keeps_public_pool(tmp_path):
    queries = [record("q1", 0, "local"), record("q2", 1, "local")]
    public = [record("a", 0, calls=[2, 1, 1, 1]), record("b", 1)]
    queries, context, audit = adaptive_cglin_context(queries, public, min_context=1)
    inputs = route_datasets(dict(queries=queries, context=context,
                                catalogue_rows=public, provenance={"adaptive_context_selection": audit}))
    analysis, landing, datasets = run_fast_datasets(inputs, tmp_path, seed=42,
                                                   bootstrap_replicates=0,
                                                   distance_threshold=0.02, tree_limit=2)
    # b is genetically identical to q1, but belongs to the other CG dataset.
    assert {(n["query_id"], n["context_id"]) for n in analysis["nearest_neighbours"]} == {("q1", "a"), ("q2", "b")}
    assert len(datasets) == 2
    assert landing.name == "profile_report.html"
    assert all(Path(d["report"]).is_file() for d in datasets)
    assert "ST39_CG39" in landing.read_text() and "ST39_CG10192" in landing.read_text()
    assert all(r["input_lineage"] == "ST39" for r in inputs["queries"])


def test_routing_preserves_fallback_context_for_unassigned_queries():
    queries = [record("typed", 0, "local"), record("unassigned", 0, "local")]
    queries[1].pop("cglin_raw")
    fallback = record("fallback", 0)
    fallback["cglin_raw"] = "0,0,107,0,1,0,0,0,0,0"
    queries, _, audit = adaptive_cglin_context(queries, [fallback])
    explicit = dict(fallback, sample_id="explicit", provided_context=True)
    result = route_datasets(dict(queries=queries, context=[fallback, explicit],
                                 catalogue_rows=[fallback],
                                 provenance={"adaptive_context_selection": audit}))
    assert result["queries"][1]["lineage"] == "ST39"
    assert result["context"][0]["lineage"] == "ST39"
    assert result["context"][1]["lineage"] == "ST39_CG39"


def test_species_case_normalization_keeps_selected_context_in_dataset(tmp_path):
    queries = [record("q", 0, "local")]
    queries[0]["species"] = "klebsiella_pneumoniae"
    public = [record("a", 0)]
    queries, context, audit = adaptive_cglin_context(queries, public, min_context=1)
    inputs = route_datasets(dict(queries=queries, context=context,
                                catalogue_rows=public, provenance={"adaptive_context_selection": audit}))
    analysis, _, datasets = run_fast_datasets(inputs, tmp_path, seed=42,
                                              bootstrap_replicates=0,
                                              distance_threshold=0.02, tree_limit=2)
    assert len(analysis["nearest_neighbours"]) == 1
    assert datasets[0]["context_count"] == 1
    assert {r["species"] for r in inputs["queries"] + inputs["context"]} == {"Klebsiella pneumoniae"}


def test_fast_uses_assembly_selection_and_retains_full_pool_evidence(tmp_path):
    queries = [record("q", 0, "local")]
    context = [record("a", 0, calls=[2, 1, 1, 1]),
               record("b", 0, calls=[2, 1, 1, 1]),
               record("c", 0, calls=[2, 2, 2, 2]),
               record("d", 0, calls=[1, 2, 2, 1]),
               record("e", 0, calls=[1, 2, 1, 2])]
    context[2]["accession"] = "PINNED_ACCESSION"
    inputs = dict(queries=queries, context=context, provenance={})
    analysis, _, _ = run_fast_datasets(
        inputs, tmp_path, seed=42, bootstrap_replicates=0, distance_threshold=0.02,
        tree_limit=2, context_size=3, include_genomes=["PINNED_ACCESSION"],
    )
    selection = analysis["shared_selection"]
    assert selection["selected_context_ids"][:2] == ["c", "a"]
    assert len(selection["selected_context_ids"]) == 5
    assert selection["budget_overrun"] == 2
    assert selection["selected_sample_ids"] == ["q"] + selection["selected_context_ids"]
    assert analysis["tree_display_selection"]["displayed_sample_ids"] == sorted(selection["selected_sample_ids"])
    assert analysis["cohorts"][0]["tree_display_limit"] == 6
    assert selection["available_context_ids"] == ["a", "b", "c", "d", "e"]
    assert selection["decisions"][0]["sample_id"] == "c"
    assert selection["decisions"][0]["reason"] == "user_requested"
    assert selection["decisions"][0]["mandatory"]
    assert analysis["context_selections"] == [dict(species="Klebsiella pneumoniae", lineage="ST39", **selection)]
    assert json.loads((tmp_path / "context_selection.json").read_text()) == selection
    assert json.loads((tmp_path / "profile_provenance.json").read_text())["shared_selection"] == selection
    # Selection changes the displayed tree, while distances still cover every record.
    with Path(analysis["paths"]["pairwise_distances"]).open() as handle:
        pairs = list(csv.DictReader(handle))
    assert len(pairs) == 15
    assert {pair["sample_id_1"] for pair in pairs} | {pair["sample_id_2"] for pair in pairs} == {"q", "a", "b", "c", "d", "e"}
