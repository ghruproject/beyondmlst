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
