from copy import deepcopy

import pytest

from chronoclade.adaptive_context import adaptive_cglin_context
from chronoclade.cglin import DEFAULT_SCHEME


def row(ident, code="0,0,107,0,0,0,0,0,0,0", **extra):
    return dict(
        sample_id=ident,
        mlst_st="39",
        lineage="ST39",
        species="Klebsiella pneumoniae",
        cglin_raw=code,
        cglin_scheme=DEFAULT_SCHEME,
        cglin_scheme_version="v1",
        cglin_status="complete",
        cglin_export_row={"Clonal Group": "39"},
        **extra,
    )


def test_narrowest_context_union_and_no_mutation():
    inputs = [row("q1"), row("q2", "0,0,107,0,1,0,0,0,0,0")]
    context = [row(f"a{i}") for i in range(20)]
    context += [row("b", "0,0,107,0,1,0,0,0,0,0"), row("outside", "0,0,107,1,0,0,0,0,0,0")]
    original = deepcopy((inputs, context))
    queries, selected, audit = adaptive_cglin_context(inputs, context)
    assert (inputs, context) == original
    assert {q["analysis_dataset"] for q in queries} == {"ST39_CG39"}
    assert len(selected) == 21
    dataset = audit["datasets"][0]
    assert dataset["public_cg_pool_count"] == 21
    assert [(s["selected_level"], s["limited_context"]) for s in dataset["subgroups"]] == [
        (7, False),
        (5, True),
    ]


def test_distinct_finer_prefix_coverage_not_hidden_by_large_group():
    inputs = [row("q1"), row("q2", "0,0,107,0,0,0,1,0,0,0")]
    context = [row(f"a{i}") for i in range(20)]
    queries, selected, audit = adaptive_cglin_context(inputs, context)
    subgroup = audit["datasets"][0]["subgroups"][0]
    assert subgroup["context_counts"]["7"] == 20
    assert subgroup["minimum_prefix_context_counts"]["7"] == 0
    assert subgroup["selected_level"] == 6
    assert len(selected) == 20


def test_two_cgs_in_same_st_are_independent_and_names_are_not_final_component():
    other = row("q2", "0,0,107,1,0,0,0,0,0,0")
    other["cglin_export_row"] = {"Clonal Group": "10192"}
    p2 = deepcopy(other)
    p2["sample_id"] = "p2"
    queries, selected, audit = adaptive_cglin_context([row("q1"), other], [row("p1"), p2])
    assert {q["analysis_dataset"] for q in queries} == {"ST39_CG39", "ST39_CG10192"}
    assert {q["analysis_dataset"] for q in selected} == {"ST39_CG39", "ST39_CG10192"}
    assert all(d["public_cg_pool_count"] == 1 for d in audit["datasets"])
    unnamed = row("q3")
    unnamed["cglin_export_row"] = {}
    assert (
        adaptive_cglin_context([unnamed], [row("p")])[0][0]["analysis_dataset"]
        == "ST39_CGprefix_0_0_107_0"
    )


def test_version_database_scope_st_and_unsupported_excluded():
    context = [row("good"), row("old"), row("digest"), row("st"), row("unsupported")]
    context[1]["cglin_scheme_version"] = "v2"
    context[2]["cglin_database_sha256"] = "b" * 64
    context[3]["mlst_st"] = "147"
    context[4]["cglin_status"] = "unsupported"
    query = row("q")
    query["cglin_database_sha256"] = "a" * 64
    assert [p["sample_id"] for p in adaptive_cglin_context([query], context)[1]] == ["good"]


def test_unversioned_requires_same_frozen_export():
    query, public = row("q"), row("p")
    query["cglin_scheme_version"] = public["cglin_scheme_version"] = "unknown"
    query["cglin_frozen_export_sha256"] = "a" * 64
    public["cglin_frozen_export_sha256"] = "b" * 64
    assert adaptive_cglin_context([query], [public])[1] == []
    public["cglin_frozen_export_sha256"] = "a" * 64
    assert len(adaptive_cglin_context([query], [public])[1]) == 1


def test_partial_codes_widen_without_inventing_unknown_components():
    query = row("q", "0,0,107,0,0,?, ?,?, ?,?")
    public = row("p")
    subgroup = adaptive_cglin_context([query], [public], min_context=1)[2]["datasets"][0][
        "subgroups"
    ][0]
    assert subgroup["selected_level"] == 5
    assert subgroup["context_counts"] == {"5": 1, "6": 0, "7": 0}
    assert not subgroup["limited_context"]


@pytest.mark.parametrize(
    "change",
    [
        dict(species="Escherichia coli"),
        dict(cglin_status="unsupported"),
        dict(cglin_raw="0,0,107,?"),
        dict(cglin_scheme="other"),
        dict(cglin_scheme_version="unknown"),
    ],
)
def test_ineligible_preserves_old_pathway(change):
    query = row("q")
    query.update(change)
    queries, context, audit = adaptive_cglin_context([query], [row("p")])
    assert "analysis_dataset" not in queries[0]
    assert context == []
    assert audit["datasets"] == []
    assert len(audit["ineligible_inputs"]) == 1


@pytest.mark.parametrize("value", [0, -1, True, 1.5])
def test_invalid_threshold(value):
    with pytest.raises(ValueError, match="positive integer"):
        adaptive_cglin_context([], [], min_context=value)


def test_full_catalogue_annotation_and_partial_cg_counts():
    from chronoclade.adaptive_context import annotate_cglin_datasets

    partial = row("partial", "0,0,107,0,?, ?,?, ?,?,?")
    inputs, selected, audit = adaptive_cglin_context([row("q")], [partial, row("p")])
    assert len(selected) == 1
    assert audit["datasets"][0]["public_cg_pool_count"] == 2
    assert annotate_cglin_datasets([partial], inputs)[0]["analysis_dataset"] == "ST39_CG39"
    assert "analysis_dataset" not in partial


def test_dataset_name_collisions_never_merge_versions():
    q1, q2 = row("q1"), row("q2")
    q2["cglin_scheme_version"] = "v2"
    p1, p2 = row("p1"), row("p2")
    p2["cglin_scheme_version"] = "v2"
    queries, selected, audit = adaptive_cglin_context([q1, q2], [p1, p2])
    names = {q["analysis_dataset"] for q in queries}
    assert len(names) == 2
    assert all(name.startswith("ST39_CG39_") for name in names)
    assert all(d["public_cg_pool_count"] == 1 for d in audit["datasets"])
    assert {q["analysis_dataset"] for q in selected} == names


def test_species_isolation_and_st_lineage_fallback():
    query = row("q")
    query.pop("mlst_st")
    context = [row("good"), row("other"), row("unknown")]
    context[1]["species"] = "Klebsiella variicola"
    context[2]["species"] = ""
    queries, selected, audit = adaptive_cglin_context([query], context)
    assert [p["sample_id"] for p in selected] == ["good"]
    assert audit["datasets"][0]["mlst_st"] == "39"
    other_query = row("q2")
    other_query["species"] = "Klebsiella variicola"
    queries, selected, audit = adaptive_cglin_context([query, other_query], context)
    assert len(audit["datasets"]) == 2
    assert len({q["analysis_dataset"] for q in queries}) == 2
    assert all(d["public_cg_pool_count"] == 1 for d in audit["datasets"])
