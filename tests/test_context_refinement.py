from copy import deepcopy

import pytest

from chronoclade.context import ContextCandidate, stratified_candidate_pool
from chronoclade.context_refinement import refine_candidate_pool


def candidate(ident, country="UK", year="2020"):
    return ContextCandidate(
        ident,
        "Klebsiella pneumoniae",
        "ST147",
        "klebsiella",
        "147",
        country=country,
        collection_date=year,
        source_genome_id=ident,
    )


def cglin(ident, code="1_2_3_4_5_6_7_8_9_10", version="v1", **extra):
    return dict(
        source_genome_id=ident,
        sample_id=ident,
        cglin_raw=code,
        cglin_scheme="scgMLST629_S",
        cglin_scheme_version=version,
        **extra,
    )


def profile(ident, alleles, version="v1", **extra):
    return dict(
        source_genome_id=ident,
        sample_id=ident,
        cgmlst_profile=alleles,
        cgmlst_scheme="test",
        cgmlst_scheme_version=version,
        **extra,
    )


def refine(rows, queries, limit=4, **kw):
    return refine_candidate_pool(
        [candidate(r["source_genome_id"]) for r in rows], rows, queries, limit=limit, seed=7, **kw
    )


def test_full_prefix_and_version_are_required():
    rows = [
        cglin("match"),
        cglin("same_final", "9_2_3_4_5_6_7_8_9_10"),
        cglin("version", version="v2"),
        cglin("partial", "1_2_3_4_5_6_?"),
    ]
    selected, audit = refine(rows, [cglin("Q")], limit=2)
    assert selected[0].sample_id == "match"
    assert audit["counts"]["priority_selected"] == 1
    assert audit["candidates"]["same_final"]["comparisons"][0]["matches"] == []
    assert (
        "cglin_incompatible_scheme_version"
        in audit["candidates"]["version"]["comparisons"][0]["reasons"]
    )
    assert audit["counts"]["background_selected"] == 1


def test_unversioned_query_falls_back_without_mutation():
    rows = [cglin("A"), cglin("B")]
    query = cglin("Q", version="unknown")
    original = deepcopy(query)
    candidates = [candidate("A"), candidate("B")]
    selected, audit = refine_candidate_pool(candidates, rows, [query], limit=1, seed=7)
    assert [c.sample_id for c in selected] == [
        c.sample_id for c in stratified_candidate_pool(candidates, limit=1, seed=7)
    ]
    assert query == original
    assert candidates[0].selection_reason == ""
    assert audit["query_status"]["Q"] == "no_comparable_evidence_same_st_fallback"


def test_exact_database_fingerprint_can_scope_unknown_versions():
    q = cglin("Q", version="unknown", cglin_database_sha256="a" * 64)
    rows = [
        cglin("A", version="unknown", cglin_database_sha256="a" * 64),
        cglin("B", version="unknown", cglin_database_sha256="b" * 64),
    ]
    _, audit = refine(rows, [q], limit=2)
    assert audit["counts"]["priority_eligible"] == 1


def test_explicit_depth_uses_entire_prefix():
    rows = [cglin("A", "1_2_3_4_5_6_9_8_9_10"), cglin("B", "9_2_3_4_5_6_9_8_9_10")]
    _, audit = refine(rows, [cglin("Q")], cglin_depth=6)
    assert audit["counts"]["priority_eligible"] == 1


@pytest.mark.parametrize("status", ["conflict", "malformed", "unsupported", "unassigned", "failed"])
def test_unresolved_codes_are_not_promoted(status):
    _, audit = refine([cglin("A")], [cglin("Q", cglin_status=status)])
    assert audit["counts"]["priority_eligible"] == 0


def test_hiercc_uses_selected_level_and_compatible_scheme():
    q = dict(
        sample_id="Q",
        hiercc_scheme="ecoli",
        hiercc_scheme_version="2026",
        hiercc_codes={"HC1100": "15", "HC10": "20"},
    )
    rows = [
        dict(source_genome_id="A", **{k: v for k, v in q.items() if k != "sample_id"}),
        dict(
            source_genome_id="B",
            hiercc_scheme="ecoli",
            hiercc_scheme_version="2026",
            hiercc_codes={"HC1100": "99", "HC10": "20"},
        ),
    ]
    _, audit = refine(rows, [q], hiercc_level="HC1100")
    assert audit["counts"]["priority_eligible"] == 1
    _, narrow = refine(rows, [q], hiercc_level="HC10")
    assert narrow["counts"]["priority_eligible"] == 2


def test_profile_distance_orders_group_matches_and_missing_calls_do_not_match():
    q = dict(
        **cglin("Q"),
        **{
            k: v
            for k, v in profile("Q", {str(i): "1" for i in range(10)}).items()
            if k not in ("sample_id", "source_genome_id")
        },
    )
    rows = []
    for ident, changes in [("far", 4), ("near", 1), ("zero", 0)]:
        calls = {str(i): "2" if i < changes else "1" for i in range(10)}
        rows.append(
            dict(
                **cglin(ident),
                **{
                    k: v
                    for k, v in profile(ident, calls).items()
                    if k not in ("sample_id", "source_genome_id")
                },
            )
        )
    rows.append(profile("missing", {str(i): "0" for i in range(10)}))
    selected, audit = refine(rows, [q])
    assert [c.sample_id for c in selected[:3]] == ["zero", "near", "far"]
    comp = audit["candidates"]["missing"]["comparisons"][0]
    assert comp["shared_called_loci"] == 0
    assert "allele_differences" not in comp


def test_sparse_profiles_require_declared_scheme_locus_universe():
    loci = [str(i) for i in range(10)]
    q = profile("Q", {key: "1" for key in loci}, cgmlst_loci=loci)
    rows = [
        profile("nine", {key: "1" for key in loci[:9]}, cgmlst_loci=loci),
        profile("eight", {key: "1" for key in loci[:8]}, cgmlst_loci=loci),
    ]
    _, audit = refine(rows, [q])
    assert audit["counts"]["priority_eligible"] == 1
    assert audit["candidates"]["nine"]["comparisons"][0]["shared_called_loci"] == 9
    assert (
        "cgmlst_insufficient_called_overlap"
        in audit["candidates"]["eight"]["comparisons"][0]["reasons"]
    )


def test_profile_incompatibility_and_locus_mismatch_are_audited():
    rows = [
        profile("version", {"a": 1}, version="v2"),
        profile("loci", {"b": 1}),
        profile("unknown", {"a": 1}, version="unknown"),
    ]
    _, audit = refine(rows, [profile("Q", {"a": 1})])
    assert audit["counts"]["priority_eligible"] == 0
    assert audit["counts"]["comparison_reasons"]["cgmlst_incompatible_locus_set"] == 1


def test_multiple_queries_take_turns_and_background_is_preserved():
    rows = [cglin(str(i)) for i in range(6)] + [
        cglin("other", "2_2_3_4_5_6_7_8_9_10"),
        cglin("wide", "3_2_3_4_5_6_7_8_9_10"),
    ]
    selected, audit = refine(rows, [cglin("Q1"), cglin("Q2", "2_2_3_4_5_6_7_8_9_10")], limit=4)
    assert selected[0].selection_reason.endswith("Q1")
    assert selected[1].sample_id == "other"
    assert audit["counts"]["reserved_background"] == 1
    assert audit["counts"]["priority_selected"] == 3
    assert audit["counts"]["background_selected"] == 1
    assert len({c.sample_id for c in selected}) == 4


def test_input_order_does_not_change_selection():
    rows = [cglin(str(i)) for i in range(10)]
    first, _ = refine(rows, [cglin("Q")])
    second, _ = refine(rows[::-1], [cglin("Q")])
    assert [c.sample_id for c in first] == [c.sample_id for c in second]


def test_empty_pool_and_empty_queries():
    selected, audit = refine([], [])
    assert selected == []
    assert audit["counts"]["eligible_same_st"] == 0
    selected, audit = refine([cglin("A")], [])
    assert len(selected) == 1
    assert audit["counts"]["priority_selected"] == 0


@pytest.mark.parametrize(
    "kwargs",
    [
        {"background_fraction": 0},
        {"background_fraction": float("nan")},
        {"min_profile_overlap": 0},
        {"hiercc_level": "1100"},
        {"cglin_depth": 11},
    ],
)
def test_invalid_configuration_is_rejected(kwargs):
    with pytest.raises(ValueError):
        refine([cglin("A")], [cglin("Q")], **kwargs)


def test_duplicate_ids_fail_instead_of_arbitrary_assignment():
    with pytest.raises(ValueError, match="Catalogue"):
        refine_candidate_pool(
            [candidate("A")], [cglin("A"), cglin("A")], [cglin("Q")], limit=1, seed=7
        )


def test_novel_exact_sequence_hashes_are_callable():
    q = profile("Q", {"a": "A" * 40, "b": "1"})
    rows = [
        profile("same", {"a": "a" * 40, "b": "1"}),
        profile("different", {"a": "b" * 40, "b": "1"}),
    ]
    _, audit = refine(rows, [q])
    assert audit["candidates"]["same"]["comparisons"][0]["allele_differences"] == 0
    assert audit["candidates"]["different"]["comparisons"][0]["allele_differences"] == 1
    assert audit["candidates"]["same"]["comparisons"][0]["shared_called_loci"] == 2


def test_native_novel_alleles_are_included_without_replacing_known_calls():
    q = profile("Q", {"b": "1"}, cgmlst_loci=["a", "b"], cgmlst_novel_alleles={"a": "a" * 40})
    rows = [profile("A", {"a": "a" * 40, "b": "1"}, cgmlst_loci=["a", "b"])]
    _, audit = refine(rows, [q])
    assert audit["counts"]["priority_eligible"] == 1
    assert audit["candidates"]["A"]["comparisons"][0]["allele_differences"] == 0


def test_conflicting_database_hashes_override_same_declared_version():
    rows = [cglin("A", cglin_database_sha256="b" * 64)]
    _, audit = refine(rows, [cglin("Q", cglin_database_sha256="a" * 64)])
    assert audit["counts"]["priority_eligible"] == 0
    rows = [profile("A", {"a": "1"}, cgmlst_database_sha256="b" * 64)]
    _, audit = refine(rows, [profile("Q", {"a": "1"}, cgmlst_database_sha256="a" * 64)])
    assert audit["counts"]["priority_eligible"] == 0


@pytest.mark.parametrize("field", ["cglin_export_sha256", "cglin_public_typing_export_sha256"])
def test_unversioned_cglin_can_match_inside_same_frozen_export(field):
    query = cglin("Q", version="unknown", **{field: "a" * 64})
    rows = [
        cglin("match", version="unknown", **{field: "a" * 64}),
        cglin("independent", version="unknown", **{field: "b" * 64}),
    ]
    selected, audit = refine(rows, [query], limit=2)
    assert selected[0].sample_id == "match"
    assert audit["counts"]["priority_eligible"] == 1
    scope = audit["candidates"]["match"]["comparisons"][0]["assignment_scopes"]["cglin"]["query"]
    assert scope == ("scgMLST629_S", "exportsha256:" + "a" * 64)
    assert (
        "cglin_incompatible_scheme_version"
        in audit["candidates"]["independent"]["comparisons"][0]["reasons"]
    )
    assert query["cglin_scheme_version"] == "unknown"
    assert "cglin_database_sha256" not in query


def test_unversioned_cglin_independent_export_never_matches():
    _, audit = refine(
        [cglin("A", version="unknown", cglin_export_sha256="b" * 64)],
        [cglin("Q", version="unknown", cglin_export_sha256="a" * 64)],
    )
    assert audit["counts"]["priority_eligible"] == 0
    assert audit["query_status"]["Q"] == "no_comparable_evidence_same_st_fallback"


def test_whole_export_scopes_mixed_original_batches(tmp_path):
    import json
    from chronoclade.cglin import load_cglin_export

    path = tmp_path / "catalogue-typing.json"
    path.write_text(
        json.dumps(
            {
                "assignments": [
                    cglin("A", version="unknown", cglin_export_sha256="a" * 64),
                    cglin("B", version="unknown", cglin_export_sha256="b" * 64),
                    cglin("newer", version="v2", cglin_export_sha256="c" * 64),
                ]
            }
        )
    )
    rows = load_cglin_export(path)
    query = dict(rows[0], sample_id="Q")
    _, audit = refine(rows, [query])
    assert audit["counts"]["priority_eligible"] == 2
    assert audit["candidates"]["B"]["comparisons"][0]["matches"] == ["cglin"]
    assert (
        "cglin_incompatible_scheme_version"
        in audit["candidates"]["newer"]["comparisons"][0]["reasons"]
    )
    second = tmp_path / "independent.json"
    second.write_text(json.dumps([cglin("C", version="unknown", cglin_export_sha256="a" * 64)]))
    independent = load_cglin_export(second)
    _, other = refine(independent, [query])
    assert other["counts"]["priority_eligible"] == 0
