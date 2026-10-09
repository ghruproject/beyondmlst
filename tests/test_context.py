import csv
from pathlib import Path


from chronoclade.context import (
    ContextCandidate,
    filter_candidates,
    parse_ska_distances,
    select_context,
    stratified_candidate_pool,
    write_candidate_table,
    write_combined_metadata,
)
from chronoclade.metadata import Sample, read_metadata


def candidate(sample_id: str, country: str, year: str) -> ContextCandidate:
    return ContextCandidate(
        sample_id=sample_id,
        species="Escherichia coli",
        lineage="ST131",
        mlst_scheme="ecoli_achtman_4",
        mlst_st="131",
        collection_date=year,
        country=country,
    )


def test_stratified_pool_is_reproducible_and_spans_strata() -> None:
    candidates = [
        candidate("UK1", "United Kingdom", "2021"),
        candidate("UK2", "United Kingdom", "2021"),
        candidate("PH1", "Philippines", "2022"),
        candidate("PH2", "Philippines", "2022"),
        candidate("IN1", "India", "2023"),
    ]

    first = stratified_candidate_pool(candidates, limit=3, seed=17)
    second = stratified_candidate_pool(candidates, limit=3, seed=17)

    assert [item.sample_id for item in first] == [item.sample_id for item in second]
    assert {item.country for item in first} == {"United Kingdom", "Philippines", "India"}


def test_filter_candidates_applies_metadata_constraints() -> None:
    candidates = [
        candidate("A", "United Kingdom:England", "2021"),
        candidate("B", "Philippines", "2024"),
    ]
    candidates[0].host = "Homo sapiens"
    candidates[0].isolation_source = "blood"

    result = filter_candidates(
        candidates,
        countries=["United Kingdom"],
        year_from=2020,
        year_to=2022,
        host="sapiens",
        isolation_source="blood",
    )

    assert [item.sample_id for item in result] == ["A"]


def test_parse_and_select_context_prioritises_focal_neighbours(tmp_path: Path) -> None:
    distance_file = tmp_path / "distances.tsv"
    distance_file.write_text(
        "Sample1\tSample2\tDistance\tMismatches (proportion)\tMatch count\tMismatch count\n"
        "F1\tC1\t2.00\t0.001\t100\t1\n"
        "F1\tC2\t30.00\t0.030\t100\t3\n"
        "F1\tC3\t40.00\t0.040\t100\t4\n"
        "F2\tC1\t35.00\t0.035\t100\t3\n"
        "F2\tC2\t3.00\t0.003\t100\t1\n"
        "F2\tC3\t20.00\t0.020\t100\t2\n",
        encoding="utf-8",
    )
    distances = parse_ska_distances(distance_file)
    candidates = [
        candidate("C1", "United Kingdom", "2021"),
        candidate("C2", "Philippines", "2022"),
        candidate("C3", "India", "2023"),
    ]

    selected, audit = select_context(
        candidates,
        ["F1", "F2"],
        distances,
        max_context=3,
        nearest_per_focal=1,
        seed=2,
    )

    assert {item.sample_id for item in selected} == {"C1", "C2", "C3"}
    assert next(item for item in selected if item.sample_id == "C1").selection_reason == (
        "nearest_to=F1"
    )
    assert next(item for item in selected if item.sample_id == "C2").selection_reason == (
        "nearest_to=F2"
    )
    assert audit["focal_neighbour_coverage"] == 1.0


def test_zero_distance_is_ranked_before_nonzero_distance() -> None:
    candidates = [
        candidate("IDENTICAL", "United Kingdom", "2021"),
        candidate("DIFFERENT", "Philippines", "2022"),
    ]
    distances = {
        ("IDENTICAL", "F1"): (0.0, 0.0),
        ("DIFFERENT", "F1"): (1.0, 0.001),
    }

    selected, _ = select_context(
        candidates,
        ["F1"],
        distances,
        max_context=2,
        nearest_per_focal=1,
        seed=2,
    )

    assert [item.sample_id for item in selected] == ["IDENTICAL", "DIFFERENT"]


def test_candidate_manifest_is_tab_separated(tmp_path: Path) -> None:
    manifest = write_candidate_table(
        tmp_path / "context_manifest.tsv",
        [candidate("SAMN1", "Philippines", "2023")],
    )

    header = manifest.read_text(encoding="utf-8").splitlines()[0]
    assert "sample_id\tspecies\tlineage" in header


def test_combined_metadata_round_trip(tmp_path: Path) -> None:
    assembly = tmp_path / "focal.fasta"
    assembly.write_text(">focal\nAAAA\n", encoding="utf-8")
    context_assembly = tmp_path / "context.fasta"
    context_assembly.write_text(">context\nAAAA\n", encoding="utf-8")
    context = candidate("SAMN1", "Philippines", "2023-02")
    context.assembly = str(context_assembly.resolve())

    focal = [
        Sample(
            "F1",
            assembly,
            "2024-01",
            "KIMS",
            "Escherichia coli",
            "ST131",
            "local",
        )
    ]
    combined = write_combined_metadata(tmp_path / "combined.csv", focal, [context])
    parsed = read_metadata(combined)

    assert [sample.origin for sample in parsed] == ["local", "context"]
    with combined.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[1]["location"] == "Philippines"
