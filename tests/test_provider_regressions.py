"""Cross-provider regression checks for resumed analysis and report provenance."""

import csv
import json
from dataclasses import replace

import pytest

from chronoclade.lineage import LineageFiles, _temporal_signal, context_evidence
from chronoclade.metadata import Sample
from chronoclade.pathogenwatch import content_hash
from chronoclade.report import _context_geography_visual
from chronoclade.temporal import file_sha256


def sample(tmp_path, name="focal", *, origin="local", collection_date="2019"):
    assembly = tmp_path / f"{name}.fasta"
    assembly.write_text(f">{name}\nACGT\n")
    return Sample(
        name, assembly, collection_date, "India", "Klebsiella pneumoniae", "ST147", origin
    )


def test_temporal_cache_invalidates_dates_clock_and_randomisation_settings(tmp_path, monkeypatch):
    files = LineageFiles.in_directory(tmp_path)
    tree = tmp_path / "tree.newick"
    tree.write_text("(focal:0.1,context:0.2);")
    clock = tmp_path / "clock.txt"
    clock.write_text("R2=0.75")
    members = [sample(tmp_path), sample(tmp_path, "context", origin="context")]
    calls = []

    def runner(**kwargs):
        calls.append([s.collection_date for s in kwargs["samples"]])
        return {
            "method": "root_to_tip",
            "tree_sha256": file_sha256(kwargs["tree"]),
            "requested_randomisations": kwargs["randomisations"],
            "seed": kwargs["seed"],
            "sequence_length": kwargs["sequence_length"],
            "generation": len(calls),
        }

    monkeypatch.setattr("chronoclade.lineage.run_date_randomisation", runner)

    def run(selected, randomisations=10):
        return _temporal_signal(
            files,
            selected,
            sequence_length=123,
            randomisations=randomisations,
            randomisation_jobs=1,
            seed=7,
            force=False,
            tree=tree,
            method="root_to_tip",
            observed_clock=clock,
        )

    first = run(members)
    assert run(members) == first and len(calls) == 1
    corrected = [replace(members[0], collection_date="2018"), members[1]]
    changed_dates = run(corrected)
    assert len(calls) == 2
    assert changed_dates["dates_sha256"] != first["dates_sha256"]
    assert changed_dates["tree_sha256"] == first["tree_sha256"]
    clock.write_text("R2=0.20")
    changed_clock = run(corrected)
    assert len(calls) == 3
    assert changed_clock["observed_clock_sha256"] != changed_dates["observed_clock_sha256"]
    assert run(corrected, randomisations=20)["requested_randomisations"] == 20
    assert len(calls) == 4
    assert run(corrected, randomisations=20)["generation"] == 4 and len(calls) == 4


def test_report_selected_context_matches_actual_members_not_entire_manifest(tmp_path, monkeypatch):
    directory = tmp_path / "analysis"
    source = tmp_path / "source"
    source.mkdir()
    catalogue = source / "context_catalogue.json"
    payload = {
        "rows": [{"source_genome_id": "included"}, {"source_genome_id": "removed"}],
        "focal_rows": [],
        "provenance": {"retrieved_at": "2026-10-09"},
    }
    digest = content_hash(payload)
    catalogue.write_text(json.dumps({**payload, "snapshot_sha256": digest}))
    rows = [
        {
            "sample_id": f"PW_{ident}",
            "species": "Klebsiella pneumoniae",
            "lineage": "ST147",
            "source": "pathogenwatch",
            "source_genome_id": ident,
            "country": "India",
            "catalogue_path": str(catalogue),
            "catalogue_sha256": digest,
        }
        for ident in ("included", "removed")
    ]
    observed = {}

    def capture(public, destination, **kwargs):
        observed["public"] = public
        observed["selected"] = kwargs["selected_source_ids"]
        return {"outputs": {}}

    monkeypatch.setattr("chronoclade.context_geography.generate_context_geography", capture)
    evidence = context_evidence(
        [sample(tmp_path), sample(tmp_path, "PW_included", origin="context")],
        rows,
        directory=directory,
    )
    assert observed["selected"] == ["included"]
    assert len(observed["public"]) == 2  # Full-catalogue composition retains the denominator.
    assert evidence["context_samples"] == 1
    assert [r["sample_id"] for r in evidence["nearest_screening_contexts"]] == ["PW_included"]
    with (directory / "context_manifest.tsv").open(newline="") as handle:
        exported = list(csv.DictReader(handle, delimiter="\t"))
    assert [r["source_genome_id"] for r in exported] == ["included"]


@pytest.mark.parametrize(
    "manifest_rows",
    [
        [],
        [
            {
                "sample_id": "unused",
                "species": "Klebsiella pneumoniae",
                "lineage": "ST147",
                "source": "pathogenwatch",
                "catalogue_path": "/missing/previous/catalogue.json",
            }
        ],
    ],
)
def test_rerun_without_current_context_removes_stale_geography_and_downloads(
    tmp_path, manifest_rows
):
    directory = tmp_path / "analysis"
    geography = directory / "context_geography"
    geography.mkdir(parents=True)
    (geography / "fragment.html").write_text("<img src='obsolete.svg'>")
    stale = [
        directory / name
        for name in ("context_manifest.tsv", "context_catalogue.json", "context_selection.json")
    ]
    for path in stale:
        path.write_text("previous run")
    assert _context_geography_visual(directory)
    evidence = context_evidence([sample(tmp_path)], manifest_rows, directory=directory)
    assert not evidence["geography_available"] and not evidence["manifest_available"]
    assert _context_geography_visual(directory) == ""
    assert not geography.exists()
    assert not any(path.exists() for path in stale)


def test_context_directory_moves_between_laptop_and_cluster(tmp_path):
    import shutil
    from chronoclade.context import ContextCandidate, write_candidate_table, write_combined_metadata
    from chronoclade.lineage import read_context_manifest
    from chronoclade.metadata import Sample, read_metadata

    original = tmp_path / "laptop"
    context = original / "context"
    context.mkdir(parents=True)
    focal_assembly = original / "focal.fasta"
    focal_assembly.write_text(">focal\nACGT\n")
    candidate_assembly = context / "public.fasta"
    candidate_assembly.write_text(">context\nACGT\n")
    catalogue = context / "context_catalogue.json"
    catalogue.write_text("{}")
    focal = Sample("F1", focal_assembly, "2020", "India", "Klebsiella pneumoniae", "ST147", "local")
    candidate = ContextCandidate(
        "PW_1",
        "Klebsiella pneumoniae",
        "ST147",
        "klebsiella",
        "147",
        collection_date="2019",
        assembly=str(candidate_assembly),
        catalogue_path=str(catalogue),
        source="pathogenwatch",
    )
    write_candidate_table(context / "context_manifest.tsv", [candidate])
    write_combined_metadata(context / "combined_metadata.csv", [focal], [candidate])
    relocated = tmp_path / "cluster"
    shutil.copytree(original, relocated)
    samples = read_metadata(relocated / "context/combined_metadata.csv")
    manifest = read_context_manifest(relocated / "context/context_manifest.tsv")
    assert samples[0].assembly == relocated / "focal.fasta"
    assert samples[1].assembly == relocated / "context/public.fasta"
    assert manifest[0]["catalogue_path"] == str(relocated / "context/context_catalogue.json")
    assert str(original) not in (context / "combined_metadata.csv").read_text()
