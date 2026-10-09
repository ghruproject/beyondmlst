import sys
from pathlib import Path

from chronoclade.lineage import LineageFiles, _run_command, _run_dated_tree


COPY_INPUT = (
    "from pathlib import Path; import sys; "
    "Path(sys.argv[2]).write_text(Path(sys.argv[1]).read_text())"
)


def test_stage_reruns_when_direct_input_changes(tmp_path: Path) -> None:
    source = tmp_path / "source.txt"
    output = tmp_path / "output.txt"
    log = tmp_path / "logs" / "copy.log"
    source.write_text("first", encoding="utf-8")
    command = [sys.executable, "-c", COPY_INPUT, str(source), str(output)]

    _run_command(command, log=log, expected=output, force=False, inputs=(source,))
    assert output.read_text(encoding="utf-8") == "first"

    source.write_text("second", encoding="utf-8")
    _run_command(command, log=log, expected=output, force=False, inputs=(source,))

    assert output.read_text(encoding="utf-8") == "second"
    assert (tmp_path / "logs" / "copy.fingerprint.json").is_file()


def test_dated_tree_keeps_the_root_used_for_temporal_testing(tmp_path: Path, monkeypatch) -> None:
    files = LineageFiles.in_directory(tmp_path)
    observed: dict[str, object] = {}

    def capture(command, **kwargs):
        observed["command"] = command
        observed["inputs"] = kwargs["inputs"]

    monkeypatch.setattr("chronoclade.lineage._run_command", capture)

    _run_dated_tree(files, sequence_length=1234, force=False)

    command = observed["command"]
    assert isinstance(command, list)
    assert str(files.clock_dir / "rerooted.newick") in command
    assert "--keep-root" in command
    assert "--reroot" not in command
    assert observed["inputs"] == (files.clock_dir / "rerooted.newick", files.metadata)


def test_corrected_stage_and_insufficient_date_finish_skip_temporal_analysis(
    tmp_path: Path, monkeypatch
) -> None:
    from types import SimpleNamespace

    from chronoclade.lineage import _run_lineage
    from chronoclade.metadata import Sample

    samples = [
        Sample(
            sample_id,
            tmp_path / f"{sample_id}.fasta",
            date,
            "UK",
            "Klebsiella pneumoniae",
            "ST147",
            "local",
        )
        for sample_id, date in (("A", "2020"), ("B", "2021"), ("C", ""))
    ]
    observed = []

    def run_case(monkeypatch, output: Path, mode: str):
        monkeypatch.setattr(
            "chronoclade.lineage.select_reference",
            lambda members: SimpleNamespace(assembly=tmp_path / "reference.fa"),
        )
        monkeypatch.setattr("chronoclade.lineage._write_lineage_inputs", lambda *args: None)
        monkeypatch.setattr("chronoclade.lineage.context_evidence", lambda *args, **kwargs: {})
        monkeypatch.setattr(
            "chronoclade.lineage._run_core_phylogeny",
            lambda files, *args: (files.tree, files.filtered_alignment, {}),
        )
        monkeypatch.setattr("chronoclade.lineage.complete_alignment_sites", lambda path: 1000)
        monkeypatch.setattr("chronoclade.lineage.alignment_length", lambda path: 1200)
        monkeypatch.setattr(
            "chronoclade.lineage.write_recombination_evidence",
            lambda **kwargs: {"inferred_importation_intervals": 0},
        )
        monkeypatch.setattr(
            "chronoclade.lineage.build_public_health_evidence",
            lambda **kwargs: {"scenario": {"code": "indeterminate"}},
        )
        monkeypatch.setattr(
            "chronoclade.lineage.build_neighbourhood_evidence", lambda **kwargs: {}
        )
        monkeypatch.setattr("chronoclade.lineage._run_location_tree", lambda *args: None)
        monkeypatch.setattr(
            "chronoclade.lineage.build_country_network", lambda *args, **kwargs: {}
        )
        def corrected_report(report, *, directory):
            observed.append((mode, report["temporal_status"], report["temporal_signal_reason"]))
            path = directory / "corrected_report.html"
            path.write_text("corrected")
            return path

        monkeypatch.setattr("chronoclade.profile_report.write_corrected_report", corrected_report)
        monkeypatch.setattr("chronoclade.lineage.write_supporting_bundle", lambda path: path)

        def unexpected(*args, **kwargs):
            raise AssertionError("Temporal analysis must not run in this stage")

        monkeypatch.setattr("chronoclade.lineage._run_observed_clock", unexpected)
        monkeypatch.setattr("chronoclade.lineage._run_dated_tree", unexpected)
        monkeypatch.setattr("chronoclade.lineage._temporal_signal", unexpected)
        item = {"slug": mode, "species": "Klebsiella pneumoniae", "lineage": "ST147"}
        return _run_lineage(
            item,
            samples,
            output=output,
            threads=1,
            randomisation_jobs=1,
            randomisations=100,
            temporal_p_value=0.05,
            seed=1,
            force=False,
            context_manifest_rows=[],
            mode=mode,
            date_randomisation_method="root_to_tip",
        )

    corrected = run_case(monkeypatch, tmp_path / "corrected", "corrected")
    assert corrected["temporal_status"] == "not_assessed"
    assert corrected["temporal_signal_reason"] == "Dating is assessed only in the finish stage."
    assert not (tmp_path / "corrected" / "corrected" / "temporal_signal.json").exists()

    finish = run_case(monkeypatch, tmp_path / "finish", "full")
    assert finish["temporal_status"] == "not_assessed"
    assert "Fewer than three usable distinct collection dates" in finish["temporal_signal_reason"]
    assert not (tmp_path / "finish" / "full" / "temporal_signal.json").exists()
    assert len(observed) == 2
    import json

    for report in (corrected, finish):
        html = Path(report["outputs"]["html_report"])
        assert html.name == "corrected_report.html"
        assert html.read_text() == "corrected"
        saved = json.loads((html.parent / "report.json").read_text())
        assert saved["outputs"]["html_report"] == str(html)


def test_unknown_locations_are_missing_traits_not_countries(tmp_path: Path) -> None:
    import csv
    from chronoclade.lineage import _write_lineage_inputs
    from chronoclade.metadata import Sample

    samples = [
        Sample(
            str(i),
            tmp_path / f"{i}.fasta",
            "2020",
            location,
            "Klebsiella pneumoniae",
            "ST147",
            "local",
        )
        for i, location in enumerate(("Italy", "Unknown", "Public_context", "?", "N/A", "Missing"))
    ]
    _write_lineage_inputs(tmp_path, samples)
    states = list(csv.DictReader((tmp_path / "states.csv").open()))
    metadata = list(csv.DictReader((tmp_path / "metadata.csv").open()))
    assert [row["location"] for row in states] == ["Italy", "?", "?", "?", "?", "?"]
    assert [row["location"] for row in metadata] == [sample.location for sample in samples]


def test_single_location_skips_model_and_clears_stale_reconstruction(
    tmp_path: Path, monkeypatch
) -> None:
    from chronoclade.lineage import _run_location_tree

    files = LineageFiles.in_directory(tmp_path)
    files.states.write_text("sample_id,location\na,Italy\nb,?\nc,Italy\n")
    directory = tmp_path / "location"
    directory.mkdir()
    for name in ("annotated_tree.nexus", "confidence.csv", "GTR.txt"):
        (directory / name).write_text("stale")

    def unexpected(*args, **kwargs):
        raise AssertionError("Single-state model must not be run")

    monkeypatch.setattr("chronoclade.lineage._run_command", unexpected)
    _run_location_tree(files, files.tree, False, seed=7)
    assert not files.location_tree.exists()
    assert not (directory / "confidence.csv").exists()


def test_location_reconstruction_is_seeded_and_repairs_missing_confidence(
    tmp_path: Path, monkeypatch
) -> None:
    from chronoclade.lineage import _run_location_tree

    files = LineageFiles.in_directory(tmp_path)
    files.states.write_text("sample_id,location\na,Italy\nb,India\n")
    observed = {}

    def capture(command, **kwargs):
        observed["command"] = command
        observed.update(kwargs)

    monkeypatch.setattr("chronoclade.lineage._run_command", capture)
    _run_location_tree(files, files.tree, False, seed=7)
    assert observed["command"][observed["command"].index("--rng-seed") + 1] == "7"
    assert "--confidence" in observed["command"]
    assert observed["force"] is True


def _native_context_fixture(tmp_path, monkeypatch):
    import json
    from chronoclade.metadata import Sample
    from chronoclade.pathogenwatch import content_hash

    source = tmp_path / "source"
    source.mkdir()
    payload = {"rows": [], "focal_rows": [], "provenance": {}}
    digest = content_hash(payload)
    catalogue = source / "context_catalogue.json"
    catalogue.write_text(json.dumps(dict(payload, snapshot_sha256=digest)))
    monkeypatch.setattr(
        "chronoclade.context_geography.generate_context_geography", lambda *args, **kwargs: {}
    )
    sample = Sample(
        "pw_public", tmp_path / "genome.fa", "2026", "UK", "Escherichia coli", "ST131", "context"
    )
    manifest = [
        {
            "sample_id": sample.sample_id,
            "species": sample.species,
            "lineage": sample.lineage,
            "source": "pathogenwatch",
            "source_genome_id": "public",
            "catalogue_path": str(catalogue),
            "catalogue_sha256": digest,
        }
    ]
    destination = tmp_path / "lineage"
    destination.mkdir()
    return source, destination, sample, manifest


def test_lineage_copies_only_hash_verified_native_public_typing(tmp_path, monkeypatch):
    import hashlib
    import json
    from chronoclade.lineage import context_evidence

    source, destination, sample, manifest = _native_context_fixture(tmp_path, monkeypatch)
    data = b'{"schema_version":1,"assignments":[{"source_genome_id":"public"}]}\n'
    (source / "native_public_typing.json").write_bytes(data)
    (source / "context_selection.json").write_text(
        json.dumps({"native_public_typing_sha256": hashlib.sha256(data).hexdigest()})
    )
    context_evidence([sample], manifest, directory=destination)
    assert (destination / "native_public_typing.json").read_bytes() == data


def test_lineage_rejects_missing_or_changed_native_typing_artifact(tmp_path, monkeypatch):
    import json
    import pytest
    from chronoclade.errors import WorkflowError
    from chronoclade.lineage import context_evidence

    source, destination, sample, manifest = _native_context_fixture(tmp_path, monkeypatch)
    (source / "context_selection.json").write_text(
        json.dumps({"native_public_typing_sha256": "a" * 64})
    )
    with pytest.raises(WorkflowError, match="typing.*missing"):
        context_evidence([sample], manifest, directory=destination)
    (source / "native_public_typing.json").write_text("changed")
    (destination / "native_public_typing.json").write_text("stale result")
    with pytest.raises(WorkflowError, match="typing hash"):
        context_evidence([sample], manifest, directory=destination)
    assert not (destination / "native_public_typing.json").exists()


def test_lineage_removes_stale_typing_when_audit_does_not_declare_it(tmp_path, monkeypatch):
    from chronoclade.lineage import context_evidence

    source, destination, sample, manifest = _native_context_fixture(tmp_path, monkeypatch)
    (source / "context_selection.json").write_text("{}")
    (source / "native_public_typing.json").write_text("unclaimed file")
    (destination / "native_public_typing.json").write_text("stale result")
    context_evidence([sample], manifest, directory=destination)
    assert not (destination / "native_public_typing.json").exists()
    (destination / "native_public_typing.json").write_text("stale result")
    context_evidence([sample], [], directory=destination)
    assert not (destination / "native_public_typing.json").exists()


def test_future_collection_dates_excluded_and_clock_uses_validated_dates(tmp_path, monkeypatch):
    from datetime import date
    from chronoclade.lineage import _dated_members, _run_observed_clock
    from chronoclade.metadata import Sample

    valid = Sample("valid", tmp_path / "valid.fa", "2020", "UK", "Klebsiella pneumoniae", "ST147", "local")
    future = Sample("future", tmp_path / "future.fa", str(date.today().year + 10),
                    "UK", "Klebsiella pneumoniae", "ST147", "context")
    assert _dated_members([valid, future]) == [valid]
    assert future.collection_date == str(date.today().year + 10)
    files = LineageFiles.in_directory(tmp_path)
    files.metadata.write_text("raw original metadata")
    clock_dates = tmp_path / "clock_dates.csv"
    clock_dates.write_text("sample_id,collection_date\nvalid,2020\nfuture,\n")
    calls = []
    monkeypatch.setattr("chronoclade.lineage._run_command",
                        lambda command, **kwargs: calls.append((command, kwargs["inputs"])))
    _run_observed_clock(files, tmp_path / "tree.nwk", 100, False)
    _run_dated_tree(files, 100, False)
    assert all(str(clock_dates) in command and clock_dates in inputs
               for command, inputs in calls)
    assert files.metadata.read_text() == "raw original metadata"
