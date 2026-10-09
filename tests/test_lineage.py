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
