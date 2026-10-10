"""Real frozen MPS vectors must generate both views without inference dependencies."""

import json
from pathlib import Path
import subprocess
import sys

import pytest
from typer.testing import CliRunner

from chronoclade.cli import app
from chronoclade.datasets import (
    AlleleMatrix,
    LocusCatalogue,
    PreparedDataset,
    samples_from_records,
    write_dataset,
)

FIXTURE = Path(__file__).resolve().parents[1] / "validation/esm2/temporal"


def test_frozen_native_vectors_generate_report_without_torch(tmp_path):
    fixture = json.loads((FIXTURE / "fixture.json").read_text())
    cat = fixture["catalogue"]
    catalogue = LocusCatalogue(
        cat["scheme_id"],
        cat["scheme_version"],
        tuple(cat["loci"]),
        cat["source"],
        complete=cat["complete"],
    )
    samples = samples_from_records(fixture["records"])
    matrix = AlleleMatrix.from_profiles(
        catalogue,
        [row["sample_id"] for row in samples],
        fixture["calls"],
    )
    manifest = write_dataset(
        PreparedDataset(samples, (matrix,), parameters=fixture["parameters"]),
        tmp_path / "prepared",
    )
    out = tmp_path / "report"
    arguments = [
        "esm2",
        "--dataset",
        str(manifest),
        "--sample-loci",
        str(FIXTURE / "sample_loci.csv"),
        "--embeddings-manifest",
        str(FIXTURE / "embeddings.json"),
        "--out",
        str(out),
    ]
    script = (
        "import sys; from typer.testing import CliRunner; from chronoclade.cli import app; "
        f"result = CliRunner().invoke(app, {arguments!r}); "
        "assert result.exit_code == 0, result.output; "
        "assert 'torch' not in sys.modules; assert 'esm' not in sys.modules"
    )
    result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    record = json.loads((out / "temporal_analysis.json").read_text())
    assert record["status"] == "complete"
    assert record["parameters"]["reference_sample_id"] == "SIM_01"
    assert (out / record["report"]["path"]).is_file()
    evidence = json.loads((out / "temporal_diagnostics.json").read_text())
    # The report keeps the two underlying diagnostics in its evidence download.
    embedding = evidence["embedding"]
    samples = {row["sample_id"]: row for row in embedding["samples"]}
    assert samples["SIM_02"]["distance"] == pytest.approx(0, abs=1e-7)
    assert samples["SIM_03"]["distance"] > 0
    assert samples["SIM_05"]["distance"] > 0
    assert len(embedding["points"]) == 4
    assert embedding["counts"]["eligible_samples"] == 5


@pytest.mark.parametrize(
    "arguments,message",
    [
        ([], "Supply a protein FASTA"),
        (["--embeddings-manifest", "saved.json"], "require --dataset"),
        (["proteins.fasta", "--dataset", "dataset.json"], "supplied together"),
        (["proteins.fasta", "--embeddings-manifest", "saved.json"], "choosing one"),
    ],
)
def test_invalid_modes_fail_before_inference(arguments, message, tmp_path):
    result = CliRunner().invoke(app, ["esm2", *arguments, "--out", str(tmp_path / "out")])
    assert result.exit_code == 2
    assert message in result.output
    assert not (tmp_path / "out").exists()
