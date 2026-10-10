"""The optional command must leave the base CLI usable without inference imports."""

import subprocess
import sys

from typer.testing import CliRunner

from chronoclade.cli import app


def test_cli_help_does_not_import_inference_runtime():
    completed = subprocess.run(
        [sys.executable, "-c", (
            "import sys; from chronoclade.cli import app; "
            "from typer.testing import CliRunner; "
            "result = CliRunner().invoke(app, ['esm2', '--help']); "
            "assert result.exit_code == 0, result.output; "
            "assert 'torch' not in sys.modules; assert 'esm' not in sys.modules"
        )],
        capture_output=True, text=True, check=False,
    )
    assert completed.returncode == 0, completed.stderr


def test_esm2_cli_reports_unavailable_explicit_checkpoint(tmp_path):
    fasta = tmp_path / "protein.fasta"
    fasta.write_text(">protein\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ\n")
    result = CliRunner().invoke(app, [
        "esm2", str(fasta), "--checkpoint", str(tmp_path / "missing.pt"),
        "--output", str(tmp_path / "out"),
    ])
    assert result.exit_code == 2
    assert "Error:" in result.output
    assert not (tmp_path / "out" / "manifest.json").exists()
