"""CLI entry point for optional protein embeddings; no model imports at startup."""

from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console


def esm2_command(
    fasta: Annotated[Path, typer.Argument(help="Protein allele FASTA (embedding foundation)")],
    output: Annotated[Path, typer.Option("--output", "--out", "-o")] = Path("esm2_results"),
    model: Annotated[str, typer.Option(help="ESM2 model: 8M or 35M")] = "8M",
    device: Annotated[str, typer.Option(help="auto, cpu, mps or cuda")] = "auto",
    checkpoint: Annotated[
        Path | None, typer.Option(help="Local checkpoint for offline inference")
    ] = None,
    allow_download: Annotated[
        bool, typer.Option(help="Allow upstream checkpoint retrieval if no local checkpoint is given")
    ] = False,
    cache_dir: Annotated[Path | None, typer.Option(help="Reusable embedding cache directory")] = None,
    token_budget: Annotated[int, typer.Option(min=3, help="Maximum padded tokens per batch")] = 4096,
    max_length: Annotated[int, typer.Option(min=1, help="Reject longer proteins; never truncate")] = 1022,
    checkpoint_sha256: Annotated[
        str | None, typer.Option(help="Expected SHA256 of the local checkpoint")
    ] = None,
) -> None:
    """Generate cached protein embeddings; genome analysis/selection is not yet implemented."""
    from chronoclade.esm2 import EmbeddingError, run_embeddings

    console = Console()
    try:
        result = run_embeddings(
            fasta.expanduser(), output.expanduser(), model=model, device=device,
            checkpoint=checkpoint.expanduser() if checkpoint else None,
            allow_download=allow_download,
            cache_dir=cache_dir.expanduser() if cache_dir else None,
            token_budget=token_budget, max_length=max_length,
            expected_checkpoint_sha256=checkpoint_sha256,
        )
    except (EmbeddingError, OSError, ValueError) as error:
        console.print(f"Error: {error}", style="bold red", markup=False)
        raise typer.Exit(code=2) from error
    console.print(f"Embedding manifest: {result.manifest_path}")
    console.print(f"Protein vectors: {result.vectors_path}")
