"""CLI entry point for optional protein embeddings; no model imports at startup."""

from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console


def esm2_command(
    fasta: Annotated[
        Path | None, typer.Argument(help="Protein allele FASTA; omit when reusing saved embeddings")
    ] = None,
    output: Annotated[Path, typer.Option("--output", "--out", "-o")] = Path("esm2_results"),
    model: Annotated[str, typer.Option(help="ESM2 model: 8M or 35M")] = "8M",
    device: Annotated[str, typer.Option(help="auto, cpu, mps or cuda")] = "auto",
    checkpoint: Annotated[
        Path | None, typer.Option(help="Local checkpoint for offline inference")
    ] = None,
    allow_download: Annotated[
        bool,
        typer.Option(help="Allow upstream checkpoint retrieval if no local checkpoint is given"),
    ] = False,
    cache_dir: Annotated[
        Path | None, typer.Option(help="Reusable embedding cache directory")
    ] = None,
    token_budget: Annotated[
        int, typer.Option(min=3, help="Maximum padded tokens per batch")
    ] = 4096,
    max_length: Annotated[
        int, typer.Option(min=1, help="Reject longer proteins; never truncate")
    ] = 1022,
    checkpoint_sha256: Annotated[
        str | None, typer.Option(help="Expected SHA256 of the local checkpoint")
    ] = None,
    dataset: Annotated[
        Path | None, typer.Option(help="Prepared dataset.json for the temporal comparison report")
    ] = None,
    sample_loci: Annotated[
        Path | None, typer.Option(help="CSV linking sample_id,locus,record_id to FASTA records")
    ] = None,
    embeddings_manifest: Annotated[
        Path | None, typer.Option(help="Reuse saved embeddings.json without loading ESM2")
    ] = None,
    reference_sample: Annotated[
        str | None, typer.Option(help="Fixed reference sample ID, chosen independently of dates")
    ] = None,
    panel_loci: Annotated[
        list[str] | None,
        typer.Option("--panel-locus", help="Repeat to specify a fixed locus panel"),
    ] = None,
) -> None:
    """Embed proteins and optionally compare embedding/date and cgMLST root-to-tip views."""
    from chronoclade.esm2 import EmbeddingError, run_embeddings

    console = Console()
    try:
        if (fasta is None) == (embeddings_manifest is None):
            raise ValueError("Supply a protein FASTA or --embeddings-manifest, choosing one")
        if (dataset is None) != (sample_loci is None):
            raise ValueError("--dataset and --sample-loci must be supplied together")
        if dataset is None and (embeddings_manifest or reference_sample or panel_loci):
            raise ValueError("Saved-vector analysis and reference/panel options require --dataset")
        output = output.expanduser()
        if fasta is not None:
            result = run_embeddings(
                fasta.expanduser(),
                output,
                model=model,
                device=device,
                checkpoint=checkpoint.expanduser() if checkpoint else None,
                allow_download=allow_download,
                cache_dir=cache_dir.expanduser() if cache_dir else None,
                token_budget=token_budget,
                max_length=max_length,
                expected_checkpoint_sha256=checkpoint_sha256,
            )
            embeddings_manifest = result.manifest_path
            console.print(f"Embedding manifest: {result.manifest_path}")
            console.print(f"Protein vectors: {result.vectors_path}")
        if dataset is not None:
            from chronoclade.esm2.temporal_workflow import run_temporal_report

            report = run_temporal_report(
                dataset.expanduser(),
                embeddings_manifest.expanduser(),
                sample_loci.expanduser(),
                output,
                reference_sample_id=reference_sample,
                panel_loci=panel_loci,
            )
            console.print(f"Temporal comparison report: {report.report_path}")
    except (EmbeddingError, OSError, ValueError) as error:
        console.print(f"Error: {error}", style="bold red", markup=False)
        raise typer.Exit(code=2) from error
