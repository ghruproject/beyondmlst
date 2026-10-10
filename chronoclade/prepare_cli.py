"""Lazy CLI entry point for independent frozen-input preparation."""

from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console


def prepare_command(
    input_path: Annotated[Path, typer.Argument(help="Frozen profile JSON or dataset.json/bundle")],
    output: Annotated[Path, typer.Option("--out", "--output", "-o")] = Path("prepared"),
    catalogues: Annotated[
        Path | None, typer.Option(help="Explicit locus catalogue JSON for frozen profile records")
    ] = None,
    input_kind: Annotated[
        str, typer.Option(help="auto, profiles or dataset; imports are entirely offline")
    ] = "auto",
    species: Annotated[
        str | None,
        typer.Option(help="Fill missing species in frozen records; conflicting values fail"),
    ] = None,
) -> None:
    """Prepare a portable dataset and audit existing typing without provider retrieval."""
    from chronoclade.prepare_stage import run_prepare

    console = Console()
    try:
        result = run_prepare(
            input_path, output, catalogues=catalogues, input_kind=input_kind, species=species
        )
    except (OSError, ValueError) as error:
        console.print(f"Error: {error}", style="bold red", markup=False)
        raise typer.Exit(code=2) from error
    console.print(f"Dataset manifest: {result.dataset_manifest}", markup=False)
    console.print(f"Preparation audit: {result.report_path}", markup=False)
