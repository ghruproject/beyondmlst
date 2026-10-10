"""Lazy CLI entry point for independent frozen-input preparation."""

from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console


def prepare_command(
    input_path: Annotated[
        str,
        typer.Argument(
            help="Profile JSON, dataset, collection UUID/URL, or accession/assembly metadata file"
        ),
    ],
    output: Annotated[Path, typer.Option("--out", "--output", "-o")] = Path("prepared"),
    catalogues: Annotated[
        Path | None, typer.Option(help="Explicit locus catalogue JSON for frozen profile records")
    ] = None,
    input_kind: Annotated[
        str, typer.Option(help="auto, profiles, dataset, collection, accessions or assemblies")
    ] = "auto",
    species: Annotated[
        str | None,
        typer.Option(help="Fill missing species in frozen records; conflicting values fail"),
    ] = None,
    metadata: Annotated[
        Path | None, typer.Option(help="Exact-ID metadata overrides for live queries")
    ] = None,
    typing_config: Annotated[
        Path | None, typer.Option(help="Configured native callers and reference databases")
    ] = None,
    query_typing: Annotated[
        Path | None, typer.Option(help="Frozen exact query typing JSON")
    ] = None,
    public_typing: Annotated[
        Path | None, typer.Option(help="Frozen provider typing assignments")
    ] = None,
    cglin_export: Annotated[
        Path | None, typer.Option(help="Authoritative frozen cgLIN export")
    ] = None,
    enrich_metadata: Annotated[
        bool,
        typer.Option(
            "--enrich-metadata/--no-enrich-metadata",
            help="Fill missing metadata from exact ENA links",
        ),
    ] = True,
) -> None:
    """Resolve input typing and publish a portable dataset with a readiness audit."""
    from chronoclade.prepare_stage import run_prepare

    console = Console()
    try:
        result = run_prepare(
            input_path,
            output,
            catalogues=catalogues,
            input_kind=input_kind,
            species=species,
            metadata=metadata,
            typing_config=typing_config,
            query_typing=query_typing,
            public_typing=public_typing,
            cglin_export=cglin_export,
            enrich_metadata=enrich_metadata,
        )
    except (OSError, ValueError) as error:
        console.print(f"Error: {error}", style="bold red", markup=False)
        raise typer.Exit(code=2) from error
    console.print(f"Dataset manifest: {result.dataset_manifest}", markup=False)
    console.print(f"Preparation audit: {result.report_path}", markup=False)
