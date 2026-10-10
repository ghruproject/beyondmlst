"""Arguments for the independent categorical-profile stage."""

from pathlib import Path
from typing import Annotated
import typer
from rich.console import Console


def cgmlst_command(
    dataset: Annotated[Path, typer.Argument(help="Prepared dataset.json")],
    output: Annotated[Path, typer.Option("--out", "--output", "-o")] = Path("cgmlst_results"),
    lin_level: Annotated[int, typer.Option(min=5, max=7)] = 5,
    hiercc_level: Annotated[str | None, typer.Option(help="Explicit E. coli HC level")] = None,
    subsamples: Annotated[int, typer.Option(min=1)] = 1,
    context_size: Annotated[int, typer.Option(min=0)] = 50,
    nearest_per_query: Annotated[int, typer.Option(min=1)] = 3,
    include: Annotated[
        list[str] | None, typer.Option(help="Pin a context ID; repeat for several")
    ] = None,
    fetch_context: Annotated[bool, typer.Option(help="Retrieve matching public profiles before analysis")] = False,
    public_typing: Annotated[Path | None, typer.Option(help="Public profile export")] = None,
    cglin_export: Annotated[Path | None, typer.Option(help="Compatible lineage export")] = None,
    catalogues: Annotated[Path | None, typer.Option(help="Explicit canonical locus catalogues")] = None,
    seed: int = 42,
    bootstrap_replicates: Annotated[int, typer.Option(min=0, max=200)] = 30,
    min_overlap: Annotated[float, typer.Option(min=0, max=1)] = 0.9,
    distance_threshold: Annotated[float, typer.Option(min=0, max=1)] = 0.02,
):
    """Analyse every matching profile in explicit frozen lineage blocks."""
    from chronoclade.cgmlst.workflow import run_cgmlst
    from chronoclade.errors import WorkflowError

    console = Console()
    try:
        source = dataset.expanduser()
        if fetch_context:
            from chronoclade.prepare_provider import discover_profile_context
            prepared = discover_profile_context(
                source, output.expanduser().parent / (output.name + "_context"),
                lin_level=lin_level, hiercc_level=hiercc_level,
                public_typing=public_typing, cglin_export=cglin_export, catalogues=catalogues,
            )
            source = prepared.dataset_manifest
        result = run_cgmlst(
            source,
            output.expanduser(),
            lin_level=lin_level,
            hiercc_level=hiercc_level,
            replicates=subsamples,
            context_size=context_size,
            nearest_per_query=nearest_per_query,
            include=include,
            seed=seed,
            bootstrap_replicates=bootstrap_replicates,
            min_overlap=min_overlap,
            distance_threshold=distance_threshold,
        )
    except (ValueError, OSError, WorkflowError) as error:
        console.print(f"Error: {error}", style="bold red", markup=False)
        raise typer.Exit(code=2) from error
    console.print(f"cgMLST manifest: {result.manifest_path}")
    console.print(f"cgMLST report: {result.report_path}")
