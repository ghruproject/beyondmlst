"""Lazy dispatch for the independent assembly-tree command."""

from pathlib import Path
from typing import Annotated
import os
import typer
from rich.console import Console


def register(app):
    app.command("tree")(tree_command)


def tree_command(
    selection: Annotated[Path, typer.Argument(help="Validated selection manifest")],
    output: Annotated[Path, typer.Option("--out", "--output", "-o")] = Path("tree_results"),
    threads: Annotated[int, typer.Option(min=1)] = 1,
    randomisations: Annotated[int, typer.Option(min=0)] = 20,
    randomisation_jobs: Annotated[int, typer.Option(min=1)] = 1,
    temporal_p_value: Annotated[float, typer.Option(min=0.000001, max=1)] = 0.05,
    seed: int = 42,
    assess_temporal: bool = True,
    force: bool = False,
    api_key_env: str = "PATHOGENWATCH_API_KEY",
    base_url: str = "https://pathogen.watch",
):
    """Build the exact selected assembly phylogeny and clustered temporal screen."""
    from chronoclade.tree_stage.workflow import run_tree
    from chronoclade.errors import WorkflowError

    try:
        result = run_tree(
            selection,
            output=output,
            threads=threads,
            randomisations=randomisations,
            randomisation_jobs=randomisation_jobs,
            temporal_p_value=temporal_p_value,
            seed=seed,
            assess_temporal=assess_temporal,
            force=force,
            api_key=_api_key(api_key_env),
            base_url=base_url,
        )
    except (ValueError, OSError, WorkflowError) as error:
        Console().print(f"Error: {error}", style="bold red", markup=False)
        raise typer.Exit(2) from error
    Console().print(f"Tree manifest: {output.resolve() / 'tree.json'}")
    Console().print(f"Temporal assessment: {result['temporal_assessment']['code']}")


def _api_key(environment):
    from chronoclade.pathogenwatch import load_api_key

    return os.environ.get(environment) or load_api_key()
