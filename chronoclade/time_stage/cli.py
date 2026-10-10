"""Lazy independent dating and ensemble-comparison dispatch."""

from pathlib import Path
from typing import Annotated
import typer
from rich.console import Console


def register(app):
    app.command("time")(time_command)
    app.command("time-compare")(comparison_command)


def time_command(
    tree: Annotated[Path, typer.Argument(help="Validated tree.json")],
    output: Annotated[Path, typer.Option("--out", "--output", "-o")] = Path("time_results"),
    force: bool = False,
    allow_unsupported: bool = False,
):
    """Date the saved tree when its temporal assessment supports dating."""
    from chronoclade.time_stage.workflow import run_time
    from chronoclade.errors import WorkflowError

    try:
        result = run_time(tree, output=output, force=force, allow_unsupported=allow_unsupported)
    except (ValueError, OSError, WorkflowError) as error:
        Console().print(f"Error: {error}", style="bold red", markup=False)
        raise typer.Exit(2) from error
    Console().print(f"Time manifest: {output.resolve() / 'time.json'}")
    Console().print(f"Dating status: {result['dating_status']}")


def comparison_command(
    ensemble: Annotated[Path, typer.Argument(help="Selection ensemble.json")],
    runs: Annotated[
        list[Path], typer.Option("--run", help="Complete time.json; repeat for each run")
    ],
    output: Annotated[Path, typer.Option("--out", "-o")] = Path("time_comparison"),
    cluster_snp_cutoff: Annotated[int, typer.Option(min=0)] = 10,
):
    """Compare independent selection dates/rates and shared-target clustering."""
    from chronoclade.time_stage.comparison import compare_time_runs
    from chronoclade.errors import WorkflowError

    try:
        compare_time_runs(ensemble, runs, output=output, cluster_snp_cutoff=cluster_snp_cutoff)
    except (ValueError, OSError, WorkflowError) as error:
        Console().print(f"Error: {error}", style="bold red", markup=False)
        raise typer.Exit(2) from error
    Console().print(f"Comparison: {output.resolve() / 'report.html'}")
