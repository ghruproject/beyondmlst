"""Command-line interface for ChronoClade."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

from chronoclade import __version__
from chronoclade.context import ContextError, prepare_context
from chronoclade.metadata import MetadataError, group_samples, read_metadata
from chronoclade.workflow import (
    WorkflowError,
    native_platform_supported,
    plan,
    run_workflow,
    tool_status,
)

app = typer.Typer(
    name="chronoclade",
    help="Recombination-aware temporal and contextual bacterial phylogenetics.",
    no_args_is_help=True,
)
console = Console()


def _fail(error: Exception) -> None:
    console.print(f"[bold red]Error:[/bold red] {error}")
    raise typer.Exit(code=2)


@app.command()
def version() -> None:
    """Print the installed version."""

    console.print(__version__)


@app.command()
def preflight() -> None:
    """Check that external workflow tools are available."""

    table = Table(title="ChronoClade preflight")
    table.add_column("Tool")
    table.add_column("Status")
    ok = True
    for tool, path in tool_status().items():
        if path:
            table.add_row(tool, f"[green]{path}[/green]")
        else:
            table.add_row(tool, "[red]missing[/red]")
            ok = False
    if native_platform_supported():
        table.add_row("native platform", "[green]supported[/green]")
    else:
        table.add_row("native platform", "[red]unsupported by the Pixi environment[/red]")
        ok = False
    console.print(table)
    if not ok:
        raise typer.Exit(code=2)


@app.command()
def validate(
    metadata: Annotated[Path, typer.Argument(help="Input metadata CSV")],
    min_samples: Annotated[
        int, typer.Option("--min-samples", min=3, help="Minimum genomes per lineage")
    ] = 10,
) -> None:
    """Validate files and show the lineage analysis plan."""

    try:
        samples = read_metadata(metadata)
        items = plan(samples, min_samples=min_samples)
    except (MetadataError, WorkflowError) as error:
        _fail(error)
    console.print(f"[green]Validated {len(samples)} samples across {len(items)} lineages.[/green]")
    _print_plan(items)


def _print_plan(items: list[dict[str, object]]) -> None:
    table = Table(title="Lineage plan")
    for heading in ("Species", "Lineage", "Samples", "Dates", "Reference", "Status"):
        table.add_column(heading)
    for item in items:
        table.add_row(
            str(item["species"]),
            str(item["lineage"]),
            str(item["sample_count"]),
            str(item["distinct_dates"]),
            str(item["reference"]),
            str(item["status"]),
        )
    console.print(table)


@app.command("prepare-context")
def prepare_context_command(
    metadata: Annotated[Path, typer.Argument(help="Focal-isolate metadata CSV")],
    scheme: Annotated[
        str,
        typer.Option("--scheme", help="MLST scheme for Pathogenwatch context"),
    ],
    output: Annotated[
        Path, typer.Option("--output", "-o", help="Context preparation directory")
    ] = Path("chronoclade_context"),
    species: Annotated[
        str | None,
        typer.Option("--species", help="Species to prepare when metadata contains several"),
    ] = None,
    lineage: Annotated[
        str | None,
        typer.Option("--lineage", help="Lineage to prepare when metadata contains several"),
    ] = None,
    st: Annotated[
        str | None,
        typer.Option("--st", help="MLST sequence type; inferred from an ST-prefixed lineage"),
    ] = None,
    cache_dir: Annotated[
        Path, typer.Option("--cache-dir", help="Source-ID/content-hash assembly cache")
    ] = Path("~/.cache/chronoclade/context"),
    country: Annotated[
        list[str] | None,
        typer.Option("--country", help="Country prefix to retain; repeat for several countries"),
    ] = None,
    year_from: Annotated[int | None, typer.Option("--year-from", min=1800, max=2200)] = None,
    year_to: Annotated[int | None, typer.Option("--year-to", min=1800, max=2200)] = None,
    host: Annotated[str | None, typer.Option("--host", help="Host substring filter")] = None,
    isolation_source: Annotated[
        str | None,
        typer.Option("--isolation-source", help="Isolation-source substring filter"),
    ] = None,
    candidate_pool: Annotated[
        int,
        typer.Option(
            "--candidate-pool",
            min=1,
            help="Maximum balanced same-ST pool downloaded before SKA screening",
        ),
    ] = 500,
    max_context: Annotated[
        int, typer.Option("--max-context", min=1, help="Maximum context genomes retained")
    ] = 150,
    nearest_per_focal: Annotated[
        int,
        typer.Option(
            "--nearest-per-focal",
            min=1,
            help="Nearest candidate set considered for each focal isolate",
        ),
    ] = 3,
    threads: Annotated[int, typer.Option("--threads", "-t", min=1)] = 4,
    catalogue: Annotated[Path | None, typer.Option(
        "--catalogue", help="Frozen Pathogenwatch JSON catalogue; no live metadata requests"
    )] = None,
    cglin_export: Annotated[Path | None, typer.Option(
        "--cglin-export", help="Validated cgLIN export for the same scheme"
    )] = None,
    focal_crosswalk: Annotated[Path | None, typer.Option(
        "--focal-crosswalk", help="Focal sample/accession to Pathogenwatch ID CSV"
    )] = None,
    refresh_catalogue: Annotated[bool, typer.Option(
        "--refresh-catalogue", help="Explicitly replace the cached public metadata snapshot"
    )] = False,
    seed: Annotated[int, typer.Option("--seed")] = 20260818,
    dry_run: Annotated[
        bool,
        typer.Option(
            "--dry-run", help="Discover and freeze the candidate pool without downloading"
        ),
    ] = False,
) -> None:
    """Fetch and select reproducible public context for one species/ST."""

    if year_from is not None and year_to is not None and year_from > year_to:
        _fail(ValueError("--year-from cannot be later than --year-to"))
    try:
        samples = read_metadata(metadata)
        groups = group_samples(samples)
        if species is None and lineage is None and len(groups) == 1:
            (selected_species, selected_lineage), focal = next(iter(groups.items()))
        else:
            if species is None or lineage is None:
                raise ContextError(
                    "Metadata contains several lineages; provide both --species and --lineage"
                )
            selected_species, selected_lineage = species, lineage
            focal = groups.get((selected_species, selected_lineage), [])
        if not focal:
            raise ContextError(
                f"No focal samples found for {selected_species} / {selected_lineage}"
            )
        selected_st = st
        if selected_st is None and selected_lineage.upper().startswith("ST"):
            selected_st = selected_lineage[2:]
        if not selected_st:
            raise ContextError("Could not infer the sequence type; provide --st")
        audit = prepare_context(
            focal,
            species=selected_species,
            lineage=selected_lineage,
            scheme=scheme,
            st=selected_st,
            output=output,
            cache_dir=cache_dir,
            countries=country,
            year_from=year_from,
            year_to=year_to,
            host=host,
            isolation_source=isolation_source,
            candidate_pool=candidate_pool,
            max_context=max_context,
            nearest_per_focal=nearest_per_focal,
            seed=seed,
            threads=threads,
            dry_run=dry_run,
            catalogue=catalogue,
            cglin_export=cglin_export,
            focal_crosswalk=focal_crosswalk,
            refresh_catalogue=refresh_catalogue,
        )
    except (ContextError, MetadataError, OSError) as error:
        _fail(error)
    if dry_run:
        console.print(
            f"[green]Prepared a dry-run pool of {audit['candidate_pool']} same-ST candidates.[/green]"
        )
        console.print(f"Audit: {output.expanduser().resolve() / 'context_selection.json'}")
        return
    console.print(
        f"[green]Selected {audit['selected_contexts']} context genomes; "
        f"focal-neighbour coverage {float(audit['focal_neighbour_coverage']):.0%}.[/green]"
    )
    console.print(f"Combined metadata: {output.expanduser().resolve() / 'combined_metadata.csv'}")


@app.command()
def run(
    metadata: Annotated[Path, typer.Argument(help="Input metadata CSV")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output directory")] = Path(
        "chronoclade_results"
    ),
    threads: Annotated[int, typer.Option("--threads", "-t", min=1)] = 4,
    lineage_jobs: Annotated[
        int,
        typer.Option(
            "--lineage-jobs",
            min=1,
            help="Maximum lineages to analyse concurrently",
        ),
    ] = 2,
    randomisation_jobs: Annotated[
        int,
        typer.Option(
            "--randomisation-jobs",
            min=1,
            help="Maximum concurrent TreeTime permutations per lineage",
        ),
    ] = 4,
    date_randomisations: Annotated[
        int,
        typer.Option(
            "--date-randomisations",
            min=0,
            help="Number of tip-date permutations; use at least 100 for analysis",
        ),
    ] = 100,
    mode: Annotated[
        str,
        typer.Option("--mode", help="Analysis mode: full or fast"),
    ] = "full",
    date_randomisation_method: Annotated[
        str,
        typer.Option(
            "--date-randomisation-method",
            help="root-to-tip (fast permutation screen) or full-tree (complete TreeTime refits)",
        ),
    ] = "root-to-tip",
    temporal_p_value: Annotated[float, typer.Option("--temporal-p-value", min=0.0, max=1.0)] = 0.05,
    min_samples: Annotated[int, typer.Option("--min-samples", min=3)] = 10,
    seed: Annotated[int, typer.Option("--seed")] = 20260818,
    context_manifest: Annotated[
        Path | None,
        typer.Option(
            "--context-manifest",
            help="Frozen context_manifest.tsv produced by prepare-context",
        ),
    ] = None,
    force: Annotated[bool, typer.Option("--force", help="Rerun completed stages")] = False,
    dry_run: Annotated[
        bool, typer.Option("--dry-run", help="Validate and print the plan only")
    ] = False,
) -> None:
    """Run lineage-specific recombination and temporal analyses."""

    try:
        normalized_method = date_randomisation_method.replace("-", "_")
        if mode not in {"full", "fast"}:
            raise WorkflowError("--mode must be full or fast")
        if normalized_method not in {"root_to_tip", "full_tree"}:
            raise WorkflowError(
                "--date-randomisation-method must be root-to-tip or full-tree"
            )
        if mode == "fast" and normalized_method != "root_to_tip":
            raise WorkflowError("--mode fast only supports root-to-tip randomisation")
        samples = read_metadata(metadata)
        items = plan(samples, min_samples=min_samples)
        if dry_run:
            _print_plan(items)
            console.print(json.dumps(items, indent=2))
            return
        summary = run_workflow(
            samples,
            output=output,
            threads=threads,
            lineage_jobs=lineage_jobs,
            randomisation_jobs=randomisation_jobs,
            randomisations=date_randomisations,
            temporal_p_value=temporal_p_value,
            min_samples=min_samples,
            seed=seed,
            force=force,
            context_manifest=context_manifest,
            mode=mode,
            date_randomisation_method=normalized_method,
        )
    except (MetadataError, WorkflowError, OSError) as error:
        _fail(error)
    supported = sum(bool(item.get("temporal_signal_supported")) for item in summary["lineages"])
    result_word = "screen" if mode == "fast" else "configured temporal-signal test"
    console.print(
        f"[green]Completed {len(summary['lineages'])} lineage records; "
        f"{supported} passed the {result_word}.[/green]"
    )
    console.print(f"Summary: {(output.expanduser().resolve() / 'summary.json')}")


if __name__ == "__main__":
    app()
