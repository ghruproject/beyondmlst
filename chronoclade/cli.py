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
from chronoclade.esm2_cli import esm2_command
from chronoclade.metadata import MetadataError, group_samples, read_metadata
from chronoclade.workflow import (
    WorkflowError,
    native_platform_supported,
    plan,
    tool_status,
)

app = typer.Typer(
    name="chronoclade",
    help="Recombination-aware temporal and contextual bacterial phylogenetics.",
    no_args_is_help=True,
)
console = Console()

# Optional inference dependencies are imported only when this command runs.
app.command("esm2")(esm2_command)


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


@app.command("setup-typing")
def setup_typing_command(
    destination: Annotated[
        Path,
        typer.Option("--destination", help="Pinned native typing tools and reference databases"),
    ] = Path("~/.local/share/chronoclade/typing"),
    species: Annotated[
        list[str] | None,
        typer.Option("--species", help="klebsiella or ecoli; repeat to install both"),
    ] = None,
    prepare_db: Annotated[
        bool,
        typer.Option(
            "--prepare-db",
            help="Also download and index publicly available cgMLST alleles; may take hours",
        ),
    ] = False,
    pasteur_secrets: Annotated[
        Path | None,
        typer.Option(
            "--pasteur-secrets",
            help="Protected Pasteur credentials file for a Klebsiella LIN database build",
        ),
    ] = None,
    enterobase_key_file: Annotated[
        Path | None,
        typer.Option(
            "--enterobase-key-file",
            help="Protected EnteroBase key file for an E. coli HierCC database build",
        ),
    ] = None,
) -> None:
    """Install Pathogenwatch typing tools via uv and report database readiness."""
    from chronoclade.typing_setup import TypingSetupError, setup_typing

    try:
        path, report = setup_typing(
            destination,
            species=species or ["klebsiella", "ecoli"],
            prepare_db=prepare_db,
            pasteur_secrets=pasteur_secrets,
            enterobase_key_file=enterobase_key_file,
        )
    except (TypingSetupError, OSError) as error:
        _fail(error)
    console.print(f"Typing setup: {path}")
    console.print(json.dumps(report["databases"], indent=2))


@app.command("type-queries")
def type_queries_command(
    metadata: Annotated[Path, typer.Argument(help="Query genome metadata CSV")],
    typing_config: Annotated[
        Path,
        typer.Option(
            "--typing-config", help="Pinned tools and reference databases from setup-typing"
        ),
    ],
    output: Annotated[Path, typer.Option("--output", "-o")] = Path("chronoclade_typing"),
) -> None:
    """Call cgMLST and Klebsiella LIN or E. coli HierCC locally."""
    from chronoclade.query_typing import QueryTypingError, type_query_assemblies

    try:
        results = type_query_assemblies(
            read_metadata(metadata), typing_config, output.expanduser().resolve()
        )
    except (QueryTypingError, MetadataError, OSError) as error:
        _fail(error)
    console.print(
        f"Typed {len(results)} queries; results: {output.expanduser().resolve() / 'query_typing.json'}"
    )


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
    catalogue: Annotated[
        Path | None,
        typer.Option(
            "--catalogue", help="Frozen Pathogenwatch JSON catalogue; no live metadata requests"
        ),
    ] = None,
    cglin_export: Annotated[
        Path | None,
        typer.Option("--cglin-export", help="Validated cgLIN export for the same scheme"),
    ] = None,
    focal_crosswalk: Annotated[
        Path | None,
        typer.Option("--focal-crosswalk", help="Focal sample/accession to Pathogenwatch ID CSV"),
    ] = None,
    typing_config: Annotated[
        Path | None,
        typer.Option(
            "--typing-config",
            help="Native Pathogenwatch tools and frozen databases; type queries and downloaded context",
        ),
    ] = None,
    query_typing: Annotated[
        Path | None,
        typer.Option(
            "--query-typing", help="Frozen type-queries results matching the focal assembly hashes"
        ),
    ] = None,
    public_typing: Annotated[
        Path | None,
        typer.Option(
            "--public-typing",
            help="Frozen cgMLST/cgLIN/HierCC assignments keyed by exact Pathogenwatch IDs",
        ),
    ] = None,
    cglin_depth: Annotated[
        int,
        typer.Option(
            "--cglin-depth",
            min=1,
            max=10,
            help="Full cgLIN prefix depth prioritised within the same-ST pool",
        ),
    ] = 7,
    hiercc_level: Annotated[
        str,
        typer.Option(
            "--hiercc-level", help="HierCC level prioritised within the same-ST pool, e.g. HC1100"
        ),
    ] = "HC1100",
    refresh_catalogue: Annotated[
        bool,
        typer.Option(
            "--refresh-catalogue", help="Explicitly replace the cached public metadata snapshot"
        ),
    ] = False,
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
            typing_config=typing_config,
            query_typing=query_typing,
            public_typing=public_typing,
            cglin_depth=cglin_depth,
            hiercc_level=hiercc_level,
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
    metadata: Annotated[
        Path | None,
        typer.Argument(
            help="Assembly/accession metadata CSV; optional with --collection or --accessions"
        ),
    ] = None,
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
        typer.Option(
            "--mode", help="Stopping stage: fast profiles, full corrected tree, finish dating"
        ),
    ] = "fast",
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
    collection: Annotated[
        str | None, typer.Option("--collection", help="Pathogenwatch collection UUID or full URL")
    ] = None,
    accessions: Annotated[
        Path | None, typer.Option("--accessions", help="Query accession list or CSV")
    ] = None,
    species: Annotated[
        str | None, typer.Option("--species", help="Declared species for new assemblies")
    ] = None,
    catalogue: Annotated[
        Path | None, typer.Option("--catalogue", help="Frozen public metadata catalogue")
    ] = None,
    public_typing: Annotated[
        Path | None, typer.Option("--public-typing", help="Frozen public cgMLST/group assignments")
    ] = None,
    query_typing: Annotated[
        Path | None, typer.Option("--query-typing", help="Verified query typing JSON")
    ] = None,
    typing_config: Annotated[
        Path | None,
        typer.Option(
            "--typing-config", help="Native caller and reference configuration for new queries"
        ),
    ] = None,
    cglin_export: Annotated[
        Path | None, typer.Option("--cglin-export", help="Frozen public cgLIN assignments")
    ] = None,
    profile_limit: Annotated[
        int,
        typer.Option(
            "--profile-limit",
            min=0,
            help="Public comparison limit for inputs without adaptive cgLIN context; 0 disables public context",
        ),
    ] = 500,
    tree_limit: Annotated[
        int, typer.Option("--tree-limit", min=2, hidden=True,
                          help="Legacy option; staged trees use the shared --context-size selection")
    ] = 80,
    lin_min_context: Annotated[
        int, typer.Option("--lin-min-context", min=1,
                          help="Public samples required per input LIN prefix before narrowing context")
    ] = 20,
    context_size: Annotated[
        int,
        typer.Option(
            "--context-size", min=0, help="Public genomes per dataset in the shared tree selection; inputs are additional"
        ),
    ] = 50,
    nearest_per_query: Annotated[int, typer.Option("--nearest-per-query", min=1)] = 3,
    include_genome: Annotated[
        list[str] | None,
        typer.Option("--include-genome", help="Require a contextual genome; repeat for more"),
    ] = None,
    bootstrap_replicates: Annotated[int, typer.Option("--profile-bootstraps", min=0, max=200)] = 10,
    group_distance: Annotated[
        float,
        typer.Option(
            "--group-distance",
            min=0,
            max=1,
            help="Exploratory complete-linkage cgMLST mismatch fraction",
        ),
    ] = 0.02,
    force: Annotated[bool, typer.Option("--force", help="Rerun completed stages")] = False,
    dry_run: Annotated[
        bool, typer.Option("--dry-run", help="Validate and print the plan only")
    ] = False,
) -> None:
    """Run cumulative profile, corrected-tree and dating stages."""
    from chronoclade.profile_inputs import ProfileInputError
    from chronoclade.staged_workflow import run_staged_workflow

    try:
        normalized_method = date_randomisation_method.replace("-", "_")
        if normalized_method not in {"root_to_tip", "full_tree"}:
            raise WorkflowError("--date-randomisation-method must be root-to-tip or full-tree")
        result = run_staged_workflow(
            metadata,
            output=output,
            mode=mode,
            collection=collection,
            accessions=accessions,
            species=species,
            catalogue=catalogue,
            public_typing=public_typing,
            query_typing=query_typing,
            typing_config=typing_config,
            cglin_export=cglin_export,
            profile_limit=profile_limit,
            tree_limit=tree_limit,
            lin_min_context=lin_min_context,
            context_size=context_size,
            nearest_per_query=nearest_per_query,
            include_genomes=include_genome,
            bootstrap_replicates=bootstrap_replicates,
            distance_threshold=group_distance,
            threads=threads,
            lineage_jobs=lineage_jobs,
            randomisation_jobs=randomisation_jobs,
            randomisations=date_randomisations,
            temporal_p_value=temporal_p_value,
            min_samples=min_samples,
            seed=seed,
            force=force,
            date_randomisation_method=normalized_method,
            context_manifest=context_manifest,
            dry_run=dry_run,
        )
    except (MetadataError, WorkflowError, ProfileInputError, ValueError, OSError) as error:
        _fail(error)
    if dry_run:
        console.print(json.dumps(result, indent=2))
    else:
        console.print(
            f"Completed through {mode}; reports: {output.expanduser().resolve() / 'index.html'}"
        )


if __name__ == "__main__":
    app()
