"""CLI entry point for optional protein embeddings; no model imports at startup."""

from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console


def esm2_command(
    fasta: Annotated[
        Path | None, typer.Argument(help="Prepared dataset.json or protein FASTA; omit with --dataset or saved embeddings")
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
    allele_catalogue: Annotated[
        Path | None, typer.Option(help="Scope-checked frozen allele DNA catalogue for genome analysis")
    ] = None,
    genome_analysis: Annotated[bool, typer.Option(help="Run the complete genome exploration/selection stage")] = False,
    lin_level: Annotated[int, typer.Option(min=5, max=7)] = 5,
    hiercc_level: Annotated[str | None, typer.Option()] = None,
    fetch_context: Annotated[bool, typer.Option(help="Retrieve matching public profiles before genome exploration")] = False,
    public_typing: Annotated[Path | None, typer.Option(help="Public profile export")] = None,
    cglin_export: Annotated[Path | None, typer.Option(help="Compatible lineage export")] = None,
    catalogues: Annotated[Path | None, typer.Option(help="Explicit canonical locus catalogues")] = None,
    selection_size: Annotated[int, typer.Option(min=0)] = 50,
    selection_runs: Annotated[int, typer.Option(min=1)] = 5,
    nearest_per_input: Annotated[int, typer.Option(min=1)] = 3,
    selection_seed: Annotated[int, typer.Option()] = 42,
    pin: Annotated[list[str] | None, typer.Option(help="Exact prepared sample ID to retain")] = None,
    neighbour_count: Annotated[int, typer.Option(min=1)] = 3,
    terminal_stop: Annotated[str, typer.Option(help="CDS terminal stop: require or allow_absent")] = "require",
    baseline: Annotated[bool, typer.Option(help="Include conventional cgMLST guide/network/date reports")] = True,
    group_threshold: Annotated[float, typer.Option(min=0, max=2, help="Exploratory cosine-distance component threshold")] = 0.0,
    reference_sample: Annotated[
        str | None, typer.Option(help="Fixed reference sample ID, chosen independently of dates")
    ] = None,
    panel_loci: Annotated[
        list[str] | None,
        typer.Option("--panel-locus", help="Repeat to specify a fixed locus panel"),
    ] = None,
) -> None:
    """Explore a prepared genome dataset, or extract optional protein embeddings."""
    from chronoclade.esm2 import EmbeddingError, run_embeddings
    from chronoclade.errors import WorkflowError

    console = Console()
    try:
        positional_dataset = fasta is not None and fasta.suffix.casefold() == ".json"
        if positional_dataset or genome_analysis or allele_catalogue:
            from chronoclade.esm2.workflow import run_esm2
            if positional_dataset:
                if dataset is not None:
                    raise ValueError("Supply the prepared dataset positionally or with --dataset, choosing one")
                dataset, fasta = fasta, None
            if dataset is None or fasta is not None:
                raise ValueError("Genome exploration requires a prepared dataset.json")
            source = dataset.expanduser()
            if fetch_context:
                from chronoclade.profile_context_provider import discover_profile_context
                prepared = discover_profile_context(
                    source, output.expanduser().parent / (output.name + "_context"),
                    lin_level=lin_level, hiercc_level=hiercc_level,
                    public_typing=public_typing, cglin_export=cglin_export, catalogues=catalogues,
                )
                source = prepared.dataset_manifest
            result = run_esm2(
                source, output.expanduser(),
                catalogue_manifest=allele_catalogue.expanduser() if allele_catalogue else None,
                embeddings_manifest=embeddings_manifest.expanduser() if embeddings_manifest else None,
                mapping_csv=sample_loci.expanduser() if sample_loci else None,
                model=model, checkpoint=checkpoint.expanduser() if checkpoint else None,
                allow_download=allow_download, device=device, cache_dir=cache_dir,
                token_budget=token_budget, max_length=max_length,
                expected_checkpoint_sha256=checkpoint_sha256,
                reference_sample_id=reference_sample, panel_loci=panel_loci,
                lin_level=lin_level, hiercc_level=hiercc_level,
                selection_size=selection_size, selection_runs=selection_runs,
                nearest_per_input=nearest_per_input, selection_seed=selection_seed,
                pinned_ids=pin or [], group_threshold=group_threshold, neighbour_count=neighbour_count,
                terminal_stop=terminal_stop, baseline=baseline,
            )
            console.print(f"Genome analysis manifest: {result.manifest_path}")
            console.print(f"Genome exploration report: {result.report_path}")
            return
        if fetch_context or public_typing or cglin_export or catalogues:
            raise ValueError("Public context options require prepared-genome exploration")
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
    except (EmbeddingError, OSError, ValueError, WorkflowError) as error:
        console.print(f"Error: {error}", style="bold red", markup=False)
        raise typer.Exit(code=2) from error
