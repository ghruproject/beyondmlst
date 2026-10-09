"""Execution stages for one bacterial lineage."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

from chronoclade.coherence import screen_alignment
from chronoclade.errors import WorkflowError
from chronoclade.evidence import build_public_health_evidence
from chronoclade.metadata import Sample, select_reference
from chronoclade.report import (
    assess_temporal_signal,
    write_fast_lineage_report,
    write_lineage_report,
    write_supporting_bundle,
)
from chronoclade.recombination import run_phipack_screen
from chronoclade.recombination_report import write_recombination_evidence
from chronoclade.temporal import (
    file_sha256,
    run_date_randomisation,
    run_full_tree_date_randomisation,
)


@dataclass(frozen=True)
class LineageFiles:
    """Canonical paths for one lineage analysis."""

    directory: Path
    inputs: Path
    metadata: Path
    states: Path
    alignment: Path
    coherence_screen: Path
    clonal_coherence_screen: Path
    ska_prefix: Path
    ska_file: Path
    iqtree_prefix: Path
    starting_tree: Path
    clonalframe_prefix: Path
    tree: Path
    importations: Path
    filtered_alignment: Path
    phipack_profile: Path
    phipack_significant_blocks: Path
    phipack_summary: Path
    fast_iqtree_prefix: Path
    fast_tree: Path
    clock_dir: Path
    clock_file: Path
    temporal_signal: Path
    time_tree: Path
    location_tree: Path

    @classmethod
    def in_directory(cls, directory: Path) -> LineageFiles:
        return cls(
            directory=directory,
            inputs=directory / "inputs.tsv",
            metadata=directory / "metadata.csv",
            states=directory / "states.csv",
            alignment=directory / "core_alignment.fasta",
            coherence_screen=directory / "lineage_coherence.tsv",
            clonal_coherence_screen=directory / "clonal_lineage_coherence.tsv",
            ska_prefix=directory / "ska",
            ska_file=directory / "ska.skf",
            iqtree_prefix=directory / "iqtree",
            starting_tree=directory / "iqtree.treefile",
            clonalframe_prefix=directory / "clonalframeml",
            tree=directory / "clonalframeml.labelled_tree.newick",
            importations=directory / "clonalframeml.importation_status.txt",
            filtered_alignment=directory / "clonalframeml.filtered.fasta",
            phipack_profile=directory / "phipack_profile.tsv",
            phipack_significant_blocks=directory / "phipack_significant_blocks.tsv",
            phipack_summary=directory / "phipack_summary.json",
            fast_iqtree_prefix=directory / "iqtree_fast",
            fast_tree=directory / "iqtree_fast.treefile",
            clock_dir=directory / "clock",
            clock_file=directory / "clock" / "molecular_clock.txt",
            temporal_signal=directory / "temporal_signal.json",
            time_tree=directory / "timetree" / "timetree.nexus",
            location_tree=directory / "location" / "annotated_tree.nexus",
        )


def _write_csv(
    path: Path,
    fieldnames: list[str],
    rows: Iterable[dict[str, object]],
    *,
    delimiter: str = ",",
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def _write_lineage_inputs(directory: Path, samples: list[Sample]) -> tuple[Path, Path, Path]:
    inputs = directory / "inputs.tsv"
    inputs.write_text(
        "".join(f"{sample.sample_id}\t{sample.assembly}\n" for sample in samples),
        encoding="utf-8",
    )
    metadata = directory / "metadata.csv"
    fields = [
        "sample_id",
        "assembly",
        "collection_date",
        "location",
        "species",
        "lineage",
        "origin",
        "is_reference",
        "patient_id",
    ]
    _write_csv(
        metadata,
        fields,
        (
            {
                **asdict(sample),
                "assembly": str(sample.assembly),
            }
            for sample in samples
        ),
    )
    states = directory / "states.csv"
    _write_csv(
        states,
        ["sample_id", "location", "origin"],
        (
            {
                "sample_id": sample.sample_id,
                "location": sample.location,
                "origin": sample.origin,
            }
            for sample in samples
        ),
    )
    return inputs, metadata, states


def _run_command(
    command: list[str],
    *,
    log: Path,
    expected: Path | tuple[Path, ...],
    force: bool,
    cwd: Path | None = None,
    inputs: tuple[Path, ...] = (),
) -> None:
    expected_paths = (expected,) if isinstance(expected, Path) else expected
    fingerprint_path = log.with_suffix(".fingerprint.json")
    fingerprint = {
        "command": command,
        "inputs": [_file_fingerprint(path) for path in inputs],
    }
    previous = None
    if fingerprint_path.is_file():
        try:
            previous = json.loads(fingerprint_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            previous = None
    if all(path.exists() for path in expected_paths) and not force and previous == fingerprint:
        return
    for path in expected_paths:
        if path.is_file():
            path.unlink()
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        handle.write("COMMAND\n" + " ".join(command) + "\n\nOUTPUT\n")
        completed = subprocess.run(
            command,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
            cwd=cwd,
        )
    missing = [path for path in expected_paths if not path.exists()]
    if completed.returncode != 0 or missing:
        missing_text = ", ".join(str(path) for path in missing)
        raise WorkflowError(
            f"Command failed or did not create expected output(s): {missing_text}. See log: {log}"
        )
    temporary_fingerprint = fingerprint_path.with_suffix(".json.tmp")
    temporary_fingerprint.write_text(json.dumps(fingerprint, indent=2) + "\n", encoding="utf-8")
    temporary_fingerprint.replace(fingerprint_path)


def _file_fingerprint(path: Path) -> dict[str, object]:
    """Hash one direct stage input so stale outputs are not silently reused."""

    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        raise WorkflowError(f"Stage input does not exist: {resolved}")
    digest = hashlib.sha256()
    with resolved.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def alignment_length(path: Path) -> int:
    """Return the number of columns in an equal-length FASTA alignment."""

    lengths: list[int] = []
    current = 0
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                    current = 0
            else:
                current += len(line)
    if current:
        lengths.append(current)
    if not lengths:
        raise WorkflowError(f"Alignment contains no sequences: {path}")
    if len(set(lengths)) != 1:
        raise WorkflowError(f"Alignment sequences are not the same length: {path}")
    return lengths[0]


def complete_alignment_sites(path: Path) -> int:
    """Count columns containing only unambiguous A/C/G/T bases."""

    length = alignment_length(path)
    incomplete = bytearray(length)
    position = 0
    seen_header = False
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seen_header and position != length:
                    raise WorkflowError(f"Alignment sequences are not the same length: {path}")
                seen_header = True
                position = 0
                continue
            if not seen_header:
                raise WorkflowError(f"Alignment sequence appears before its FASTA header: {path}")
            for base in line.upper():
                if position >= length:
                    raise WorkflowError(f"Alignment sequences are not the same length: {path}")
                if base not in "ACGT":
                    incomplete[position] = 1
                position += 1
    complete = length - sum(incomplete)
    if complete == 0:
        raise WorkflowError(f"Alignment contains no complete A/C/G/T sites: {path}")
    return complete


def read_context_manifest(path: Path | None) -> list[dict[str, str]]:
    """Read a context-selection manifest for report provenance."""

    if path is None:
        return []
    path = path.expanduser().resolve()
    if not path.is_file():
        raise WorkflowError(f"Context manifest does not exist: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "species", "lineage"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise WorkflowError(
                "Context manifest must be tab-separated and contain sample_id, species and lineage"
            )
        rows = [dict(row) for row in reader]
        for row in rows:
            for key in ("assembly", "catalogue_path"):
                if row.get(key):
                    row[key] = str((path.parent / row[key]).resolve())
        return rows


def context_evidence(
    members: list[Sample], manifest_rows: list[dict[str, str]], *, directory: Path
) -> dict[str, object]:
    """Summarise contextual sampling and preserve the lineage manifest subset."""

    directory.mkdir(parents=True, exist_ok=True)
    origins: dict[str, int] = {}
    for member in members:
        key = member.origin.strip().lower()
        origins[key] = origins.get(key, 0) + 1
    contexts = [member for member in members if member.origin.strip().lower() == "context"]
    local = [member for member in members if member.origin.strip().lower() == "local"]
    context_ids = {member.sample_id for member in contexts}
    relevant_rows = [
        row
        for row in manifest_rows
        if row.get("species") == members[0].species
        and row.get("lineage") == members[0].lineage
        and row.get("sample_id") in context_ids
    ]
    nearest: list[dict[str, object]] = []
    for row in relevant_rows:
        raw_distance = row.get("min_ska_distance", "")
        try:
            distance: float | None = float(raw_distance) if raw_distance else None
        except ValueError:
            distance = None
        nearest.append(
            {
                "sample_id": row.get("sample_id", ""),
                "country": row.get("country", ""),
                "collection_date": row.get("collection_date", ""),
                "nearest_focal": row.get("nearest_focal", ""),
                "min_ska_distance": distance,
                "selection_reason": row.get("selection_reason", ""),
            }
        )
    nearest.sort(
        key=lambda row: (
            row["min_ska_distance"] is None,
            row["min_ska_distance"] if row["min_ska_distance"] is not None else float("inf"),
            str(row["sample_id"]),
        )
    )

    subset_path: Path | None = None
    if relevant_rows:
        subset_path = directory / "context_manifest.tsv"
        fieldnames = list(relevant_rows[0])
        _write_csv(subset_path, fieldnames, relevant_rows, delimiter="\t")
    geography = None
    catalogue_paths = {
        row.get("catalogue_path", "")
        for row in relevant_rows
        if row.get("source") == "pathogenwatch" and row.get("catalogue_path")
    }
    if len(catalogue_paths) > 1:
        raise WorkflowError("Context manifest refers to multiple frozen catalogues")
    if not catalogue_paths:
        shutil.rmtree(directory / "context_geography", ignore_errors=True)
        for name in ("context_catalogue.json", "context_selection.json"):
            (directory / name).unlink(missing_ok=True)
        if not relevant_rows:
            (directory / "context_manifest.tsv").unlink(missing_ok=True)
    if catalogue_paths:
        from chronoclade.context_geography import generate_context_geography
        from chronoclade.pathogenwatch import content_hash

        frozen_path = Path(next(iter(catalogue_paths)))
        if not frozen_path.is_file():
            raise WorkflowError("Frozen Pathogenwatch catalogue required by manifest is missing")
        payload = json.loads(frozen_path.read_text(encoding="utf-8"))
        expected = payload.pop("snapshot_sha256", "")
        if (
            not expected
            or content_hash(payload) != expected
            or any(
                row.get("catalogue_sha256") != expected
                for row in relevant_rows
                if row.get("source") == "pathogenwatch"
            )
        ):
            raise WorkflowError("Frozen Pathogenwatch catalogue hash does not match manifest")
        provenance = payload.get("provenance", {})
        geography = generate_context_geography(
            payload["rows"],
            directory / "context_geography",
            selected_source_ids=[
                row["source_genome_id"] for row in relevant_rows if row.get("source_genome_id")
            ],
            focal_rows=payload.get("focal_rows", []),
            scope={
                "description": f"Public same-ST Pathogenwatch catalogue: "
                f"{members[0].species} {members[0].lineage}",
                "snapshot": provenance.get("retrieved_at", "unknown"),
                "filters": "QC pass; accession-deduplicated; undated records retained",
            },
        )
        (directory / "context_catalogue.json").write_text(
            json.dumps({**payload, "snapshot_sha256": expected}, indent=2) + "\n", encoding="utf-8"
        )
        selection = frozen_path.parent / "context_selection.json"
        if selection.is_file():
            (directory / "context_selection.json").write_bytes(selection.read_bytes())
    return {
        "geography_available": geography is not None,
        "local_samples": len(local),
        "context_samples": len(contexts),
        "retrospective_samples": origins.get("retrospective", 0),
        "origin_counts": origins,
        "context_locations": sorted({member.location for member in contexts}),
        "manifest_available": bool(relevant_rows),
        "manifest_path": str(subset_path) if subset_path else None,
        "nearest_screening_contexts": nearest[:10],
        "interpretation": (
            "Contextual genomes are available for topology and relatedness assessment."
            if contexts
            else "No contextual genomes were supplied; introductions cannot be distinguished "
            "from local circulation."
        ),
    }


def _run_core_phylogeny(
    files: LineageFiles,
    members: list[Sample],
    reference: Sample,
    threads: int,
    force: bool,
    mode: str,
) -> tuple[Path, Path, dict[str, object] | None]:
    mapping_stages = [
        (
            [
                "ska",
                "build",
                "-f",
                str(files.inputs),
                "-o",
                str(files.ska_prefix),
                "--threads",
                str(threads),
            ],
            "ska_build.log",
            files.ska_file,
            (files.inputs, *(sample.assembly for sample in members)),
        ),
        (
            [
                "ska",
                "map",
                str(reference.assembly),
                str(files.ska_file),
                "-o",
                str(files.alignment),
                "--ambig-mask",
                "--repeat-mask",
                "--threads",
                str(threads),
            ],
            "ska_map.log",
            files.alignment,
            (reference.assembly, files.ska_file),
        ),
    ]
    for command, log_name, expected, inputs in mapping_stages:
        _run_command(
            command,
            log=files.directory / "logs" / log_name,
            expected=expected,
            force=force,
            cwd=files.directory,
            inputs=inputs,
        )

    coherence = screen_alignment(
        files.alignment,
        members,
        output=files.coherence_screen,
    )
    if coherence["flagged_samples"]:
        names = ", ".join(str(name) for name in coherence["flagged_samples"])
        raise WorkflowError(
            "Lineage-coherence screen found extreme raw-distance outlier(s): "
            f"{names}. Verify their accessions, species and lineage assignment, or split the "
            f"input lineage. Evidence: {files.coherence_screen}"
        )

    if mode == "fast":
        phipack = run_phipack_screen(
            alignment=files.alignment,
            reference=reference.assembly,
            profile_output=files.phipack_profile,
            significant_blocks_output=files.phipack_significant_blocks,
            summary_output=files.phipack_summary,
            threads=threads,
            force=force,
        )
        _run_command(
            [
                "iqtree",
                "-s",
                str(files.alignment),
                "-m",
                "GTR+G",
                "--fast",
                "-T",
                str(threads),
                "-pre",
                str(files.fast_iqtree_prefix),
                "-redo",
            ],
            log=files.directory / "logs" / "iqtree_fast.log",
            expected=files.fast_tree,
            force=force,
            cwd=files.directory,
            inputs=(files.alignment,),
        )
        return files.fast_tree, files.alignment, phipack

    phylogeny_stages = [
        (
            [
                "iqtree",
                "-s",
                str(files.alignment),
                "-m",
                "GTR+G",
                "-T",
                str(threads),
                "-pre",
                str(files.iqtree_prefix),
                "-redo",
            ],
            "iqtree.log",
            files.starting_tree,
            (files.alignment,),
        ),
        (
            [
                "ClonalFrameML",
                str(files.starting_tree),
                str(files.alignment),
                str(files.clonalframe_prefix),
                "-ignore_incomplete_sites",
                "true",
                "-output_filtered",
                "true",
                "-num_threads",
                str(threads),
            ],
            "clonalframeml.log",
            (files.tree, files.importations, files.filtered_alignment),
            (files.starting_tree, files.alignment),
        ),
    ]
    for command, log_name, expected, inputs in phylogeny_stages:
        _run_command(
            command,
            log=files.directory / "logs" / log_name,
            expected=expected,
            force=force,
            cwd=files.directory,
            inputs=inputs,
        )

    clonal_coherence = screen_alignment(
        files.filtered_alignment,
        members,
        output=files.clonal_coherence_screen,
    )
    if clonal_coherence["flagged_samples"]:
        names = ", ".join(str(name) for name in clonal_coherence["flagged_samples"])
        raise WorkflowError(
            "Lineage-coherence screen found extreme recombination-filtered distance "
            f"outlier(s): {names}. Verify their accessions and lineage assignment before "
            f"temporal analysis. Evidence: {files.clonal_coherence_screen}"
        )
    return files.tree, files.filtered_alignment, None


def _run_observed_clock(files: LineageFiles, tree: Path, sequence_length: int, force: bool) -> None:
    _run_command(
        [
            "treetime",
            "clock",
            "--tree",
            str(tree),
            "--dates",
            str(files.metadata),
            "--name-column",
            "sample_id",
            "--date-column",
            "collection_date",
            "--sequence-length",
            str(sequence_length),
            "--reroot",
            "least-squares",
            "--allow-negative-rate",
            "--clock-filter",
            "0",
            "--plot-rtt",
            "root_to_tip_regression.svg",
            "--outdir",
            str(files.clock_dir),
        ],
        log=files.directory / "logs" / "clock.log",
        expected=files.clock_file,
        force=force,
        inputs=(tree, files.metadata),
    )


def _temporal_signal(
    files: LineageFiles,
    members: list[Sample],
    sequence_length: int,
    randomisations: int,
    randomisation_jobs: int,
    seed: int,
    force: bool,
    tree: Path,
    method: str,
    observed_clock: Path,
) -> dict[str, object]:
    dates_sha256 = hashlib.sha256(
        json.dumps(
            [(sample.sample_id, sample.collection_date) for sample in members],
            separators=(",", ":"),
        ).encode()
    ).hexdigest()
    observed_sha256 = file_sha256(observed_clock)
    if not force and files.temporal_signal.exists():
        previous = json.loads(files.temporal_signal.read_text(encoding="utf-8"))
        if (
            previous.get("method", "root_to_tip") == method
            and int(previous.get("requested_randomisations", -1)) == randomisations
            and int(previous.get("seed", -1)) == seed
            and int(previous.get("sequence_length", -1)) == sequence_length
            and previous.get("tree_sha256") == file_sha256(tree)
            and previous.get("dates_sha256") == dates_sha256
            and previous.get("observed_clock_sha256") == observed_sha256
        ):
            return previous
    runner = run_full_tree_date_randomisation if method == "full_tree" else run_date_randomisation
    result = runner(
        tree=tree,
        sequence_length=sequence_length,
        samples=members,
        observed_clock=observed_clock,
        randomisations=randomisations,
        randomisation_jobs=randomisation_jobs,
        seed=seed,
        output=files.temporal_signal,
    )
    result.update(dates_sha256=dates_sha256, observed_clock_sha256=observed_sha256)
    files.temporal_signal.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def _run_dated_tree(
    files: LineageFiles,
    sequence_length: int,
    force: bool,
    tree: Path | None = None,
    *,
    reroot: bool = False,
) -> None:
    rooted_tree = files.clock_dir / "rerooted.newick"
    input_tree = tree if tree is not None else files.tree
    command = [
        "treetime",
        "--tree",
        str(input_tree if reroot else rooted_tree),
        "--dates",
        str(files.metadata),
        "--name-column",
        "sample_id",
        "--date-column",
        "collection_date",
        "--sequence-length",
        str(sequence_length),
        "--confidence",
        "--time-marginal",
        "only-final",
        "--covariation",
        "--clock-filter",
        "0",
        "--plot-tree",
        "timetree.svg",
        "--plot-rtt",
        "root_to_tip_regression.svg",
        "--outdir",
        str(files.directory / "timetree"),
    ]
    if reroot:
        command.extend(["--reroot", "least-squares"])
    else:
        command.append("--keep-root")
    _run_command(
        command,
        log=files.directory / "logs" / "timetree.log",
        expected=files.time_tree,
        force=force,
        inputs=((input_tree if reroot else rooted_tree), files.metadata),
    )


def _run_location_tree(files: LineageFiles, tree: Path, force: bool) -> None:
    _run_command(
        [
            "treetime",
            "mugration",
            "--tree",
            str(tree),
            "--states",
            str(files.states),
            "--name-column",
            "sample_id",
            "--attribute",
            "location",
            "--confidence",
            "--outdir",
            str(files.directory / "location"),
        ],
        log=files.directory / "logs" / "location.log",
        expected=files.location_tree,
        force=force,
        inputs=(tree, files.states),
    )


def _report_record(
    item: dict[str, object],
    files: LineageFiles,
    reference: Sample,
    sequence_length: int,
    threads: int,
    randomisation_jobs: int,
    temporal: dict[str, object],
    assessment: dict[str, object],
    context_summary: dict[str, object],
    public_health: dict[str, object],
    recombination_masking: dict[str, object],
    time_tree_available: bool,
) -> dict[str, object]:
    def optional(path: Path) -> str | None:
        return str(path) if time_tree_available else None

    return {
        **item,
        "analysis_mode": "full",
        "date_randomisation_method": temporal.get("method", "root_to_tip"),
        "reference_path": str(reference.assembly),
        "alignment_length": alignment_length(files.alignment),
        "complete_alignment_sites": sequence_length,
        "resources": {"threads": threads, "date_randomisation_jobs": randomisation_jobs},
        "temporal_status": assessment["code"],
        "temporal_signal_reason": assessment["reason"],
        "temporal_signal_supported": bool(assessment["supported"]),
        "temporal_signal": temporal,
        "context": context_summary,
        "public_health_status": public_health["scenario"]["code"],
        "public_health_label": public_health["scenario"]["label"],
        "public_health_confidence": public_health["scenario"]["confidence"],
        "public_health": public_health,
        "recombination_masking": recombination_masking,
        "outputs": {
            "alignment": str(files.alignment),
            "lineage_coherence": str(files.coherence_screen),
            "clonal_lineage_coherence": str(files.clonal_coherence_screen),
            "starting_tree": str(files.starting_tree),
            "tree": str(files.tree),
            "recombination_importations": str(files.importations),
            "recombination_intervals": str(files.directory / "recombination_intervals.tsv"),
            "recombination_summary": str(files.directory / "recombination_summary.json"),
            "recombination_profile": str(files.directory / "recombination_genome_profile.csv"),
            "recombination_map": str(files.directory / "recombination_map.svg"),
            "recombination_map_png": str(files.directory / "recombination_map.png"),
            "filtered_alignment": str(files.filtered_alignment),
            "root_to_tip_plot": str(files.clock_dir / "root_to_tip_regression.svg"),
            "root_to_tip_png": str(files.directory / "root_to_tip.png"),
            "root_to_tip_data": str(files.clock_dir / "rtt.csv"),
            "date_randomisation_plot": str(files.directory / "date_randomisation.svg"),
            "date_randomisation_png": str(files.directory / "date_randomisation.png"),
            "date_randomisation_csv": str(files.directory / "date_randomisation.csv"),
            "timetree": optional(files.time_tree),
            "timetree_plot": optional(files.directory / "timetree" / "timetree.svg"),
            "timetree_png": optional(files.directory / "timetree.png"),
            "timetree_confidence_plot": optional(files.directory / "timetree_with_confidence.svg"),
            "location_tree": str(files.location_tree),
            "public_health_evidence": str(files.directory / "public_health_evidence.json"),
            "clonal_pairwise_distances": str(files.directory / "clonal_pairwise_distances.tsv"),
            "clonal_snp_matrix": str(files.directory / "clonal_snp_matrix.tsv"),
            "pairwise_callable_sites": str(files.directory / "pairwise_callable_sites.tsv"),
            "clonal_snp_heatmap": str(files.directory / "clonal_snp_heatmap.svg"),
            "clonal_snp_heatmap_png": str(files.directory / "clonal_snp_heatmap.png"),
            "timetree_confidence": optional(files.directory / "timetree_confidence.csv"),
            "node_dates": optional(files.directory / "node_dates.csv"),
            "supporting_results": str(files.directory / "supporting_results.zip"),
            "html_report": str(files.directory / "report.html"),
        },
    }


def _run_lineage(
    item: dict[str, object],
    members: list[Sample],
    *,
    output: Path,
    threads: int,
    randomisation_jobs: int,
    randomisations: int,
    temporal_p_value: float,
    seed: int,
    force: bool,
    context_manifest_rows: list[dict[str, str]],
    mode: str,
    date_randomisation_method: str,
) -> dict[str, object]:
    """Run one ready lineage; independent lineages may call this concurrently."""

    files = LineageFiles.in_directory(output / str(item["slug"]))
    files.directory.mkdir(parents=True, exist_ok=True)
    reference = select_reference(members)
    _write_lineage_inputs(files.directory, members)
    context_summary = context_evidence(members, context_manifest_rows, directory=files.directory)
    tree, filtered_alignment, recombination = _run_core_phylogeny(
        files, members, reference, threads, force, mode
    )
    sequence_length = complete_alignment_sites(filtered_alignment)
    _run_observed_clock(files, tree, sequence_length, force)
    observed_clock = files.clock_file
    if date_randomisation_method == "full_tree":
        _run_dated_tree(files, sequence_length, force, tree, reroot=True)
        observed_clock = files.directory / "timetree" / "molecular_clock.txt"
    temporal = _temporal_signal(
        files,
        members,
        sequence_length,
        randomisations,
        randomisation_jobs,
        seed,
        force,
        tree,
        date_randomisation_method,
        observed_clock,
    )
    assessment = assess_temporal_signal(temporal, p_value_threshold=temporal_p_value)
    if mode == "fast":
        report = {
            **item,
            "analysis_mode": "fast",
            "date_randomisation_method": "root_to_tip",
            "reference_path": str(reference.assembly),
            "alignment_length": alignment_length(files.alignment),
            "complete_alignment_sites": sequence_length,
            "temporal_status": assessment["code"],
            "temporal_signal_reason": assessment["reason"],
            "temporal_signal_supported": bool(assessment["supported"]),
            "temporal_signal": temporal,
            "recombination": recombination,
            "outputs": {
                "alignment": str(files.alignment),
                "lineage_coherence": str(files.coherence_screen),
                "phipack_profile": str(files.phipack_profile),
                "phipack_significant_blocks": str(files.phipack_significant_blocks),
                "screening_tree": str(files.fast_tree),
                "root_to_tip_plot": str(files.clock_dir / "root_to_tip_regression.svg"),
                "date_randomisation_csv": str(files.directory / "date_randomisation.csv"),
                "html_report": str(files.directory / "report.html"),
            },
        }
        (files.directory / "report.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8"
        )
        write_fast_lineage_report(
            report, directory=files.directory, p_value_threshold=temporal_p_value
        )
        write_supporting_bundle(files.directory)
        return report
    rooted_tree = files.clock_dir / "rerooted.newick"
    recombination_masking = write_recombination_evidence(
        alignment=files.alignment,
        filtered_alignment=files.filtered_alignment,
        importations=files.importations,
        output_directory=files.directory,
        reference=reference.assembly,
    )
    public_health = build_public_health_evidence(
        filtered_alignment=filtered_alignment,
        rooted_tree=rooted_tree if rooted_tree.is_file() else files.tree,
        samples=members,
        output=files.directory,
        temporal_assessment=assessment,
    )
    temporal_supported = bool(assessment["supported"])
    if temporal_supported:
        if date_randomisation_method != "full_tree":
            _run_dated_tree(files, sequence_length, force, tree)
    time_tree_available = temporal_supported and files.time_tree.exists()
    _run_location_tree(files, files.time_tree if time_tree_available else files.tree, force)
    report = _report_record(
        item,
        files,
        reference,
        sequence_length,
        threads,
        randomisation_jobs,
        temporal,
        assessment,
        context_summary,
        public_health,
        recombination_masking,
        time_tree_available,
    )
    (files.directory / "report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    write_lineage_report(report, directory=files.directory, p_value_threshold=temporal_p_value)
    write_supporting_bundle(files.directory)
    return report
