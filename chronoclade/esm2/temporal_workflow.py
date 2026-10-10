"""Combine two descriptive temporal views without fitting an embedding clock."""

from dataclasses import dataclass
from pathlib import Path

from chronoclade import __version__
from chronoclade.datasets import load_dataset
from chronoclade.temporal_diagnostics import cgmlst_root_to_tip

from .proteins import file_sha256
from .sample_distances import analyze_embedding_dates
from .storage import write_json
from .temporal_report import write_temporal_report


@dataclass(frozen=True)
class TemporalReportResult:
    report_path: Path
    manifest_path: Path


def run_temporal_report(
    dataset_manifest: Path,
    embeddings_manifest: Path,
    mapping_csv: Path,
    output: Path,
    *,
    reference_sample_id: str | None = None,
    panel_loci: list[str] | None = None,
) -> TemporalReportResult:
    """Read frozen evidence and publish an ESM2 report with its cgMLST baseline.

    This operation does not import a model, retrieve sequences, or infer ancestry
    from embeddings. A saved vector manifest can be analysed on a CPU-only host.
    """
    dataset = load_dataset(dataset_manifest)
    embedding = analyze_embedding_dates(
        dataset,
        embeddings_manifest,
        mapping_csv,
        reference_sample_id=reference_sample_id,
        panel_loci=panel_loci,
    )
    reference = embedding["reference"]["sample_id"]
    guide = cgmlst_root_to_tip(dataset, reference_sample_id=reference)
    report = write_temporal_report(embedding, guide, output)
    manifest = {
        "schema": "chronoclade.esm2.temporal-report",
        "schema_version": 1,
        "software_version": __version__,
        "status": "complete",
        "experimental": True,
        "scope": "descriptive-temporal-comparison",
        "dataset_id": dataset.dataset_id,
        "inputs": {
            name: {"path": str(Path(path).resolve()), "sha256": file_sha256(path)}
            for name, path in (
                ("dataset", dataset_manifest),
                ("embeddings", embeddings_manifest),
                ("sample_loci", mapping_csv),
            )
        },
        "parameters": {
            "reference_sample_id": reference,
            "requested_panel_loci": panel_loci,
            "reference_selection": embedding["reference"],
        },
        "report": {"path": report.relative_to(output).as_posix(), "sha256": file_sha256(report)},
    }
    manifest_path = output / "temporal_analysis.json"
    write_json(manifest_path, manifest)
    return TemporalReportResult(report, manifest_path)
