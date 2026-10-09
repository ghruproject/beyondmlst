"""Execution of the ChronoClade workflow."""

from __future__ import annotations

import json
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path

from chronoclade.errors import WorkflowError
from chronoclade.evidence import EvidenceError
from chronoclade.lineage import (
    _run_lineage,
    complete_alignment_sites as complete_alignment_sites,
    context_evidence as context_evidence,
    read_context_manifest,
)
from chronoclade.metadata import (
    Sample,
    group_samples,
    select_reference,
    slugify_lineage,
)
from chronoclade.report import write_summary_report

TOOLS_BY_MODE = {
    "full": ("ska", "iqtree", "ClonalFrameML", "treetime"),
    "corrected": ("ska", "iqtree", "ClonalFrameML", "treetime"),
    "fast": ("ska", "iqtree", "Profile", "treetime"),
}
REQUIRED_TOOLS = tuple(dict.fromkeys(tool for tools in TOOLS_BY_MODE.values() for tool in tools))


@dataclass(frozen=True)
class ResourcePlan:
    """Bounded concurrency derived from the user's global CPU budget."""

    total_threads: int
    lineage_workers: int
    threads_per_lineage: int
    randomisation_workers_per_lineage: int


def allocate_resources(
    *,
    total_threads: int,
    lineage_jobs: int,
    randomisation_jobs: int,
    ready_lineages: int,
) -> ResourcePlan:
    """Allocate lineage and permutation workers without oversubscribing CPUs."""

    if min(total_threads, lineage_jobs, randomisation_jobs) < 1:
        raise ValueError("thread and job counts must be at least 1")
    if ready_lineages < 0:
        raise ValueError("ready_lineages cannot be negative")

    lineage_workers = min(lineage_jobs, ready_lineages, total_threads)
    if lineage_workers == 0:
        return ResourcePlan(total_threads, 0, total_threads, min(randomisation_jobs, total_threads))

    threads_per_lineage = max(1, total_threads // lineage_workers)
    return ResourcePlan(
        total_threads=total_threads,
        lineage_workers=lineage_workers,
        threads_per_lineage=threads_per_lineage,
        randomisation_workers_per_lineage=min(randomisation_jobs, threads_per_lineage),
    )


def tool_status() -> dict[str, str | None]:
    return {tool: shutil.which(tool) for tool in REQUIRED_TOOLS}


def native_platform_supported() -> bool:
    """Return whether the locked Pixi tools support this operating system."""

    return sys.platform.startswith(("linux", "darwin"))


def check_tools(mode: str = "full") -> None:
    if not native_platform_supported():
        raise WorkflowError(
            "Native execution is supported on macOS and Linux. This platform is not in the "
            "locked Pixi environment."
        )
    missing = [tool for tool in TOOLS_BY_MODE[mode] if shutil.which(tool) is None]
    if missing:
        raise WorkflowError(
            "Missing workflow tools: "
            + ", ".join(missing)
            + ". Run this command through `pixi run`."
        )


def plan(samples: list[Sample], *, min_samples: int) -> list[dict[str, object]]:
    """Return the deterministic lineage plan without running external tools."""

    result: list[dict[str, object]] = []
    for (species, lineage), members in sorted(group_samples(samples).items()):
        reference = select_reference(members)
        result.append(
            {
                "species": species,
                "lineage": lineage,
                "slug": slugify_lineage(species, lineage),
                "sample_count": len(members),
                "distinct_dates": len({sample.collection_date for sample in members}),
                "first_collection_date": min(sample.collection_date for sample in members),
                "last_collection_date": max(sample.collection_date for sample in members),
                "locations": sorted({sample.location for sample in members}),
                "reference": reference.sample_id,
                "status": "ready" if len(members) >= min_samples else "too_few_samples",
            }
        )
    return result


def run_workflow(
    samples: list[Sample],
    *,
    output: Path,
    threads: int,
    lineage_jobs: int,
    randomisation_jobs: int,
    randomisations: int,
    temporal_p_value: float,
    min_samples: int,
    seed: int,
    force: bool,
    context_manifest: Path | None = None,
    mode: str = "full",
    date_randomisation_method: str = "root_to_tip",
) -> dict[str, object]:
    """Run every lineage through alignment, recombination and dating stages."""

    if mode not in TOOLS_BY_MODE:
        raise WorkflowError("mode must be 'full', 'corrected' or 'fast'")
    if date_randomisation_method not in {"root_to_tip", "full_tree"}:
        raise WorkflowError("date randomisation method must be 'root_to_tip' or 'full_tree'")
    if mode == "fast" and date_randomisation_method != "root_to_tip":
        raise WorkflowError("fast mode only supports root-to-tip date randomisation")
    check_tools() if mode == "full" else check_tools(mode)
    output = output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    lineage_plan = plan(samples, min_samples=min_samples)
    groups = group_samples(samples)
    manifest_rows = read_context_manifest(context_manifest)
    ready_items = [item for item in lineage_plan if item["status"] == "ready"]
    resources = allocate_resources(
        total_threads=threads,
        lineage_jobs=lineage_jobs,
        randomisation_jobs=randomisation_jobs,
        ready_lineages=len(ready_items),
    )

    def run_ready_item(indexed_item: tuple[int, dict[str, object]]) -> dict[str, object]:
        index, item = indexed_item
        key = (str(item["species"]), str(item["lineage"]))
        try:
            return _run_lineage(
                item,
                groups[key],
                output=output,
                threads=resources.threads_per_lineage,
                randomisation_jobs=resources.randomisation_workers_per_lineage,
                randomisations=randomisations,
                temporal_p_value=temporal_p_value,
                seed=seed + index,
                force=force,
                context_manifest_rows=manifest_rows,
                mode=mode,
                date_randomisation_method=date_randomisation_method,
            )
        except (EvidenceError, OSError, WorkflowError) as error:
            raise WorkflowError(f"Lineage {key[0]} / {key[1]} failed: {error}") from error

    completed_by_slug: dict[str, dict[str, object]] = {}
    if ready_items:
        indexed_items = list(enumerate(ready_items))
        with ThreadPoolExecutor(max_workers=resources.lineage_workers) as executor:
            for report in executor.map(run_ready_item, indexed_items):
                completed_by_slug[str(report["slug"])] = report

    reports = [
        completed_by_slug[str(item["slug"])] if item["status"] == "ready" else item
        for item in lineage_plan
    ]

    summary = {
        "workflow": "chronoclade",
        "analysis_mode": mode,
        "date_randomisation_method": date_randomisation_method,
        "resources": asdict(resources),
        "lineages": reports,
        "context_manifest": str(context_manifest.expanduser().resolve())
        if context_manifest
        else None,
        "guardrail": (
            "Location-state reconstructions are exploratory and do not by themselves establish "
            "direct transmission or a definitive number of introductions."
        ),
    }
    summary["report"] = str(write_summary_report(summary, output=output))
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary
