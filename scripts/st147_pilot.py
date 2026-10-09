#!/usr/bin/env python3
"""Bounded ST147 pilot orchestration, with an entirely offline replay mode."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import subprocess
import sys
import time
from dataclasses import fields
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(path: Path, delimiter="\t") -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def replay(context: Path, output: Path) -> dict:
    """Replay stored selection distances and annotated catalogue without any network/tool calls."""
    from chronoclade.context import ContextCandidate, parse_ska_distances, select_context
    from chronoclade.context_geography import generate_context_geography
    from chronoclade.pathogenwatch import content_hash

    audit = json.loads((context / "context_selection.json").read_text())
    catalogue = json.loads((context / "context_catalogue.json").read_text())
    payload = {key: value for key, value in catalogue.items() if key != "snapshot_sha256"}
    if content_hash(payload) != catalogue["snapshot_sha256"]:
        raise ValueError("Annotated catalogue changed; re-run context preparation")
    if catalogue["snapshot_sha256"] != audit["catalogue_sha256"]:
        raise ValueError("Selection uses another catalogue snapshot")
    ledger = audit["download_ledger"]
    successful = {row["source_genome_id"] for row in ledger if row["status"] != "failed"}
    names = {f.name for f in fields(ContextCandidate)}
    candidates = [
        ContextCandidate(**{k: v for k, v in row.items() if k in names})
        for row in table(context / "candidate_pool.tsv")
        if row["source_genome_id"] in successful
    ]
    focal_ids = [
        row["sample_id"]
        for row in table(context / "combined_metadata.csv", ",")
        if row.get("origin") != "context"
    ]
    settings = audit["filters"]
    selected, _ = select_context(
        candidates,
        focal_ids,
        parse_ska_distances(context / "context_distances.tsv"),
        max_context=settings["max_context"],
        nearest_per_focal=settings["nearest_per_focal"],
        seed=audit["seed"],
    )
    actual = [(r.source_genome_id, r.selection_reason) for r in selected]
    expected = [
        (r["source_genome_id"], r["selection_reason"])
        for r in table(context / "context_manifest.tsv")
    ]
    if actual != expected:
        raise ValueError("Offline selection IDs/reasons differ from saved manifest")
    saved_geography = json.loads(
        (context / "context_geography/country_composition_audit.json").read_text()
    )
    regenerated = generate_context_geography(
        catalogue["rows"],
        output / "context_geography",
        focal_rows=catalogue["focal_rows"],
        selected_source_ids=[r[0] for r in actual],
        scope=saved_geography["scope"],
    )
    if regenerated["rows"] != saved_geography["rows"]:
        raise ValueError("Offline country/group counts differ from saved source table")
    result = {
        "offline": True,
        "selection_identical": True,
        "country_counts_identical": True,
        "selected_ids_and_reasons": actual,
        "catalogue_sha256": catalogue["snapshot_sha256"],
        "rows": len(regenerated["rows"]),
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "replay_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage", choices=("prepare", "analyse", "replay", "figures", "live-smoke"), required=True
    )
    parser.add_argument("--output", type=Path, default=Path("validation/st147_pathogenwatch/run"))
    parser.add_argument("--context", type=Path)
    parser.add_argument("--frozen-fixture", type=Path)
    parser.add_argument("--focal", type=Path)
    parser.add_argument("--catalogue", type=Path)
    parser.add_argument("--cglin-export", type=Path)
    parser.add_argument("--focal-crosswalk", type=Path)
    parser.add_argument("--candidate-pool", type=int, default=8)
    parser.add_argument("--max-context", type=int, default=4)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20261009)
    parser.add_argument("--date-randomisations", type=int, default=100)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--public-focal-demonstration",
        action="store_true",
        help="Label analysis reports as demonstrations using public focal stand-ins",
    )
    parser.add_argument(
        "--confirm-live", action="store_true", help="Explicit opt-in to contact upstream"
    )
    args = parser.parse_args(argv)
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    context = (args.context or args.output / "context").resolve()
    if args.stage == "figures":
        if not args.frozen_fixture:
            parser.error("figures requires --frozen-fixture")
        import gzip
        from chronoclade.pathogenwatch import content_hash
        from chronoclade.context_geography import generate_context_geography

        fixture = json.loads(gzip.decompress(args.frozen_fixture.read_bytes()))
        payload = {k: v for k, v in fixture.items() if k != "fixture_sha256"}
        if content_hash(payload) != fixture["fixture_sha256"]:
            raise ValueError("Compact frozen fixture hash mismatch")
        generate_context_geography(
            fixture["rows"],
            args.output / "context_geography",
            focal_rows=fixture["focal_rows"],
            selected_source_ids=fixture.get("selected_source_ids", ()),
            scope=fixture["scope"],
        )
        return 0
    if args.stage == "replay":
        print(json.dumps(replay(context, args.output / "offline_replay"), indent=2))
        return 0
    if args.stage == "live-smoke":
        if not args.confirm_live:
            parser.error("live-smoke requires --confirm-live")
        from chronoclade.pathogenwatch import PathogenwatchClient

        client = PathogenwatchClient()
        search = client.search("573", "147")
        # Public metadata only; small capability-contract probe, never a complete freeze claim.
        first = search["genomes"][0] if search["genomes"] else None
        result = {
            "utc": datetime.now(timezone.utc).isoformat(),
            "public_only": True,
            "search_count": search["expected_count"],
            "catalogue_frozen": False,
            "first_uuid": first["uuid"] if first else None,
        }
        if first:
            detail = client.request_json(
                "GET", "/api/genomes/details", params={"id": first["uuid"]}
            )
            if detail.get("uuid") != first["uuid"]:
                raise ValueError("Live details UUID contract changed")
            result["details_uuid_verified"] = True
        (args.output / "live_smoke.json").write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2))
        return 0
    prefix = [sys.executable, "-m", "chronoclade.cli"]
    inputs = {}
    if args.stage == "prepare":
        if not args.focal or not args.catalogue or not args.cglin_export:
            parser.error("prepare requires --focal, --catalogue and --cglin-export")
        command = prefix + [
            "prepare-context",
            str(args.focal.resolve()),
            "--scheme",
            "klebsiella",
            "--context-source",
            "pathogenwatch",
            "--catalogue",
            str(args.catalogue.resolve()),
            "--cglin-export",
            str(args.cglin_export.resolve()),
            "--output",
            str(context),
            "--candidate-pool",
            str(args.candidate_pool),
            "--max-context",
            str(args.max_context),
            "--nearest-per-focal",
            "1",
            "--seed",
            str(args.seed),
            "--threads",
            str(args.threads),
            "--cache-dir",
            str(args.output / "assembly_cache"),
        ]
        if args.focal_crosswalk:
            command += ["--focal-crosswalk", str(args.focal_crosswalk.resolve())]
        for field in ("focal", "catalogue", "cglin_export", "focal_crosswalk"):
            path = getattr(args, field)
            if path:
                inputs[field] = {"path": str(path.resolve()), "sha256": sha256(path)}
    else:
        command = prefix + [
            "run",
            str(context / "combined_metadata.csv"),
            "--context-manifest",
            str(context / "context_manifest.tsv"),
            "--output",
            str(args.output / "analysis"),
            "--threads",
            str(args.threads),
            "--lineage-jobs",
            "1",
            "--randomisation-jobs",
            "1",
            "--min-samples",
            "4",
            "--date-randomisations",
            str(args.date_randomisations),
            "--seed",
            str(args.seed),
        ]
        for name in ("combined_metadata.csv", "context_manifest.tsv"):
            path = context / name
            inputs[name] = {"path": str(path), "sha256": sha256(path)}
    if args.dry_run:
        command.append("--dry-run")
    started = time.monotonic()
    result = subprocess.run(command, check=False)
    if (
        args.stage == "analyse"
        and result.returncode == 0
        and not args.dry_run
        and args.public_focal_demonstration
    ):
        label_public_demonstration(args.output / "analysis")
    provenance = {
        "stage": args.stage,
        "utc": datetime.now(timezone.utc).isoformat(),
        "command": command,
        "inputs": inputs,
        "returncode": result.returncode,
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "threads": args.threads,
        "seed": args.seed,
        "analysis_gates": "ChronoClade defaults; temporal p=0.05; divergence checks retained",
    }
    (args.output / f"{args.stage}_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    return result.returncode


def label_public_demonstration(analysis: Path) -> None:
    """Mark only explicitly requested public stand-in runs, without altering their results."""
    from chronoclade.report import (
        write_fast_lineage_report,
        write_lineage_report,
        write_supporting_bundle,
    )

    for path in analysis.glob("*/report.json"):
        report = json.loads(path.read_text())
        report["demonstration"] = {
            "label": "Public-data demonstration",
            "description": (
                "The focal samples in this pilot are public genomes chosen to demonstrate "
                "the workflow. They come from different countries and are not a defined "
                "local outbreak or patient cohort. Automated introduction/persistence "
                "labels are illustrative and must not be treated as an epidemiological finding."
            ),
        }
        path.write_text(json.dumps(report, indent=2) + "\n")
        writer = (
            write_fast_lineage_report
            if report.get("analysis_mode") == "fast"
            else write_lineage_report
        )
        writer(report, directory=path.parent, p_value_threshold=0.05)
        write_supporting_bundle(path.parent)


if __name__ == "__main__":
    raise SystemExit(main())
