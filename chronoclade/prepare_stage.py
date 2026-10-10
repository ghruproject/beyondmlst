"""Independent, offline preparation of frozen profiles and validated bundles.

This stage publishes typed evidence and audits readiness. It does not call the
legacy resolver, retrieve public context, or imply that missing typing was run.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from html import escape
import json
from pathlib import Path
import shutil
import tempfile

from chronoclade import __version__
from chronoclade.artifacts import file_sha256, write_json
from chronoclade.datasets import (
    DatasetError,
    LocusCatalogue,
    SCHEMA_NAME,
    from_profile_records,
    load_dataset,
    write_dataset,
)
from chronoclade.datasets.model import EVIDENCE_TABLES, check_json
from chronoclade.report_components.styles import report_styles
from chronoclade.typing_scopes import typing_scope


@dataclass(frozen=True)
class PrepareResult:
    dataset_manifest: Path
    report_path: Path
    audit_path: Path
    stage_manifest: Path


def _json(path: Path):
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        check_json(value, location=path.name)
        return value
    except (OSError, ValueError) as error:
        raise DatasetError(f"Cannot read preparation JSON {path}: {error}") from error


def _catalogues(path: Path) -> tuple[LocusCatalogue, ...]:
    raw = _json(path)
    if isinstance(raw, dict):
        if set(raw) != {"catalogues"}:
            raise DatasetError("Catalogue JSON object requires only a catalogues list")
        raw = raw["catalogues"]
    if not isinstance(raw, list) or not raw:
        raise DatasetError("Catalogue JSON requires a nonempty catalogue list")
    result = []
    for index, value in enumerate(raw):
        if not isinstance(value, dict) or not isinstance(value.get("loci"), list):
            raise DatasetError(f"catalogues[{index}] requires an explicit loci list")
        if "complete" not in value:
            raise DatasetError(f"catalogues[{index}] requires explicit complete true/false")
        try:
            catalogue = LocusCatalogue(**dict(value, loci=tuple(value["loci"])))
            catalogue.validate()
        except (TypeError, ValueError) as error:
            raise DatasetError(f"catalogues[{index}]: {error}") from error
        result.append(catalogue)
    return tuple(result)


def _frozen_dataset(raw, catalogues, species):
    evidence = {}
    if isinstance(raw, dict):
        allowed = {"records", "parameters", *EVIDENCE_TABLES}
        if set(raw) - allowed or "records" not in raw:
            raise DatasetError(
                "Frozen input requires records and optional evidence tables/parameters; "
                "legacy queries/context resolver outputs need explicit adaptation"
            )
        evidence = {key: value for key, value in raw.items() if key != "records"}
        raw = raw["records"]
    if not isinstance(raw, list) or not raw or any(not isinstance(row, dict) for row in raw):
        raise DatasetError("Frozen input requires a nonempty list of canonical record objects")
    for name in EVIDENCE_TABLES:
        if name in evidence and not isinstance(evidence[name], list):
            raise DatasetError(f"Frozen {name} must be an explicit list")
    if "parameters" in evidence and not isinstance(evidence["parameters"], dict):
        raise DatasetError("Frozen parameters must be an object")
    records = []
    for index, row in enumerate(raw):
        if not isinstance(row.get("sample_id"), str):
            raise DatasetError(f"records[{index}] requires a stable sample_id")
        role = row.get("role")
        origin = row.get("origin")
        if role is None and origin is None:
            role = "input"
        record = dict(row, role=role)
        if species:
            existing = row.get("species")
            if existing and (
                not isinstance(existing, str) or existing.casefold() != species.casefold()
            ):
                raise DatasetError(f"Sample {row['sample_id']!r} species differs from --species")
            record["species"] = species
        records.append(record)
    try:
        return from_profile_records(records, catalogues, **evidence)
    except (KeyError, TypeError, AttributeError) as error:
        raise DatasetError(f"Malformed frozen profile evidence: {error}") from error


def _portable_assemblies(dataset, source_root, staging):
    """Copy explicitly linked local assemblies; leave external identities unresolved."""
    rows, artifacts = [], {}
    for row in dataset.samples:
        row = dict(row)
        reference = row.get("assembly_reference")
        if reference:
            if "://" in reference:
                rows.append(row)
                continue
            path = Path(reference).expanduser()
            path = path if path.is_absolute() else source_root / path
            if path.is_file():
                digest = file_sha256(path)
                relative = f"assemblies/{digest}{path.suffix}"
                destination = staging / relative
                destination.parent.mkdir(exist_ok=True)
                if not destination.exists():
                    shutil.copyfile(path, destination)
                row["assembly_reference"] = relative
                artifacts[relative] = {
                    "path": relative,
                    "sha256": digest,
                    "bytes": destination.stat().st_size,
                }
            elif Path(reference).is_absolute() or "/" in reference or Path(reference).suffix:
                raise DatasetError(
                    f"Sample {row['sample_id']!r} assembly file is missing: {reference}"
                )
        rows.append(row)
    return replace(dataset, samples=tuple(rows)), list(artifacts.values())


def _validate_previous_stage(root):
    """Validate extra assembly/source/report artifacts when importing our output."""
    path = root / "prepare.json"
    if not path.exists():
        return
    stage = _json(path)
    if (
        not isinstance(stage, dict)
        or stage.get("schema") != "chronoclade.prepare-result"
        or type(stage.get("schema_version")) is not int
        or stage["schema_version"] != 1
        or stage.get("status") != "complete"
        or not isinstance(stage.get("artifacts"), list)
    ):
        raise DatasetError("Existing prepare result requires a complete supported manifest")
    seen = set()
    for item in stage["artifacts"]:
        if not isinstance(item, dict):
            raise DatasetError("Existing prepare artifact requires a descriptor")
        raw = item.get("path")
        if (
            not isinstance(raw, str)
            or not raw
            or "\\" in raw
            or Path(raw).is_absolute()
            or any(part in {".", "..", ""} for part in raw.split("/"))
        ):
            raise DatasetError("Existing prepare artifact requires a safe relative path")
        artifact = (root / raw).resolve()
        if not artifact.is_relative_to(root) or raw in seen or not artifact.is_file():
            raise DatasetError(
                "Existing prepare artifact is missing, duplicated or outside the bundle"
            )
        seen.add(raw)
        if (
            type(item.get("bytes")) is not int
            or artifact.stat().st_size != item["bytes"]
            or file_sha256(artifact) != item.get("sha256")
        ):
            raise DatasetError(f"Existing prepare artifact checksum/size mismatch: {raw}")
    if not {"dataset.json", "readiness.json", "report.html"}.issubset(seen):
        raise DatasetError("Existing prepare result is missing its required artifacts")


def _expected_lineage(species):
    species = (species or "").casefold()
    if "klebsiella" in species:
        return "cglin"
    if species in {"escherichia coli", "e. coli", "ecoli"}:
        return "hiercc"
    return None


def _lineage_ready(evidence, kind):
    """Retain the provider's complete cgLIN status without accepting provisional codes."""
    status = evidence.get("resolution_status")
    if status not in ({"resolved", "complete"} if kind == "cglin" else {"resolved"}):
        return False
    frozen = evidence.get("evidence") or {}
    if not isinstance(frozen, dict):
        raise DatasetError(f"Sample {evidence['sample_id']!r} lineage evidence must be an object")
    provisional = frozen.get(kind + "_provisional")
    if provisional is True or str(provisional).casefold() in {"true", "1", "yes"}:
        return False
    scope = typing_scope(
        dict(
            frozen,
            **{
                kind + "_scheme": evidence.get("scheme_id"),
                kind + "_scheme_version": evidence.get("scheme_version"),
            },
        ),
        kind,
    )
    return bool(evidence.get("kind") == kind and evidence.get("assignment") and scope)


def _readiness(dataset, manifest):
    profiles, lineages = {}, {}
    for matrix in dataset.profiles:
        for ident, codes in zip(matrix.sample_ids, matrix.codes, strict=True):
            profiles.setdefault(ident, []).append(
                {
                    "scheme_id": matrix.catalogue.scheme_id,
                    "scheme_version": matrix.catalogue.scheme_version,
                    "database_version": matrix.catalogue.database_version,
                    "database_sha256": matrix.catalogue.database_sha256,
                    "catalogue_complete": matrix.catalogue.complete,
                    "loci": len(matrix.catalogue.loci),
                    "called_loci": int((codes != 0).sum()),
                    "missing_loci": int((codes == 0).sum()),
                }
            )
    for row in dataset.lineages:
        lineages.setdefault(row["sample_id"], []).append(row)
    samples = []
    for row in dataset.samples:
        ident = row["sample_id"]
        profile_rows = profiles.get(ident, [])
        kind = _expected_lineage(row.get("species"))
        lineage_rows = lineages.get(ident, [])
        lineage_available = (
            any(_lineage_ready(evidence, kind) for evidence in lineage_rows) if kind else False
        )
        reasons = []
        if not any(evidence["called_loci"] for evidence in profile_rows):
            reasons.append("cgMLST profile unavailable; no profile calling was performed")
        if any(not evidence["catalogue_complete"] for evidence in profile_rows):
            reasons.append("Locus catalogue is explicitly incomplete")
        if kind and not lineage_available:
            reasons.append(f"Resolved {kind} assignment with a known scope is unavailable")
        if not kind:
            reasons.append("Lineage readiness is not assessed for this species")
        reference = row.get("assembly_reference")
        samples.append(
            {
                "sample_id": ident,
                "label": row["label"],
                "role": row["role"],
                "profiles": profile_rows,
                "required_lineage_kind": kind,
                "lineage_ready": bool(lineage_available),
                "profile_ready": bool(any(evidence["called_loci"] for evidence in profile_rows)),
                "readiness": "ready" if not reasons else "incomplete",
                "reasons": reasons,
                "assembly_status": "bundled"
                if reference and (manifest.parent / reference).is_file()
                else "external_reference"
                if reference
                else "unavailable",
                "missing_metadata": [
                    field
                    for field in ("country", "region", "nuts2", "host", "isolation_source")
                    if not row.get(field)
                ]
                + (
                    ["collection_date"]
                    if row.get("date_precision") in {None, "missing", "invalid"}
                    else []
                ),
            }
        )
    return {
        "schema": "chronoclade.prepare-readiness",
        "schema_version": 1,
        "status": "complete",
        "dataset_manifest": "dataset.json",
        "dataset_sha256": file_sha256(manifest),
        "dataset_id": json.loads(manifest.read_text())["dataset_id"],
        "sample_ids": list(dataset.sample_ids),
        "input_ids": [row["sample_id"] for row in dataset.samples if row["role"] == "input"],
        "context_ids": [row["sample_id"] for row in dataset.samples if row["role"] == "context"],
        "readiness": "ready"
        if all(row["readiness"] == "ready" for row in samples)
        else "incomplete",
        "samples": samples,
        "operations": {
            "typing": "not_run_frozen_import",
            "lineage_assignment": "not_run_frozen_import",
            "provider_retrieval": "not_run_offline",
            "context_discovery": "not_run",
            "typing_tool_installation": "not_checked_not_required_for_import",
            "reference_database": "not_checked_not_required_for_import",
        },
    }


def _report(audit):
    def text(value):
        return escape(str(value), quote=True)

    ready = sum(row["readiness"] == "ready" for row in audit["samples"])
    rows = []
    for row in audit["samples"]:
        calls = (
            "; ".join(
                f"{profile['called_loci']} / {profile['loci']} called ({profile['scheme_id']})"
                for profile in row["profiles"]
            )
            or "Unavailable"
        )
        reasons = "; ".join(row["reasons"]) or "Frozen profile and lineage evidence available"
        rows.append(
            "<tr>"
            + "".join(
                f"<td>{text(value)}</td>"
                for value in (
                    row["label"],
                    row["role"],
                    calls,
                    row["readiness"],
                    reasons,
                    ", ".join(row["missing_metadata"]) or "None",
                    row["assembly_status"],
                )
            )
            + "</tr>"
        )
    return (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        "<title>ChronoClade preparation audit</title><style>"
        + report_styles()
        + '</style></head><body><main class="shell"><header class="identity">'
        '<div class="identity-mark">ChronoClade · PREPARATION AUDIT</div>'
        '<div class="identity-copy"><h1>Prepared dataset</h1>'
        "<p>Frozen evidence is saved in a validated portable bundle. Missing profiles, "
        "lineage assignments and metadata remain visible.</p></div></header>"
        '<section class="stage"><div class="stage-body"><h2>Typing readiness</h2>'
        "<p>Existing typing was imported. Profile calling, lineage assignment, provider "
        "retrieval and reference database availability were not assessed or run.</p>"
        f'<div class="measure-strip"><div><span>Samples</span><b>{len(rows)}</b></div>'
        f"<div><span>Ready frozen typing</span><b>{ready}</b></div>"
        f"<div><span>Incomplete typing</span><b>{len(rows) - ready}</b></div></div>"
        "<p>Missing metadata does not exclude a sample. Imported context rows are "
        "retained as supplied; this stage discovers no public context.</p>"
        '<div class="table-scroll" tabindex="0" role="region" aria-label="Sample readiness">'
        "<table><thead><tr>"
        + "".join(
            f'<th scope="col">{name}</th>'
            for name in (
                "Sample",
                "Role",
                "Profile calls",
                "Typing",
                "Reason",
                "Missing metadata",
                "Assembly",
            )
        )
        + "</tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table></div>"
        '<details class="evidence-files" open><summary>Saved evidence</summary>'
        '<ul><li><a href="dataset.json">Dataset manifest</a></li>'
        '<li><a href="readiness.json">Machine-readable readiness audit</a></li>'
        '<li><a href="samples.csv">Samples and metadata</a></li>'
        '<li><a href="lineages.csv">Complete imported lineage evidence</a></li>'
        '<li><a href="provenance.csv">Field provenance</a></li></ul></details>'
        "</div></section></main></body></html>"
    )


def run_prepare(
    input_path: str | Path,
    output: str | Path,
    *,
    catalogues: str | Path | None = None,
    input_kind: str = "auto",
    species: str | None = None,
) -> PrepareResult:
    """Import frozen canonical JSON or an existing bundle into a new stage output.

    Explicit catalogues are mandatory for canonical records, even when all calls
    are missing. Existing outputs are immutable. Local assembly references are
    copied and checksummed, without provider acquisition or context selection.
    """
    source, target = Path(input_path).expanduser().resolve(), Path(output).expanduser().absolute()
    if target.exists():
        raise DatasetError(f"Preparation output already exists: {target}")
    if input_kind not in {"auto", "profiles", "dataset"}:
        raise DatasetError(
            "prepare supports profiles or dataset inputs; collection/accession resolution "
            "and new assembly typing await a query-only provider stage"
        )
    if species is not None and (not isinstance(species, str) or not species.strip()):
        raise DatasetError("--species must be nonempty text")
    if source.is_dir():
        source = source / "dataset.json"
    raw = _json(source)
    kind = input_kind
    if kind == "auto":
        kind = (
            "dataset" if isinstance(raw, dict) and raw.get("schema") == SCHEMA_NAME else "profiles"
        )
    catalogue_path = Path(catalogues).expanduser().resolve() if catalogues else None
    if kind == "dataset":
        if catalogue_path or species:
            raise DatasetError(
                "Existing dataset import preserves species/catalogues; omit overrides"
            )
        _validate_previous_stage(source.parent)
        dataset = load_dataset(source, mmap_mode=None)
    else:
        if catalogue_path is None:
            raise DatasetError(
                "Frozen profiles require --catalogues with an explicit locus catalogue"
            )
        dataset = _frozen_dataset(raw, _catalogues(catalogue_path), species)
    target.parent.mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(prefix=f".{target.name}.prepare-", dir=target.parent))
    try:
        assembly_staging = scratch / "assembly_files"
        assembly_staging.mkdir()
        dataset, assemblies = _portable_assemblies(dataset, source.parent, assembly_staging)
        source_evidence = {"path": "sources/input.json", "sha256": file_sha256(source)}
        parameters = dict(dataset.parameters)
        parameters["prepare"] = {
            "input_kind": kind,
            "source": source_evidence,
            "source_dataset_id": dataset.dataset_id,
            "species": species,
            "assemblies": assemblies,
        }
        if catalogue_path:
            parameters["prepare"]["catalogues"] = {
                "path": "sources/catalogues.json",
                "sha256": file_sha256(catalogue_path),
            }
        dataset = replace(dataset, parameters=parameters, dataset_id=None)
        bundle = scratch / "bundle"
        manifest = write_dataset(dataset, bundle)
        sources = bundle / "sources"
        sources.mkdir()
        shutil.copyfile(source, sources / "input.json")
        if catalogue_path:
            shutil.copyfile(catalogue_path, sources / "catalogues.json")
        if (assembly_staging / "assemblies").exists():
            shutil.move(str(assembly_staging / "assemblies"), bundle / "assemblies")
        audit = _readiness(dataset, manifest)
        write_json(bundle / "readiness.json", audit)
        (bundle / "report.html").write_text(_report(audit), encoding="utf-8")
        # Validate the actual saved dataset, then seal all stage-owned artifacts.
        load_dataset(manifest)
        artifacts = []
        for path in sorted(bundle.rglob("*")):
            if path.is_file():
                artifacts.append(
                    {
                        "path": path.relative_to(bundle).as_posix(),
                        "sha256": file_sha256(path),
                        "bytes": path.stat().st_size,
                    }
                )
        write_json(
            bundle / "prepare.json",
            {
                "schema": "chronoclade.prepare-result",
                "schema_version": 1,
                "software": {"name": "chronoclade", "version": __version__},
                "status": "complete",
                "input_kind": kind,
                "dataset_manifest": "dataset.json",
                "report": "report.html",
                "readiness_audit": "readiness.json",
                "readiness": audit["readiness"],
                "sample_ids": list(dataset.sample_ids),
                "artifacts": artifacts,
            },
        )
        bundle.rename(target)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    return PrepareResult(
        target / "dataset.json",
        target / "report.html",
        target / "readiness.json",
        target / "prepare.json",
    )
