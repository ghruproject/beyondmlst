"""Run Pathogenwatch's native cgMLST and lineage assigners on query assemblies.

Prepared databases are supplied explicitly; this module never builds databases,
handles API credentials, or substitutes a neighbour's unresolved lineage code.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import shutil
import subprocess
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from chronoclade.cglin import DEFAULT_SCHEME, normalise_assignment
from chronoclade.metadata import SAFE_IDENTIFIER, Sample


class QueryTypingError(ValueError):
    """A typing configuration or result cannot safely be used."""


def normalise_typing_species(value: Any) -> str:
    """Return a supported taxon key while retaining user labels in output rows."""
    if not isinstance(value, str):
        raise QueryTypingError("Query typing requires a supported species identity")
    key = " ".join(value.replace("_", " ").casefold().split())
    if key not in {"klebsiella pneumoniae", "escherichia coli"}:
        raise QueryTypingError(f"Unsupported native typing species: {value}")
    return key


def _same_species(left: Any, right: Any) -> bool:
    try:
        return normalise_typing_species(left) == normalise_typing_species(right)
    except QueryTypingError:
        return False


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _database_hash(path: Path) -> str:
    if path.is_file():
        return _sha256(path)
    files = sorted(item for item in path.rglob("*") if item.is_file())
    if not files:
        raise QueryTypingError(f"Empty typing database: {path}")
    digest = hashlib.sha256()
    for item in files:
        digest.update(item.relative_to(path).as_posix().encode())
        digest.update(b"\0")
        digest.update(_sha256(item).encode())
        digest.update(b"\n")
    return digest.hexdigest()


def _path(value: Any, base: Path) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise QueryTypingError("Typing database paths must be nonempty strings")
    path = Path(value).expanduser()
    path = (base / path).resolve() if not path.is_absolute() else path.resolve()
    if not path.exists():
        raise QueryTypingError(f"Typing resource does not exist: {path}")
    return path


def _command(value: Any, base: Path) -> list[str]:
    if (
        not isinstance(value, list)
        or not value
        or any(not isinstance(x, str) or not x for x in value)
    ):
        raise QueryTypingError("Tool command must be a nonempty JSON array of arguments")
    command = list(value)
    executable = shutil.which(command[0])
    if executable is None:
        candidate = base / command[0]
        if candidate.is_file():
            executable = str(candidate.resolve())
    if executable is None:
        raise QueryTypingError(f"Typing executable is unavailable: {command[0]}")
    command[0] = executable
    # Script paths can be relative to the configuration, not the process cwd.
    for index, arg in enumerate(command[1:], 1):
        if (base / arg).is_file():
            command[index] = str((base / arg).resolve())
    return command


def load_typing_config(path: Path) -> dict[str, Any]:
    """Validate native tool and frozen database paths, resolving relative paths."""
    path = path.expanduser().resolve()
    try:
        config = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        raise QueryTypingError("Cannot read typing configuration JSON") from exc
    if isinstance(config, dict) and config.get("ready") is False:
        raise QueryTypingError(
            "Typing databases are not ready. Run setup-typing and supply the prepared database files described in its setup instructions."
        )
    organisms = config.get("organisms") if isinstance(config, dict) else None
    if not isinstance(organisms, dict) or not organisms:
        raise QueryTypingError("Typing config requires a nonempty organisms mapping")
    species_keys = set()
    for species, organism in organisms.items():
        if not isinstance(species, str) or not isinstance(organism, dict):
            raise QueryTypingError("Invalid organism typing configuration")
        species_key = normalise_typing_species(species)
        if species_key in species_keys:
            raise QueryTypingError("Typing config contains duplicate aliases of the same species")
        species_keys.add(species_key)
        for name in ("cgmlst", "assignment"):
            tool = organism.get(name)
            if not isinstance(tool, dict):
                raise QueryTypingError(f"{species}: missing {name} configuration")
            for field in ("scheme", "scheme_version", "tool_version"):
                if (
                    not isinstance(tool.get(field), str)
                    or not tool[field].strip()
                    or tool[field] == "unknown"
                ):
                    raise QueryTypingError(f"{species}: {name} requires pinned {field}")
            tool["command"] = _command(tool.get("command"), path.parent)
        caller, assignment = organism["cgmlst"], organism["assignment"]
        caller["index_dir"] = _path(caller.get("index_dir"), path.parent)
        if not caller["index_dir"].is_dir():
            raise QueryTypingError("cgMLST index_dir must be a directory")
        assignment["loci_file"] = _path(assignment.get("loci_file"), path.parent)
        assignment["canonical_loci"] = _canonical_loci(assignment["loci_file"])
        kind = assignment.get("kind")
        expected_kind = "plincer" if species_key == "klebsiella pneumoniae" else "hclink"
        if kind != expected_kind:
            raise QueryTypingError("Native typing assigner does not match its configured species")

        if kind == "plincer":
            if caller["scheme"] != "klebsiella_1" or assignment["scheme"] != DEFAULT_SCHEME:
                raise QueryTypingError("plincer requires klebsiella_1 and scgMLST629_S schemes")
            for field in ("scheme_file", "profiles_file", "alleles_db"):
                assignment[field] = _path(assignment.get(field), path.parent)
                if not assignment[field].is_file():
                    raise QueryTypingError(f"plincer {field} must be a file")
        elif kind == "hclink":
            if caller["scheme"] != "ecoli_1":
                raise QueryTypingError("E. coli hclink requires the ecoli_1 cgMLST scheme")
            assignment["reference_db"] = _path(assignment.get("reference_db"), path.parent)
            try:
                caller_date = date.fromisoformat(caller["scheme_version"])
                assignment_date = date.fromisoformat(assignment["scheme_version"])
                metadata = json.loads((assignment["reference_db"] / "metadata.json").read_text())
                actual_date = datetime.fromisoformat(metadata["datestamp"]).date()
            except (ValueError, KeyError, OSError, TypeError) as exc:
                raise QueryTypingError(
                    "hclink requires ISO database dates and metadata.json datestamp"
                ) from exc
            if caller_date > assignment_date or actual_date != assignment_date:
                raise QueryTypingError(
                    "hclink database must match its declared date and be at least as recent as cgMLST"
                )
        else:
            raise QueryTypingError("Assignment kind must be plincer or hclink")
    return config


def _canonical_loci(path: Path) -> list[str]:
    try:
        raw = json.loads(path.read_text())
        genes = raw.get("genes") if isinstance(raw, dict) else raw
    except (OSError, ValueError) as exc:
        raise QueryTypingError("Cannot read assignment database ordered loci JSON") from exc
    if (
        not isinstance(genes, list)
        or not genes
        or any(not isinstance(gene, str) or not gene for gene in genes)
    ):
        raise QueryTypingError("Assignment database requires a nonempty ordered loci list")
    genes = [gene.replace(" ", "_") for gene in genes]
    if len(set(genes)) != len(genes):
        raise QueryTypingError("Assignment database locus IDs must be unique")
    return genes


def _check_locus_order(genes: list[str], canonical: list[str]) -> None:
    if [gene.replace(" ", "_") for gene in genes] != canonical:
        raise QueryTypingError(
            "cgMLST locus order does not match the assignment database; rebuild compatible databases before typing"
        )


def _provenance(tool: Mapping[str, Any], database_fields: Sequence[str]) -> dict[str, Any]:
    command = tool["command"]
    return {
        "version": tool["tool_version"],
        "command": command,
        "tool_files_sha256": {arg: _sha256(Path(arg)) for arg in command if Path(arg).is_file()},
        "databases": {
            field: {"path": str(tool[field]), "sha256": _database_hash(tool[field])}
            for field in database_fields
        },
        "scheme": tool["scheme"],
        "scheme_version": tool["scheme_version"],
    }


def _run(command: list[str], data: bytes, destination: Path) -> dict[str, Any]:
    try:
        result = subprocess.run(command, input=data, capture_output=True, timeout=3600, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise QueryTypingError(
            f"Native typing tool could not complete: {Path(command[0]).name}"
        ) from exc
    if result.returncode:
        # Never expose subprocess stderr, which may include configured credentials.
        raise QueryTypingError(
            f"Native typing tool {Path(command[0]).name} failed with exit code {result.returncode}"
        )
    try:
        parsed = json.loads(result.stdout)
    except (UnicodeDecodeError, ValueError) as exc:
        raise QueryTypingError("Native typing tool did not return a JSON object") from exc
    if not isinstance(parsed, dict):
        raise QueryTypingError("Native typing tool did not return a JSON object")
    destination.write_text(json.dumps(parsed, indent=2) + "\n")
    return parsed


def parse_cgmlst_result(raw: Mapping[str, Any], scheme: str, version: str) -> dict[str, Any]:
    """Keep known alleles, missing calls, novel hashes and duplicate hits distinct."""
    genes = raw.get("genes")
    code = raw.get("code")
    if raw.get("scheme") != scheme or not isinstance(genes, list) or not genes:
        raise QueryTypingError("cgMLST output has an incompatible scheme or no ordered loci")
    if any(not isinstance(gene, str) or not gene for gene in genes) or len(set(genes)) != len(
        genes
    ):
        raise QueryTypingError("cgMLST output locus IDs must be unique nonempty strings")
    if (
        not isinstance(code, str)
        or len(code.split("_")) != len(genes)
        or raw.get("schemeSize", len(genes)) != len(genes)
    ):
        raise QueryTypingError("cgMLST code length does not match its ordered loci")
    profile, novel, missing, ambiguous = {}, {}, [], []
    hits = raw.get("alleles", {})
    for gene, allele in zip(genes, code.split("_"), strict=True):
        gene_hits = hits.get(gene, [])
        if isinstance(gene_hits, list) and len(gene_hits) > 1:
            ambiguous.append(gene)
            continue
        if not allele:
            missing.append(gene)
        elif allele.isdigit() and int(allele) > 0:
            profile[gene] = str(int(allele))
        elif re.fullmatch(r"[0-9a-fA-F]{40}", allele):
            novel[gene] = allele.lower()
        else:
            raise QueryTypingError("cgMLST output contains an invalid allele token")
    return {
        "cgmlst_scheme": scheme,
        "cgmlst_scheme_version": version,
        "cgmlst_loci": genes,
        "cgmlst_profile": profile,
        "cgmlst_novel_alleles": novel,
        "cgmlst_missing_loci": missing,
        "cgmlst_ambiguous_loci": ambiguous,
        "cgmlst_called_fraction": (len(profile) + len(novel)) / len(genes),
        "cgmlst_known_fraction": len(profile) / len(genes),
        "cgmlst_code": code,
    }


def _assignment_result(raw: dict[str, Any], config: Mapping[str, Any]) -> dict[str, Any]:
    if config["kind"] == "plincer":
        code = raw.get("LINcode", [])
        if not isinstance(code, list):
            raise QueryTypingError("plincer LINcode must be an array")
        # Multiple best profiles can disagree: retain only their common prefix.
        codes = [code] + [match.get("LINcode", []) for match in raw.get("matches", [])]
        conservative = []
        for index, part in enumerate(code):
            if str(part) == "*" or any(
                len(other) <= index or str(other[index]) != str(part) for other in codes
            ):
                break
            conservative.append(str(part))
        if len(conservative) < 10 and conservative:
            conservative += ["*"] * (10 - len(conservative))
        return normalise_assignment(
            {
                "LINcode": "_".join(conservative),
                "cgST": raw.get("cgST", ""),
                "cglin_source": "pathogenwatch_plincer",
            },
            scheme=config["scheme"],
            scheme_version=config["scheme_version"],
        )
    codes = raw.get("hierCC")
    if not isinstance(codes, list):
        raise QueryTypingError("hclink output requires a hierCC array")
    result = {}
    for pair in codes:
        if (
            not isinstance(pair, list)
            or len(pair) != 2
            or not re.fullmatch(r"(?:d|HC)\d+", str(pair[0]))
            or (str(pair[1]) != "" and not str(pair[1]).isdigit())
        ):
            raise QueryTypingError("Invalid hclink HierCC assignment")
        if str(pair[1]) == "":
            continue
        level = "HC" + re.sub(r"^(?:d|HC)", "", str(pair[0]))
        if level in result:
            raise QueryTypingError("Duplicate hclink HierCC level")
        result[level] = str(pair[1])
    return {
        "hiercc_codes": result,
        "hiercc_scheme": config["scheme"],
        "hiercc_scheme_version": config["scheme_version"],
        "hiercc_source": "pathogenwatch_hclink",
        "hiercc_status": "predicted" if result else "unassigned",
        "hiercc_distance": raw.get("hierccDistance", raw.get("hierCC_distance")),
    }


def type_query_assemblies(
    samples: Sequence[Sample],
    config_path: Path,
    output: Path,
    *,
    existing_assignments: Mapping[str, dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """Type each assembly using pinned native tools; write query_typing.json."""
    config = load_typing_config(config_path)
    output.mkdir(parents=True, exist_ok=True)
    organisms = {
        normalise_typing_species(species): settings
        for species, settings in config["organisms"].items()
    }
    results, provenance = [], {}
    seen = set()
    for sample in samples:
        if not SAFE_IDENTIFIER.fullmatch(sample.sample_id) or sample.sample_id in seen:
            raise QueryTypingError("Query sample IDs must be unique safe identifiers")
        seen.add(sample.sample_id)
        assembly_hash = _sha256(sample.assembly)
        existing = (existing_assignments or {}).get(sample.sample_id)
        species_key = normalise_typing_species(sample.species)
        organism = organisms.get(species_key)
        if organism is None:
            raise QueryTypingError(f"No native typing configuration for {sample.species}")
        caller, assigner = organism["cgmlst"], organism["assignment"]
        if species_key not in provenance:
            fields = (
                ("scheme_file", "profiles_file", "alleles_db", "loci_file")
                if assigner["kind"] == "plincer"
                else ("reference_db", "loci_file")
            )
            provenance[species_key] = {
                "cgmlst": _provenance(caller, ("index_dir",)),
                "assignment": _provenance(assigner, fields),
            }
        if (
            existing
            and existing.get("assembly_sha256") == assembly_hash
            and _same_species(existing.get("species"), sample.species)
            and existing.get("typing_provenance") == provenance[species_key]
        ):
            results.append(
                dict(
                    existing, sample_id=sample.sample_id, species=sample.species, typing_reused=True
                )
            )
            continue
        directory = output / sample.sample_id
        directory.mkdir(exist_ok=True)
        opener = gzip.open if sample.assembly.suffix == ".gz" else open
        with opener(sample.assembly, "rb") as handle:
            assembly = handle.read()
        raw = _run(
            caller["command"]
            + [f"--scheme={caller['scheme']}", f"--indexDir={caller['index_dir']}"],
            assembly,
            directory / "cgmlst.json",
        )
        parsed = parse_cgmlst_result(raw, caller["scheme"], caller["scheme_version"])
        _check_locus_order(raw["genes"], assigner["canonical_loci"])
        # Upstream assigns numeric code even for duplicate loci. Blank ambiguous
        # loci before classification so they cannot drive lineage assignments.
        safe_raw = dict(raw)
        safe_raw["code"] = "_".join(
            "" if gene in parsed["cgmlst_ambiguous_loci"] else allele
            for gene, allele in zip(raw["genes"], raw["code"].split("_"), strict=True)
        )
        if parsed["cgmlst_ambiguous_loci"]:
            safe_raw["st"] = "unassigned"
        if assigner["kind"] == "plincer":
            argv = [
                "classify",
                "-",
                "--scheme-file",
                str(assigner["scheme_file"]),
                "--profiles-file",
                str(assigner["profiles_file"]),
                "--alleles-db",
                str(assigner["alleles_db"]),
            ]
        else:
            argv = ["assign", "-", "--exact", "--reference-db", str(assigner["reference_db"])]
        if not parsed["cgmlst_profile"] and not parsed["cgmlst_novel_alleles"]:
            assignment = {"typing_status": "unassigned"}
        else:
            assignment_raw = _run(
                assigner["command"] + argv,
                json.dumps(safe_raw).encode(),
                directory / "assignment.json",
            )
            assignment = _assignment_result(assignment_raw, assigner)
            assignment["typing_status"] = (
                "typed"
                if assignment.get("cglin_resolved_depth", 0) or assignment.get("hiercc_codes")
                else "unassigned"
            )
            assignment["typing_assignment_evidence"] = assignment_raw

        results.append(
            {
                "sample_id": sample.sample_id,
                "species": sample.species,
                "assembly_sha256": assembly_hash,
                "typing_status": "typed",
                **parsed,
                **assignment,
                "typing_provenance": provenance[species_key],
                "cgmlst_database_sha256": provenance[species_key]["cgmlst"]["databases"][
                    "index_dir"
                ]["sha256"],
                (
                    "cglin_database_sha256"
                    if assigner["kind"] == "plincer"
                    else "hiercc_database_sha256"
                ): hashlib.sha256(
                    json.dumps(
                        {
                            field: info["sha256"]
                            for field, info in provenance[species_key]["assignment"][
                                "databases"
                            ].items()
                        },
                        sort_keys=True,
                    ).encode()
                ).hexdigest(),
            }
        )
    envelope = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "assignments": results,
    }
    (output / "query_typing.json").write_text(json.dumps(envelope, indent=2) + "\n")
    return results


def load_query_typing(path: Path) -> list[dict[str, Any]]:
    """Load native typing output; assembly identity is checked by its consumer."""
    try:
        envelope = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        raise QueryTypingError("Cannot read query typing JSON") from exc
    if (
        not isinstance(envelope, dict)
        or envelope.get("schema_version") != 1
        or not isinstance(envelope.get("assignments"), list)
    ):
        raise QueryTypingError("Unsupported query typing schema")
    seen = set()
    for row in envelope["assignments"]:
        if (
            not isinstance(row, dict)
            or not isinstance(row.get("sample_id"), str)
            or not SAFE_IDENTIFIER.fullmatch(row["sample_id"])
            or row["sample_id"] in seen
        ):
            raise QueryTypingError("Query typing IDs must be unique safe identifiers")
        seen.add(row["sample_id"])
        if not re.fullmatch(r"[0-9a-f]{64}", str(row.get("assembly_sha256", ""))):
            raise QueryTypingError("Query typing requires an assembly SHA256")
        if "cgmlst_profile" in row and (
            not isinstance(row["cgmlst_profile"], dict)
            or not isinstance(row.get("cgmlst_loci"), list)
        ):
            raise QueryTypingError("Query cgMLST profile requires its declared locus universe")
        if not isinstance(row.get("species"), str) or not row["species"].strip():
            raise QueryTypingError("Query typing requires species identity")
        species_key = normalise_typing_species(row["species"])
        if "cgmlst_profile" in row:
            loci, profile = row["cgmlst_loci"], row["cgmlst_profile"]
            if (
                not loci
                or any(not isinstance(x, str) or not x for x in loci)
                or len(set(loci)) != len(loci)
            ):
                raise QueryTypingError("Invalid query typing locus universe")
            if not set(profile).issubset(loci) or any(
                not str(v).isdigit() or int(v) < 1 for v in profile.values()
            ):
                raise QueryTypingError("Invalid known query cgMLST alleles")
            for field in ("cgmlst_scheme", "cgmlst_scheme_version"):
                if not isinstance(row.get(field), str) or row[field] in {"", "unknown"}:
                    raise QueryTypingError("Query profile requires a pinned scheme")
        expected_scheme = "klebsiella_1" if species_key == "klebsiella pneumoniae" else "ecoli_1"
        if "cgmlst_scheme" in row and row["cgmlst_scheme"] != expected_scheme:
            raise QueryTypingError("Query cgMLST scheme does not match its species")
        if "cgmlst_novel_alleles" in row:
            novel = row["cgmlst_novel_alleles"]
            if (
                not isinstance(novel, dict)
                or not set(novel).issubset(row.get("cgmlst_loci", []))
                or any(not re.fullmatch(r"[0-9a-f]{40}", str(v)) for v in novel.values())
            ):
                raise QueryTypingError("Invalid novel query cgMLST allele checksums")
            if set(novel).intersection(row.get("cgmlst_profile", {})):
                raise QueryTypingError("Query locus cannot be both known and novel")
        if "cglin_raw" in row:
            if species_key != "klebsiella pneumoniae":
                raise QueryTypingError("cgLIN query assignment requires Klebsiella pneumoniae")
            if (
                row.get("cglin_scheme") != DEFAULT_SCHEME
                or not isinstance(row.get("cglin_scheme_version"), str)
                or row["cglin_scheme_version"] in {"", "unknown"}
            ):
                raise QueryTypingError("Query cgLIN requires a pinned compatible scheme")
            assignment = normalise_assignment(row)
            if assignment["cglin_code_status"] in {"malformed", "unsupported"}:
                raise QueryTypingError("Invalid query cgLIN assignment")
        if "hiercc_codes" in row:
            if species_key != "escherichia coli":
                raise QueryTypingError("HierCC query assignment requires Escherichia coli")
            codes = row["hiercc_codes"]
            if not isinstance(codes, dict) or any(
                not re.fullmatch(r"HC\d+", str(k)) or not str(v).isdigit() for k, v in codes.items()
            ):
                raise QueryTypingError("Invalid query HierCC assignment")
            if any(
                not isinstance(row.get(k), str) or row[k] in {"", "unknown"}
                for k in ("hiercc_scheme", "hiercc_scheme_version")
            ):
                raise QueryTypingError("Query HierCC requires a pinned scheme")
    return envelope["assignments"]
