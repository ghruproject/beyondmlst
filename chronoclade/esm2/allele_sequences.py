"""Offline, scope-checked allele DNA translation with explicit sample mappings.

Allele categories are lookup keys only. Neither proteins nor missing DNA are
inferred from an allele identifier, related allele or neighbouring sample.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
import hashlib
import io
import json
from pathlib import Path
import re

from Bio.Data import CodonTable
from Bio.Seq import Seq

from chronoclade import __version__
from chronoclade.artifacts import file_sha256, write_json
from chronoclade.datasets import PreparedDataset
from .proteins import EmbeddingError

CATALOGUE_SCHEMA = "chronoclade.esm2.allele-sequence-catalogue"
MAPPING_SCHEMA = "chronoclade.esm2.allele-sequence-mapping"
_SCOPE_FIELDS = ("scheme_id", "scheme_version", "database_version", "database_sha256")


class SequenceMappingError(EmbeddingError):
    """Sequence preparation failed; a post-validation audit may have been saved."""

    def __init__(self, message, *, manifest_path=None):
        super().__init__(message)
        self.manifest_path = manifest_path


@dataclass(frozen=True)
class SequenceMappingResult:
    manifest_path: Path
    fasta_path: Path
    mapping_path: Path
    manifest: dict


def _digest(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _text(value):
    return isinstance(value, str) and bool(value) and value == value.strip()


def _read_catalogue(path, dataset, matrix):
    """Require frozen input bytes and exact scheme/database/species scope."""
    path = Path(path).resolve()
    try:
        raw = path.read_bytes()
        manifest = json.loads(raw)
        if (
            manifest.get("schema") != CATALOGUE_SCHEMA
            or manifest.get("schema_version") != 1
            or manifest.get("status") != "complete"
        ):
            raise SequenceMappingError("Expected a complete allele-sequence catalogue v1")
        scope = manifest["scope"]
        if not isinstance(scope, dict) or set(scope) != {"species", *_SCOPE_FIELDS}:
            raise SequenceMappingError("Sequence catalogue requires an explicit species and scope")
        for key in ("species", "scheme_id", "scheme_version"):
            if not _text(scope[key]):
                raise SequenceMappingError(f"Sequence catalogue requires nonempty {key}")
        if scope["database_version"] is not None and not _text(scope["database_version"]):
            raise SequenceMappingError("Catalogue database_version must be nonempty text or null")
        fingerprint = scope["database_sha256"]
        if not isinstance(fingerprint, str) or not re.fullmatch(r"[0-9a-f]{64}", fingerprint):
            raise SequenceMappingError(
                "Sequence catalogue requires a known database SHA256 fingerprint"
            )
        for key in _SCOPE_FIELDS:
            if scope[key] != getattr(matrix.catalogue, key):
                raise SequenceMappingError(
                    f"Sequence catalogue {key} is incompatible with prepared data"
                )
        # No taxonomic aliases are guessed. Unknown or mixed species fail closed,
        # even for samples with no allele calls, which remain in the same cohort.
        if any(row.get("species") != scope["species"] for row in dataset.samples):
            raise SequenceMappingError(
                "Sequence catalogue species is incompatible with prepared samples"
            )
        descriptor = manifest["artifact"]
        relative = Path(descriptor["path"])
        if relative.is_absolute() or ".." in relative.parts or not relative.parts:
            raise SequenceMappingError("Catalogue DNA CSV path must be a local relative path")
        dna_path = (path.parent / relative).resolve()
        if not dna_path.is_relative_to(path.parent):
            raise SequenceMappingError(
                "Catalogue DNA CSV must remain within its manifest directory"
            )
        csv_bytes = dna_path.read_bytes()
        actual_hash = hashlib.sha256(csv_bytes).hexdigest()
        if descriptor["sha256"] != actual_hash:
            raise SequenceMappingError("Catalogue DNA CSV SHA256 mismatch")
        reader = csv.DictReader(io.StringIO(csv_bytes.decode("utf-8"), newline=""))
        if reader.fieldnames != ["locus", "allele", "dna"]:
            raise SequenceMappingError(
                "Catalogue CSV requires exactly locus,allele,dna columns in order"
            )
        sequences = {}
        loci = set(matrix.catalogue.loci)
        for number, row in enumerate(reader, 2):
            if None in row or any(not _text(row.get(key)) for key in reader.fieldnames):
                raise SequenceMappingError(
                    f"Catalogue CSV row {number} has missing or malformed values"
                )
            locus, allele = row["locus"], row["allele"]
            if locus not in loci:
                raise SequenceMappingError(
                    f"Catalogue locus {locus!r} is outside the prepared scheme"
                )
            if allele.casefold() in {"0", "-", "?", "none", "null", "unknown", "na", "n/a"}:
                raise SequenceMappingError(
                    f"Catalogue row {number} uses a missing marker as an allele"
                )
            key = (locus, allele)
            if key in sequences:
                raise SequenceMappingError(f"Catalogue repeats locus/allele {key!r}")
            sequences[key] = row["dna"]
        return (
            sequences,
            scope,
            {
                "manifest_path": str(path),
                "manifest_sha256": hashlib.sha256(raw).hexdigest(),
                "dna_csv_path": str(dna_path),
                "dna_csv_sha256": actual_hash,
                "catalogue_records": len(sequences),
            },
        )
    except SequenceMappingError:
        raise
    except (
        OSError,
        UnicodeError,
        ValueError,
        KeyError,
        TypeError,
        AttributeError,
        csv.Error,
    ) as exc:
        raise SequenceMappingError(f"Invalid frozen allele-sequence catalogue: {exc}") from exc


def _translate(dna, *, terminal_stop, max_length):
    dna = dna.upper()
    evidence = {"dna_sha256": _digest(dna), "dna_length": len(dna)}
    if set(dna) - set("ACGT"):
        return None, "ambiguous_or_unsupported_dna", evidence
    if len(dna) % 3:
        return None, "out_of_frame", evidence
    table = CodonTable.unambiguous_dna_by_id[11]
    if dna[:3] not in table.start_codons:
        return None, "invalid_start_codon", evidence
    has_stop = dna[-3:] in table.stop_codons
    evidence.update(terminal_stop_removed=has_stop, initiation_to_methionine=dna[:3] != "ATG")
    if terminal_stop == "require" and not has_stop:
        return None, "missing_terminal_stop", evidence
    coding = dna[:-3] if has_stop else dna
    codons = [coding[position : position + 3] for position in range(0, len(coding), 3)]
    if any(codon in table.stop_codons for codon in codons):
        return None, "internal_stop", evidence
    # The first codon uses explicit table-11 CDS initiation semantics, including
    # alternative initiators such as GTG. This is recorded, not a guessed repair.
    protein = "M" + str(Seq(coding[3:]).translate(table=11))
    evidence["protein_length"] = len(protein)
    if len(protein) > max_length:
        return None, "protein_exceeds_max_length", evidence
    return protein, None, evidence


def prepare_allele_sequences(
    dataset: PreparedDataset,
    catalogue_manifest: Path,
    out: Path,
    *,
    require_complete: bool = False,
    invalid_policy: str = "exclude",
    terminal_stop: str = "require",
    genetic_code: int = 11,
    max_length: int = 1022,
) -> SequenceMappingResult:
    """Translate called alleles and publish a checked, auditable offline bundle.

    Exactly one prepared profile matrix and a matching known database fingerprint
    are required. ``terminal_stop='allow_absent'`` explicitly permits a stop-less
    coding sequence; any present terminal stop is removed. Internal stops, gaps,
    ambiguous DNA, invalid initiators, frames and overlength proteins are excluded,
    never repaired. ``invalid_policy='error'`` or ``require_complete=True`` saves a
    failed audit and raises SequenceMappingError with its ``manifest_path``.

    The output directory must be new or empty. Every mapped allele has its own
    FASTA record ID; the existing embedding reader deduplicates identical proteins
    while retaining those record IDs. All samples and all locus cells are audited.
    """
    dataset.validate()
    if len(dataset.profiles) != 1:
        raise SequenceMappingError("Sequence mapping requires exactly one profile scheme matrix")
    if type(require_complete) is not bool:
        raise SequenceMappingError("require_complete must be boolean")
    if invalid_policy not in {"error", "exclude"}:
        raise SequenceMappingError("invalid_policy must be 'error' or 'exclude'")
    if terminal_stop not in {"require", "allow_absent"}:
        raise SequenceMappingError("terminal_stop must be 'require' or 'allow_absent'")
    if type(genetic_code) is not int or genetic_code != 11:
        raise SequenceMappingError("Only explicit bacterial genetic_code=11 is supported")
    if type(max_length) is not int or not 1 <= max_length <= 1022:
        raise SequenceMappingError(
            "max_length must be between 1 and 1022; truncation is unsupported"
        )
    matrix = dataset.profiles[0]
    sequences, scope, input_evidence = _read_catalogue(catalogue_manifest, dataset, matrix)
    profiles = {ident: matrix.profile(ident) for ident in matrix.sample_ids}
    called = sorted(
        {
            (locus, allele)
            for profile in profiles.values()
            for locus, allele in profile.items()
            if allele is not None
        }
    )
    allele_audit, proteins, by_allele = [], {}, {}
    for locus, allele in called:
        record_id = "allele_" + _digest(
            json.dumps(
                {**scope, "locus": locus, "allele": allele}, sort_keys=True, separators=(",", ":")
            )
        )
        record = {"locus": locus, "allele": allele, "record_id": record_id}
        dna = sequences.get((locus, allele))
        if dna is None:
            record.update(status="missing_sequence", reason="allele_not_in_frozen_catalogue")
        else:
            protein, reason, evidence = _translate(
                dna, terminal_stop=terminal_stop, max_length=max_length
            )
            record.update(evidence)
            if reason:
                record.update(status="invalid_cds", reason=reason)
            else:
                record.update(status="mapped", reason=None, protein_id=_digest(protein))
                proteins[record_id] = protein
        allele_audit.append(record)
        by_allele[locus, allele] = record
    cells, samples, mapping = [], [], []
    for sample in dataset.sample_ids:
        mapped = 0
        for locus in matrix.catalogue.loci:
            allele = profiles.get(sample, {}).get(locus)
            cell = {"sample_id": sample, "locus": locus, "allele": allele}
            if sample not in profiles:
                cell.update(status="missing_profile", reason="sample_not_in_profile_matrix")
            elif allele is None:
                cell.update(status="missing_call", reason="no_allele_call")
            else:
                record = by_allele[locus, allele]
                cell.update(status=record["status"], reason=record["reason"])
                if record["status"] == "mapped":
                    cell.update(record_id=record["record_id"], protein_id=record["protein_id"])
                    mapping.append({key: cell[key] for key in ("sample_id", "locus", "record_id")})
                    mapped += 1
            cells.append(cell)
        samples.append(
            {
                "sample_id": sample,
                "mapped_loci": mapped,
                "total_loci": len(matrix.catalogue.loci),
                "complete": mapped == len(matrix.catalogue.loci),
            }
        )
    unavailable = len(cells) - len(mapping)
    invalid = sum(record["status"] == "invalid_cds" for record in allele_audit)
    failed = (require_complete and unavailable > 0) or (invalid_policy == "error" and invalid > 0)
    out = Path(out).resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise SequenceMappingError("Sequence mapping output directory must be new or empty")
    out.mkdir(parents=True, exist_ok=True)
    fasta_path, mapping_path = out / "proteins.fasta", out / "sample_loci.csv"
    with fasta_path.open("x", encoding="utf-8") as stream:
        for record_id, protein in proteins.items():
            stream.write(f">{record_id}\n{protein}\n")
    with mapping_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=("sample_id", "locus", "record_id"))
        writer.writeheader()
        writer.writerows(mapping)
    manifest = {
        "schema": MAPPING_SCHEMA,
        "schema_version": 1,
        "software_version": __version__,
        "status": "failed" if failed else "complete",
        "experimental": True,
        "coverage": "complete" if not unavailable else "partial",
        "dataset_id": dataset.dataset_id,
        "prepared_input": {
            "samples_sha256": _digest(
                json.dumps(dataset.samples, sort_keys=True, separators=(",", ":"), allow_nan=False)
            ),
            "profiles_sha256": _digest(
                json.dumps(profiles, sort_keys=True, separators=(",", ":"), allow_nan=False)
            ),
            "loci_sha256": _digest(json.dumps(matrix.catalogue.loci, separators=(",", ":"))),
        },
        "scope": scope,
        "input": input_evidence,
        "parameters": {
            "require_complete": require_complete,
            "invalid_policy": invalid_policy,
            "terminal_stop": terminal_stop,
            "genetic_code": genetic_code,
            "max_length": max_length,
            "dna_normalisation": "uppercase",
            "truncation": False,
            "imputation": False,
            "unreferenced_catalogue_alleles": "not-translated",
        },
        "counts": {
            "samples": len(samples),
            "sample_locus_cells": len(cells),
            "mapped_cells": len(mapping),
            "unavailable_cells": unavailable,
            "called_alleles": len(called),
            "mapped_alleles": len(proteins),
            "unique_proteins": len(set(proteins.values())),
            "invalid_alleles": invalid,
            "missing_sequence_alleles": sum(
                r["status"] == "missing_sequence" for r in allele_audit
            ),
        },
        "samples": samples,
        "alleles": allele_audit,
        "cells": cells,
        "artifacts": {
            "proteins": {"path": fasta_path.name, "sha256": file_sha256(fasta_path)},
            "sample_loci": {"path": mapping_path.name, "sha256": file_sha256(mapping_path)},
        },
    }
    manifest_path = out / "sequence_mapping.json"
    write_json(manifest_path, manifest)
    if failed:
        raise SequenceMappingError(
            "Incomplete or invalid sequence mapping; see saved audit", manifest_path=manifest_path
        )
    return SequenceMappingResult(manifest_path, fasta_path, mapping_path, manifest)
