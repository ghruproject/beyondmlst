"""Experimental date diagnostics from a fixed panel of mapped protein loci.

The mean is taken over locus-specific cosine distances to one actual sample.
Neither averaged genome vectors nor inferred ancestral embeddings are used.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
import re
from zipfile import BadZipFile

import numpy as np

from chronoclade.datasets import PreparedDataset
from chronoclade.temporal_diagnostics import date_regression, sample_date_interval
from .proteins import EmbeddingError, file_sha256
from .runtime import POOLING, model_spec


def _load_embeddings(path: Path):
    """Read only a complete, checksummed artifact with auditable model identity."""
    path = Path(path)
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
        if (
            manifest.get("schema") != "chronoclade.esm2.protein-embeddings"
            or manifest.get("schema_version") != 1
            or manifest.get("status") != "complete"
        ):
            raise EmbeddingError("Expected a complete ESM2 protein-embedding manifest v1")
        provenance = manifest["provenance"]
        spec = model_spec(provenance["model"])
        requested = model_spec(manifest["parameters"]["requested_model"])
        if provenance["model"] != spec["name"] or requested != spec:
            raise EmbeddingError("Actual model provenance does not match the requested model")
        if (
            not re.fullmatch(r"[0-9a-f]{64}", provenance.get("model_sha256", ""))
            or provenance.get("precision") != "float32"
            or provenance.get("pooling") != POOLING
            or provenance.get("actual_device") not in {"cpu", "mps", "cuda"}
            or provenance.get("actual_device") != manifest.get("actual_device")
            or provenance.get("dimension", spec["dimension"]) != spec["dimension"]
            or provenance.get("representation_layer", spec["layers"]) != spec["layers"]
        ):
            raise EmbeddingError("Invalid actual model, pooling or device provenance")
        descriptor = manifest["artifacts"]["embeddings"]
        relative = Path(descriptor["path"])
        if relative.is_absolute() or ".." in relative.parts or not relative.parts:
            raise EmbeddingError("Embedding vector path must be contained and relative")
        vectors_path = (path.parent / relative).resolve()
        if not vectors_path.is_relative_to(path.parent.resolve()):
            raise EmbeddingError("Embedding vector path escapes its manifest directory")
        if file_sha256(vectors_path) != descriptor["sha256"]:
            raise EmbeddingError("Embedding vector artifact SHA256 mismatch")
        ids = manifest["protein_ids"]
        if (
            not ids
            or any(not isinstance(p, str) or not re.fullmatch(r"[0-9a-f]{64}", p) for p in ids)
            or len(ids) != len(set(ids))
        ):
            raise EmbeddingError("Manifest requires unique protein SHA256 IDs")
        known_proteins = set(ids)
        shape = (len(ids), spec["dimension"])
        with np.load(vectors_path, allow_pickle=False) as arrays:
            saved_ids = arrays["protein_ids"]
            vectors = arrays["embeddings"]
        if (
            saved_ids.ndim != 1
            or saved_ids.tolist() != ids
            or vectors.shape != shape
            or descriptor.get("shape") != list(shape)
            or descriptor.get("dtype") != "float32"
            or vectors.dtype != np.float32
            or not np.isfinite(vectors).all()
        ):
            raise EmbeddingError(
                "Embedding artifact IDs, shape, precision or finite values invalid"
            )
        squared_norms = np.einsum("ij,ij->i", vectors, vectors, dtype=np.float64)
        if not np.isfinite(squared_norms).all() or np.any(squared_norms == 0):
            raise EmbeddingError("Embedding artifact contains an invalid zero vector norm")
        record_proteins = {}
        for record in manifest["mapping"]:
            ident, protein = record["record_id"], record["protein_id"]
            if not isinstance(ident, str) or not ident or ident in record_proteins:
                raise EmbeddingError("Manifest contains invalid or repeated record IDs")
            if protein not in known_proteins:
                raise EmbeddingError("Manifest record mapping refers to an unknown protein ID")
            record_proteins[ident] = protein
        if not record_proteins:
            raise EmbeddingError("Manifest has no retained protein record mapping")
        return manifest, vectors, dict(zip(ids, range(len(ids)), strict=True)), record_proteins
    except EmbeddingError:
        raise
    except (
        OSError,
        UnicodeError,
        ValueError,
        KeyError,
        TypeError,
        AttributeError,
        EOFError,
        BadZipFile,
    ) as exc:
        raise EmbeddingError(f"Invalid saved embedding artifact: {exc}") from exc


def _read_mapping(path, dataset, matrix, record_proteins):
    """Require explicit sample/locus/record identity and consistent allele translation."""
    ids, loci = set(dataset.sample_ids), set(matrix.catalogue.loci)
    matrix_ids = set(matrix.sample_ids)
    profiles = {ident: matrix.profile(ident) for ident in matrix.sample_ids}
    cells, allele_proteins, rows = {}, {}, 0
    try:
        with Path(path).open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            if (
                reader.fieldnames is None
                or len(set(reader.fieldnames)) != len(reader.fieldnames)
                or set(reader.fieldnames) != {"sample_id", "locus", "record_id"}
            ):
                raise EmbeddingError(
                    "Mapping CSV requires exactly sample_id,locus,record_id columns"
                )
            for row in reader:
                rows += 1
                sample, locus, record = (row.get(k) for k in ("sample_id", "locus", "record_id"))
                if None in row or any(not x or x != x.strip() for x in (sample, locus, record)):
                    raise EmbeddingError(f"Mapping row {rows + 1} has missing or malformed IDs")
                if sample not in ids:
                    raise EmbeddingError(f"Mapping refers to unknown sample {sample!r}")
                if locus not in loci:
                    raise EmbeddingError(f"Mapping locus {locus!r} is outside the selected scheme")
                if record not in record_proteins:
                    raise EmbeddingError(f"Mapping refers to unknown retained record {record!r}")
                if (sample, locus) in cells:
                    raise EmbeddingError(
                        f"Duplicate mapping cell for sample {sample!r}, locus {locus!r}"
                    )
                if sample not in matrix_ids or profiles[sample][locus] is None:
                    raise EmbeddingError(
                        f"Mapped sample {sample!r}, locus {locus!r} has no allele call"
                    )
                protein = record_proteins[record]
                allele_key = (locus, profiles[sample][locus])
                if allele_key in allele_proteins and allele_proteins[allele_key] != protein:
                    raise EmbeddingError(
                        f"Locus {locus!r} has inconsistent protein IDs for one allele"
                    )
                allele_proteins[allele_key] = protein
                cells[sample, locus] = protein
        if not cells:
            raise EmbeddingError("Mapping CSV contains no mapped sample/locus cells")
        return cells, profiles, rows
    except EmbeddingError:
        raise
    except (OSError, UnicodeError, csv.Error, TypeError) as exc:
        raise EmbeddingError(f"Invalid sample/locus mapping CSV: {exc}") from exc


def analyze_embedding_dates(
    dataset: PreparedDataset,
    embeddings_manifest: Path,
    mapping_csv: Path,
    *,
    reference_sample_id=None,
    panel_loci=None,
) -> dict:
    """Describe mapped sample distances against canonical collection-date intervals.

    The default panel is the intersection across all samples with any mappings;
    wholly unmapped samples are retained in the audit. An explicit panel excludes
    incomplete samples. Reference choice and panel selection never use dates.
    """
    dataset.validate()
    if len(dataset.profiles) != 1:
        raise EmbeddingError("ESM2 date diagnostics require exactly one profile scheme matrix")
    matrix = dataset.profiles[0]
    manifest, vectors, positions, record_proteins = _load_embeddings(embeddings_manifest)
    cells, profiles, mapping_rows = _read_mapping(mapping_csv, dataset, matrix, record_proteins)
    mapped = {ident: set() for ident in dataset.sample_ids}
    for ident, locus in cells:
        mapped[ident].add(locus)
    mapped_ids = sorted(ident for ident, loci in mapped.items() if loci)
    if panel_loci is None:
        panel = sorted(set.intersection(*(mapped[ident] for ident in mapped_ids)))
        selection = "mapped-intersection"
    else:
        if isinstance(panel_loci, str):
            raise EmbeddingError("Fixed panel must be a sequence of locus IDs")
        try:
            panel = list(panel_loci)
        except TypeError as exc:
            raise EmbeddingError("Fixed panel must be a sequence of locus IDs") from exc
        if any(not isinstance(locus, str) for locus in panel) or len(set(panel)) != len(panel):
            raise EmbeddingError("Fixed panel requires unique locus IDs")
        if set(panel) - set(matrix.catalogue.loci):
            raise EmbeddingError("Fixed panel contains loci outside the selected scheme")
        panel.sort()
        selection = "explicit-fixed"
    if not panel:
        raise EmbeddingError("A fixed shared panel of at least one mapped locus is required")
    eligible = sorted(ident for ident in mapped_ids if set(panel) <= mapped[ident])
    if not eligible:
        raise EmbeddingError("No sample has complete protein mapping for the fixed panel")
    metadata_by_id = {row["sample_id"]: row for row in dataset.samples}
    species_groups = {
        str(metadata_by_id[ident].get("species") or "").casefold().replace("_", " ")
        for ident in eligible
    }
    if len(species_groups) > 1:
        raise EmbeddingError(
            "ESM2 date diagnostics require one species group among eligible samples; "
            "prepare separate species cohorts"
        )
    if reference_sample_id is not None and reference_sample_id not in eligible:
        raise EmbeddingError(
            "Explicit reference sample must have complete mapping for the fixed panel"
        )
    reference = reference_sample_id if reference_sample_id is not None else eligible[0]
    mapped_proteins = set(cells.values())
    units = {}
    for protein in sorted(mapped_proteins):
        vector = vectors[positions[protein]].astype(np.float64)
        norm = np.linalg.norm(vector)
        if not np.isfinite(norm) or norm == 0:
            raise EmbeddingError(f"Mapped protein {protein!r} has an invalid zero vector norm")
        units[protein] = vector / norm
    samples, points, exclusions = [], [], []
    for metadata in sorted(dataset.samples, key=lambda row: row["sample_id"]):
        ident = metadata["sample_id"]
        date = sample_date_interval(metadata)
        reason = None
        distance = None
        locus_distances = {}
        if not mapped[ident]:
            reason = (
                "missing_profile"
                if ident not in profiles or not any(profiles[ident].values())
                else "no_mapped_proteins"
            )
        elif ident not in eligible:
            reason = "incomplete_fixed_panel"
        else:
            locus_distances = {
                locus: 0.0
                if cells[ident, locus] == cells[reference, locus]
                else float(
                    np.clip(
                        1.0 - np.dot(units[cells[ident, locus]], units[cells[reference, locus]]),
                        0.0,
                        2.0,
                    )
                )
                for locus in panel
            }
            distance = float(np.mean(list(locus_distances.values())))
            points.append({**metadata, "distance": distance})
        sample = {
            **metadata,
            "distance": distance,
            "locus_distances": locus_distances,
            "status": "excluded" if reason else "included",
            "reason": reason,
            "mapped_loci": sorted(mapped[ident]),
            "missing_panel_loci": sorted(set(panel) - mapped[ident]),
            "date_status": date["status"],
            "date_interval": date,
        }
        samples.append(sample)
        if reason:
            exclusions.append(
                {
                    "sample_id": ident,
                    "reason": reason,
                    "missing_panel_loci": sample["missing_panel_loci"],
                }
            )
    distance_units = "mean locus cosine distance"
    regression = date_regression(
        points,
        distance_field="distance",
        distance_units=distance_units,
        label="Experimental ESM2 distance–date association",
    )
    exclusions.extend(regression["exclusions"])
    catalogue = matrix.catalogue
    return {
        "schema": "chronoclade.esm2.temporal-diagnostics",
        "schema_version": 1,
        "mode": "esm2-distance",
        "experimental": True,
        "status": regression["status"],
        "units": {"distance": distance_units, "slope": distance_units + "/year"},
        "reference": {
            "sample_id": reference,
            "selection": "explicit"
            if reference_sample_id is not None
            else "lexicographic-eligible-id",
            "date_independent": True,
        },
        "panel": {"loci": panel, "count": len(panel), "selection": selection},
        "scheme": {
            "scheme_id": catalogue.scheme_id,
            "scheme_version": catalogue.scheme_version,
            "database_version": catalogue.database_version,
            "database_sha256": catalogue.database_sha256,
        },
        "provenance": {
            "embeddings_manifest": str(Path(embeddings_manifest).resolve()),
            "mapping_csv": str(Path(mapping_csv).resolve()),
            "manifest_sha256": file_sha256(embeddings_manifest),
            "mapping_sha256": file_sha256(mapping_csv),
            "model": manifest["provenance"],
            "vector_artifact": manifest["artifacts"]["embeddings"],
            "dataset_id": dataset.dataset_id,
        },
        "counts": {
            "total_samples": len(samples),
            "mapped_samples": len(mapped_ids),
            "eligible_samples": len(eligible),
            "dated_points": len(regression["points"]),
            "excluded_samples": len(samples) - len(eligible),
            "mapping_rows": mapping_rows,
            "unique_mapped_proteins": len(mapped_proteins),
        },
        "samples": samples,
        "points": regression["points"],
        "regression": regression,
        "exclusions": exclusions,
        "limitations": [
            "Experimental descriptive embedding-distance association.",
            "Partial collection dates use interval midpoints and retain their bounds.",
            "No substitution rate, ancestral embedding or MRCA estimate is inferred.",
        ],
    }
