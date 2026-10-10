"""Locus-preserving genome distances from validated, frozen protein vectors.

Missing loci are exclusions, never zero vectors. No embedding tree is inferred.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np

from .proteins import EmbeddingError
from .sample_distances import _load_embeddings, _read_mapping, analyze_embedding_dates

MAX_GENOMES = 1500
DISTANCE_DEFINITION = "mean cosine distance over the same fixed protein locus panel"


def nearest_ties(matrix, ids, count=3):
    """Retain every tie at the requested neighbour boundary, excluding self."""
    result = []
    for i, ident in enumerate(ids):
        others = sorted((float(matrix[i, j]), other) for j, other in enumerate(ids) if i != j)
        cutoff = others[min(count, len(others)) - 1][0] if others else None
        chosen = [(distance, other) for distance, other in others if distance <= cutoff + 1e-12]
        result.append({"sample_id": ident, "requested_count": count,
                       "boundary_distance": cutoff,
                       "neighbours": [{"sample_id": other, "distance": distance}
                                      for distance, other in chosen]})
    return result


def _components(distances, ids, threshold):
    """Explicit connected components, not calibrated biological clades."""
    groups, unseen = {}, set(range(len(ids)))
    while unseen:
        pending = [min(unseen)]
        members = set()
        while pending:
            i = pending.pop()
            if i in members:
                continue
            members.add(i)
            unseen.discard(i)
            pending.extend(j for j in sorted(unseen) if distances[i, j] <= threshold + 1e-12)
        label = f"embedding-group-{len(set(groups.values())) + 1:03d}"
        groups.update({ids[i]: label for i in members})
    return groups


def analyse_genomes(dataset, embeddings_manifest: Path, mapping_csv: Path, *,
                    panel_loci=None, reference_sample_id=None, neighbour_count=3,
                    group_threshold=0.0):
    """Use a fixed equal-weight panel; ordination retains locus identity.

    The PCoA geometry is sqrt(2 * mean cosine), the Euclidean chord distance of
    concatenated unit vectors divided by sqrt(panel size). Its axes are not allele
    changes. The dense experimental backend has an explicit 1,500-genome limit.
    """
    if type(neighbour_count) is not int or neighbour_count < 1:
        raise EmbeddingError("neighbour_count must be positive")
    if not isinstance(group_threshold, (int, float)) or not np.isfinite(group_threshold) or not 0 <= group_threshold <= 2:
        raise EmbeddingError("group_threshold must be finite and between zero and two")
    diagnostic = analyze_embedding_dates(dataset, embeddings_manifest, mapping_csv,
                                          panel_loci=panel_loci,
                                          reference_sample_id=reference_sample_id)
    manifest, vectors, positions, record_proteins = _load_embeddings(embeddings_manifest)
    matrix = dataset.profiles[0]
    cells, profiles, _ = _read_mapping(mapping_csv, dataset, matrix, record_proteins)
    ids = [row["sample_id"] for row in diagnostic["samples"] if row["status"] == "included"]
    if len(ids) > MAX_GENOMES:
        raise EmbeddingError(f"Experimental dense ESM2 genome analysis supports at most {MAX_GENOMES} genomes; no samples were subsampled")
    panel = diagnostic["panel"]["loci"]
    all_loci = list(matrix.catalogue.loci)
    # Retain per-locus protein IDs and a complete genome/locus missing mask. This
    # is a compact representation, not an average that discards locus identity.
    protein_ids = np.asarray([[cells.get((ident, locus), "") for locus in all_loci]
                              for ident in dataset.sample_ids])
    mask = protein_ids != ""
    gram = np.zeros((len(ids), len(ids)), dtype=np.float64)
    variable_loci = []
    for locus in panel:
        proteins = [cells[ident, locus] for ident in ids]
        if len(set(proteins)) > 1:
            variable_loci.append(locus)
        unit = np.asarray([vectors[positions[protein]] for protein in proteins], dtype=np.float64)
        unit /= np.linalg.norm(unit, axis=1)[:, None]
        gram += unit @ unit.T / len(panel)
    distances = np.clip(1.0 - gram, 0.0, 2.0)
    np.fill_diagonal(distances, 0.0)
    # Identical protein panels are exactly zero even after floating-point dot products.
    for i in range(len(ids)):
        for j in range(i):
            if all(cells[ids[i], locus] == cells[ids[j], locus] for locus in panel):
                distances[i, j] = distances[j, i] = 0.0
    centered = gram - gram.mean(axis=0)[None, :] - gram.mean(axis=1)[:, None] + gram.mean()
    values, axes = np.linalg.eigh(centered)
    order = np.argsort(values)[::-1]
    values, axes = values[order], axes[:, order]
    tolerance = max(1.0, float(np.max(np.abs(values)))) * 1e-10
    positive = np.where(values > tolerance, values, 0.0)
    coordinates = np.zeros((len(ids), 2), dtype=np.float64)
    for axis in range(min(2, len(ids))):
        coordinates[:, axis] = axes[:, axis] * np.sqrt(positive[axis])
        # Stable sign convention for saved figures; no interpretation of direction.
        if coordinates[np.argmax(np.abs(coordinates[:, axis])), axis] < 0:
            coordinates[:, axis] *= -1
    fractions = (positive[:2] / positive.sum()).tolist() if positive.sum() else [0.0, 0.0]
    fractions += [0.0] * (2 - len(fractions))
    pairs, allele_distances = [], np.zeros_like(distances)
    synonymous = []
    for i, a in enumerate(ids):
        for j in range(i + 1, len(ids)):
            b = ids[j]
            changed = [locus for locus in panel if profiles[a][locus] != profiles[b][locus]]
            invisible = [locus for locus in changed if cells[a, locus] == cells[b, locus]]
            allele_distances[i, j] = allele_distances[j, i] = len(changed) / len(panel)
            pairs.append({"sample_id_1": a, "sample_id_2": b,
                          "distance": float(distances[i, j]), "panel_loci": len(panel),
                          "allele_differences": len(changed), "allele_mismatch_fraction": len(changed) / len(panel),
                          "protein_identical_changed_alleles": len(invisible)})
            if invisible:
                synonymous.append({"sample_id_1": a, "sample_id_2": b, "loci": invisible})
    neighbours = nearest_ties(distances, ids, neighbour_count)
    baseline = nearest_ties(allele_distances, ids, neighbour_count)
    comparison = []
    for embedding, allele in zip(neighbours, baseline, strict=True):
        a = {r["sample_id"] for r in embedding["neighbours"]}
        b = {r["sample_id"] for r in allele["neighbours"]}
        comparison.append({"sample_id": embedding["sample_id"], "embedding_ids": sorted(a),
                           "allele_ids": sorted(b), "shared_ids": sorted(a & b),
                           "jaccard": len(a & b) / len(a | b) if a | b else None})
    summary = {
        "schema": "chronoclade.esm2.genome-distances", "schema_version": 1,
        "status": "complete", "experimental": True,
        "sample_ids": ids, "dataset_id": dataset.dataset_id,
        "distance_definition": DISTANCE_DEFINITION, "distance_units": "mean locus cosine distance",
        "panel": {**diagnostic["panel"], "weighting": "equal weight per locus",
                  "variable_protein_loci": variable_loci,
                  "omitted_scheme_loci": sorted(set(all_loci) - set(panel))},
        "counts": diagnostic["counts"], "samples": diagnostic["samples"],
        "model": manifest["provenance"], "provenance": diagnostic["provenance"],
        "ordination": {"method": "exact centered locus-preserving unit-vector Gram eigendecomposition",
                       "distance_definition": "sqrt(2 * mean locus cosine distance)",
                       "axis_units": "protein chord distance", "axis_variance_fraction": fractions,
                       "negative_eigenvalues": int(np.sum(values < -tolerance))},
        "groups": {"method": "connected components at an explicit embedding-distance threshold",
                   "threshold": float(group_threshold), "calibrated": False,
                   "by_sample": _components(distances, ids, group_threshold)},
        "neighbours": neighbours,
        "biological_comparison": {"scope": "same frozen genomes and fixed callable loci",
                                  "allele_neighbour_agreement": comparison,
                                  "synonymous_protein_identity": synonymous,
                                  "corrected_snp_validation": "not_supplied",
                                  "scientific_default_validated": False},
        "exclusions": [r for r in diagnostic["samples"] if r["status"] == "excluded"],
        "limitations": ["Synonymous DNA changes can be invisible to protein embeddings.",
                        "Groups are exploratory connected components, not calibrated clades.",
                        "No embedding ancestry, substitution rate or molecular clock is inferred.",
                        f"Dense experimental analysis is limited to {MAX_GENOMES} genomes."],
    }
    return summary, {"distances": distances, "coordinates": coordinates,
                     "sample_ids": np.asarray(ids), "all_sample_ids": np.asarray(dataset.sample_ids),
                     "loci": np.asarray(all_loci), "protein_ids": protein_ids, "mapped_mask": mask,
                     "panel_loci": np.asarray(panel)}, pairs, diagnostic


def add_dna_comparison(summary, pairs, dataset, sequences, neighbour_count=3):
    """Add exact same-length DNA comparisons from the frozen catalogue.

    Length-changing locus pairs are recorded as unavailable rather than aligned
    with an unrecorded heuristic. Distances have one fixed panel denominator.
    """
    ids, panel = summary["sample_ids"], summary["panel"]["loci"]
    matrix = dataset.profiles[0]
    profiles = {ident: matrix.profile(ident) for ident in ids}
    index = {ident: i for i, ident in enumerate(ids)}
    distances = np.zeros((len(ids), len(ids)))
    unavailable = []
    for row in pairs:
        a, b = row["sample_id_1"], row["sample_id_2"]
        locus_distances = []
        failures = []
        differences, sites = 0, 0
        for locus in panel:
            left = sequences[locus, profiles[a][locus]].upper()
            right = sequences[locus, profiles[b][locus]].upper()
            if len(left) != len(right):
                failures.append(locus)
                continue
            mismatch = sum(x != y for x, y in zip(left, right, strict=True))
            locus_distances.append(mismatch / len(left))
            differences += mismatch
            sites += len(left)
        value = float(np.mean(locus_distances)) if not failures else None
        row.update(dna_mean_locus_mismatch_fraction=value,
                   dna_nucleotide_differences=differences if not failures else None,
                   dna_compared_sites=sites if not failures else None,
                   dna_incompatible_length_loci=failures)
        distances[index[a], index[b]] = distances[index[b], index[a]] = value if value is not None else np.nan
        if failures:
            unavailable.append({"sample_id_1": a, "sample_id_2": b, "loci": failures,
                                "reason": "unequal_cds_lengths_no_alignment_requested"})
    comparison = {"status": "partial" if unavailable else "complete",
                  "distance_definition": "mean per-locus DNA mismatch fraction on same-length CDS",
                  "fixed_panel_loci": panel, "incompatible_pairs": unavailable}
    if not unavailable:
        neighbours = nearest_ties(distances, ids, neighbour_count)
        comparison["nearest_neighbours"] = neighbours
        comparison["embedding_neighbour_agreement"] = []
        for dna, embedding in zip(neighbours, summary["neighbours"], strict=True):
            a = {r["sample_id"] for r in dna["neighbours"]}
            b = {r["sample_id"] for r in embedding["neighbours"]}
            comparison["embedding_neighbour_agreement"].append(
                {"sample_id": dna["sample_id"], "shared_ids": sorted(a & b),
                 "jaccard": len(a & b) / len(a | b) if a | b else None})
    summary["biological_comparison"]["dna_comparison"] = comparison
