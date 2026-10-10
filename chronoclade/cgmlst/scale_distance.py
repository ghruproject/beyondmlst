"""Batched categorical distance computation with disk-backed exact evidence."""

from collections import Counter, defaultdict
from pathlib import Path
import json
import numpy as np

from chronoclade.artifacts import file_sha256
from chronoclade.context_refinement import _allele, _profile, _loci
from chronoclade.profile_analysis import _pair
from chronoclade.typing_scopes import typing_scope
from chronoclade.matrix_distances import BinaryDistanceEvidence


def prepare_records(records):
    """Use precisely the legacy declared/observed denominator rules."""
    records = [dict(row) for row in sorted(records, key=lambda row: row["sample_id"])]
    universes = defaultdict(set)

    def scope(row):
        return (
            str(row.get("species", "")).casefold().replace("_", " "),
            typing_scope(row, "cgmlst"),
            row.get("cgmlst_database_sha256", ""),
        )

    for row in records:
        if row.get("cgmlst_locus_universe_complete") is False:
            universes[scope(row)].update(row.get("cgmlst_loci") or [])
            universes[scope(row)].update((_profile(row) or {}).keys())
    for row in records:
        if row.get("cgmlst_locus_universe_complete") is False:
            row["cgmlst_loci"] = sorted(universes[scope(row)])
    return records


def encode_calls(rows, loci):
    """Encode each locus independently; zero exclusively means missing."""
    encoded = np.zeros((len(rows), len(loci)), dtype=np.float64)
    profiles = [_profile(row) for row in rows]
    for column, locus in enumerate(loci):
        labels = {}
        for i, profile in enumerate(profiles):
            value = _allele(profile.get(locus))
            if value is not None:
                encoded[i, column] = labels.setdefault(value, len(labels) + 1)
    return encoded


def categorical_blocks(calls, batch_size=128):
    """Compiled SciPy Hamming kernels; missing calls never become matches.

    Hamming on encoded calls counts mismatches including one-missing pairs.
    Hamming on presence masks counts exactly those one-missing positions; their
    subtraction gives called mismatches. The presence XOR and row call totals
    give shared counts without any numeric interpretation of allele labels.
    """
    from scipy.spatial.distance import cdist

    n, loci = calls.shape
    present = np.asarray(calls != 0, dtype=np.float64)
    counts = present.sum(axis=1).astype(np.int64)
    complete = bool(np.all(counts == loci))
    for i in range(0, n, batch_size):
        a = slice(i, min(n, i + batch_size))
        for j in range(i, n, batch_size):
            b = slice(j, min(n, j + batch_size))
            differences = np.rint(cdist(calls[a], calls[b], metric="hamming") * loci).astype(
                np.int64
            )
            if complete:
                shared = np.full(differences.shape, loci, dtype=np.int64)
            else:
                xor = np.rint(cdist(present[a], present[b], metric="hamming") * loci).astype(
                    np.int64
                )
                shared = (counts[a, None] + counts[None, b] - xor) // 2
                differences -= xor
            yield a, b, differences.astype(np.uint32), shared.astype(np.uint32)


def write_distances(records, output, *, min_overlap=0.9, batch_size=128):
    """Write all pair evidence and deterministic complete-comparability cohorts."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    n = len(records)
    arrays = {
        key: np.lib.format.open_memmap(
            output / f"{key}.npy",
            mode="w+",
            dtype="float64" if key == "distance" else "uint32",
            shape=(n, n),
        )
        for key in ("distance", "allele_differences", "shared_called_loci", "scheme_loci")
    }
    arrays["distance"][:] = np.nan
    np.fill_diagonal(arrays["distance"], 0)
    eligible, exclusions, groups = [], [], defaultdict(list)
    for i, row in enumerate(records):
        evidence, reason = _pair(row, row, min_overlap)
        if reason:
            exclusions.append(dict(sample_id=row["sample_id"], reason=reason))
        else:
            eligible.append(i)
        # All profiles, including self-ineligible ones, retain exact cross-pair
        # counts where their namespace is defined. Eligible coverage is separate.
        profile = _profile(row)
        loci = tuple(sorted(_loci(row, profile))) if profile is not None else ()
        species = str(row.get("species", "")).casefold().replace("_", " ")
        signature = (species, typing_scope(row, "cgmlst"), loci)
        groups[signature].append(i)
    incompatible = Counter()
    signatures = list(groups)
    audit_buckets = {}
    for signature, indices in groups.items():
        buckets = defaultdict(list)
        for index in indices:
            row = records[index]
            buckets[str(row.get("cgmlst_database_sha256") or "").strip().lower()].append(index)
        audit_buckets[signature] = list(buckets.values())
    for k, left in enumerate(signatures):
        for right in signatures[k + 1 :]:
            for left_indices in audit_buckets[left]:
                for right_indices in audit_buckets[right]:
                    a, b = records[left_indices[0]], records[right_indices[0]]
                    _, reason = _pair(a, b, min_overlap)
                    incompatible[reason or "different_comparison_namespace"] += len(
                        left_indices
                    ) * len(right_indices)
    for signature, indices in groups.items():
        rows = [records[i] for i in indices]
        if signature[1] is None or not signature[2] or any(_profile(row) is None for row in rows):
            _, reason = _pair(rows[0], rows[0], min_overlap)
            incompatible[reason or "profile_unavailable"] += len(rows) * (len(rows) - 1) // 2
            continue
        loci = signature[2]
        calls = encode_calls(rows, loci)
        hashes = np.asarray(
            [str(row.get("cgmlst_database_sha256") or "").strip().lower() for row in rows]
        )
        for a, b, raw, shared in categorical_blocks(calls, batch_size):
            ii, jj = np.asarray(indices[a]), np.asarray(indices[b])
            distance = np.full(raw.shape, np.nan)
            compatible = (
                (hashes[a, None] == hashes[None, b])
                | (hashes[a, None] == "")
                | (hashes[None, b] == "")
            )
            overlap = (shared > 0) & (shared / len(loci) >= min_overlap)
            valid = compatible & overlap
            np.divide(raw, shared, out=distance, where=valid)
            unique = (
                np.ones(raw.shape, dtype=bool)
                if a.start != b.start
                else np.triu(np.ones(raw.shape, dtype=bool), 1)
            )
            incompatible["incompatible_or_unversioned_scheme"] += int(np.sum(unique & ~compatible))
            incompatible["insufficient_called_overlap"] += int(
                np.sum(unique & compatible & ~overlap)
            )
            for key, values in (
                ("distance", distance),
                ("allele_differences", raw),
                ("shared_called_loci", shared),
                ("scheme_loci", np.full(raw.shape, len(loci), dtype=np.uint32)),
            ):
                arrays[key][np.ix_(ii, jj)] = values
                arrays[key][np.ix_(jj, ii)] = values.T
    np.fill_diagonal(arrays["distance"], 0)
    for array in arrays.values():
        array.flush()
    cohorts = []
    for index in eligible:
        for cohort in cohorts:
            if np.isfinite(arrays["distance"][index, cohort]).all():
                cohort.append(index)
                break
        else:
            cohorts.append([index])
    count = sum(int(np.isfinite(arrays["distance"][i, :i]).sum()) for i in range(n))
    manifest = dict(
        schema="chronoclade.binary_distances",
        schema_version=1,
        method="cgmlst",
        distance_definition="fraction of mismatching jointly called cgMLST loci",
        sample_ids=[row["sample_id"] for row in records],
        scope=dict(
            min_overlap=min_overlap,
            sample_count=n,
            comparable_pairs=count,
            unavailable_pairs_by_reason={k: v for k, v in incompatible.items() if v},
            cohort_sample_ids=[[records[i]["sample_id"] for i in indexes] for indexes in cohorts],
        ),
        paths={
            key: dict(path=f"{key}.npy", sha256=file_sha256(output / f"{key}.npy"))
            for key in arrays
        },
    )
    path = output / "distance_evidence.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return BinaryDistanceEvidence(path, validate=False), cohorts, exclusions


def cohort_matrix(evidence, indexes, path):
    """Save a cohort square matrix and SciPy condensed vector without RAM copies."""
    n = len(indexes)
    matrix = np.lib.format.open_memmap(path, mode="w+", dtype="float64", shape=(n, n))
    source = evidence.arrays["distance"]
    for i, source_index in enumerate(indexes):
        matrix[i] = source[source_index, indexes]
    matrix.flush()
    condensed = np.lib.format.open_memmap(
        Path(path).with_name(Path(path).stem + "_condensed.npy"),
        mode="w+",
        dtype="float64",
        shape=(n * (n - 1) // 2,),
    )
    cursor = 0
    for i in range(n - 1):
        width = n - i - 1
        condensed[cursor : cursor + width] = matrix[i, i + 1 :]
        cursor += width
    condensed.flush()
    return matrix, condensed
