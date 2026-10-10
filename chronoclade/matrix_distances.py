"""Validated, memory-mapped categorical distance evidence and lazy pair maps."""

from __future__ import annotations

import json
from pathlib import Path
import numpy as np
from chronoclade.artifacts import file_sha256


class MatrixDistanceMap:
    """Pair lookup without materialising quadratic Python dictionaries."""

    def __init__(self, evidence, raw=False):
        self.evidence, self.raw = evidence, raw

    def __contains__(self, pair):
        return self.evidence.get(*pair) is not None

    def __getitem__(self, pair):
        item = self.evidence.get(*pair)
        if item is None:
            raise KeyError(pair)
        return item["allele_differences" if self.raw else "distance"]

    def get(self, pair, default=None):
        try:
            return self[pair]
        except KeyError:
            return default


class BinaryDistanceEvidence:
    """Load checksummed manifests and check arrays in bounded row batches.

    Missing or incompatible pairs have NaN distance; integer arrays retain callable
    counts for low-overlap pairs, which are never made selectable.
    """

    def __init__(self, manifest_path, *, validate=True):
        self.manifest_path = Path(manifest_path).resolve()
        data = json.loads(self.manifest_path.read_text())
        if not isinstance(data, dict):
            raise ValueError("Binary distance manifest must be an object")
        if (
            data.get("schema") != "chronoclade.binary_distances"
            or type(data.get("schema_version")) is not int
            or data.get("schema_version") != 1
        ):
            raise ValueError("Unsupported binary distance evidence schema")
        identifiers = data.get("sample_ids")
        if not isinstance(identifiers, list):
            raise ValueError("Binary distance sample_ids must be a list")
        self.sample_ids = tuple(identifiers)
        if (
            not self.sample_ids
            or any(not isinstance(x, str) or not x for x in self.sample_ids)
            or len(set(self.sample_ids)) != len(self.sample_ids)
        ):
            raise ValueError("Binary distance sample IDs must be unique and nonempty")
        self.method = data.get("method")
        self.distance_definition = data.get("distance_definition")
        if (
            self.method != "cgmlst"
            or self.distance_definition != "fraction of mismatching jointly called cgMLST loci"
        ):
            raise ValueError("Unsupported binary distance method or definition")
        self.scope = data.get("scope", {})
        if not isinstance(self.scope, dict):
            raise ValueError("Binary distance scope must be an object")
        if type(self.scope.get("sample_count")) is not int or self.scope.get("sample_count") != len(
            self.sample_ids
        ):
            raise ValueError("Binary distance scope sample_count disagrees with IDs")
        overlap = self.scope.get("min_overlap")
        if type(overlap) not in (int, float) or not 0 < overlap <= 1:
            raise ValueError("Binary distance scope needs a callable-overlap fraction in (0,1]")
        cohorts = self.scope.get("cohort_sample_ids")
        if not isinstance(cohorts, list) or any(
            not isinstance(group, list) or not group for group in cohorts
        ):
            raise ValueError("Binary distance scope requires complete-comparability cohorts")
        members = [name for group in cohorts for name in group]
        if any(not isinstance(name, str) or not name for name in members):
            raise ValueError("Binary distance cohort identifiers must be nonempty strings")
        if len(set(members)) != len(members) or not set(members).issubset(self.sample_ids):
            raise ValueError("Binary distance cohort identifiers repeat or escape the sample scope")
        maximum_pairs = len(self.sample_ids) * (len(self.sample_ids) - 1) // 2
        comparable = self.scope.get("comparable_pairs")
        reasons = self.scope.get("unavailable_pairs_by_reason")
        if type(comparable) is not int or not 0 <= comparable <= maximum_pairs:
            raise ValueError("Binary distance scope comparable_pairs must be a valid pair count")
        if not isinstance(reasons, dict) or any(
            not isinstance(reason, str) or type(count) is not int or count < 0
            for reason, count in reasons.items()
        ):
            raise ValueError(
                "Binary distance unavailable-pair audit must have nonnegative integer counts"
            )
        if sum(reasons.values()) + comparable != maximum_pairs:
            raise ValueError(
                "Binary distance pair audit does not account for the complete sample scope"
            )
        self.paths = data.get("paths")
        if not isinstance(self.paths, dict):
            raise ValueError("Binary distance paths must be an object")
        self.descriptor = data
        self._index = {name: i for i, name in enumerate(self.sample_ids)}
        self.arrays = {}
        for kind in ("distance", "allele_differences", "shared_called_loci", "scheme_loci"):
            spec = self.paths.get(kind)
            if (
                not isinstance(spec, dict)
                or not isinstance(spec.get("path"), str)
                or not spec["path"]
            ):
                raise ValueError(f"Binary distance file reference is invalid: {kind}")
            digest = spec.get("sha256")
            if (
                not isinstance(digest, str)
                or len(digest) != 64
                or any(c not in "0123456789abcdef" for c in digest)
            ):
                raise ValueError(f"Binary distance checksum reference is invalid: {kind}")
            path = (self.manifest_path.parent / spec["path"]).resolve()
            if not path.is_relative_to(self.manifest_path.parent):
                raise ValueError("Binary evidence path escapes its manifest directory")
            if validate and file_sha256(path) != spec["sha256"]:
                raise ValueError(f"Binary distance checksum mismatch: {kind}")
            array = np.load(path, mmap_mode="r", allow_pickle=False)
            if not isinstance(array, np.ndarray):
                raise ValueError(f"Binary distance file must contain one NumPy array: {kind}")
            if array.shape != (len(self.sample_ids),) * 2:
                raise ValueError(f"Binary distance shape mismatch: {kind}")
            expected = np.dtype("float64" if kind == "distance" else "uint32")
            if array.dtype != expected:
                raise ValueError(f"Binary distance dtype mismatch: {kind}")
            self.arrays[kind] = array
        if validate:
            self._validate()

    def _validate(self):
        n = len(self.sample_ids)
        comparable = 0
        for start in range(0, n, 128):
            stop = min(start + 128, n)
            d, raw, shared, scheme = (
                self.arrays[k][start:stop]
                for k in ("distance", "allele_differences", "shared_called_loci", "scheme_loci")
            )
            valid = np.isfinite(d)
            if np.any(np.isinf(d)) or np.any(d[valid] < 0) or np.any(d[valid] > 1):
                raise ValueError("Invalid binary distance values")
            if np.any(raw > shared) or np.any(shared > scheme):
                raise ValueError("Invalid binary callable counts")
            offdiag = valid.copy()
            offdiag[np.arange(stop - start), np.arange(start, stop)] = False
            if np.any(shared[offdiag] == 0) or not np.allclose(
                d[offdiag], raw[offdiag] / shared[offdiag], rtol=0, atol=1e-12
            ):
                raise ValueError("Binary distances disagree with mismatch/callable counts")
            if np.any(shared[offdiag] / scheme[offdiag] < self.scope["min_overlap"]):
                raise ValueError(
                    "Finite binary distances do not meet the declared callable overlap"
                )
            for offset, row in enumerate(valid):
                comparable += int(row[: start + offset].sum())
            for key, array in self.arrays.items():
                if not np.array_equal(array[start:stop], array[:, start:stop].T, equal_nan=True):
                    raise ValueError(f"Asymmetric binary distance evidence: {key}")
        if comparable != self.scope.get("comparable_pairs"):
            raise ValueError("Binary distance scope comparable_pairs disagrees with matrices")
        for group in self.scope["cohort_sample_ids"]:
            indices = [self._index[name] for name in group]
            for index in indices:
                shared = int(self.arrays["shared_called_loci"][index, index])
                scheme = int(self.arrays["scheme_loci"][index, index])
                if not scheme or not shared or shared / scheme < self.scope["min_overlap"]:
                    raise ValueError("Binary distance cohort contains self-ineligible profiles")
                if not np.isfinite(self.arrays["distance"][index, indices]).all():
                    raise ValueError("Binary distance cohort contains unavailable pairs")
        if np.any(np.diag(self.arrays["allele_differences"]) != 0):
            raise ValueError("Binary raw allele differences must have a zero diagonal")
        if not np.all(np.diag(self.arrays["distance"]) == 0):
            raise ValueError("Binary distance diagonal must be zero")

    def value_map(self, raw=False):
        return MatrixDistanceMap(self, raw)

    def get(self, left, right):
        if left not in self._index or right not in self._index:
            return None
        i, j = self._index[left], self._index[right]
        value = float(self.arrays["distance"][i, j])
        if not np.isfinite(value):
            return None
        shared = int(self.arrays["shared_called_loci"][i, j])
        scheme = int(self.arrays["scheme_loci"][i, j])
        if i == j and (not shared or not scheme or shared / scheme < self.scope["min_overlap"]):
            return None
        return dict(
            sample_a=left,
            sample_b=right,
            distance=value,
            allele_differences=int(self.arrays["allele_differences"][i, j]),
            shared_called_loci=shared,
            scheme_loci=scheme,
            call_overlap=shared / scheme if scheme else 0.0,
        )

    def row(self, sample_id):
        return self.arrays["distance"][self._index[sample_id]]

    def iter_pairs(self, sample_ids=None):
        ids = self.sample_ids if sample_ids is None else tuple(sample_ids)
        for i, left in enumerate(ids):
            for right in ids[i + 1 :]:
                item = self.get(left, right)
                if item is not None:
                    yield item

    def export_pairs(self, output, *, sample_ids=None, pairs=None):
        """Export explicitly requested rows, including auditable unavailable pairs."""
        import csv
        from itertools import combinations

        if (sample_ids is None) == (pairs is None):
            raise ValueError("Choose explicit sample_ids or explicit pairs for distance export")
        requested = combinations(tuple(sample_ids), 2) if sample_ids is not None else pairs
        output = Path(output)
        fields = [
            "sample_id_1",
            "sample_id_2",
            "status",
            "distance",
            "allele_differences",
            "shared_called_loci",
            "scheme_loci",
            "call_overlap",
        ]
        with output.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for left, right in requested:
                if left not in self._index or right not in self._index:
                    raise ValueError("Requested export identifiers are outside binary evidence")
                item = self.get(left, right)
                i, j = self._index[left], self._index[right]
                shared = int(self.arrays["shared_called_loci"][i, j])
                scheme = int(self.arrays["scheme_loci"][i, j])
                writer.writerow(
                    dict(
                        sample_id_1=left,
                        sample_id_2=right,
                        status="comparable" if item else "unavailable",
                        distance=item["distance"] if item else "",
                        allele_differences=int(self.arrays["allele_differences"][i, j]),
                        shared_called_loci=shared,
                        scheme_loci=scheme,
                        call_overlap=shared / scheme if scheme else "",
                    )
                )
        return output
