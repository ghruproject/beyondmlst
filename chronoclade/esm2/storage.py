"""Atomic, checksum-validated numerical artifacts without pickle deserialisation."""

import hashlib
import json
import os
from pathlib import Path
import tempfile
from zipfile import BadZipFile

import numpy as np


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def provenance_key(provenance):
    return hashlib.sha256(canonical_json(provenance).encode()).hexdigest()


def vector_hash(vector):
    return hashlib.sha256(np.asarray(vector, dtype="<f4").tobytes()).hexdigest()


def _atomic(path, writer):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=f".{path.name}.",
                                         suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
            writer(handle)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def write_json(path, value):
    _atomic(path, lambda handle: handle.write((canonical_json(value) + "\n").encode()))


def write_npz(path, **arrays):
    _atomic(path, lambda handle: np.savez_compressed(handle, **arrays))


class VectorCache:
    """Content-addressed vectors scoped by full extraction/runtime provenance."""

    def __init__(self, directory, provenance):
        self.provenance = provenance
        self.directory = Path(directory) / provenance_key(provenance)
        self.invalid = 0

    def read(self, protein_id, dimension):
        path = self.directory / f"{protein_id}.npz"
        if not path.exists():
            return None
        try:
            with np.load(path, allow_pickle=False) as arrays:
                vector = arrays["vector"]
                info = json.loads(str(arrays["metadata"].item()))
            if (info["protein_id"] != protein_id or info["provenance"] != self.provenance
                    or vector.shape != (dimension,) or vector.dtype != np.float32
                    or not np.isfinite(vector).all() or info["sha256"] != vector_hash(vector)):
                raise ValueError("Cache identity, shape, precision or checksum mismatch")
            return vector
        except (OSError, ValueError, KeyError, EOFError, TypeError, BadZipFile):
            self.invalid += 1
            return None

    def write(self, protein_id, vector):
        info = {"protein_id": protein_id, "provenance": self.provenance,
                "sha256": vector_hash(vector)}
        write_npz(self.directory / f"{protein_id}.npz", vector=np.asarray(vector, dtype=np.float32),
                  metadata=np.asarray(canonical_json(info)))
