"""Immutable prepared bundles with checksums and atomic completion publication."""

from __future__ import annotations

import csv
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import tempfile

import numpy as np

from chronoclade import __version__
from .model import (
    AlleleMatrix,
    DatasetError,
    EVIDENCE_TABLES,
    LocusCatalogue,
    PreparedDataset,
    SAMPLE_COLUMNS,
    SCHEMA_NAME,
    SCHEMA_VERSION,
    check_json,
)


def _json(value):
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, allow_nan=False, separators=(",", ":")
    )


def _digest(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path, value):
    path.write_text(_json(value) + "\n", encoding="utf-8")


def _artifact(root, path, format_name):
    return {
        "path": path.relative_to(root).as_posix(),
        "sha256": _digest(path),
        "bytes": path.stat().st_size,
        "format": format_name,
    }


def _sample_array(rows):
    # Unicode fields are typed, column-addressable and pickle-free. Null masks
    # distinguish absent values from empty strings. Extra metadata stays JSON.
    extras = [
        _json({key: value for key, value in row.items() if key not in SAMPLE_COLUMNS})
        for row in rows
    ]
    present = [_json([key for key in SAMPLE_COLUMNS if key in row]) for row in rows]
    dtype = []
    for key in SAMPLE_COLUMNS:
        length = max(1, *(len(row.get(key) or "") for row in rows))
        dtype.extend([(key, f"U{length}"), (f"{key}__null", "?")])
    dtype.extend(
        [
            ("__extra_json", f"U{max(map(len, extras), default=1)}"),
            ("__present_json", f"U{max(map(len, present), default=1)}"),
        ]
    )
    array = np.empty(len(rows), dtype=dtype)
    for key in SAMPLE_COLUMNS:
        array[key] = [row.get(key) or "" for row in rows]
        array[f"{key}__null"] = [row.get(key) is None for row in rows]
    array["__extra_json"], array["__present_json"] = extras, present
    return array


def _read_samples(array):
    expected = {name for key in SAMPLE_COLUMNS for name in (key, f"{key}__null")}
    if array.ndim != 1 or set(array.dtype.names or ()) != expected | {
        "__extra_json",
        "__present_json",
    }:
        raise DatasetError("Samples binary table has an invalid schema")
    for key in SAMPLE_COLUMNS:
        if array.dtype[key].kind != "U" or array.dtype[f"{key}__null"].kind != "b":
            raise DatasetError("Samples binary table has invalid column types")
    rows = []
    for item in array:
        row = json.loads(str(item["__extra_json"]))
        present = json.loads(str(item["__present_json"]))
        if (
            not isinstance(row, dict)
            or set(row) & set(SAMPLE_COLUMNS)
            or not isinstance(present, list)
        ):
            raise DatasetError("Samples binary table has invalid extra metadata")
        if len(set(present)) != len(present) or set(present) - set(SAMPLE_COLUMNS):
            raise DatasetError("Samples binary table has invalid present columns")
        row.update({key: None if item[f"{key}__null"] else str(item[key]) for key in present})
        rows.append(row)
    return tuple(rows)


def _write_csv(path, rows, columns):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: _json(value) if isinstance(value, (dict, list, tuple)) else value
                    for key, value in row.items()
                }
            )


def _write_bundle(dataset, root):
    tables = {}
    sample_path = root / "samples.npy"
    np.save(sample_path, _sample_array(dataset.samples), allow_pickle=False)
    sample_csv = root / "samples.csv"
    columns = list(SAMPLE_COLUMNS) + sorted(
        set().union(*(row.keys() for row in dataset.samples)) - set(SAMPLE_COLUMNS)
    )
    _write_csv(sample_csv, dataset.samples, columns)
    tables["samples"] = {
        "rows": len(dataset.samples),
        "data": _artifact(root, sample_path, "npy"),
        "inspection": _artifact(root, sample_csv, "csv"),
    }
    for name in EVIDENCE_TABLES:
        rows = getattr(dataset, name)
        data, inspection = root / f"{name}.jsonl", root / f"{name}.csv"
        with data.open("w", encoding="utf-8") as stream:
            for row in rows:
                stream.write(_json(row) + "\n")
        columns = ["sample_id"] + sorted(set().union(*(row.keys() for row in rows)) - {"sample_id"})
        _write_csv(inspection, rows, columns)
        tables[name] = {
            "rows": len(rows),
            "data": _artifact(root, data, "jsonl"),
            "inspection": _artifact(root, inspection, "csv"),
        }
    profiles = []
    profile_dir = root / "profiles"
    profile_dir.mkdir()
    for index, matrix in enumerate(dataset.profiles):
        directory = profile_dir / f"{index:04d}"
        directory.mkdir()
        catalogue, vocabulary, codes = (
            directory / name for name in ("catalogue.json", "categories.json", "alleles.npy")
        )
        _write_json(catalogue, asdict(matrix.catalogue))
        _write_json(vocabulary, matrix.categories)
        np.save(codes, matrix.codes, allow_pickle=False)
        export = directory / "profiles.csv"
        with export.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["sample_id", *matrix.catalogue.loci])
            for ident, row_codes in zip(matrix.sample_ids, matrix.codes, strict=True):
                writer.writerow(
                    [
                        ident,
                        *(matrix.categories[int(code) - 1] if code else "0" for code in row_codes),
                    ]
                )
        profiles.append(
            {
                "sample_ids": list(matrix.sample_ids),
                "catalogue": _artifact(root, catalogue, "json"),
                "categories": _artifact(root, vocabulary, "json"),
                "codes": _artifact(root, codes, "npy"),
                "inspection": _artifact(root, export, "csv"),
            }
        )
    manifest = {
        "schema": SCHEMA_NAME,
        "schema_version": SCHEMA_VERSION,
        "software": {"name": "chronoclade", "version": __version__},
        "status": "complete",
        "sample_ids": list(dataset.sample_ids),
        "parameters": dataset.parameters,
        "tables": tables,
        "profiles": profiles,
    }
    manifest["dataset_id"] = "sha256:" + hashlib.sha256(_json(manifest).encode()).hexdigest()
    manifest["created_at"] = datetime.now(timezone.utc).isoformat()
    _write_json(root / "dataset.json", manifest)


def write_dataset(dataset: PreparedDataset, directory: str | Path) -> Path:
    """Publish a new immutable bundle and return its authoritative manifest path.

    A sibling staging directory is renamed only after every artifact and the
    completed manifest have been written and flushed. Existing targets are never
    replaced. Failed writes remove staging data, leaving no completed bundle.
    """
    dataset.validate()
    target = Path(directory).expanduser().absolute()
    if target.exists():
        raise DatasetError(f"Dataset output already exists: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{target.name}.pending-", dir=target.parent))
    try:
        _write_bundle(dataset, staging)
        for path in staging.rglob("*"):
            if path.is_file():
                with path.open("rb") as stream:
                    os.fsync(stream.fileno())
        # Readers can only find the completed directory after this rename.
        staging.rename(target)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return target / "dataset.json"


def _read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def _verified_artifact(root, spec, expected_format, seen):
    if not isinstance(spec, dict) or spec.get("format") != expected_format:
        raise DatasetError("Dataset artifact has an invalid descriptor or format")
    raw = spec.get("path")
    if not isinstance(raw, str) or not raw or "\\" in raw:
        raise DatasetError("Dataset artifacts require relative POSIX paths")
    relative = PurePosixPath(raw)
    if relative.is_absolute() or any(part in {".", ".."} for part in raw.split("/")):
        raise DatasetError("Dataset artifact path escapes the bundle")
    path = (root / raw).resolve()
    if not path.is_relative_to(root) or path in seen or not path.is_file():
        raise DatasetError("Dataset artifact is missing, duplicated, or outside the bundle")
    seen.add(path)
    checksum, size = spec.get("sha256"), spec.get("bytes")
    if not isinstance(checksum, str) or not re.fullmatch(r"[0-9a-f]{64}", checksum):
        raise DatasetError("Dataset artifact requires a SHA256 checksum")
    if (
        type(size) is not int
        or size < 0
        or path.stat().st_size != size
        or _digest(path) != checksum
    ):
        raise DatasetError(f"Dataset artifact checksum/size mismatch: {raw}")
    return path


def load_dataset(manifest_path: str | Path, *, mmap_mode="r") -> PreparedDataset:
    """Validate all checksums and completion before exposing a prepared dataset.

    Matrices default to read-only memory maps. Loading performs no network access
    and never unpickles arrays. Unknown schema versions fail explicitly.
    """
    if mmap_mode not in {None, "r"}:
        raise DatasetError(
            "Prepared datasets support only read-only memory maps or in-memory loading"
        )
    path = Path(manifest_path).expanduser().resolve()
    try:
        manifest = _read_json(path)
        return _load_bundle(path.parent, manifest, mmap_mode)
    except DatasetError:
        raise
    except (OSError, ValueError, TypeError, KeyError, AttributeError) as error:
        raise DatasetError(f"Cannot load prepared dataset: {path}: {error}") from error


def _load_bundle(root, manifest, mmap_mode):
    if not isinstance(manifest, dict) or manifest.get("schema") != SCHEMA_NAME:
        raise DatasetError("Unsupported prepared dataset schema")
    if (
        type(manifest.get("schema_version")) is not int
        or manifest["schema_version"] != SCHEMA_VERSION
    ):
        raise DatasetError("Unsupported prepared dataset schema version")
    if manifest.get("status") != "complete":
        raise DatasetError("Prepared dataset is not complete")
    software = manifest.get("software")
    if (
        not isinstance(software, dict)
        or software.get("name") != "chronoclade"
        or not software.get("version")
    ):
        raise DatasetError("Dataset manifest requires software version")
    identity = {
        key: value for key, value in manifest.items() if key not in {"dataset_id", "created_at"}
    }
    expected = "sha256:" + hashlib.sha256(_json(identity).encode()).hexdigest()
    if manifest.get("dataset_id") != expected:
        raise DatasetError("Dataset manifest identity does not match its contents")
    check_json(manifest, location="manifest")
    seen, rows = set(), {}
    tables = manifest["tables"]
    if not isinstance(tables, dict) or set(tables) != {"samples", *EVIDENCE_TABLES}:
        raise DatasetError("Dataset manifest requires all evidence tables")
    for name, table in tables.items():
        if type(table.get("rows")) is not int or table["rows"] < 0:
            raise DatasetError("Dataset table requires a valid row count")
        data = _verified_artifact(
            root, table["data"], "npy" if name == "samples" else "jsonl", seen
        )
        _verified_artifact(root, table["inspection"], "csv", seen)
        if name == "samples":
            rows[name] = _read_samples(np.load(data, allow_pickle=False, mmap_mode="r"))
        else:
            with data.open(encoding="utf-8") as stream:
                rows[name] = tuple(json.loads(line) for line in stream)
        if len(rows[name]) != table["rows"]:
            raise DatasetError(f"Dataset table row count mismatch: {name}")
    profiles = []
    if not isinstance(manifest["profiles"], list):
        raise DatasetError("Dataset profiles descriptor must be a list")
    for entry in manifest["profiles"]:
        catalogue = _read_json(_verified_artifact(root, entry["catalogue"], "json", seen))
        catalogue["loci"] = tuple(catalogue["loci"])
        categories = _read_json(_verified_artifact(root, entry["categories"], "json", seen))
        if not isinstance(categories, list) or not isinstance(entry["sample_ids"], list):
            raise DatasetError("Dataset profiles require explicit sample IDs and vocabulary")
        codes = np.load(
            _verified_artifact(root, entry["codes"], "npy", seen),
            allow_pickle=False,
            mmap_mode=mmap_mode,
        )
        _verified_artifact(root, entry["inspection"], "csv", seen)
        profiles.append(
            AlleleMatrix(
                LocusCatalogue(**catalogue), tuple(entry["sample_ids"]), tuple(categories), codes
            )
        )
    dataset = PreparedDataset(
        **rows,
        profiles=tuple(profiles),
        parameters=manifest["parameters"],
        dataset_id=manifest["dataset_id"],
    )
    dataset.validate()
    if list(dataset.sample_ids) != manifest["sample_ids"]:
        raise DatasetError("Dataset manifest sample IDs do not match the samples table")
    return dataset
