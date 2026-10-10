"""Protein-level ESM2 extraction; no genome aggregation or selection claims."""

from dataclasses import dataclass
import hashlib
from pathlib import Path
import platform
import time

import numpy as np

from chronoclade import __version__
from .proteins import EmbeddingError, file_sha256, read_proteins, token_batches
from .runtime import ESMBackend, model_spec
from .storage import VectorCache, write_json, write_npz


@dataclass(frozen=True)
class EmbeddingResult:
    manifest_path: Path
    vectors_path: Path
    manifest: dict


def _extract(inputs, backend, token_budget, cache_dir):
    provenance = {**backend.provenance, "token_budget": token_budget}
    cache = VectorCache(cache_dir, provenance) if cache_dir is not None else None
    vectors, missing = {}, []
    for protein in inputs.proteins:
        cached = cache.read(protein.id, backend.dimension) if cache else None
        if cached is None:
            missing.append(protein)
        else:
            vectors[protein.id] = cached
    reused = len(vectors)
    backend.synchronize()
    started = time.perf_counter()
    batch_count, padded_tokens = 0, 0
    for batch in token_batches(missing, token_budget):
        result = np.asarray(backend.embed_batch(batch))
        if result.shape != (len(batch), backend.dimension):
            raise EmbeddingError("Backend returned vectors with the wrong dimensions")
        if not np.isfinite(result).all():
            raise EmbeddingError("Backend returned non-finite vectors")
        result = result.astype(np.float32)
        for protein, vector in zip(batch, result, strict=True):
            vectors[protein.id] = vector
        batch_count += 1
        padded_tokens += len(batch) * (max(len(p.sequence) for p in batch) + 2)
    backend.synchronize()
    inference_seconds = time.perf_counter() - started if missing else 0.0
    # Only publish cache vectors after the complete extraction succeeds. A device
    # fallback cannot mix partial GPU vectors into the CPU provenance namespace.
    if cache:
        for protein in missing:
            cache.write(protein.id, vectors[protein.id])
    stats = {"cache_hits": reused, "embedded_proteins": len(missing),
             "invalid_cache_entries": cache.invalid if cache else 0,
             "batches": batch_count, "padded_tokens": padded_tokens,
             "inference_seconds": inference_seconds,
             "residues_embedded": sum(len(p.sequence) for p in missing)}
    return np.stack([vectors[protein.id] for protein in inputs.proteins]), provenance, stats


def run_embeddings(fasta: Path, out: Path, *, model="8M", checkpoint=None,
                   allow_download=False, device="auto", token_budget=4096,
                   max_length=1022, cache_dir=None, invalid_policy="error",
                   expected_checkpoint_sha256=None, _backend=None) -> EmbeddingResult:
    """Embed unique proteins, saving vectors, all source mappings and audited exclusions.

    Local checkpoint loading is offline. ``allow_download`` must be explicitly set
    to acquire standard weights. ``_backend`` is a private test/benchmark seam.
    """
    started = time.perf_counter()
    spec = model_spec(model)
    if device not in {"auto", "cpu", "mps", "cuda"}:
        raise EmbeddingError("device must be auto, cpu, mps or cuda")
    inputs = read_proteins(fasta, max_length=max_length, invalid_policy=invalid_policy)
    # Validate every protein's resource limit before optional dependencies or weights load.
    list(token_batches(inputs.proteins, token_budget))
    try:
        backend = _backend or ESMBackend(model=model, checkpoint=checkpoint,
                                        allow_download=allow_download, device=device,
                                        expected_checkpoint_sha256=expected_checkpoint_sha256)
        if backend.dimension != spec["dimension"]:
            raise EmbeddingError("Backend dimension does not match the requested model")
        try:
            vectors, provenance, stats = _extract(inputs, backend, token_budget, cache_dir)
        except RuntimeError as exc:
            if device != "auto" or backend.device == "cpu":
                raise
            backend.fallback_to_cpu(str(exc))
            vectors, provenance, stats = _extract(inputs, backend, token_budget, cache_dir)
    except EmbeddingError:
        raise
    except (RuntimeError, OSError) as exc:
        raise EmbeddingError(f"ESM2 extraction failed: {exc}") from exc
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    # An immutable artifact name keeps an existing manifest valid if publication
    # of the next run is interrupted before its final manifest replacement.
    identity = hashlib.sha256(vectors.tobytes() + "".join(p.id for p in inputs.proteins).encode())
    vectors_path = out / f"embeddings.{identity.hexdigest()[:20]}.npz"
    write_npz(vectors_path, protein_ids=np.asarray([p.id for p in inputs.proteins]),
              embeddings=vectors, lengths=np.asarray([len(p.sequence) for p in inputs.proteins]))
    lengths = [len(p.sequence) for p in inputs.proteins]
    manifest = {
        "schema": "chronoclade.esm2.protein-embeddings", "schema_version": 1,
        "software_version": __version__, "status": "complete", "experimental": True,
        "scope": "unique-protein-embeddings",
        "input": {"path": str(Path(fasta).resolve()), "sha256": inputs.input_sha256,
                  "format": "protein-fasta", "records": inputs.records},
        "parameters": {"requested_model": model, "requested_device": device,
                       "token_budget": token_budget, "max_length": max_length,
                       "invalid_policy": invalid_policy, "truncation": False,
                       "sequence_normalisation": "uppercase",
                       "protein_id_method": "sha256-normalised-amino-acids"},
        "provenance": provenance,
        "actual_device": backend.device,
        "device_fallback_reason": backend.fallback_reason,
        "protein_ids": [p.id for p in inputs.proteins],
        "mapping": list(inputs.mapping), "exclusions": list(inputs.exclusions),
        "counts": {"input_records": inputs.records, "retained_records": len(inputs.mapping),
                   "unique_proteins": len(inputs.proteins), "excluded_records": len(inputs.exclusions)},
        "sequence_lengths": {"minimum": min(lengths), "maximum": max(lengths),
                             "mean": float(np.mean(lengths)), "median": float(np.median(lengths))},
        "artifacts": {"embeddings": {"path": vectors_path.name,
                                     "sha256": file_sha256(vectors_path),
                                     "shape": list(vectors.shape), "dtype": "float32"}},
        "runtime": {"python": platform.python_version(), "platform": platform.platform(),
                    "model_load_seconds": backend.load_seconds if _backend is None else 0.0, **stats,
                    "accelerator_memory": backend.accelerator_memory(),
                    "elapsed_seconds": time.perf_counter() - started},
    }
    manifest_path = out / "embeddings.json"
    write_json(manifest_path, manifest)
    return EmbeddingResult(manifest_path, vectors_path, manifest)
