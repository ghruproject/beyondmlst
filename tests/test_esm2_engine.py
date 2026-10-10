import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from chronoclade.esm2 import EmbeddingError, run_embeddings
from chronoclade.esm2.proteins import Protein, file_sha256, read_proteins, token_batches
from chronoclade.esm2.runtime import choose_device, mean_pool
from chronoclade.esm2.storage import VectorCache, provenance_key, write_npz


class FakeBackend:
    dimension = 320
    load_seconds = 0.01
    fallback_reason = None

    def __init__(self, *, device="cpu", hash_value="a" * 64):
        self.device = device
        self.hash_value = hash_value
        self.batches = []

    @property
    def provenance(self):
        return {"model": "esm2_t6_8M_UR50D", "model_sha256": self.hash_value,
                "precision": "float32", "actual_device": self.device, "pooling": "mean-v1"}

    def synchronize(self):
        pass

    def accelerator_memory(self):
        return None

    def embed_batch(self, proteins):
        self.batches.append([p.id for p in proteins])
        return np.stack([np.arange(self.dimension, dtype=np.float32)
                         + sum(map(ord, p.sequence)) for p in proteins])

    def fallback_to_cpu(self, reason):
        self.device = "cpu"
        self.fallback_reason = reason


def fasta_file(tmp_path, text=">sample1 locus=abc\nMKT\n>sample2 locus=abc\nMKT\n>sample3\nMKTA\n"):
    path = tmp_path / "proteins.faa"
    path.write_text(text)
    return path


def test_deduplication_retains_all_record_headers_and_npz_order(tmp_path):
    fasta = fasta_file(tmp_path)
    backend = FakeBackend()
    result = run_embeddings(fasta, tmp_path / "out", _backend=backend, token_budget=12)
    manifest = result.manifest
    assert manifest["counts"] == {"input_records": 3, "retained_records": 3,
                                  "unique_proteins": 2, "excluded_records": 0}
    assert sum(map(len, backend.batches)) == 2
    assert manifest["mapping"][0]["protein_id"] == manifest["mapping"][1]["protein_id"]
    assert manifest["mapping"][0]["header"] == "sample1 locus=abc"
    assert manifest["input"]["sha256"] == file_sha256(fasta)
    assert manifest["scope"] == "unique-protein-embeddings"
    assert manifest["artifacts"]["embeddings"]["sha256"] == file_sha256(result.vectors_path)
    with np.load(result.vectors_path, allow_pickle=False) as vectors:
        assert list(vectors["protein_ids"]) == manifest["protein_ids"]
        assert vectors["embeddings"].shape == (2, 320)
        assert vectors["embeddings"].dtype == np.float32
        np.testing.assert_array_equal(vectors["lengths"], [3, 4])


@pytest.mark.parametrize("sequence", ["M*K", "MK-", "MK T", "", "M" * 1023])
def test_invalid_or_overlength_sequences_never_truncate(tmp_path, sequence):
    fasta = fasta_file(tmp_path, f">bad\n{sequence}\n>good\nmkt\n")
    with pytest.raises(EmbeddingError, match="no truncation"):
        read_proteins(fasta)
    inputs = read_proteins(fasta, invalid_policy="exclude")
    assert inputs.records == 2
    assert inputs.exclusions[0]["record_id"] == "bad"
    assert inputs.proteins[0].sequence == "MKT"
    assert inputs.mapping[0]["record_id"] == "good"


@pytest.mark.parametrize("text,match", [("MKT\n", "before FASTA"),
                                      (">\nMKT", "Empty FASTA"),
                                      (">same\nMKT\n>same abc\nMKTA", "Duplicate"),
                                      ("", "no records")])
def test_fasta_structure_errors(tmp_path, text, match):
    with pytest.raises(EmbeddingError, match=match):
        read_proteins(fasta_file(tmp_path, text))


def test_padded_token_budget_includes_special_tokens_and_every_protein():
    proteins = [Protein(str(i), "M" * length) for i, length in enumerate([8, 2, 7, 3])]
    batches = list(token_batches(proteins, 18))
    assert sorted(p.id for batch in batches for p in batch) == ["0", "1", "2", "3"]
    assert all(len(batch) * (max(len(p.sequence) for p in batch) + 2) <= 18
               for batch in batches)
    with pytest.raises(EmbeddingError, match="needs 10 tokens"):
        list(token_batches(proteins, 9))


def test_mean_pool_excludes_bos_eos_and_padding():
    # A short protein's EOS and padding are inside the long protein's residue span.
    representations = np.array([[[1000, 1000], [2, 4], [4, 8], [2000, 2000],
                                 [3000, 3000], [4000, 4000]],
                                [[1000, 1000], [1, 3], [3, 5], [5, 7],
                                 [7, 9], [2000, 2000]]], dtype=np.float32)
    pooled = mean_pool(representations, [2, 4])
    np.testing.assert_array_equal(pooled, [[3, 6], [4, 6]])


def test_valid_cache_reuses_vectors_and_changed_model_hash_invalidates(tmp_path):
    fasta = fasta_file(tmp_path)
    cache = tmp_path / "cache"
    first = run_embeddings(fasta, tmp_path / "first", cache_dir=cache, _backend=FakeBackend())
    backend = FakeBackend()
    second = run_embeddings(fasta, tmp_path / "second", cache_dir=cache, _backend=backend)
    assert not backend.batches
    assert second.manifest["runtime"]["cache_hits"] == 2
    assert second.manifest["runtime"]["embedded_proteins"] == 0
    assert first.vectors_path.read_bytes() == second.vectors_path.read_bytes()
    changed = FakeBackend(hash_value="b" * 64)
    third = run_embeddings(fasta, tmp_path / "third", cache_dir=cache, _backend=changed)
    assert sum(map(len, changed.batches)) == 2
    assert third.manifest["runtime"]["cache_hits"] == 0


@pytest.mark.parametrize("mutation", ["bad-zip", "wrong-shape", "wrong-provenance",
                                     "changed-vector", "non-finite"])
def test_corrupted_cache_recomputed(tmp_path, mutation):
    fasta = fasta_file(tmp_path)
    cache_dir = tmp_path / "cache"
    result = run_embeddings(fasta, tmp_path / "first", cache_dir=cache_dir, _backend=FakeBackend())
    cache_file = next((cache_dir / provenance_key(result.manifest["provenance"])).glob("*.npz"))
    if mutation == "bad-zip":
        cache_file.write_bytes(b"PK\x03\x04not-a-zip")
    else:
        with np.load(cache_file, allow_pickle=False) as arrays:
            vector = arrays["vector"].copy()
            info = json.loads(str(arrays["metadata"].item()))
        if mutation == "wrong-shape":
            vector = vector[:10]
        elif mutation == "wrong-provenance":
            info["provenance"]["model_sha256"] = "c" * 64
        elif mutation == "changed-vector":
            vector[0] += 1
        else:
            vector[0] = np.nan
        write_npz(cache_file, vector=vector, metadata=np.asarray(json.dumps(info)))
    backend = FakeBackend()
    rerun = run_embeddings(fasta, tmp_path / "second", cache_dir=cache_dir, _backend=backend)
    assert rerun.manifest["runtime"]["cache_hits"] == 1
    assert rerun.manifest["runtime"]["invalid_cache_entries"] == 1
    assert sum(map(len, backend.batches)) == 1


def test_actual_device_changes_cache_namespace(tmp_path):
    cache = tmp_path / "cache"
    fasta = fasta_file(tmp_path)
    run_embeddings(fasta, tmp_path / "cpu", cache_dir=cache, _backend=FakeBackend())
    backend = FakeBackend(device="mps")
    result = run_embeddings(fasta, tmp_path / "mps", cache_dir=cache, _backend=backend)
    assert result.manifest["actual_device"] == "mps"
    assert result.manifest["runtime"]["cache_hits"] == 0


@pytest.mark.parametrize("kind", ["dimension", "non-finite"])
def test_bad_backend_vectors_never_publish_completed_output(tmp_path, kind):
    backend = FakeBackend()
    backend.embed_batch = lambda batch: (np.zeros((len(batch), 12)) if kind == "dimension"
                                        else np.full((len(batch), 320), np.inf))
    with pytest.raises(EmbeddingError, match="wrong dimensions|non-finite"):
        run_embeddings(fasta_file(tmp_path), tmp_path / "out", _backend=backend)
    assert not (tmp_path / "out" / "embeddings.json").exists()


def test_auto_fallback_restarts_and_records_cpu_without_gpu_cache(tmp_path):
    backend = FakeBackend(device="mps")
    original = backend.embed_batch

    def fail_on_mps(batch):
        if backend.device == "mps":
            raise RuntimeError("unsupported operator")
        return original(batch)

    backend.embed_batch = fail_on_mps
    result = run_embeddings(fasta_file(tmp_path), tmp_path / "out", device="auto",
                            cache_dir=tmp_path / "cache", _backend=backend)
    assert result.manifest["actual_device"] == "cpu"
    assert result.manifest["device_fallback_reason"] == "unsupported operator"
    assert result.manifest["provenance"]["actual_device"] == "cpu"
    assert len(list((tmp_path / "cache").iterdir())) == 1


def test_explicit_device_does_not_fallback(tmp_path):
    backend = FakeBackend(device="mps")

    def fail(batch):
        raise RuntimeError("unsupported operator")

    backend.embed_batch = fail
    with pytest.raises(EmbeddingError, match="unsupported operator"):
        run_embeddings(fasta_file(tmp_path), tmp_path / "out", device="mps", _backend=backend)
    assert backend.device == "mps"


@pytest.mark.parametrize("requested,available,expected", [
    ("auto", (False, False), "cpu"), ("auto", (False, True), "mps"),
    ("auto", (True, True), "cuda"), ("cpu", (True, True), "cpu"),
    ("cuda", (False, True), None), ("mps", (True, False), None),
    ("gpu", (True, True), None),
])
def test_device_choice_is_audited_and_unavailable_explicit_device_fails(requested, available,
                                                                     expected):
    torch = SimpleNamespace(cuda=SimpleNamespace(is_available=lambda: available[0]),
                            backends=SimpleNamespace(mps=SimpleNamespace(
                                is_available=lambda: available[1])))
    if expected is None:
        with pytest.raises(EmbeddingError):
            choose_device(torch, requested)
    else:
        assert choose_device(torch, requested) == expected


def test_failed_manifest_replacement_preserves_previous_complete_run(tmp_path, monkeypatch):
    from chronoclade.esm2 import storage

    fasta = fasta_file(tmp_path)
    out = tmp_path / "out"
    first = run_embeddings(fasta, out, _backend=FakeBackend())
    old_manifest = first.manifest_path.read_bytes()
    fasta.write_text(">changed\nAAAC\n")
    real_replace = storage.os.replace

    def interrupted(source, destination):
        if Path(destination).name == "embeddings.json":
            raise OSError("simulated disk failure")
        return real_replace(source, destination)

    monkeypatch.setattr(storage.os, "replace", interrupted)
    with pytest.raises(OSError, match="disk failure"):
        run_embeddings(fasta, out, _backend=FakeBackend())
    assert first.manifest_path.read_bytes() == old_manifest
    assert file_sha256(first.vectors_path) == first.manifest["artifacts"]["embeddings"]["sha256"]
    assert not list(out.glob(".*.tmp"))


def test_cache_rejects_wrong_id_and_pickle_without_loading_objects(tmp_path):
    cache = VectorCache(tmp_path, {"model": "test"})
    cache.directory.mkdir(parents=True)
    np.savez(cache.directory / "id.npz", vector=np.asarray([object()], dtype=object),
             metadata=np.asarray("{}"))
    assert cache.read("id", 320) is None
    assert cache.invalid == 1
