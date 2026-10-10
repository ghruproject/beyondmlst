from contextlib import nullcontext
import sys
from types import SimpleNamespace

import pytest

from chronoclade.esm2 import EmbeddingError
from chronoclade.esm2.proteins import file_sha256
from chronoclade.esm2.runtime import ESMBackend, model_spec


@pytest.fixture
def fake_runtime(tmp_path, monkeypatch):
    checkpoint = tmp_path / "renamed-local-weights.pt"
    checkpoint.write_bytes(b"test checkpoint contents")
    calls = []
    cfg = SimpleNamespace(encoder_layers=6, encoder_embed_dim=320)

    def forbidden_download(*args, **kwargs):
        raise AssertionError("Offline loading must not invoke any downloader")

    def load(path, *, map_location, weights_only):
        calls.append((path, map_location, weights_only))
        return {"cfg": {"model": cfg}, "model": {}}

    model = SimpleNamespace()
    model.eval = lambda: model
    model.float = lambda: model
    model.to = lambda device: model
    alphabet = SimpleNamespace(prepend_bos=True, append_eos=True,
                               get_batch_converter=lambda: "converter")
    torch = SimpleNamespace(load=load, __version__="test", get_num_threads=lambda: 1,
                            cuda=SimpleNamespace(is_available=lambda: False),
                            backends=SimpleNamespace(mps=SimpleNamespace(is_available=lambda: False)),
                            hub=SimpleNamespace(download_url_to_file=forbidden_download,
                                                get_dir=forbidden_download),
                            inference_mode=nullcontext)

    def core(name, data, regression_data):
        calls.append((name, regression_data))
        return model, alphabet

    esm = SimpleNamespace(pretrained=SimpleNamespace(load_model_and_alphabet_core=core,
                                                    load_model_and_alphabet_local=forbidden_download,
                                                    load_model_and_alphabet_hub=forbidden_download))
    monkeypatch.setitem(sys.modules, "torch", torch)
    monkeypatch.setitem(sys.modules, "esm", esm)
    monkeypatch.setattr("chronoclade.esm2.runtime.metadata.version", lambda name: "test")
    return checkpoint, cfg, calls


def test_explicit_local_checkpoint_load_never_uses_hub_or_regression_sidecar(fake_runtime):
    checkpoint, _, calls = fake_runtime
    backend = ESMBackend(checkpoint=checkpoint, device="cpu")
    assert calls[0] == (checkpoint, "cpu", False)
    assert calls[1] == ("esm2_t6_8M_UR50D", None)
    assert backend.provenance["model_sha256"] == file_sha256(checkpoint)
    assert backend.provenance["actual_device"] == "cpu"


def test_wrong_model_architecture_fails_before_core_loader(fake_runtime):
    checkpoint, cfg, calls = fake_runtime
    cfg.encoder_layers, cfg.encoder_embed_dim = 12, 480
    with pytest.raises(EmbeddingError, match="architecture"):
        ESMBackend(checkpoint=checkpoint, model="8M")
    assert len(calls) == 1


def test_wrong_checkpoint_hash_fails_before_deserialisation(fake_runtime):
    checkpoint, _, calls = fake_runtime
    with pytest.raises(EmbeddingError, match="SHA256"):
        ESMBackend(checkpoint=checkpoint, expected_checkpoint_sha256="0" * 64)
    assert calls == []


def test_download_requires_explicit_opt_in():
    with pytest.raises(EmbeddingError, match="explicitly enable"):
        ESMBackend()


def test_nonexistent_checkpoint_fails_before_optional_imports(tmp_path):
    with pytest.raises(EmbeddingError, match="does not exist"):
        ESMBackend(checkpoint=tmp_path / "absent.pt")


def test_small_model_names():
    assert model_spec("8m")["dimension"] == 320
    assert model_spec("esm2_t12_35M_UR50D")["dimension"] == 480
    with pytest.raises(EmbeddingError, match="Supported models"):
        model_spec("650M")
    with pytest.raises(EmbeddingError, match="Supported models"):
        model_spec(None)
