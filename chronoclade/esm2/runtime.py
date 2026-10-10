"""Lazy fair-esm runtime; loading a local checkpoint never invokes the hub loader."""

from importlib import metadata
from pathlib import Path
import platform
import time
import warnings

from .proteins import EmbeddingError, file_sha256

MODELS = {
    "8M": {"name": "esm2_t6_8M_UR50D", "layers": 6, "dimension": 320},
    "35M": {"name": "esm2_t12_35M_UR50D", "layers": 12, "dimension": 480},
}
POOLING = "residue-mean-excluding-bos-eos-pad-v1"


def model_spec(model: str) -> dict:
    if not isinstance(model, str):
        raise EmbeddingError("Supported models are 8M and 35M")
    for alias, spec in MODELS.items():
        if model.lower() in {alias.lower(), spec["name"].lower()}:
            return dict(spec)
    raise EmbeddingError("Supported models are 8M and 35M")


def choose_device(torch, requested: str) -> str:
    if requested not in {"auto", "cpu", "mps", "cuda"}:
        raise EmbeddingError("device must be auto, cpu, mps or cuda")
    available = {"cpu": True, "cuda": torch.cuda.is_available(),
                 "mps": bool(getattr(torch.backends, "mps", None)
                             and torch.backends.mps.is_available())}
    if requested == "auto":
        return next(name for name in ("cuda", "mps", "cpu") if available[name])
    if not available[requested]:
        raise EmbeddingError(f"Requested device {requested!r} is unavailable")
    return requested


def mean_pool(representations, lengths):
    """Slice exact residue spans: special tokens and padding never enter the mean."""
    return [representations[index, 1:length + 1].mean(axis=0)
            for index, length in enumerate(lengths)]


class ESMBackend:
    """Small ESM2 models with fixed float32 extraction and explicit weight provenance."""

    def __init__(self, *, model="8M", checkpoint=None, allow_download=False, device="auto",
                 expected_checkpoint_sha256=None):
        started = time.perf_counter()
        spec = model_spec(model)
        if checkpoint is None and not allow_download:
            raise EmbeddingError("Supply a local checkpoint or explicitly enable model download")
        if checkpoint is not None and not Path(checkpoint).is_file():
            raise EmbeddingError(f"Checkpoint does not exist: {checkpoint}")
        try:
            import torch
            import esm
        except ImportError as exc:
            raise EmbeddingError("ESM2 requires the optional torch and fair-esm dependencies") from exc
        self.torch = torch
        self.requested_device = device
        self.device = choose_device(torch, device)
        checkpoint_source = "explicit-local" if checkpoint is not None else "opt-in-hub"
        if checkpoint is None:
            checkpoint = Path(torch.hub.get_dir()) / "checkpoints" / f"{spec['name']}.pt"
            checkpoint.parent.mkdir(parents=True, exist_ok=True)
            if not checkpoint.is_file():
                # This branch is reachable only after explicit allow_download=True.
                url = f"https://dl.fbaipublicfiles.com/fair-esm/models/{spec['name']}.pt"
                torch.hub.download_url_to_file(url, str(checkpoint), progress=True)
        checkpoint = Path(checkpoint)
        digest = file_sha256(checkpoint)
        if expected_checkpoint_sha256 and digest != expected_checkpoint_sha256.lower():
            raise EmbeddingError("Checkpoint SHA256 does not match the expected hash")
        try:
            # Official ESM2 files contain configuration objects as well as tensor weights.
            # Only checkpoints acquired from a trusted source should be opened.
            data = torch.load(checkpoint, map_location="cpu", weights_only=False)
            cfg = data["cfg"]["model"]
            if (cfg.encoder_layers, cfg.encoder_embed_dim) != (spec["layers"], spec["dimension"]):
                raise EmbeddingError("Checkpoint architecture does not match the requested model")
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="Regression weights not found.*")
                loaded_model, alphabet = esm.pretrained.load_model_and_alphabet_core(
                    spec["name"], data, regression_data=None
                )
        except EmbeddingError:
            raise
        except Exception as exc:
            raise EmbeddingError(f"Could not load ESM2 checkpoint {checkpoint}: {exc}") from exc
        self.model = loaded_model.eval().float()
        self.alphabet = alphabet
        if not alphabet.prepend_bos or not alphabet.append_eos:
            raise EmbeddingError("ESM2 alphabet must supply BOS and EOS tokens")
        self.dimension = spec["dimension"]
        self.layer = spec["layers"]
        self.converter = alphabet.get_batch_converter()
        self.fallback_reason = None
        try:
            self.model.to(self.device)
        except RuntimeError as exc:
            if device != "auto" or self.device == "cpu":
                raise EmbeddingError(f"Cannot load model on {self.device}: {exc}") from exc
            self.fallback_to_cpu(str(exc))
        self.synchronize()
        self.load_seconds = time.perf_counter() - started
        self._provenance = {
            "model": spec["name"], "model_sha256": digest, "dimension": self.dimension,
            "representation_layer": self.layer, "pooling": POOLING,
            "precision": "float32", "torch": torch.__version__,
            "fair_esm": metadata.version("fair-esm"),
            "python": platform.python_version(), "platform": platform.platform(),
            "cpu_threads": torch.get_num_threads(),
            "checkpoint_source": checkpoint_source,
        }

    @property
    def provenance(self):
        device_name = platform.machine()
        if self.device == "cuda":
            device_name = self.torch.cuda.get_device_name()
        return {**self._provenance, "actual_device": self.device, "device_name": device_name}

    def synchronize(self):
        if self.device == "cuda":
            self.torch.cuda.synchronize()
        elif self.device == "mps":
            self.torch.mps.synchronize()

    def fallback_to_cpu(self, reason):
        self.model.to("cpu")
        self.device = "cpu"
        self.fallback_reason = reason

    def embed_batch(self, proteins):
        _, _, tokens = self.converter([(protein.id, protein.sequence) for protein in proteins])
        if tokens.shape[1] != max(len(p.sequence) for p in proteins) + 2:
            raise EmbeddingError("Tokenizer changed sequence lengths; refusing truncated extraction")
        with self.torch.inference_mode():
            result = self.model(tokens.to(self.device), repr_layers=[self.layer],
                                return_contacts=False)["representations"][self.layer]
            pooled = mean_pool(result, [len(p.sequence) for p in proteins])
            vectors = self.torch.stack(pooled).float().cpu().numpy()
        return vectors

    def accelerator_memory(self):
        if self.device == "cuda":
            return {"allocated_bytes": self.torch.cuda.memory_allocated(),
                    "peak_allocated_bytes": self.torch.cuda.max_memory_allocated()}
        if self.device == "mps":
            return {"allocated_bytes": self.torch.mps.current_allocated_memory(),
                    "driver_allocated_bytes": self.torch.mps.driver_allocated_memory()}
        return None
