"""Optional experimental protein embeddings. Torch and fair-esm are loaded on demand."""

from .engine import EmbeddingResult, run_embeddings
from .proteins import EmbeddingError

__all__ = ["EmbeddingError", "EmbeddingResult", "run_embeddings"]
