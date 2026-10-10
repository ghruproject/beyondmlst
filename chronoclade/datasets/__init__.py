"""Provider-independent, versioned prepared dataset persistence.

Use explicit LocusCatalogue objects with from_profile_records for frozen legacy
records, or construct PreparedDataset and AlleleMatrix directly. write_dataset
publishes a new immutable bundle; load_dataset validates its completed manifest.
"""

from .model import (
    AlleleMatrix,
    DatasetError,
    LocusCatalogue,
    PreparedDataset,
    SCHEMA_NAME,
    SCHEMA_VERSION,
)
from .adapters import from_profile_records, samples_from_records
from .storage import load_dataset, write_dataset

__all__ = [
    "AlleleMatrix",
    "DatasetError",
    "LocusCatalogue",
    "PreparedDataset",
    "SCHEMA_NAME",
    "SCHEMA_VERSION",
    "from_profile_records",
    "load_dataset",
    "samples_from_records",
    "write_dataset",
]
